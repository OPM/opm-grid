/*
  Copyright 2026 SINTEF Digital.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <config.h>

#define BOOST_TEST_MODULE ConformingBuilderTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/refinement/conforming/ConformingBlockBuilder.hpp>
#include <opm/grid/cpgrid/refinement/conforming/GrdeclRefinement.hpp>
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <dune/common/parallel/mpihelper.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <vector>

namespace
{

struct MPIFixture
{
    MPIFixture()
    {
        auto& argv = boost::unit_test::framework::master_test_suite().argv;
        auto& argc = boost::unit_test::framework::master_test_suite().argc;
        Dune::MPIHelper::instance(argc, argv);
    }
};

struct TestGrdecl
{
    std::array<int,3> dims;
    std::vector<double> coord;
    std::vector<double> zcorn;
    std::vector<int> actnum;

    grdecl raw() const
    {
        grdecl g;
        g.dims[0] = dims[0]; g.dims[1] = dims[1]; g.dims[2] = dims[2];
        g.coord = coord.data();
        g.zcorn = zcorn.data();
        g.actnum = actnum.empty() ? nullptr : actnum.data();
        return g;
    }
};

int cellOf(int doubled) { return doubled / 2; }
int sideOf(int doubled) { return doubled % 2; }

TestGrdecl makeVerticalPillarGrid(const std::array<int,3>& dims,
                                  const std::function<double(int,int,int)>& depth)
{
    TestGrdecl g;
    g.dims = dims;
    const auto& [nx, ny, nz] = dims;

    g.coord.resize(6*(nx + 1)*(ny + 1));
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            double* p = &g.coord[6*(static_cast<std::size_t>(j)*(nx + 1) + i)];
            p[0] = i; p[1] = j; p[2] = 0.0;
            p[3] = i; p[4] = j; p[5] = 100.0;
        }
    }
    g.zcorn.resize(8*static_cast<std::size_t>(nx)*ny*nz);
    for (int k_ = 0; k_ < 2*nz; ++k_) {
        for (int j_ = 0; j_ < 2*ny; ++j_) {
            for (int i_ = 0; i_ < 2*nx; ++i_) {
                g.zcorn[static_cast<std::size_t>(i_) + 2*static_cast<std::size_t>(nx)*j_
                        + 4*static_cast<std::size_t>(nx)*ny*k_] = depth(i_, j_, k_);
            }
        }
    }
    return g;
}

void useConformingBuilder(Dune::CpGrid& grid, std::unique_ptr<Opm::Refinement::Builder> builder)
{
    grid.setLgrBackend(Opm::Refinement::Backend::Conforming);
    grid.setRefinementBuilder(std::move(builder));
}

double totalVolume(const Dune::CpGrid& grid)
{
    double vol = 0.0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        vol += element.geometry().volume();
    }
    return vol;
}

// Every leaf face's corners wind counter-clockwise around its stored normal: the
// polygon's area vector (sum of corner cross products) points along it.
void checkFacesWindAroundNormal(const Dune::CpGrid& grid)
{
    const auto& leaf = *grid.currentData().back();
    int inverted = 0;
    for (int face = 0; face < leaf.numFaces(); ++face) {
        const auto nodes = leaf.faceToPoint(face);
        const auto count = static_cast<int>(nodes.size());
        if (count < 3) {
            continue;
        }
        Dune::FieldVector<double,3> area(0.0);
        for (int n = 0; n < count; ++n) {
            const auto a = Dune::cpgrid::Entity<3>(leaf, nodes[n], true).geometry().center();
            const auto b = Dune::cpgrid::Entity<3>(leaf, nodes[(n + 1) % count], true).geometry().center();
            area[0] += a[1]*b[2] - a[2]*b[1];
            area[1] += a[2]*b[0] - a[0]*b[2];
            area[2] += a[0]*b[1] - a[1]*b[0];
        }
        if (area.dot(leaf.faceNormals(face)) <= 0.0) {
            ++inverted;
        }
    }
    BOOST_CHECK_EQUAL(inverted, 0);
}

// Structural validity of a leaf with a faulted box boundary (phase 1): one
// refined level, volume conserved, every cell closed (sum of area*outward-normal
// is zero), and the box connected to coarse neighbours.
void checkValidFaultedLeaf(const Dune::CpGrid& grid, double volumeBefore)
{
    BOOST_CHECK_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
    int refinedToCoarse = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        Dune::FieldVector<double,3> closure(0.0);
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            auto n = is.centerUnitOuterNormal();
            n *= is.geometry().volume();
            closure += n;
            if (is.neighbor() && is.inside().hasFather() != is.outside().hasFather()) {
                ++refinedToCoarse;
            }
        }
        BOOST_CHECK_SMALL(closure.two_norm(), 1e-9);
    }
    BOOST_CHECK_GT(refinedToCoarse, 0);
    checkFacesWindAroundNormal(grid);
}

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_AUTO_TEST_CASE(endToEndSingleBox)
{
    // Distorted vertical-pillar 4x3x3 grid.
    auto depth = [](int i_, int j_, int k_) {
        const double x = cellOf(i_) + sideOf(i_);
        const double y = cellOf(j_) + sideOf(j_);
        return 2.0*(cellOf(k_) + sideOf(k_)) + 0.3*std::sin(0.9*x) + 0.2*std::cos(0.6*y);
    };
    auto parent = makeVerticalPillarGrid({4, 3, 3}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);
    const int cellsBefore = grid.size(0);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // Box of 2x1x1 parents, refined 2x2x2.
    grid.addLgrsUpdateLeafView({{2,2,2}}, {{1,1,1}}, {{3,2,2}}, {"LGR1"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_EQUAL(grid.size(0), cellsBefore - 2 + 2*8);
    BOOST_CHECK_EQUAL(grid.size(0, 0), cellsBefore);    // level 0 view
    BOOST_CHECK_EQUAL(grid.size(1, 0), 16);             // refined level view
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
    BOOST_CHECK_EQUAL(grid.getLgrNameToLevel().at("LGR1"), 1);

    // Parent-child structure through the leaf view.
    std::map<int, double> volumePerParent;
    int refinedCount = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        if (element.hasFather()) {
            ++refinedCount;
            const auto father = element.father();
            BOOST_CHECK_EQUAL(father.level(), 0);
            volumePerParent[father.index()] += element.geometry().volume();
            BOOST_CHECK_CLOSE(element.geometryInFather().volume(), 1.0/8.0, 1e-9);
            // Leaf global_cell_ points at the parent's Cartesian index.
            BOOST_CHECK_EQUAL(grid.globalCell()[element.index()],
                              grid.currentData()[0]->globalCell()[father.index()]);
        }
    }
    BOOST_CHECK_EQUAL(refinedCount, 16);
    BOOST_CHECK_EQUAL(volumePerParent.size(), 2u);
    for (const auto& [parentIdx, childVolume] : volumePerParent) {
        const Dune::cpgrid::Entity<0> parentCell(*grid.currentData()[0], parentIdx, true);
        BOOST_CHECK_CLOSE(childVolume, parentCell.geometry().volume(), 1e-8);
    }

    // Intersection consistency over the whole leaf: every interior
    // intersection must be seen from both sides with matching geometry.
    std::map<std::pair<int,int>, int> pairCount;
    double interiorArea = 0.0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& intersection : Dune::intersections(grid.leafGridView(), element)) {
            if (intersection.neighbor()) {
                const int inside = intersection.inside().index();
                const int outside = intersection.outside().index();
                BOOST_REQUIRE(inside != outside);
                pairCount[{std::min(inside, outside), std::max(inside, outside)}] += 1;
                interiorArea += intersection.geometry().volume();
            }
        }
    }
    for (const auto& [cells, count] : pairCount) {
        BOOST_CHECK_EQUAL(count, 2); // seen once from each side
    }

    // The coarse i-neighbor of the box (cell (0,1,1)) sees the refined
    // mosaic: 5 unrefined faces + 4 refined ones.
    const auto& dims = grid.logicalCartesianSize();
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        if (element.hasFather()) {
            continue;
        }
        const int cart = grid.globalCell()[element.index()];
        const int i = cart % dims[0];
        const int j = (cart / dims[0]) % dims[1];
        const int k = cart / (dims[0]*dims[1]);
        if (i == 0 && j == 1 && k == 1) {
            int intersections = 0;
            int refinedNeighbors = 0;
            for (const auto& intersection : Dune::intersections(grid.leafGridView(), element)) {
                ++intersections;
                if (intersection.neighbor() && intersection.outside().hasFather()) {
                    ++refinedNeighbors;
                }
            }
            BOOST_CHECK_EQUAL(intersections, 9);
            BOOST_CHECK_EQUAL(refinedNeighbors, 4);
        }
    }

    // Global ids: unique over leaf cells and over leaf points.
    const auto& ids = grid.globalIdSet();
    std::set<std::int64_t> cellIds;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        cellIds.insert(ids.id(element));
    }
    BOOST_CHECK_EQUAL(cellIds.size(), static_cast<std::size_t>(grid.size(0)));

    std::set<std::int64_t> pointIds;
    for (const auto& vertex : Dune::vertices(grid.leafGridView())) {
        pointIds.insert(ids.id(vertex));
    }
    BOOST_CHECK_EQUAL(pointIds.size(), static_cast<std::size_t>(grid.size(3)));
}

BOOST_AUTO_TEST_CASE(stableCellIdDisambiguatesRefinedSiblings)
{
    // global_cell_ is NOT unique on the leaf -- refined siblings share the
    // parent's Cartesian index -- so it cannot key per-cell output gathered
    // across ranks. stableCellId() must give every leaf cell a unique id:
    // coarse cells keep their Cartesian index; refined cells encode
    // (parent Cartesian, child index) under a tag bit.
    auto parent = makeVerticalPillarGrid({4, 3, 3}, [](int, int, int k_) {
        return 2.0*(cellOf(k_) + sideOf(k_));
    });

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));
    grid.addLgrsUpdateLeafView({{2,2,2}}, {{1,1,1}}, {{3,2,2}}, {"LGR1"}); // 2 parents -> 16 refined

    constexpr std::int64_t tag = std::int64_t(1) << 62;
    constexpr std::int64_t childMask = (std::int64_t(1) << 20) - 1;

    const auto sid = grid.stableCellId();
    BOOST_REQUIRE_EQUAL(sid.size(), static_cast<std::size_t>(grid.size(0)));

    // global_cell_ really does collide for siblings, but stableCellId does not.
    const auto& gc = grid.globalCell();
    BOOST_CHECK_LT(std::set<int>(gc.begin(), gc.end()).size(), gc.size());
    BOOST_CHECK_EQUAL(std::set<std::int64_t>(sid.begin(), sid.end()).size(), sid.size());

    std::map<int, std::set<int>> childrenOfParent; // father Cartesian -> child indices seen
    int refined = 0, coarse = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        const int c = element.index();
        if (element.hasFather()) {
            ++refined;
            const int fatherCart = grid.currentData()[0]->globalCell()[element.father().index()];
            BOOST_CHECK((sid[c] & tag) != 0);                       // tagged refined
            BOOST_CHECK_EQUAL((sid[c] & ~tag) >> 20, fatherCart);    // decodes to parent Cartesian
            childrenOfParent[fatherCart].insert(static_cast<int>(sid[c] & childMask));
        }
        else {
            ++coarse;
            BOOST_CHECK_EQUAL(sid[c] & tag, 0);                      // untagged coarse
            BOOST_CHECK_EQUAL(sid[c], static_cast<std::int64_t>(gc[c]));
        }
    }
    BOOST_CHECK_EQUAL(refined, 16);
    BOOST_CHECK_GT(coarse, 0);
    // Each of the two refined parents has all 8 distinct child indices 0..7.
    BOOST_CHECK_EQUAL(childrenOfParent.size(), 2u);
    for (const auto& [father, kids] : childrenOfParent) {
        BOOST_CHECK_EQUAL(kids.size(), 8u);
        BOOST_CHECK_EQUAL(*kids.begin(), 0);
        BOOST_CHECK_EQUAL(*kids.rbegin(), 7);
    }

    // Construction-stable: recomputing yields the identical ids.
    BOOST_CHECK(grid.stableCellId() == sid);
}

BOOST_AUTO_TEST_CASE(twoSeparatedBoxesAndGuards)
{
    auto parent = makeVerticalPillarGrid({6, 3, 3}, [](int, int, int k_) {
        return 2.0*(cellOf(k_) + sideOf(k_));
    });

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // Overlapping boxes are rejected (validation).
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView({{2,2,2}, {2,2,2}},
                                                 {{0,0,0}, {1,0,0}},
                                                 {{3,2,2}, {4,2,2}},
                                                 {"A", "B"}),
                      std::invalid_argument);
    BOOST_CHECK_EQUAL(grid.maxLevel(), 0); // grid unchanged

    // Two separated boxes work.
    grid.addLgrsUpdateLeafView({{2,2,2}, {3,3,3}},
                               {{0,0,0}, {4,1,1}},
                               {{1,1,1}, {6,3,3}},
                               {"A", "B"});
    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    BOOST_CHECK_EQUAL(grid.getLgrNameToLevel().at("A"), 1);
    BOOST_CHECK_EQUAL(grid.getLgrNameToLevel().at("B"), 2);
    // 54 - 1 - 8 parents + 8 + 8*27 children
    BOOST_CHECK_EQUAL(grid.size(0), 54 - 1 - 8 + 8 + 8*27);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);

    std::set<std::int64_t> cellIds;
    const auto& ids = grid.globalIdSet();
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        cellIds.insert(ids.id(element));
    }
    BOOST_CHECK_EQUAL(cellIds.size(), static_cast<std::size_t>(grid.size(0)));
}

// Two diagonally-adjacent boxes share only an edge (the CARFIN deck
// configuration). The refined corners along the shared edge that are
// interior to a parent cell must merge into single leaf vertices.
BOOST_AUTO_TEST_CASE(edgeSharingBoxesMergeCorners)
{
    auto parent = makeVerticalPillarGrid({4, 4, 2}, [](int, int, int k_) {
        return 2.0*(cellOf(k_) + sideOf(k_));
    });

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // Box A = i,j in [0,2); box B = i,j in [2,4); both k in [0,2).
    // They meet along the vertical edge (i=2, j=2). Factors 2x2x2.
    grid.addLgrsUpdateLeafView({{2,2,2}, {2,2,2}},
                               {{0,0,0}, {2,2,0}},
                               {{2,2,2}, {4,4,2}},
                               {"A", "B"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
    // 32 - 2*8 parents + 2*64 children
    BOOST_CHECK_EQUAL(grid.size(0), 32 - 2*8 + 2*64);

    // No two distinct leaf vertices may sit at the same coordinate: the
    // shared-edge interior corners must have merged.
    std::set<std::array<double,3>> coords;
    int vertexCount = 0;
    for (const auto& vertex : Dune::vertices(grid.leafGridView())) {
        const auto& c = vertex.geometry().center();
        coords.insert({c[0], c[1], c[2]});
        ++vertexCount;
    }
    BOOST_CHECK_EQUAL(coords.size(), static_cast<std::size_t>(vertexCount));
    BOOST_CHECK_EQUAL(vertexCount, grid.size(3));

    // Two-sided intersection symmetry over the whole leaf.
    std::map<std::pair<int,int>, int> pairCount;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& intersection : Dune::intersections(grid.leafGridView(), element)) {
            if (intersection.neighbor()) {
                const int in = intersection.inside().index();
                const int out = intersection.outside().index();
                pairCount[{std::min(in, out), std::max(in, out)}] += 1;
            }
        }
    }
    for (const auto& [cells, count] : pairCount) {
        BOOST_CHECK_EQUAL(count, 2);
    }

    // Unique global ids for cells and points.
    const auto& ids = grid.globalIdSet();
    std::set<std::int64_t> cellIds, pointIds;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        cellIds.insert(ids.id(element));
    }
    for (const auto& vertex : Dune::vertices(grid.leafGridView())) {
        pointIds.insert(ids.id(vertex));
    }
    BOOST_CHECK_EQUAL(cellIds.size(), static_cast<std::size_t>(grid.size(0)));
    BOOST_CHECK_EQUAL(pointIds.size(), static_cast<std::size_t>(grid.size(3)));
}

// Face-sharing boxes: the boundary faces of the two refined blocks on the
// shared plane pair into interior faces. The out-of-face subdivision may
// differ between the two boxes; the in-face subdivisions must match.
BOOST_AUTO_TEST_CASE(faceSharingBoxes)
{
    auto parent = makeVerticalPillarGrid({4, 2, 2}, [](int, int, int k_) {
        return 2.0*(cellOf(k_) + sideOf(k_));
    });

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // Box A = i in [0,2), box B = i in [2,4); both full in j,k -> share the
    // i=2 face. Different in-plane-normal factor (rx 2 vs 3), matching
    // in-face factors (ry=rz=2).
    grid.addLgrsUpdateLeafView({{2,2,2}, {3,2,2}},
                               {{0,0,0}, {2,0,0}},
                               {{2,2,2}, {4,2,2}},
                               {"A", "B"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
    // Box A (8 parents x 2x2x2) + box B (8 parents x 3x2x2); no coarse left.
    BOOST_CHECK_EQUAL(grid.size(0), 8*8 + 8*12);

    // No coincident-distinct vertices (shared-face corners merged).
    std::set<std::array<double,3>> coords;
    int vertexCount = 0;
    for (const auto& vertex : Dune::vertices(grid.leafGridView())) {
        const auto& c = vertex.geometry().center();
        coords.insert({c[0], c[1], c[2]});
        ++vertexCount;
    }
    BOOST_CHECK_EQUAL(coords.size(), static_cast<std::size_t>(vertexCount));

    // Two-sided intersection symmetry, and the A/B blocks are actually
    // connected (interior faces pairing an A-child with a B-child exist).
    std::map<std::pair<int,int>, int> pairCount;
    int abConnections = 0;
    const auto& dims = grid.logicalCartesianSize();
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& intersection : Dune::intersections(grid.leafGridView(), element)) {
            if (!intersection.neighbor()) {
                continue;
            }
            const int in = intersection.inside().index();
            const int out = intersection.outside().index();
            pairCount[{std::min(in, out), std::max(in, out)}] += 1;
            if (intersection.inside().hasFather() && intersection.outside().hasFather()) {
                const int ci = grid.globalCell()[in] % dims[0];
                const int co = grid.globalCell()[out] % dims[0];
                // One side in box A's i-range [0,2), the other in B's [2,4).
                if ((ci < 2) != (co < 2)) {
                    ++abConnections;
                }
            }
        }
    }
    for (const auto& [cells, count] : pairCount) {
        BOOST_CHECK_EQUAL(count, 2);
    }
    // The shared i=2 face carries ry*rz = 2*2 = 4 connections per parent
    // pair, 2x2 parent pairs in (j,k) -> 16 connections, counted twice.
    BOOST_CHECK_EQUAL(abConnections, 16 * 2);

    const auto& ids = grid.globalIdSet();
    std::set<std::int64_t> cellIds, pointIds;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        cellIds.insert(ids.id(element));
    }
    for (const auto& vertex : Dune::vertices(grid.leafGridView())) {
        pointIds.insert(ids.id(vertex));
    }
    BOOST_CHECK_EQUAL(cellIds.size(), static_cast<std::size_t>(grid.size(0)));
    BOOST_CHECK_EQUAL(pointIds.size(), static_cast<std::size_t>(grid.size(3)));
}

// Non-matching in-face subdivisions are rejected.
BOOST_AUTO_TEST_CASE(faceSharingNonMatchingSubdivisionsThrow)
{
    auto parent = makeVerticalPillarGrid({4, 2, 2}, [](int, int, int k_) {
        return 2.0*(cellOf(k_) + sideOf(k_));
    });

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // ry differs (2 vs 3) on the shared i=2 face.
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView({{2,2,2}, {2,3,2}},
                                                 {{0,0,0}, {2,0,0}},
                                                 {{2,2,2}, {4,2,2}},
                                                 {"A", "B"}),
                      std::logic_error);
    BOOST_CHECK_EQUAL(grid.maxLevel(), 0);
}

// Worst-cell closure: max over leaf cells of |sum of area-weighted outward
// normals|. Zero iff every cell is closed (a necessary condition for a conformal
// leaf with no unmatched/hanging faces).
double maxClosure(const Dune::CpGrid& grid)
{
    double worst = 0.0;
    for (const auto& e : Dune::elements(grid.leafGridView())) {
        Dune::FieldVector<double,3> closure(0.0);
        for (const auto& is : Dune::intersections(grid.leafGridView(), e)) {
            auto n = is.centerUnitOuterNormal();
            n *= is.geometry().volume();
            closure += n;
        }
        worst = std::max(worst, closure.two_norm());
    }
    return worst;
}

// Every interior face is seen from both sides (a hanging/unmatched face would be
// seen once). Returns the count of leaf cells touched, for a sanity assertion.
void checkEveryInteriorFaceTwoSided(const Dune::CpGrid& grid)
{
    std::map<std::pair<int,int>, int> pairCount;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (!is.neighbor()) {
                continue;
            }
            const int in = is.inside().index();
            const int out = is.outside().index();
            pairCount[{std::min(in, out), std::max(in, out)}] += 1;
        }
    }
    for (const auto& [cells, count] : pairCount) {
        BOOST_CHECK_EQUAL(count, 2);
    }
    checkFacesWindAroundNormal(grid);
}

// ---------------------------------------------------------------------------
// A2 compatible sub-face mosaic (LGR_GAPS A2): two boxes meeting on a shared
// face with DIFFERENT but compatible in-face subdivisions (one a multiple of the
// other). The coarser side's interface cells become >6-face hexes tiled by the
// finer side's sub-faces. Implemented in assembleLeafGrid; these lock it in.
// ---------------------------------------------------------------------------

// STACKED pair (shared horizontal face): in-face directions I,J. Box LO is x4 in
// BOTH I and J, box HI is x2 (compatible 4-vs-2 in both). This mirrors the
// TLGR_VSTACK_HCOMPAT flow deck (TOP 4x4x2 over BOT 2x2x2). A horizontal
// interface never carries a fault, so this is the pure mosaic case.
BOOST_AUTO_TEST_CASE(stackedCompatibleInFaceMosaicBuilds)
{
    auto parent = makeVerticalPillarGrid({2, 2, 2}, [](int, int, int k_) {
        return 2.0*(cellOf(k_) + sideOf(k_));
    });
    Dune::CpGrid grid;
    auto raw = parent.raw();
    grid.processEclipseFormat(raw, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // Lower slab k in [0,1): x4 in I,J, x2 in K. Upper slab k in [1,2): x2 all.
    grid.addLgrsUpdateLeafView(
        {{4,4,2}, {2,2,2}}, {{0,0,0}, {0,0,1}}, {{2,2,1}, {2,2,2}}, {"LO", "HI"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    // LO: 4 parents x (4*4*2)=32 -> 128. HI: 4 parents x (2*2*2)=8 -> 32.
    BOOST_CHECK_EQUAL(grid.size(0), 128 + 32);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
    BOOST_CHECK_SMALL(maxClosure(grid), 1e-9);          // every cell closed
    checkEveryInteriorFaceTwoSided(grid);               // no hanging faces

    // The interface couples the two boxes: x4-vs-x2 in I and J -> each coarser
    // bottom cell on the HI side is met by 2x2 = 4 finer faces. Count
    // finer->coarser connections across the k=0/k=1 parent boundary.
    int crossConnections = 0;
    const auto& dims = grid.logicalCartesianSize();
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (!is.neighbor() || !is.inside().hasFather() || !is.outside().hasFather()) {
                continue;
            }
            const int kin  = grid.globalCell()[is.inside().index()]  / (dims[0]*dims[1]);
            const int kout = grid.globalCell()[is.outside().index()] / (dims[0]*dims[1]);
            if (kin != kout) {
                ++crossConnections;
            }
        }
    }
    // LO top layer is 8x8 (I,J) = 64 sub-faces, each a connection, counted from
    // both sides.
    BOOST_CHECK_EQUAL(crossConnections, 64 * 2);

    // No coincident-distinct vertices (shared-face corners merged across boxes).
    std::set<std::array<double,3>> coords;
    int vertexCount = 0;
    for (const auto& vertex : Dune::vertices(grid.leafGridView())) {
        const auto& c = vertex.geometry().center();
        coords.insert({c[0], c[1], c[2]});
        ++vertexCount;
    }
    BOOST_CHECK_EQUAL(coords.size(), static_cast<std::size_t>(vertexCount));
}

// SIDE-BY-SIDE pair (shared vertical face), no fault: in-face directions J,K.
// Box A is x4 in K, box B is x2 in K (compatible); equal in J. Confirms the
// mosaic also works on a vertical interface when there is no fault.
BOOST_AUTO_TEST_CASE(sideBySideCompatibleInFaceMosaicBuilds)
{
    auto parent = makeVerticalPillarGrid({2, 2, 2}, [](int, int, int k_) {
        return 2.0*(cellOf(k_) + sideOf(k_));
    });
    Dune::CpGrid grid;
    auto raw = parent.raw();
    grid.processEclipseFormat(raw, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // Box A i in [0,1): x4 in K, x2 in I,J. Box B i in [1,2): x2 all. Shared i=1
    // vertical face; in-face J equal, in-face K differs 4 vs 2.
    grid.addLgrsUpdateLeafView(
        {{2,2,4}, {2,2,2}}, {{0,0,0}, {1,0,0}}, {{1,2,2}, {2,2,2}}, {"A", "B"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    // A: 4 parents x (2*2*4)=16 -> 64. B: 4 parents x 8 -> 32.
    BOOST_CHECK_EQUAL(grid.size(0), 64 + 32);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
    BOOST_CHECK_SMALL(maxClosure(grid), 1e-9);
    checkEveryInteriorFaceTwoSided(grid);
}

// ---------------------------------------------------------------------------
// A5 (LGR_GAPS A5): two refined boxes meeting side-by-side across a FAULT. The
// shared vertical face coincides with the fault throw, so the two refined sides'
// sub-faces are vertically staggered. The faulted-boundary machinery now maps a
// "refined away" neighbour to the neighbour box's child cell (instead of skipping
// it), assembling the staggered interface conformally. First cut: equal
// subdivisions both sides.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(faceSharingAcrossFaultEqualSubdiv)
{
    auto depth = [](int i_, int, int k_) {
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;     // fault between i=1 and i=2
    };
    auto parent = makeVerticalPillarGrid({4, 2, 2}, depth);
    Dune::CpGrid grid;
    auto raw = parent.raw();
    grid.processEclipseFormat(raw, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // Box A i in [0,2), box B i in [2,4); shared i=2 face IS the fault. Equal
    // subdivisions both sides.
    grid.addLgrsUpdateLeafView(
        {{2,2,2}, {2,2,2}}, {{0,0,0}, {2,0,0}}, {{2,2,2}, {4,2,2}}, {"A", "B"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
    BOOST_CHECK_SMALL(maxClosure(grid), 1e-9);          // every cell closed
    checkEveryInteriorFaceTwoSided(grid);               // no hanging faces

    // The two boxes are actually connected across the fault: there exist interior
    // faces pairing an A-child (i<2) with a B-child (i>=2).
    int abConnections = 0;
    const auto& dims = grid.logicalCartesianSize();
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (!is.neighbor() || !is.inside().hasFather() || !is.outside().hasFather()) {
                continue;
            }
            const int ii = grid.globalCell()[is.inside().index()]  % dims[0];
            const int oi = grid.globalCell()[is.outside().index()] % dims[0];
            if ((ii < 2) != (oi < 2)) {
                ++abConnections;
            }
        }
    }
    BOOST_CHECK_GT(abConnections, 0);
}

// A5 + A2 combined: two refined boxes across a FAULT with DIFFERENT but compatible
// in-face subdivisions. Box A (i<2) is x4 in K, box B (i>=2) is x2 in K (equal in
// J); the shared i=2 vertical face is the fault throw. The finer side (A) emits
// every interface sub-face -- carrying the throw -- and references box B's
// coarser child that contains it, so the leaf is conformal.
BOOST_AUTO_TEST_CASE(faceSharingAcrossFaultCompatibleSubdiv)
{
    auto depth = [](int i_, int, int k_) {
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;     // fault between i=1 and i=2
    };
    auto parent = makeVerticalPillarGrid({4, 2, 2}, depth);
    Dune::CpGrid grid;
    auto raw = parent.raw();
    grid.processEclipseFormat(raw, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // Box A i in [0,2): x4 in K, x2 in I,J. Box B i in [2,4): x2 all. Shared i=2
    // face IS the fault; in-face J equal, in-face K differs 4 vs 2 (compatible).
    grid.addLgrsUpdateLeafView(
        {{2,2,4}, {2,2,2}}, {{0,0,0}, {2,0,0}}, {{2,2,2}, {4,2,2}}, {"A", "B"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
    BOOST_CHECK_SMALL(maxClosure(grid), 1e-9);          // every cell closed
    checkEveryInteriorFaceTwoSided(grid);               // no hanging faces

    int abConnections = 0;
    const auto& dims = grid.logicalCartesianSize();
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (!is.neighbor() || !is.inside().hasFather() || !is.outside().hasFather()) {
                continue;
            }
            const int ii = grid.globalCell()[is.inside().index()]  % dims[0];
            const int oi = grid.globalCell()[is.outside().index()] % dims[0];
            if ((ii < 2) != (oi < 2)) {
                ++abConnections;
            }
        }
    }
    BOOST_CHECK_GT(abConnections, 0);
}

// Nested LGR (a box whose parent is another LGR). The level grids are built
// over their parent (Phase A/B of docs/NESTED_LGR_PLAN.md), but the recursive
// leaf stitching (Phase C) is not implemented yet, so the build is expected to
// reach that boundary and throw a precise, staged message. Locks in: (a) parent
// ordering is honoured, (b) the staged failure is the *leaf* one (i.e. the
// nested level grids were successfully assembled first).
BOOST_AUTO_TEST_CASE(nestedRefinementReachesLeafBoundary)
{
    auto parent = makeVerticalPillarGrid({4, 2, 2}, [](int, int, int k_) {
        return 2.0*(cellOf(k_) + sideOf(k_));
    });

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // LGR1: 2x1x1 GLOBAL parents -> local refined dims 4x2x2.
    // NEST1 {0,0,0}-{2,2,2} touches LGR1's boundary (startIJK has zeros), so it
    // is rejected by the containment check (only fully-interior nesting yet).
    const auto isNestedLeafBoundary = [](const std::logic_error& e) {
        return std::string(e.what()).find("touches the boundary")
               != std::string::npos;
    };
    BOOST_CHECK_EXCEPTION(
        grid.addLgrsUpdateLeafView({{2,2,2}, {2,2,2}},
                                   {{1,0,0}, {0,0,0}},
                                   {{3,1,1}, {2,2,2}},
                                   {"LGR1", "NEST1"},
                                   {"GLOBAL", "LGR1"}),
        std::logic_error, isNestedLeafBoundary);
}

// The builder appends level grids in request order and resolves each box's
// parent by name, so a child must follow its parent. A child listed before its
// parent must be rejected with an actionable message (the simulator's
// topological sort is what guarantees this never happens in practice).
BOOST_AUTO_TEST_CASE(nestedChildBeforeParentThrows)
{
    auto parent = makeVerticalPillarGrid({4, 2, 2}, [](int, int, int k_) {
        return 2.0*(cellOf(k_) + sideOf(k_));
    });

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    const auto isUnbuiltParent = [](const std::logic_error& e) {
        return std::string(e.what()).find("has not been built yet")
               != std::string::npos;
    };
    BOOST_CHECK_EXCEPTION(
        grid.addLgrsUpdateLeafView({{2,2,2}, {2,2,2}},
                                   {{0,0,0}, {1,0,0}},
                                   {{2,2,2}, {3,1,1}},
                                   {"NEST1", "LGR1"},
                                   {"LGR1", "GLOBAL"}),
        std::logic_error, isUnbuiltParent);
}

// A CARFIN whose parent is another LGR (nested), with the child box strictly
// interior to the parent LGR - it touches none of the parent's refined boundary
// faces. This is the simplest nesting case: the child's leaf stitching is
// entirely within the parent LGR's refined region, with no interaction with
// GLOBAL. Until the leaf assembler handles nesting (docs/NESTED_LGR_PLAN.md
// Phase C) the nested level grids are built (Phases A/B) and then a precise
// "not implemented" is thrown.
//
// Phase-C acceptance: replace BOOST_CHECK_EXCEPTION below with the build call,
// then assert (mirroring endToEndSingleBox):
//   - grid.maxLevel() == 2 and getLgrNameToLevel() == {LGR1:1, NEST1:2};
//   - total leaf volume == volume before refinement (conservation);
//   - every NEST1 leaf cell's father is an LGR1 cell whose father is GLOBAL
//     (the child -> parent-LGR -> GLOBAL chain), geometryInFather == 1/8;
//   - leaf intersections are conformal (each interior face seen from both sides
//     with matching area).
BOOST_AUTO_TEST_CASE(nestedFullyContainedBuilds)
{
    // GLOBAL 4x3x3 with distorted vertical pillars.
    auto depth = [](int i_, int j_, int k_) {
        const double x = cellOf(i_) + sideOf(i_);
        const double y = cellOf(j_) + sideOf(j_);
        return 2.0*(cellOf(k_) + sideOf(k_)) + 0.3*std::sin(0.9*x) + 0.2*std::cos(0.6*y);
    };
    auto parent = makeVerticalPillarGrid({4, 3, 3}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // LGR1: GLOBAL cells i:1-2, j:1, k:1 (2x1x1) refined 3x3x3 -> LGR1-local 6x3x3 = 54.
    // NEST1: LGR1-local i:1-4, j:1, k:1 (4 LGR1 cells) refined 2x2x2 = 32 -> strictly
    //        interior to LGR1 (touches no LGR1 boundary face), parent grid "LGR1".
    grid.addLgrsUpdateLeafView({{3,3,3}, {2,2,2}},
                               {{1,1,1}, {1,1,1}},
                               {{3,2,2}, {5,2,2}},
                               {"LGR1", "NEST1"},
                               {"GLOBAL", "LGR1"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    BOOST_CHECK_EQUAL(grid.getLgrNameToLevel().at("LGR1"), 1);
    BOOST_CHECK_EQUAL(grid.getLgrNameToLevel().at("NEST1"), 2);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);

    // Father chain: NEST1 cell (level 2) -> LGR1 cell (level 1) -> GLOBAL (level 0).
    int nestCount = 0, lgr1Count = 0, coarseCount = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        const int lvl = element.level();
        if (lvl == 2) {
            ++nestCount;
            const auto father = element.father();
            BOOST_CHECK_EQUAL(father.level(), 1);
            BOOST_REQUIRE(father.hasFather());
            BOOST_CHECK_EQUAL(father.father().level(), 0);
            BOOST_CHECK_CLOSE(element.geometryInFather().volume(), 1.0/8.0, 1e-9);
        }
        else if (lvl == 1) {
            ++lgr1Count;
            BOOST_CHECK_EQUAL(element.father().level(), 0);
            BOOST_CHECK_CLOSE(element.geometryInFather().volume(), 1.0/27.0, 1e-9);
        }
        else {
            ++coarseCount;
            BOOST_CHECK(!element.hasFather());
        }
    }
    BOOST_CHECK_EQUAL(nestCount, 32);          // 4 LGR1 cells x 2x2x2
    BOOST_CHECK_EQUAL(lgr1Count, 54 - 4);      // 6x3x3 minus the 4 NEST refined away
    BOOST_CHECK_EQUAL(coarseCount, 4*3*3 - 2); // GLOBAL minus the 2 LGR1 refined away
    BOOST_CHECK_EQUAL(grid.size(0), nestCount + lgr1Count + coarseCount);

    // Conformal leaf: every interior face seen from both sides with matching area.
    std::map<std::pair<int,int>, int> pairCount;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& intersection : Dune::intersections(grid.leafGridView(), element)) {
            if (intersection.neighbor()) {
                const int inside = intersection.inside().index();
                const int outside = intersection.outside().index();
                BOOST_REQUIRE(inside != outside);
                pairCount[{std::min(inside, outside), std::max(inside, outside)}] += 1;
            }
        }
    }
    for (const auto& [cells, count] : pairCount) {
        BOOST_CHECK_EQUAL(count, 2);
    }

    // Global ids unique over leaf cells and points (delegated across 3 levels).
    const auto& ids = grid.globalIdSet();
    std::set<std::int64_t> cellIds;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        cellIds.insert(ids.id(element));
    }
    BOOST_CHECK_EQUAL(cellIds.size(), static_cast<std::size_t>(grid.size(0)));
    std::set<std::int64_t> pointIds;
    for (const auto& vertex : Dune::vertices(grid.leafGridView())) {
        pointIds.insert(ids.id(vertex));
    }
    BOOST_CHECK_EQUAL(pointIds.size(), static_cast<std::size_t>(grid.size(3)));
}

BOOST_AUTO_TEST_CASE(faultInsideBoxEndToEnd)
{
    // Fault between i=1 and i=2 columns, throw 0.6; box covers it.
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;
    };
    auto parent = makeVerticalPillarGrid({4, 2, 2}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // The box {1..3} has the fault *inside* (between i=1 and i=2) and
    // unfaulted boundaries (i=0|1 and i=3|4 columns are unfaulted... note
    // the fault plane is i=2, the box boundary planes are i=1 and i=3,
    // both unfaulted since the throw is constant on each side).
    grid.addLgrsUpdateLeafView({{2,2,2}}, {{1,0,0}}, {{3,2,2}}, {"LGR1"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);

    // Both-sided intersection consistency across the refined fault.
    std::map<std::pair<int,int>, int> pairCount;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& intersection : Dune::intersections(grid.leafGridView(), element)) {
            if (intersection.neighbor()) {
                const int inside = intersection.inside().index();
                const int outside = intersection.outside().index();
                pairCount[{std::min(inside, outside), std::max(inside, outside)}] += 1;
            }
        }
    }
    for (const auto& [cells, count] : pairCount) {
        BOOST_CHECK_EQUAL(count, 2);
    }
}

BOOST_AUTO_TEST_CASE(faultAtBoxBoundaryBuilds)
{
    // Fault between i=1 and i=2 (throw 0.6); box boundary right on the fault
    // plane. The boundary parents each connect to two coarse cells, so the
    // boundary is rebuilt as split faces. Phase 1: assert the leaf is a valid
    // polyhedral grid (builds, volume conserved, every interior face two-sided,
    // refined boundary cells connected to the coarse neighbours).
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;
    };
    auto parent = makeVerticalPillarGrid({4, 2, 2}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    BOOST_REQUIRE_NO_THROW(grid.addLgrsUpdateLeafView({{2,2,2}}, {{2,0,0}}, {{4,2,2}}, {"LGR1"}));
    BOOST_CHECK_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);

    // Geometric validity: every leaf cell is closed, i.e. the sum of
    // area * outward-unit-normal over all its faces is zero (divergence
    // theorem). This validates the synthetic split faces' areas AND normal
    // orientation, independent of how the interface is subdivided. Also: every
    // interior face is two-sided, and the box connects to coarse neighbours.
    std::map<std::pair<int,int>, int> pairCount;
    int refinedToCoarse = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        Dune::FieldVector<double,3> closure(0.0);
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            auto n = is.centerUnitOuterNormal();
            n *= is.geometry().volume();
            closure += n;
            if (is.neighbor()) {
                const int a = is.inside().index();
                const int b = is.outside().index();
                BOOST_REQUIRE(a != b);
                pairCount[{std::min(a, b), std::max(a, b)}] += 1;
                if (is.inside().hasFather() != is.outside().hasFather()) {
                    ++refinedToCoarse;
                }
            }
        }
        BOOST_CHECK_SMALL(closure.two_norm(), 1e-9);
    }
    for (const auto& [cells, count] : pairCount) {
        BOOST_CHECK_EQUAL(count % 2, 0); // each shared face seen from both sides
    }
    BOOST_CHECK_GT(refinedToCoarse, 0);  // the box connects to coarse neighbours
}

// Note on the boundary throw branches: with vertical pillars + a ZCORN throw
// (the helper here), a faulted box boundary always produces *split* faces (a
// parent connecting to two coarse cells) — the 'fault-split face' branch. The
// sibling 'partial (faulted or degenerate) face' branch needs the geometric
// offset-pillar fault that the helper cannot express; it is exercised by the
// flow deck opm-tests/flow_diagnostic_test/SIMPLE_2PH_W_FAULT_LGR.DATA.

BOOST_AUTO_TEST_CASE(faultAtBoxBoundaryJDirectionBuilds)
{
    // Fault and box boundary in the J direction — the other lateral axis.
    auto depth = [](int i_, int j_, int k_) {
        (void)i_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(j_) >= 2) ? base + 0.6 : base;
    };
    auto parent = makeVerticalPillarGrid({2, 4, 2}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    BOOST_REQUIRE_NO_THROW(grid.addLgrsUpdateLeafView({{2,2,2}}, {{0,2,0}}, {{2,4,2}}, {"LGR1"}));
    checkValidFaultedLeaf(grid, volumeBefore);
}

BOOST_AUTO_TEST_CASE(faultAtBoxBoundaryUpThrownSideBuilds)
{
    // Box on the *up-thrown* side: fault between i=1 and i=2 (cells i>=2 thrown
    // down), box i0-1 so its right boundary (i=2) faces the thrown neighbour.
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;
    };
    auto parent = makeVerticalPillarGrid({4, 2, 2}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    BOOST_REQUIRE_NO_THROW(grid.addLgrsUpdateLeafView({{2,2,2}}, {{0,0,0}}, {{2,2,2}}, {"LGR1"}));
    checkValidFaultedLeaf(grid, volumeBefore);
}

BOOST_AUTO_TEST_CASE(faultLargeThrowAtBoxBoundaryBuilds)
{
    // Throw larger than one layer (3.0 > layer thickness 2.0), three layers, so
    // a boundary parent spans two coarse neighbours — a wider overlap pattern.
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 3.0 : base;
    };
    auto parent = makeVerticalPillarGrid({4, 2, 3}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    BOOST_REQUIRE_NO_THROW(grid.addLgrsUpdateLeafView({{2,2,2}}, {{2,0,0}}, {{4,2,3}}, {"LGR1"}));
    checkValidFaultedLeaf(grid, volumeBefore);
}

namespace
{

// Face area between cells whose parent is at i = iIn, layers [k0, k1), and cells
// whose parent is at i = iOut: the connection across a fault at an i-face.
double areaAcrossIFace(const Dune::CpGrid& grid, int iIn, int iOut, int k0, int k1)
{
    const auto& dims = grid.logicalCartesianSize();
    const auto ijk = [&dims](int cart) {
        return std::array<int,3>{ cart % dims[0], (cart / dims[0]) % dims[1],
                                  cart / (dims[0]*dims[1]) };
    };
    double area = 0.0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        const auto in = ijk(grid.globalCell()[element.index()]);
        if (in[0] != iIn || in[2] < k0 || in[2] >= k1) {
            continue;
        }
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (is.neighbor() && ijk(grid.globalCell()[is.outside().index()])[0] == iOut) {
                area += is.geometry().volume();
            }
        }
    }
    return area;
}

} // anonymous namespace

// A throw taller than the box: the cells across the fault lie outside the box's
// own layers, and the box must still reach them with the level-zero face area.
BOOST_AUTO_TEST_CASE(faultThrowTallerThanBoxKeepsConnections)
{
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 5.0 : base;
    };
    auto parent = makeVerticalPillarGrid({4, 2, 8}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);
    const double areaBefore = areaAcrossIFace(grid, 2, 1, 3, 5);
    BOOST_REQUIRE_GT(areaBefore, 0.0);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    BOOST_REQUIRE_NO_THROW(grid.addLgrsUpdateLeafView({{2,2,2}}, {{2,0,3}}, {{4,2,5}}, {"LGR1"}));
    checkValidFaultedLeaf(grid, volumeBefore);
    BOOST_CHECK_CLOSE(areaAcrossIFace(grid, 2, 1, 3, 5), areaBefore, 1e-8);
}

// Two boxes that meet only through a throw, with in-face subdivisions that do
// not nest: refused, rather than sealing the interface between them.
BOOST_AUTO_TEST_CASE(boxesMeetingOnlyAcrossThrowMustNest)
{
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 5.0 : base;
    };
    auto parent = makeVerticalPillarGrid({4, 2, 8}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView({{2,3,2}, {2,2,3}},
                                                 {{0,0,6}, {2,0,3}},
                                                 {{2,2,8}, {4,2,5}},
                                                 {"A", "B"}),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(faultNotOnBoxBoundaryBuilds)
{
    // Control: a fault exists (between i=0 and i=1) but the box (i=2..3) sits
    // entirely on the down-thrown side, so its boundary toward the coarse
    // neighbour (i=1) is a full, matching face. This must build today and
    // conserve volume — it isolates "fault on the boundary" from "fault merely
    // present in the grid".
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 1) ? base + 0.6 : base;
    };
    auto parent = makeVerticalPillarGrid({4, 2, 2}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    BOOST_CHECK_NO_THROW(grid.addLgrsUpdateLeafView({{2,2,2}}, {{2,0,0}}, {{4,2,2}}, {"LGR1"}));
    BOOST_CHECK_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
}

BOOST_AUTO_TEST_CASE(processGrdeclSplitsFaultedCoarseFineColumnPair)
{
    // Kernel proof for faults-at-box-boundary (docs/PLAN.md Track 1): the planned
    // fix builds the connecting faces at a faulted box boundary by running
    // process_grdecl on a mini-grdecl (coarse neighbour + refined box boundary)
    // and reading the split connections. This checks that mechanism in isolation.
    //
    // Column 0 (coarse): one cell z[0,4]; its second logical cell is pinched
    //   (z[4,4]) so process_grdecl removes it -> one coarse cell.
    // Column 1 (refined): two cells z[1,3] and z[3,5] -- the same span shifted by
    //   a 1.0 throw. The coarse cell must connect to BOTH refined cells.
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        if (cellOf(i_) == 0) {                 // coarse column: cell0 [0,4], cell1 pinched at 4
            return (k_ == 0) ? 0.0 : 4.0;
        }
        static const double fine[4] = {1.0, 3.0, 3.0, 5.0};   // refined, throw +1
        return fine[k_];
    };
    auto g = makeVerticalPillarGrid({2, 1, 2}, depth);
    auto raw = g.raw();

    struct processed_grid out;
    const int ok = process_grdecl(/*pinchActive=*/1, /*edge_conformal=*/0, 1e-6, &raw, nullptr, &out);
    BOOST_REQUIRE(ok);

    // Pinched coarse cell removed: one coarse + two refined.
    BOOST_CHECK_EQUAL(out.number_of_cells, 3);

    // The split connection lives in the interior I-faces between the two columns
    // (face_neighbors holds compressed active cell ids; -1 is a domain boundary).
    // The coarse cell (0) must connect to both refined cells (1 and 2).
    int interiorIFaces = 0;
    std::set<int> coarseConnectedTo;
    for (unsigned f = 0; f < out.number_of_faces; ++f) {
        const int a = out.face_neighbors[2*f];
        const int b = out.face_neighbors[2*f + 1];
        if (out.face_tag[f] == I_FACE && a >= 0 && b >= 0) {
            ++interiorIFaces;
            BOOST_CHECK((a == 0) || (b == 0));        // one side is the coarse cell
            coarseConnectedTo.insert(a == 0 ? b : a);
        }
    }
    BOOST_CHECK_EQUAL(interiorIFaces, 2);
    BOOST_CHECK((coarseConnectedTo == std::set<int>{1, 2}));

    free_processed_grid(&out);
}

// ---------------------------------------------------------------------------
// Skew-pillar variants of the face-sharing cases. The corner pool identifies
// shared box-boundary corners by exact coordinate, relying on both boxes
// resampling with bitwise-identical results (docs/CORNER-POOL-EXACTNESS.md).
// Inclined pillars and binary-inexact depth increments make that a real test.
// A merge miss leaves two vertices an ulp apart (caught by the pair-distance
// check) and, for equal factors, seals the interface (caught by the connection
// counts and the assembler's unmatched-owner throw).
// ---------------------------------------------------------------------------

namespace
{

// Every pillar inclined with its own tilt, so sub-pillar interpolation and
// zcorn resampling produce full-precision doubles.
TestGrdecl makeSkewPillarGrid(const std::array<int,3>& dims,
                              const std::function<double(int,int,int)>& depth)
{
    TestGrdecl g = makeVerticalPillarGrid(dims, depth);
    const auto& [nx, ny, nz] = dims;
    static_cast<void>(nz);
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            double* p = &g.coord[6*(static_cast<std::size_t>(j)*(nx + 1) + i)];
            p[3] = i + 0.4 + 0.05*i - 0.03*j;
            p[4] = j - 0.3 + 0.04*j + 0.02*i;
        }
    }
    return g;
}

// Layers tilted in i and j with binary-inexact increments; node-based, so
// faces conform bitwise across cells.
double tiltedDepth(int i_, int j_, int k_)
{
    return 2.0*((k_ + 1)/2) + 0.1*((i_ + 1)/2) + 0.07*((j_ + 1)/2);
}

// Smallest distance between distinct leaf vertices. Complements the exact
// duplicate-coordinate check: a corner-pool merge miss produces two vertices
// an ulp apart, which no exact-set check can see.
double minVertexPairDistance(const Dune::CpGrid& grid)
{
    std::vector<Dune::FieldVector<double,3>> pts;
    for (const auto& vertex : Dune::vertices(grid.leafGridView())) {
        pts.push_back(vertex.geometry().center());
    }
    double minDist = std::numeric_limits<double>::max();
    for (std::size_t a = 0; a < pts.size(); ++a) {
        for (std::size_t b = a + 1; b < pts.size(); ++b) {
            auto d = pts[a];
            d -= pts[b];
            minDist = std::min(minDist, d.two_norm());
        }
    }
    return minDist;
}

// Refined-to-refined leaf connections whose parents lie on opposite sides of
// splitAt along axis (i.e. across the box-box interface), counted twice.
int crossBoxConnections(const Dune::CpGrid& grid, int axis, int splitAt)
{
    const auto& dims = grid.logicalCartesianSize();
    const auto parentAlong = [&](int leafCell) {
        const int cart = grid.globalCell()[leafCell];
        const std::array<int,3> ijk = { cart % dims[0],
                                        (cart / dims[0]) % dims[1],
                                        cart / (dims[0]*dims[1]) };
        return ijk[axis];
    };
    int count = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (!is.neighbor() || !is.inside().hasFather() || !is.outside().hasFather()) {
                continue;
            }
            const int pin = parentAlong(is.inside().index());
            const int pout = parentAlong(is.outside().index());
            if ((pin < splitAt) != (pout < splitAt)) {
                ++count;
            }
        }
    }
    return count;
}

} // anonymous namespace

// Equal in-face factors on inclined pillars: 1:1 face pairing by merged
// corner set. A merge miss would seal the i=2 interface.
BOOST_AUTO_TEST_CASE(skewPillarFaceSharingBoxes)
{
    auto parent = makeSkewPillarGrid({4, 2, 2}, tiltedDepth);
    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    grid.addLgrsUpdateLeafView({{2,2,2}, {3,2,2}},
                               {{0,0,0}, {2,0,0}},
                               {{2,2,2}, {4,2,2}},
                               {"A", "B"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    BOOST_CHECK_EQUAL(grid.size(0), 8*8 + 8*12);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-6);
    BOOST_CHECK_SMALL(maxClosure(grid), 1e-9);
    checkEveryInteriorFaceTwoSided(grid);
    // ry*rz = 4 connections per parent pair, 2x2 pairs, counted twice.
    BOOST_CHECK_EQUAL(crossBoxConnections(grid, 0, 2), 16*2);
    // Shared corners merged to one vertex; nothing an ulp apart either.
    BOOST_CHECK_GT(minVertexPairDistance(grid), 1e-3);
}

// Compatible (one side uniformly finer) in-face factors on inclined pillars:
// the finer side's even-index fractions (2/4, ...) must reproduce the coarser
// side's (1/2, ...) bitwise for the shared corners to merge.
BOOST_AUTO_TEST_CASE(skewPillarCompatibleSubdivisionMosaic)
{
    auto parent = makeSkewPillarGrid({4, 2, 2}, tiltedDepth);
    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    grid.addLgrsUpdateLeafView({{2,4,4}, {2,2,2}},
                               {{0,0,0}, {2,0,0}},
                               {{2,2,2}, {4,2,2}},
                               {"A", "B"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    BOOST_CHECK_EQUAL(grid.size(0), 8*32 + 8*8);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-6);
    // The coarser interface cell is closed by the finer side's sub-face
    // triangulations; on inclined pillars the shared face is curved, so the
    // cell closes only to the face's non-planarity, not to roundoff.
    BOOST_CHECK_SMALL(maxClosure(grid), 1e-7);
    checkEveryInteriorFaceTwoSided(grid);
    // Finer side emits 4x4 sub-faces per parent pair, 2x2 pairs, twice.
    BOOST_CHECK_EQUAL(crossBoxConnections(grid, 0, 2), 64*2);
    BOOST_CHECK_GT(minVertexPairDistance(grid), 1e-3);
}

// Stacked boxes sharing a k-face on inclined pillars: the shared corners
// collapse in the vertical fraction onto the tilted layer surface.
BOOST_AUTO_TEST_CASE(skewPillarStackedBoxesShareKFace)
{
    auto parent = makeSkewPillarGrid({2, 2, 4}, tiltedDepth);
    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    grid.addLgrsUpdateLeafView({{2,2,2}, {2,2,3}},
                               {{0,0,0}, {0,0,2}},
                               {{2,2,2}, {2,2,4}},
                               {"LOW", "HIGH"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    BOOST_CHECK_EQUAL(grid.size(0), 8*8 + 8*12);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-6);
    BOOST_CHECK_SMALL(maxClosure(grid), 1e-9);
    checkEveryInteriorFaceTwoSided(grid);
    // rx*ry = 4 connections per parent pair, 2x2 pairs, counted twice.
    BOOST_CHECK_EQUAL(crossBoxConnections(grid, 2, 2), 16*2);
    BOOST_CHECK_GT(minVertexPairDistance(grid), 1e-3);
}

// Boxes stacked in a column get bitwise the same refined pillars; with
// pillarsFromBoxLayer (the reference) only a box starting at the column top does.
BOOST_AUTO_TEST_CASE(skewPillarRefinedPillarsFromColumn)
{
    const auto parent = makeSkewPillarGrid({2, 2, 4}, tiltedDepth);
    const auto coordOf = [&parent](int k0, int k1, bool fromBox) {
        Opm::Refinement::BlockRefinement req;
        req.name = "B";
        req.cellsPerDim = {2, 2, 1};
        req.startIJK = {0, 0, k0};
        req.endIJK = {2, 2, k1};
        req.pillarsFromBoxLayer = fromBox;
        return Opm::Refinement::refineBlock(parent.dims, parent.coord.data(),
                                            parent.zcorn.data(), parent.actnum.data(),
                                            req).coord;
    };
    BOOST_CHECK(coordOf(0, 2, false) == coordOf(2, 4, false));
    BOOST_CHECK(coordOf(0, 2, true) == coordOf(0, 2, false));
}
