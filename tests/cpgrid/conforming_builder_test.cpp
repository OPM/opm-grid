/*
  Copyright 2026 Equinor ASA.

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
#include <opm/grid/cpgrid/refinement/conforming/GrdeclRefinement.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <tests/cpgrid/lgr/ConformingTestHelpers.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Opm::ConformingTest;

namespace
{

// Cartesian index triplet of a leaf cell; a refined cell has its level-zero parent's.
std::array<int,3> ijkOf(const Dune::CpGrid& grid, int cell)
{
    std::array<int,3> ijk;
    grid.getIJK(cell, ijk);
    return ijk;
}

// Every leaf face's corners wind counter-clockwise around its stored normal: the
// face's area vector (CpGrid::faceAreaNormalEcl) points along it.
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
        const auto area = grid.faceAreaNormalEcl(face);
        if (area.dot(leaf.faceNormals(face)) <= 0.0) {
            ++inverted;
        }
    }
    BOOST_CHECK_EQUAL(inverted, 0);
}

// Structural validity of a leaf refined by one box: one refined level, volume
// conserved, every cell closed (sum of area*outward-normal is zero), and the box
// connected to coarse neighbours.
void checkLeafClosedAndConnectionsToCoarse(const Dune::CpGrid& grid, double volumeBefore)
{
    BOOST_CHECK_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
    BOOST_CHECK_SMALL(maxClosure(grid), 1e-9);
    int refinedCoarseConnections = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (is.neighbor() && is.inside().hasFather() != is.outside().hasFather()) {
                ++refinedCoarseConnections; // each face counted twice, once from each side
            }
        }
    }
    BOOST_CHECK_GT(refinedCoarseConnections, 0);
    checkFacesWindAroundNormal(grid);
}

// Every interior face is seen from both sides (a hanging or unmatched face would be
// seen once), joins at most two cells, and winds around its normal.
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
            BOOST_REQUIRE(in != out);
            BOOST_CHECK_LT(grid.currentLeafData().faceToCellSize(is.id()), 3); // at most 2 cells
            pairCount[{std::min(in, out), std::max(in, out)}] += 1;
        }
    }
    for (const auto& [cells, count] : pairCount) {
        BOOST_CHECK_EQUAL(count, 2);
    }
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
    auto parent = makeGrid({4, 3, 3}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);
    const int cellsBefore = grid.size(0);

    useConformingBuilder(grid, parent);

    // Box of 2x1x1 parents, refined 2x2x2.
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}},
                               /* startIJK_vec = */      {{1,1,1}},
                               /* endIJK_vec = */        {{3,2,2}},
                               /* lgr_name_vec = */      {"LGR1"});

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

    // Intersection consistency over the whole leaf.
    checkEveryInteriorFaceTwoSided(grid);

    // The coarse i-neighbor of the box (cell (0,1,1)) sees the refined
    // mosaic: 5 unrefined faces + 4 refined ones.
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        if (element.hasFather()) {
            continue;
        }
        const auto [i, j, k] = ijkOf(grid, element.index());
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

    // Unique global ids for cells and points.
    checkUniqueGlobalIds(grid);
}

BOOST_AUTO_TEST_CASE(stableCellIdDisambiguatesRefinedSiblings)
{
    // global_cell_ is NOT unique on the leaf -- refined siblings share the
    // parent's Cartesian index -- so it cannot key per-cell output gathered
    // across ranks. stableCellId() must give every leaf cell a unique id:
    // coarse cells keep their Cartesian index; refined cells encode
    // (parent Cartesian index, index in parent cell) under a tag bit.
    auto parent = makeGrid({4, 3, 3}, flat);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);

    useConformingBuilder(grid, parent);
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}},
                               /* startIJK_vec = */      {{1,1,1}},
                               /* endIJK_vec = */        {{3,2,2}},
                               /* lgr_name_vec = */      {"LGR1"}); // 2 parents -> 16 refined

    constexpr std::int64_t tag = std::int64_t(1) << 62;
    constexpr std::int64_t childMask = (std::int64_t(1) << 20) - 1;

    const auto sid = grid.stableCellId();
    BOOST_REQUIRE_EQUAL(sid.size(), static_cast<std::size_t>(grid.size(0)));

    // global_cell_ really does collide for siblings, but stableCellId does not.
    const auto& gc = grid.globalCell();
    BOOST_CHECK_LT(std::set<int>(gc.begin(), gc.end()).size(), gc.size());
    BOOST_CHECK_EQUAL(std::set<std::int64_t>(sid.begin(), sid.end()).size(), sid.size());

    std::map<int, std::set<int>> childrenOfParent; // father Cartesian index -> indices in parent cell seen
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
    BOOST_CHECK_EQUAL(coarse, 4*3*3 - 2); // level zero cell count - 2 parent cells
    // Each of the two refined parents has all 8 distinct indices in parent cell, 0..7.
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
    auto parent = makeGrid({6, 3, 3}, flat);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, parent);

    // Overlapping boxes are rejected (validation).
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}, {2,2,2}},
                                                 /* startIJK_vec = */      {{0,0,0}, {1,0,0}},
                                                 /* endIJK_vec = */        {{3,2,2}, {4,2,2}},
                                                 /* lgr_name_vec = */      {"A", "B"}),
                      std::invalid_argument);
    BOOST_CHECK_EQUAL(grid.maxLevel(), 0); // grid unchanged

    // Two separated boxes work.
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}, {3,3,3}},
                               /* startIJK_vec = */      {{0,0,0}, {4,1,1}},
                               /* endIJK_vec = */        {{1,1,1}, {6,3,3}},
                               /* lgr_name_vec = */      {"A", "B"});
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

// A box with no host left is an empty LGR: its level grid has no cells and
// the leaf keeps every other cell. The hosts are either inactive or collapsed
// to zero thickness, which is how MINPV leaves the cells it removes.
BOOST_AUTO_TEST_CASE(boxWithoutHostsIsEmptyLgr)
{
    const auto isHost = [](int i, int j) { return (i < 2) && (j < 2); };
    for (const bool collapsed : {false, true}) {
        auto parent = makeGrid({4, 4, 1}, [&](int i_, int j_, int k_) {
            const bool host = isHost(cellOf(i_), cellOf(j_));
            return (collapsed && host) ? 0.0 : 10.0*sideOf(k_);
        });
        if (!collapsed) {
            parent.actnum.assign(16, 1);
            for (int j = 0; j < 4; ++j) {
                for (int i = 0; i < 4; ++i) {
                    parent.actnum[i + 4*j] = isHost(i, j) ? 0 : 1;
                }
            }
        }
        auto rawParent = parent.raw();

        // The empty box alone.
        {
            Dune::CpGrid grid;
            grid.processEclipseFormat(rawParent, false);
            BOOST_REQUIRE_EQUAL(grid.size(0), 12);
            useConformingBuilder(grid, parent);

            grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{3,3,1}},
                                       /* startIJK_vec = */      {{0,0,0}},
                                       /* endIJK_vec = */        {{2,2,1}},
                                       /* lgr_name_vec = */      {"LGR1"});

            BOOST_REQUIRE_EQUAL(grid.maxLevel(), 1);
            BOOST_CHECK_EQUAL(grid.getLgrNameToLevel().at("LGR1"), 1);
            BOOST_CHECK_EQUAL(grid.size(1, 0), 0);  // level 1 grid view
            BOOST_CHECK_EQUAL(grid.size(0), 12);    // leaf grid view
        }

        // Next to a box that refines.
        {
            Dune::CpGrid grid;
            grid.processEclipseFormat(rawParent, false);
            const double volumeBefore = totalVolume(grid);
            useConformingBuilder(grid, parent);

            grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{3,3,1}, {2,2,2}},
                                       /* startIJK_vec = */      {{0,0,0}, {3,3,0}},
                                       /* endIJK_vec = */        {{2,2,1}, {4,4,1}},
                                       /* lgr_name_vec = */      {"EMPTY", "B"});

            BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
            BOOST_CHECK_EQUAL(grid.size(1, 0), 0); // level 1 grid view
            BOOST_CHECK_EQUAL(grid.size(2, 0), 8); // level 2 grid view
            BOOST_CHECK_EQUAL(grid.size(0), 12 - 1 + 8); // leaf grid view
            BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
        }
    }
}

// Two diagonally-adjacent boxes share only an edge (the CARFIN deck
// configuration). The refined corners along the shared edge that are
// interior to a parent cell must merge into single leaf vertices.
BOOST_AUTO_TEST_CASE(edgeSharingBoxesMergeCorners)
{
    auto parent = makeGrid({4, 4, 2}, flat);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, parent);

    // Box A = i,j in [0,2); box B = i,j in [2,4); both k in [0,2).
    // They meet along the vertical edge (i=2, j=2). Factors 2x2x2.
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}, {2,2,2}},
                               /* startIJK_vec = */      {{0,0,0}, {2,2,0}},
                               /* endIJK_vec = */        {{2,2,2}, {4,4,2}},
                               /* lgr_name_vec = */      {"A", "B"});

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

    // Leaf vertices: level-zero corners outside both boxes (each box covers 27; the 3 on
    // the shared edge belong to both) + both level grids' corners - the 5 corners on the
    // shared edge, which are in both level grids.
    const int expectedVertexCount = grid.size(/* level = */ 0, /* codim = */ 3) - (2*3*3*3 - 3)
                                  + grid.size(/* level = */ 1, /* codim = */ 3)
                                  + grid.size(/* level = */ 2, /* codim = */ 3) - 5;
    BOOST_CHECK_EQUAL(grid.size(/* codim = */ 3), expectedVertexCount);

    // Two-sided intersection symmetry over the whole leaf.
    checkEveryInteriorFaceTwoSided(grid);

    // Unique global ids for cells and points.
    checkUniqueGlobalIds(grid);
}

// Face-sharing boxes: the boundary faces of the two refined blocks on the
// shared plane pair into interior faces. The out-of-face subdivision may
// differ between the two boxes; the in-face subdivisions must match.
BOOST_AUTO_TEST_CASE(faceSharingBoxes)
{
    auto parent = makeGrid({4, 2, 2}, flat);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, parent);

    // Box A = i in [0,2), box B = i in [2,4); both full in j,k -> share the
    // i=2 face. Different in-plane-normal factor (rx 2 vs 3), matching
    // in-face factors (ry=rz=2).
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}, {3,2,2}},
                               /* startIJK_vec = */      {{0,0,0}, {2,0,0}},
                               /* endIJK_vec = */        {{2,2,2}, {4,2,2}},
                               /* lgr_name_vec = */      {"A", "B"});

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

    // Leaf vertices: no level-zero corner is outside the boxes; both level grids'
    // corners, less the 5 x 5 on the shared face x = 2, which are in both.
    const int expectedVertexCount = grid.size(/* level = */ 1, /* codim = */ 3)
                                  + grid.size(/* level = */ 2, /* codim = */ 3) - 5*5;
    BOOST_CHECK_EQUAL(grid.size(/* codim = */ 3), expectedVertexCount);

    // Two-sided intersection symmetry, and the A/B blocks are actually
    // connected (interior faces pairing an A-child with a B-child exist).
    std::map<std::pair<int,int>, int> pairCount;
    int abConnections = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& intersection : Dune::intersections(grid.leafGridView(), element)) {
            if (!intersection.neighbor()) {
                continue;
            }
            const int in = intersection.inside().index();
            const int out = intersection.outside().index();
            pairCount[{std::min(in, out), std::max(in, out)}] += 1;
            BOOST_CHECK_LT(grid.currentLeafData().faceToCellSize(intersection.id()), 3); // at most 2 cells
            if (intersection.inside().hasFather() && intersection.outside().hasFather()) {
                const int ci = ijkOf(grid, in)[0];
                const int co = ijkOf(grid, out)[0];
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

    // Unique global ids for cells and points.
    checkUniqueGlobalIds(grid);
}

// Non-matching in-face subdivisions are rejected.
BOOST_AUTO_TEST_CASE(faceSharingNonMatchingSubdivisionsThrow)
{
    auto parent = makeGrid({4, 2, 2}, flat);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);

    useConformingBuilder(grid, parent);

    // ry differs (2 vs 3) on the shared i=2 face.
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}, {2,3,2}},
                                                 /* startIJK_vec = */      {{0,0,0}, {2,0,0}},
                                                 /* endIJK_vec = */        {{2,2,2}, {4,2,2}},
                                                 /* lgr_name_vec = */      {"A", "B"}),
                      std::logic_error);
    BOOST_CHECK_EQUAL(grid.maxLevel(), 0);
}

// Face-sharing boxes, each finer in a different in-face direction, are rejected.
BOOST_AUTO_TEST_CASE(faceSharingMixedNestingThrows)
{
    auto parent = makeGrid(/* dims = */ {2, 1, 1}, flat);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, /* remove_ij_boundary = */ false);

    useConformingBuilder(grid, parent);

    // On the shared i=1 face, A is finer in K (6 vs 3) and B is finer in J (4 vs 2).
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{1,2,6}, {1,4,3}},
                                                 /* startIJK_vec = */      {{0,0,0}, {1,0,0}},
                                                 /* endIJK_vec = */        {{1,1,1}, {2,1,1}},
                                                 /* lgr_name_vec = */      {"A", "B"}),
                      std::logic_error);
    BOOST_CHECK_EQUAL(grid.maxLevel(), 0);
}

// ---------------------------------------------------------------------------
// Graded boxes next to other boxes. A 2x2x2 grid of unit columns, layers 2 thick.
// ---------------------------------------------------------------------------
namespace
{

// Two parent cells, each split 0.3 / 0.7.
Opm::Refinement::AxisSubdivision twoParentsSplit03()
{
    return { {0, 0, 1, 1}, {0.0, 0.3, 0.0, 0.3}, {0.3, 1.0, 0.3, 1.0} };
}

// Columns as opm-common builds them from N*FIN counts and H*FIN widths: cumulative
// shares over their sum within each parent cell.
Opm::Refinement::AxisSubdivision columnsFromWidths(const std::vector<int>& counts,
                                                   const std::vector<double>& widths)
{
    Opm::Refinement::AxisSubdivision sub;
    std::size_t column = 0;
    for (std::size_t parent = 0; parent < counts.size(); ++parent) {
        double total = 0.0;
        for (int c = 0; c < counts[parent]; ++c) {
            total += widths[column + c];
        }
        double cumulative = 0.0;
        for (int c = 0; c < counts[parent]; ++c, ++column) {
            sub.parentOffset.push_back(static_cast<int>(parent));
            sub.fracLo.push_back(cumulative / total);
            cumulative += widths[column];
            sub.fracHi.push_back((c + 1 == counts[parent]) ? 1.0 : cumulative / total);
        }
    }
    return sub;
}

// Refines a fresh 2x2x2 grid; returns it for the checks.
std::unique_ptr<Dune::CpGrid> refineGradedPair(const std::vector<Opm::Refinement::BlockRefinement>& requests,
                                               double* volumeBefore = nullptr)
{
    auto parent = makeGrid({2, 2, 2}, flat);
    auto grid = std::make_unique<Dune::CpGrid>();
    auto rawParent = parent.raw();
    grid->processEclipseFormat(rawParent, false);
    if (volumeBefore) {
        *volumeBefore = totalVolume(*grid);
    }
    useConformingBuilder(*grid, parent);
    try {
        grid->addLgrsUpdateLeafView(requests);
    }
    catch (const std::invalid_argument& e) {
        BOOST_CHECK_MESSAGE(std::string(e.what()).find("Touching graded boxes are not supported")
                            != std::string::npos, "unexpected refusal: " << e.what());
        BOOST_CHECK_EQUAL(grid->maxLevel(), 0);     // refused before any level grid
        return nullptr;
    }
    return grid;
}

} // anonymous namespace

// A graded box sharing a face with another box is refused before any level grid is
// built, also when the graded direction's cellsPerDim is 0 or differs from the
// other box's.
BOOST_AUTO_TEST_CASE(gradedBoxSharingAFaceIsRefused)
{
    // A is the column I = 0, B the column I = 1: they share the face x = 1.
    auto a = box("A", {2,2,2}, {0,0,0}, {1,2,2});
    a.subdivision[1] = twoParentsSplit03();
    const auto b = box("B", {2,2,2}, {1,0,0}, {2,2,2});
    BOOST_CHECK(!refineGradedPair({a, b}));

    auto unset = a;
    unset.cellsPerDim = {2,0,2};
    BOOST_CHECK(!refineGradedPair({unset, b}));

    // Two and three columns in the two J parents: 5 columns, 5/2 = 2 per parent.
    auto mixed = a;
    mixed.subdivision[1] = { {0, 0, 1, 1, 1}, {0.0, 0.3, 0.0, 0.2, 0.5}, {0.3, 1.0, 0.2, 0.5, 1.0} };
    BOOST_CHECK(!refineGradedPair({mixed, box("B", {2,3,2}, {1,0,0}, {2,2,2})}));
}

// A graded box sharing only an edge with another box is refused: the shared points
// would not merge (columns of 0.3 + 0.3 + 0.3 are not exactly thirds) or would not
// match (0.3 / 0.7 against halves).
BOOST_AUTO_TEST_CASE(gradedBoxSharingAnEdgeIsRefused)
{
    // A is the column (0,0), B the column (1,1): they share the vertical edge on the
    // pillar (1,1).
    auto thirds = box("A", {2,2,3}, {0,0,0}, {1,1,2});
    thirds.subdivision[2] = columnsFromWidths({3, 3}, {0.3, 0.3, 0.3, 0.3, 0.3, 0.3});
    BOOST_CHECK(!refineGradedPair({thirds, box("B", {2,2,3}, {1,1,0}, {2,2,2})}));

    auto split = box("A", {2,2,2}, {0,0,0}, {1,1,2});
    split.subdivision[2] = twoParentsSplit03();
    BOOST_CHECK(!refineGradedPair({split, box("B", {2,2,2}, {1,1,0}, {2,2,2})}));
}

// A graded box meeting another box at a corner only is built: the corner is a corner
// of level zero.
BOOST_AUTO_TEST_CASE(gradedBoxMeetingAtACornerBuilds)
{
    auto a = box("A", {2,2,2}, {0,0,0}, {1,1,1});
    for (int d = 0; d < 3; ++d) {
        a.subdivision[d] = { {0, 0}, {0.0, 0.3}, {0.3, 1.0} };
    }
    double volumeBefore = 0.0;
    const auto grid = refineGradedPair({a, box("B", {2,2,2}, {1,1,1}, {2,2,2})},
                                       &volumeBefore);
    BOOST_REQUIRE(grid);
    BOOST_CHECK_EQUAL(grid->maxLevel(), 2);
    BOOST_CHECK_EQUAL(grid->size(0), 8 - 2 + 2*8);
    BOOST_CHECK_CLOSE(totalVolume(*grid), volumeBefore, 1e-10);
    BOOST_CHECK_SMALL(maxClosure(*grid), 1e-12);
}

// ---------------------------------------------------------------------------
// Compatible sub-face mosaic: two boxes meeting on a shared
// face with DIFFERENT but compatible in-face subdivisions (one a multiple of the
// other). The coarser side's interface cells become >6-face hexes tiled by the
// finer side's sub-faces. Implemented in LeafGridAssembler; these lock it in.
// ---------------------------------------------------------------------------

// STACKED pair (shared horizontal face): in-face directions I,J. Box LO is x4 in
// BOTH I and J, box HI is x2 (compatible 4-vs-2 in both). This mirrors the
// TLGR_VSTACK_HCOMPAT flow deck (TOP 4x4x2 over BOT 2x2x2). A horizontal
// interface never carries a fault, so this is the pure mosaic case.
BOOST_AUTO_TEST_CASE(stackedCompatibleInFaceMosaicBuilds)
{
    auto parent = makeGrid({2, 2, 2}, flat);
    Dune::CpGrid grid;
    auto raw = parent.raw();
    grid.processEclipseFormat(raw, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, parent);

    // Upper slab k in [0,1): x4 in I,J, x2 in K. Lower slab k in [1,2): x2 all.
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{4,4,2}, {2,2,2}},
                               /* startIJK_vec = */      {{0,0,0}, {0,0,1}},
                               /* endIJK_vec = */        {{2,2,1}, {2,2,2}},
                               /* lgr_name_vec = */      {"LO", "HI"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    // LO: 4 parents x (4*4*2)=32 -> 128. HI: 4 parents x (2*2*2)=8 -> 32.
    BOOST_CHECK_EQUAL(grid.size(0), 128 + 32);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
    BOOST_CHECK_SMALL(maxClosure(grid), 1e-9);        // every cell closed
    checkEveryInteriorFaceTwoSided(grid);               // no hanging faces

    // The interface couples the two boxes: x4-vs-x2 in I and J -> each coarser
    // top cell of HI is met by 2x2 = 4 finer faces. Count finer->coarser
    // connections across the k=0/k=1 parent boundary.
    int crossConnections = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (!is.neighbor() || !is.inside().hasFather() || !is.outside().hasFather()) {
                continue;
            }
            // kin and kout are the k of the level-zero parents. Both boxes refine their
            // whole slab, so the parents differ in k only where the two boxes meet: the
            // plane at depth 2 between parent layers k = 0 and k = 1.
            const int kin  = ijkOf(grid, is.inside().index())[2];
            const int kout = ijkOf(grid, is.outside().index())[2];
            if (kin != kout) {
                ++crossConnections;
            }
        }
    }
    // LO's bottom sub-layer is 8x8 (I,J) = 64 sub-faces, each a connection,
    // counted from both sides.
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

    // Leaf vertices: no level-zero corner is outside the boxes; both level grids'
    // corners, less the 5 x 5 of HI's top face at depth 2, which are also corners of
    // LO's bottom face.
    const int expectedVertexCount = grid.size(/* level = */ 1, /* codim = */ 3)
                                  + grid.size(/* level = */ 2, /* codim = */ 3) - 5*5;
    BOOST_CHECK_EQUAL(grid.size(/* codim = */ 3), expectedVertexCount);
}

// SIDE-BY-SIDE pair (shared vertical face), no fault: in-face directions J,K.
// Box A is x4 in K, box B is x2 in K (compatible); equal in J. Confirms the
// mosaic also works on a vertical interface when there is no fault.
BOOST_AUTO_TEST_CASE(sideBySideCompatibleInFaceMosaicBuilds)
{
    auto parent = makeGrid({2, 2, 2}, flat);
    Dune::CpGrid grid;
    auto raw = parent.raw();
    grid.processEclipseFormat(raw, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, parent);

    // Box A i in [0,1): x4 in K, x2 in I,J. Box B i in [1,2): x2 all. Shared i=1
    // vertical face; in-face J equal, in-face K differs 4 vs 2.
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,4}, {2,2,2}},
                               /* startIJK_vec = */      {{0,0,0}, {1,0,0}},
                               /* endIJK_vec = */        {{1,2,2}, {2,2,2}},
                               /* lgr_name_vec = */      {"A", "B"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    // A: 4 parents x (2*2*4)=16 -> 64. B: 4 parents x 8 -> 32.
    BOOST_CHECK_EQUAL(grid.size(0), 64 + 32);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
    BOOST_CHECK_SMALL(maxClosure(grid), 1e-9);
    checkEveryInteriorFaceTwoSided(grid);
}

// ---------------------------------------------------------------------------
// Two refined boxes meeting side-by-side across a FAULT. The
// shared vertical face coincides with the fault throw, so the two refined sides'
// sub-faces are vertically staggered. The faulted-boundary machinery maps a
// "refined away" neighbour to the neighbour box's child cell, assembling the
// staggered interface conformally. Equal subdivisions both sides.
// ---------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(faceSharingAcrossFaultEqualSubdiv)
{
    auto depth = [](int i_, int, int k_) {
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;     // fault between i=1 and i=2
    };
    auto parent = makeGrid({4, 2, 2}, depth);
    Dune::CpGrid grid;
    auto raw = parent.raw();
    grid.processEclipseFormat(raw, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, parent);

    // Box A i in [0,2), box B i in [2,4); shared i=2 face IS the fault. Equal
    // subdivisions both sides.
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}, {2,2,2}},
                               /* startIJK_vec = */      {{0,0,0}, {2,0,0}},
                               /* endIJK_vec = */        {{2,2,2}, {4,2,2}},
                               /* lgr_name_vec = */      {"A", "B"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
    BOOST_CHECK_SMALL(maxClosure(grid), 1e-9);        // every cell closed
    checkEveryInteriorFaceTwoSided(grid);               // no hanging faces

    // The two boxes are actually connected across the fault: there exist interior
    // faces pairing an A-child (i<2) with a B-child (i>=2).
    int abConnections = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (!is.neighbor() || !is.inside().hasFather() || !is.outside().hasFather()) {
                continue;
            }
            const int ii = ijkOf(grid, is.inside().index())[0];
            const int oi = ijkOf(grid, is.outside().index())[0];
            if ((ii < 2) != (oi < 2)) {
                ++abConnections;
            }
        }
    }
    // On x = 2 the A sub-layers (0, 1, 2, 3, 4) and the B sub-layers (0.6, 1.6, 2.6,
    // 3.6, 4.6) overlap in 7 pieces between depths 0.6 and 4, with 4 children in y:
    // 28 faces, counted twice (once from each side).
    BOOST_CHECK_EQUAL(abConnections, 7*4*2);
}

// Two refined boxes across a FAULT with DIFFERENT but compatible in-face
// subdivisions. Box A (i<2) is x4 in K, box B (i>=2) is x2 in K (equal in J); the
// shared i=2 vertical face is the fault throw. The finer side (A) emits every
// interface sub-face -- carrying the throw -- and references box B's coarser child
// that contains it, so the leaf is conformal.
BOOST_AUTO_TEST_CASE(faceSharingAcrossFaultCompatibleSubdiv)
{
    auto depth = [](int i_, int, int k_) {
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;     // fault between i=1 and i=2
    };
    auto parent = makeGrid({4, 2, 2}, depth);
    Dune::CpGrid grid;
    auto raw = parent.raw();
    grid.processEclipseFormat(raw, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, parent);

    // Box A i in [0,2): x4 in K, x2 in I,J. Box B i in [2,4): x2 all. Shared i=2
    // face IS the fault; in-face J equal, in-face K differs 4 vs 2 (compatible).
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,4}, {2,2,2}},
                               /* startIJK_vec = */      {{0,0,0}, {2,0,0}},
                               /* endIJK_vec = */        {{2,2,2}, {4,2,2}},
                               /* lgr_name_vec = */      {"A", "B"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
    BOOST_CHECK_SMALL(maxClosure(grid), 1e-9);        // every cell closed
    checkEveryInteriorFaceTwoSided(grid);               // no hanging faces

    int abConnections = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (!is.neighbor() || !is.inside().hasFather() || !is.outside().hasFather()) {
                continue;
            }
            const int ii = ijkOf(grid, is.inside().index())[0];
            const int oi = ijkOf(grid, is.outside().index())[0];
            if ((ii < 2) != (oi < 2)) {
                ++abConnections;
            }
        }
    }
    // On x = 2 the A sub-layers (0, 0.5, ..., 4) and the B sub-layers (0.6, 1.6, 2.6,
    // 3.6, 4.6) overlap in 10 pieces between depths 0.6 and 4, with 4 children in y:
    // 40 faces, counted twice (once from each side).
    BOOST_CHECK_EQUAL(abConnections, 10*4*2);
}

// Nested LGR (a box whose parent is another LGR). A nested box must lie inside
// its parent LGR without touching its boundary; one that touches it is refused
// by the leaf assembly, after the nested level grids were built.
BOOST_AUTO_TEST_CASE(nestedRefinementReachesLeafBoundary)
{
    auto parent = makeGrid({4, 2, 2}, flat);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);

    useConformingBuilder(grid, parent);

    // LGR1: 2x1x1 GLOBAL parents -> local refined dims 4x2x2.
    // NEST1 {0,0,0}-{2,2,2} touches LGR1's boundary (startIJK has zeros), so it
    // is rejected by the containment check (only fully-interior nesting yet).
    const auto isNestedLeafBoundary = [](const std::logic_error& e) {
        return std::string(e.what()).find("touches the boundary")
               != std::string::npos;
    };
    BOOST_CHECK_EXCEPTION(
        grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */        {{2,2,2}, {2,2,2}},
                                   /* startIJK_vec = */             {{1,0,0}, {0,0,0}},
                                   /* endIJK_vec = */               {{3,1,1}, {2,2,2}},
                                   /* lgr_name_vec = */             {"LGR1", "NEST1"},
                                   /* lgr_parent_grid_name_vec = */ {"GLOBAL", "LGR1"}),
        std::logic_error, isNestedLeafBoundary);
}

// The builder appends level grids in request order and resolves each box's
// parent by name, so a child must follow its parent. A child listed before its
// parent must be rejected with an actionable message (the simulator's
// topological sort is what guarantees this never happens in practice).
BOOST_AUTO_TEST_CASE(nestedChildBeforeParentThrows)
{
    auto parent = makeGrid({4, 2, 2}, flat);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);

    useConformingBuilder(grid, parent);

    const auto isUnbuiltParent = [](const std::logic_error& e) {
        return std::string(e.what()).find("has not been built yet")
               != std::string::npos;
    };
    BOOST_CHECK_EXCEPTION(
        grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */        {{2,2,2}, {2,2,2}},
                                   /* startIJK_vec = */             {{0,0,0}, {1,0,0}},
                                   /* endIJK_vec = */               {{2,2,2}, {3,1,1}},
                                   /* lgr_name_vec = */             {"NEST1", "LGR1"},
                                   /* lgr_parent_grid_name_vec = */ {"LGR1", "GLOBAL"}),
        std::logic_error, isUnbuiltParent);
}

// A CARFIN whose parent is another LGR (nested), with the child box strictly
// interior to the parent LGR - it touches none of the parent's refined boundary
// faces. This is the simplest nesting case: the child's leaf stitching is
// entirely within the parent LGR's refined region, with no interaction with
// GLOBAL. Checked (mirroring endToEndSingleBox):
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
    auto parent = makeGrid({4, 3, 3}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, parent);

    // LGR1: GLOBAL cells i:1-2, j:1, k:1 (2x1x1) refined 3x3x3 -> LGR1-local 6x3x3 = 54.
    // NEST1: LGR1-local i:1-4, j:1, k:1 (4 LGR1 cells) refined 2x2x2 = 32 -> strictly
    //        interior to LGR1 (touches no LGR1 boundary face), parent grid "LGR1".
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */        {{3,3,3}, {2,2,2}},
                               /* startIJK_vec = */             {{1,1,1}, {1,1,1}},
                               /* endIJK_vec = */               {{3,2,2}, {5,2,2}},
                               /* lgr_name_vec = */             {"LGR1", "NEST1"},
                               /* lgr_parent_grid_name_vec = */ {"GLOBAL", "LGR1"});

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
    checkEveryInteriorFaceTwoSided(grid);

    // Global ids unique over leaf cells and points (delegated across 3 levels).
    checkUniqueGlobalIds(grid);

    // NEST1 refines two LGR1 cells of one level-zero cell, which stableCellId() cannot tell apart.
    BOOST_CHECK_THROW(grid.stableCellId(), std::logic_error);
}

BOOST_AUTO_TEST_CASE(faultInsideBoxEndToEnd)
{
    // Fault between i=1 and i=2 columns, throw 0.6; box covers it.
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;
    };
    auto parent = makeGrid({4, 2, 2}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, parent);

    // The box {1..3} has the fault *inside* (between i=1 and i=2) and
    // unfaulted boundaries (i=0|1 and i=3|4 columns are unfaulted... note
    // the fault plane is i=2, the box boundary planes are i=1 and i=3,
    // both unfaulted since the throw is constant on each side).
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}},
                               /* startIJK_vec = */      {{1,0,0}},
                               /* endIJK_vec = */        {{3,2,2}},
                               /* lgr_name_vec = */      {"LGR1"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);

    // Both-sided intersection consistency across the refined fault.
    checkEveryInteriorFaceTwoSided(grid);

    // Across the fault inside the box (parents i = 1 and i = 2) the sub-layers 0 to 4 and
    // 0.6 to 4.6 overlap in 7 pieces between depths 0.6 and 4, with 4 children in y:
    // 28 faces, counted twice (once from each side).
    int acrossFault = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (!is.neighbor()) {
                continue;
            }
            const int iIn = ijkOf(grid, is.inside().index())[0];
            const int iOut = ijkOf(grid, is.outside().index())[0];
            if (((iIn == 1) && (iOut == 2)) || ((iIn == 2) && (iOut == 1))) {
                ++acrossFault;
            }
        }
    }
    BOOST_CHECK_EQUAL(acrossFault, 7*4*2);
}

BOOST_AUTO_TEST_CASE(faultAtBoxBoundaryBuilds)
{
    // Fault between i=1 and i=2 (throw 0.6); box boundary right on the fault
    // plane. The boundary parents each connect to two coarse cells, so the
    // boundary is rebuilt as split faces. Assert the leaf is a valid
    // polyhedral grid (builds, volume conserved, every interior face two-sided,
    // refined boundary cells connected to the coarse neighbours).
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;
    };
    auto parent = makeGrid({4, 2, 2}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, parent);

    BOOST_REQUIRE_NO_THROW(grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}},
                                                      /* startIJK_vec = */      {{2,0,0}},
                                                      /* endIJK_vec = */        {{4,2,2}},
                                                      /* lgr_name_vec = */      {"LGR1"}));
    BOOST_CHECK_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);

    // Geometric validity: every leaf cell is closed, i.e. the sum of
    // area * outward-unit-normal over all its faces is zero (divergence
    // theorem). This validates the synthetic split faces' areas AND normal
    // orientation, independent of how the interface is subdivided. Also: every
    // interior face is two-sided, and the box connects to coarse neighbours.
    std::map<std::pair<int,int>, int> pairCount;
    int refinedCoarseConnections = 0;
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
                    ++refinedCoarseConnections;
                }
            }
        }
        BOOST_CHECK_SMALL(closure.two_norm(), 1e-9);
    }
    for (const auto& [cells, count] : pairCount) {
        BOOST_CHECK_EQUAL(count % 2, 0); // each shared face seen from both sides
    }
    // On x = 2 each of the 4 (j) box columns has 5 faces toward column i = 1: sub-layers
    // 0 and 2 face one coarse cell, sub-layer 1 two, and sub-layer 3 one (the domain
    // below depth 4). 20 faces, counted twice (once from each side).
    BOOST_CHECK_EQUAL(refinedCoarseConnections, 5*4*2);

    // Split-face corners exist only on the leaf; they still get unique ids.
    checkUniqueGlobalIds(grid);
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
    auto parent = makeGrid({2, 4, 2}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, parent);

    BOOST_REQUIRE_NO_THROW(grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}},
                                                      /* startIJK_vec = */      {{0,2,0}},
                                                      /* endIJK_vec = */        {{2,4,2}},
                                                      /* lgr_name_vec = */      {"LGR1"}));
    checkLeafClosedAndConnectionsToCoarse(grid, volumeBefore);
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
    auto parent = makeGrid({4, 2, 2}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, parent);

    BOOST_REQUIRE_NO_THROW(grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}},
                                                      /* startIJK_vec = */      {{0,0,0}},
                                                      /* endIJK_vec = */        {{2,2,2}},
                                                      /* lgr_name_vec = */      {"LGR1"}));
    checkLeafClosedAndConnectionsToCoarse(grid, volumeBefore);
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
    auto parent = makeGrid({4, 2, 3}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, parent);

    BOOST_REQUIRE_NO_THROW(grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}},
                                                      /* startIJK_vec = */      {{2,0,0}},
                                                      /* endIJK_vec = */        {{4,2,3}},
                                                      /* lgr_name_vec = */      {"LGR1"}));
    checkLeafClosedAndConnectionsToCoarse(grid, volumeBefore);
}

// A throw of half a layer puts the thrown depths on the box's sub-layer boundary, equal
// only up to rounding (feet converted to metres, dipping layers).
BOOST_AUTO_TEST_CASE(faultThrowOnSubLayerBoundaryBuilds)
{
    auto depth = [](int i_, int, int k_) {
        const double ft = 0.3048;
        const double z = ft*(8325.0 + 6.0*(cellOf(i_) + sideOf(i_)) + 20.0*(cellOf(k_) + sideOf(k_)));
        return (cellOf(i_) >= 2) ? z + 10.0*ft : z;
    };
    auto parent = makeGrid(/* dims = */ {4, 2, 4}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, parent);

    BOOST_REQUIRE_NO_THROW(grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{3,3,2}},
                                                      /* startIJK_vec = */      {{0,0,1}},
                                                      /* endIJK_vec = */        {{2,2,3}},
                                                      /* lgr_name_vec = */      {"LGR1"}));
    checkLeafClosedAndConnectionsToCoarse(grid, volumeBefore);
}

namespace
{

// Face area between cells whose parent is at i = iIn, layers [k0, k1), and cells
// whose parent is at i = iOut: the connection across a fault at an i-face.
double areaAcrossIFace(const Dune::CpGrid& grid, int iIn, int iOut, int k0, int k1)
{
    double area = 0.0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        const auto in = ijkOf(grid, element.index());
        if (in[0] != iIn || in[2] < k0 || in[2] >= k1) {
            continue;
        }
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (is.neighbor() && ijkOf(grid, is.outside().index())[0] == iOut) {
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
    auto parent = makeGrid({4, 2, 8}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);
    const double areaBefore = areaAcrossIFace(grid, /* iIn = */ 2, /* iOut = */ 1,
                                              /* k0 = */ 3, /* k1 = */ 5);
    BOOST_REQUIRE_GT(areaBefore, 0.0);

    useConformingBuilder(grid, parent);

    BOOST_REQUIRE_NO_THROW(grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}},
                                                      /* startIJK_vec = */      {{2,0,3}},
                                                      /* endIJK_vec = */        {{4,2,5}},
                                                      /* lgr_name_vec = */      {"LGR1"}));
    checkLeafClosedAndConnectionsToCoarse(grid, volumeBefore);
    BOOST_CHECK_CLOSE(areaAcrossIFace(grid, /* iIn = */ 2, /* iOut = */ 1,
                                      /* k0 = */ 3, /* k1 = */ 5), areaBefore, 1e-8);
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
    auto parent = makeGrid({4, 2, 8}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);

    useConformingBuilder(grid, parent);

    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,3,2}, {2,2,3}},
                                                 /* startIJK_vec = */      {{0,0,6}, {2,0,3}},
                                                 /* endIJK_vec = */        {{2,2,8}, {4,2,5}},
                                                 /* lgr_name_vec = */      {"A", "B"}),
                      std::invalid_argument);
}

// The leaf assembly refuses after the level grids are built; the grid must be left
// unrefined and accept a valid request afterwards.
BOOST_AUTO_TEST_CASE(failedRefinementLeavesGridUnchanged)
{
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 5.0 : base;
    };
    auto parent = makeGrid(/* dims = */ {4, 2, 8}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const int cellsBefore = grid.leafGridView().size(0);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, parent);

    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,3,2}, {2,2,3}},
                                                 /* startIJK_vec = */      {{0,0,6}, {2,0,3}},
                                                 /* endIJK_vec = */        {{2,2,8}, {4,2,5}},
                                                 /* lgr_name_vec = */      {"A", "B"}),
                      std::invalid_argument);
    BOOST_CHECK_EQUAL(grid.currentData().size(), 1u);
    BOOST_CHECK_EQUAL(grid.maxLevel(), 0);
    BOOST_CHECK_EQUAL(grid.leafGridView().size(0), cellsBefore);
    BOOST_CHECK(grid.currentData().front()->getParentToChildren().empty());

    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,3,2}},
                               /* startIJK_vec = */      {{0,0,6}},
                               /* endIJK_vec = */        {{2,2,8}},
                               /* lgr_name_vec = */      {"A"});
    BOOST_CHECK_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
}

BOOST_AUTO_TEST_CASE(faultNotOnBoxBoundaryBuilds)
{
    // Control: a fault exists (between i=0 and i=1) but the box (i=2..3) sits
    // entirely on the down-thrown side, so its boundary toward the coarse
    // neighbour (i=1) is a full, matching face. This must build and
    // conserve volume — it isolates "fault on the boundary" from "fault merely
    // present in the grid".
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 1) ? base + 0.6 : base;
    };
    auto parent = makeGrid({4, 2, 2}, depth);

    Dune::CpGrid grid;
    auto rawParent = parent.raw();
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);

    useConformingBuilder(grid, parent);

    BOOST_CHECK_NO_THROW(grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}},
                                                    /* startIJK_vec = */      {{2,0,0}},
                                                    /* endIJK_vec = */        {{4,2,2}},
                                                    /* lgr_name_vec = */      {"LGR1"}));
    BOOST_CHECK_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-8);
}

// If a cell inside an LGR box is connected across the fault to a cell in another layer
// that is outside the box, the refinement throws, naming the two cells. A box over every
// layer the throw joins is built.
BOOST_AUTO_TEST_CASE(throwInsideBoxJoiningCellsOutsideIsRefused)
{
    // Columns i >= 2 are 0.6 deeper; layers 2 thick.
    auto parent = makeGrid({4, 3, 2}, [](int i_, int, int k_) {
        return 2.0*(cellOf(k_) + sideOf(k_)) + ((cellOf(i_) >= 2) ? 0.6 : 0.0);
    });
    auto rawParent = parent.raw();
    const auto refine = [&](Dune::CpGrid& grid, std::array<int,3> start, std::array<int,3> end) {
        grid.processEclipseFormat(rawParent, false);
        useConformingBuilder(grid, parent);
        grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{3,1,3}},
                                   /* startIJK_vec = */      {start},
                                   /* endIJK_vec = */        {end},
                                   /* lgr_name_vec = */      {"LGR1"});
    };
    const auto namesTheCells = [](const std::string& box, const std::string& outside) {
        return [box, outside](const std::invalid_argument& e) {
            const std::string what = e.what();
            return (what.find("joins its cell " + box + " to the cell " + outside) != std::string::npos)
                && (what.find("must cover every layer") != std::string::npos);
        };
    };
    // Lower layer only: the box cell (3,1,1) (from 1) faces (2,1,2) above the box; upper
    // layer only: the box cell (2,1,2) faces (3,1,1) below it.
    {
        Dune::CpGrid grid;
        BOOST_CHECK_EXCEPTION(refine(grid, {1,0,0}, {3,3,1}), std::invalid_argument,
                              namesTheCells("(3,1,1)", "(2,1,2)"));
    }
    {
        Dune::CpGrid grid;
        BOOST_CHECK_EXCEPTION(refine(grid, {1,0,1}, {3,3,2}), std::invalid_argument,
                              namesTheCells("(2,1,2)", "(3,1,1)"));
    }
    // Every layer: the faces across the fault keep the level-zero area, 3 x (4 - 0.6).
    Dune::CpGrid grid;
    grid.processEclipseFormat(rawParent, false);
    const double volumeBefore = totalVolume(grid);
    BOOST_REQUIRE_CLOSE(areaAcrossIFace(grid, /* iIn = */ 2, /* iOut = */ 1,
                                        /* k0 = */ 0, /* k1 = */ 2), 3.0*3.4, 1e-10);
    Dune::CpGrid refined;
    refine(refined, {1,0,0}, {3,3,2});
    checkLeafClosedAndConnectionsToCoarse(refined, volumeBefore);
    BOOST_CHECK_CLOSE(areaAcrossIFace(refined, /* iIn = */ 2, /* iOut = */ 1,
                                      /* k0 = */ 0, /* k1 = */ 2), 3.0*3.4, 1e-10);
}

BOOST_AUTO_TEST_CASE(processGrdeclSplitsFaultedCoarseFineColumnPair)
{
    // The faces at a faulted box boundary are built by running process_grdecl on
    // a small grid (coarse neighbour + refined box boundary) and reading the
    // split connections. This checks that mechanism in isolation.
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
    auto g = makeGrid({2, 1, 2}, depth);
    auto raw = g.raw();

    struct processed_grid out;
    const int ok = process_grdecl(/*pinchActive=*/1, /*edge_conformal=*/0, /*tol=*/1e-6, &raw,
                                  /*is_aquifer_cell=*/nullptr, &out);
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
// resampling with bitwise-identical results.
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
    TestGrdecl g = makeGrid(dims, depth);
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
    const auto parentAlong = [&](int leafCell) {
        return ijkOf(grid, leafCell)[axis];
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

    useConformingBuilder(grid, parent);

    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}, {3,2,2}},
                               /* startIJK_vec = */      {{0,0,0}, {2,0,0}},
                               /* endIJK_vec = */        {{2,2,2}, {4,2,2}},
                               /* lgr_name_vec = */      {"A", "B"});

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

    useConformingBuilder(grid, parent);

    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,4,4}, {2,2,2}},
                               /* startIJK_vec = */      {{0,0,0}, {2,0,0}},
                               /* endIJK_vec = */        {{2,2,2}, {4,2,2}},
                               /* lgr_name_vec = */      {"A", "B"});

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

    useConformingBuilder(grid, parent);

    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}, {2,2,3}},
                               /* startIJK_vec = */      {{0,0,0}, {0,0,2}},
                               /* endIJK_vec = */        {{2,2,2}, {2,2,4}},
                               /* lgr_name_vec = */      {"LOW", "HIGH"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    BOOST_CHECK_EQUAL(grid.size(0), 8*8 + 8*12);
    BOOST_CHECK_CLOSE(totalVolume(grid), volumeBefore, 1e-6);
    BOOST_CHECK_SMALL(maxClosure(grid), 1e-9);
    checkEveryInteriorFaceTwoSided(grid);
    // rx*ry = 4 connections per parent pair, 2x2 pairs, counted twice.
    BOOST_CHECK_EQUAL(crossBoxConnections(grid, 2, 2), 16*2);
    BOOST_CHECK_GT(minVertexPairDistance(grid), 1e-3);
}

// Boxes stacked in a column get bitwise the same refined pillars, also where the
// layers' thickness changes with depth: the pillars come from the parent pillar lines.
BOOST_AUTO_TEST_CASE(skewPillarStackedBoxesShareRefinedPillars)
{
    const auto coordOf = [](const auto& parent, int k0, int k1) {
        const auto req = box("B", {2, 2, 1}, {0, 0, k0}, {2, 2, k1});
        return Opm::Refinement::refineBlock(parent.dims, parent.coord.data(),
                                            parent.zcorn.data(), nullptr, req).coord;
    };
    const auto parent = makeSkewPillarGrid({2, 2, 4}, tiltedDepth);
    BOOST_CHECK(coordOf(parent, 0, 2) == coordOf(parent, 2, 4));

    // Layers that thicken with depth, more on some pillars.
    const auto thickening = makeSkewPillarGrid({2, 2, 4}, [](int i_, int j_, int k_) {
        const int k = (k_ + 1)/2;
        return 2.0*k + 0.3*((i_ + 1)/2 + (j_ + 1)/2)*k*k;
    });
    BOOST_CHECK(coordOf(thickening, 0, 2) == coordOf(thickening, 2, 4));
}
