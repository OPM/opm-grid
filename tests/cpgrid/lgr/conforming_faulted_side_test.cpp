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

#define BOOST_TEST_MODULE ConformingFaultedSideTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/refinement/conforming/ConformingBlockBuilder.hpp>
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>
#include <opm/grid/cpgrid/refinement/RefinementRequest.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <dune/common/parallel/mpihelper.hh>

#include <algorithm>
#include <array>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

// A refined box whose side a fault crosses: the side is rebuilt by processing the box and
// the cells across it as one corner-point grid. These cases check that the rebuilt side
// gives each cell pair one face, uses the box's own corners and keeps every cell closed.

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

using Pillar = std::function<std::array<double,6>(int, int)>;
using Depth = std::function<double(int, int, int)>;

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
int nodeOf(int doubled) { return cellOf(doubled) + sideOf(doubled); }

std::array<double,6> vertical(int i, int j)
{
    return { double(i), double(j), 0.0, double(i), double(j), 100.0 };
}

// The direction of a pillar changes with j: pillars are not parallel.
std::array<double,6> fanning(int i, int j)
{
    return { double(i), double(j), 0.0, i + 0.5*j, 1.3*j, 100.0 };
}

// J runs toward -y: a left-handed grid.
std::array<double,6> leftHanded(int i, int j)
{
    return { double(i), -double(j), 0.0, double(i), -double(j), 100.0 };
}

TestGrdecl makeGrid(const std::array<int,3>& dims, const Pillar& pillar, const Depth& depth)
{
    TestGrdecl g;
    g.dims = dims;
    const auto& [nx, ny, nz] = dims;
    g.coord.resize(6*(nx + 1)*(ny + 1));
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            const auto p = pillar(i, j);
            std::copy(p.begin(), p.end(), &g.coord[6*(static_cast<std::size_t>(j)*(nx + 1) + i)]);
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

Opm::Refinement::BlockRefinement box(const std::string& name, std::array<int,3> cellsPerDim,
                                     std::array<int,3> start, std::array<int,3> end)
{
    Opm::Refinement::BlockRefinement request;
    request.name = name;
    request.cellsPerDim = cellsPerDim;
    request.startIJK = start;
    request.endIJK = end;
    return request;
}

void refine(Dune::CpGrid& grid, const TestGrdecl& g, std::vector<Opm::Refinement::BlockRefinement> requests,
            bool edgeConformal = false, const std::vector<int>* builderActnum = nullptr)
{
    grid.setLgrBackend(Opm::Refinement::Backend::Conforming);
    grid.setRefinementBuilder(std::make_shared<Opm::Refinement::ConformingBlockBuilder>(
        g.dims, g.coord, g.zcorn, builderActnum ? *builderActnum : g.actnum, edgeConformal));
    grid.addLgrsUpdateLeafView(std::move(requests));
}

double totalVolume(const Dune::CpGrid& grid)
{
    double volume = 0.0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        volume += element.geometry().volume();
    }
    return volume;
}

// Largest |sum of area * outward normal| of a leaf cell, relative to its surface.
double worstClosure(const Dune::CpGrid& grid)
{
    double worst = 0.0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        Dune::FieldVector<double,3> sum(0.0);
        double area = 0.0;
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            auto n = is.centerUnitOuterNormal();
            n *= is.geometry().volume();
            sum += n;
            area += is.geometry().volume();
        }
        worst = std::max(worst, sum.two_norm() / area);
    }
    return worst;
}

// Faces between each pair of leaf cells, each face counted once.
std::map<std::pair<int,int>, int> facesPerPair(const Dune::CpGrid& grid)
{
    std::map<std::pair<int,int>, int> pairs;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (is.neighbor() && (is.inside().index() < is.outside().index())) {
                ++pairs[{is.inside().index(), is.outside().index()}];
            }
        }
    }
    return pairs;
}

int maxFacesPerPair(const Dune::CpGrid& grid)
{
    int most = 0;
    for (const auto& [cells, count] : facesPerPair(grid)) {
        most = std::max(most, count);
    }
    return most;
}

// Cells i >= faultI are `throwBy` deeper; layers 2 thick.
Depth faulted(int faultI, double throwBy)
{
    return [faultI, throwBy](int i_, int, int k_) {
        return 2.0*nodeOf(k_) + ((cellOf(i_) >= faultI) ? throwBy : 0.0);
    };
}

// Layer boundaries that are not on one line per column.
Depth wedge(int faultI, double throwBy)
{
    return [faultI, throwBy](int i_, int j_, int k_) {
        const int k = nodeOf(k_);
        return 2.0*k + 0.7*nodeOf(j_)*k*k + ((cellOf(i_) >= faultI) ? throwBy : 0.0);
    };
}

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(MPIFixture);

// The helper grid cuts the coarse neighbour at the box's sub-layers; the pieces of one
// cell pair are one face.
BOOST_AUTO_TEST_CASE(oneFacePerCellPair)
{
    const auto g = makeGrid({4, 2, 2}, vertical, faulted(2, 0.6));
    Dune::CpGrid grid;
    grid.processEclipseFormat(g.raw(), false);
    const double volume = totalVolume(grid);
    refine(grid, g, {box("LGR1", {2,2,2}, {2,0,0}, {4,2,2})});

    BOOST_CHECK_EQUAL(maxFacesPerPair(grid), 1);
    BOOST_CHECK_SMALL(worstClosure(grid), 1e-12);
    BOOST_CHECK_CLOSE(totalVolume(grid), volume, 1e-10);

    // A refined child spanning two of the coarse cell's sub-layers keeps the whole area.
    double refinedToCoarse = 0.0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (is.neighbor() && element.hasFather() && !is.outside().hasFather()) {
                refinedToCoarse += is.geometry().volume();
            }
        }
    }
    // The fault plane x = 2 overlaps the coarse cells over z in [0.6, 4] along y in [0, 2].
    BOOST_CHECK_CLOSE(refinedToCoarse, (4.0 - 0.6)*2.0, 1e-10);
}

// With non-parallel pillars, a faulted low side uses the box's own sub-pillars: the leaf
// has the same corners as with vertical pillars and the same depths.
BOOST_AUTO_TEST_CASE(lowSideUsesTheBoxPillars)
{
    const auto request = box("LGR1", {2,2,2}, {2,0,0}, {4,2,3});
    const auto vertices = [&](const Pillar& pillar) {
        const auto g = makeGrid({4, 2, 3}, pillar, wedge(2, 0.6));
        Dune::CpGrid grid;
        grid.processEclipseFormat(g.raw(), false);
        refine(grid, g, {request});
        BOOST_CHECK_EQUAL(maxFacesPerPair(grid), 1);
        return grid.size(3);
    };
    BOOST_CHECK_EQUAL(vertices(fanning), vertices(vertical));
}

// Corners made where the fault cuts the side have ids of their own.
BOOST_AUTO_TEST_CASE(faultedSideCornersHaveIds)
{
    const auto g = makeGrid({4, 2, 2}, vertical, faulted(2, 0.6));
    Dune::CpGrid grid;
    grid.processEclipseFormat(g.raw(), false);
    refine(grid, g, {box("LGR1", {2,2,2}, {2,0,0}, {4,2,2})});

    BOOST_REQUIRE_GT(static_cast<std::size_t>(grid.size(3)),
                     grid.currentData().back()->cornerHistorySize());
    std::set<Dune::CpGrid::GlobalIdSet::IdType> ids;
    for (const auto& vertex : Dune::vertices(grid.leafGridView())) {
        ids.insert(grid.localIdSet().id(vertex));
    }
    BOOST_CHECK_EQUAL(ids.size(), static_cast<std::size_t>(grid.size(3)));
}

// On the edge of an edge-conformal grid, a box side a J-fault reaches is rebuilt as
// faces toward the domain: no refined cell loses its side.
BOOST_AUTO_TEST_CASE(faultedSideOnTheGridEdge)
{
    const auto depth = [](int, int j_, int k_) {
        return 2.0*nodeOf(k_) + ((cellOf(j_) >= 1) ? 0.6 : 0.0);
    };
    const auto g = makeGrid({3, 2, 2}, vertical, depth);
    Dune::CpGrid grid;
    grid.processEclipseFormat(g.raw(), false, false, /*edge_conformal=*/true);
    const double volume = totalVolume(grid);
    refine(grid, g, {box("LGR1", {2,2,2}, {0,0,0}, {2,2,2})}, /*edgeConformal=*/true);

    BOOST_CHECK_SMALL(worstClosure(grid), 1e-12);
    BOOST_CHECK_CLOSE(totalVolume(grid), volume, 1e-10);
}

// Two boxes share a face; a fault crosses only part of it. Both sides are rebuilt from the
// fault, so the shared faces are made once.
BOOST_AUTO_TEST_CASE(sharedFacePartlyFaulted)
{
    const auto depth = [](int i_, int j_, int k_) {
        return 2.0*nodeOf(k_) + (((cellOf(i_) >= 2) && (cellOf(j_) == 0)) ? 0.6 : 0.0);
    };
    const auto g = makeGrid({4, 2, 2}, vertical, depth);
    Dune::CpGrid grid;
    grid.processEclipseFormat(g.raw(), false);
    const double volume = totalVolume(grid);
    BOOST_REQUIRE_NO_THROW(refine(grid, g, {box("A", {2,2,2}, {0,0,0}, {2,2,2}),
                                            box("B", {2,2,2}, {2,1,0}, {4,2,2})}));

    BOOST_CHECK_EQUAL(maxFacesPerPair(grid), 1);
    BOOST_CHECK_SMALL(worstClosure(grid), 1e-12);
    BOOST_CHECK_CLOSE(totalVolume(grid), volume, 1e-10);
    int boxToBox = 0;
    for (const auto& [cells, count] : facesPerPair(grid)) {
        const auto inBox = [&](int cell, int level) {
            const auto element = Dune::cpgrid::Entity<0>(*grid.currentData().back(), cell, true);
            return element.level() == level;
        };
        if ((inBox(cells.first, 1) && inBox(cells.second, 2))
            || (inBox(cells.first, 2) && inBox(cells.second, 1))) {
            ++boxToBox;
        }
    }
    BOOST_CHECK_GT(boxToBox, 0);
}

// A cell across the faulted side that the grid does not have leaves the refined cells
// facing it with boundary faces, not holes.
BOOST_AUTO_TEST_CASE(missingNeighbourAcrossTheFault)
{
    auto g = makeGrid({4, 1, 2}, vertical, faulted(2, 0.6));
    const std::vector<int> allActive(8, 1);
    g.actnum = allActive;
    g.actnum[1 + 4*1] = 0;                  // (1,0,1), across the box's low side
    Dune::CpGrid grid;
    grid.processEclipseFormat(g.raw(), false);
    refine(grid, g, {box("LGR1", {2,1,2}, {2,0,0}, {4,1,2})}, false, &allActive);

    BOOST_CHECK_SMALL(worstClosure(grid), 1e-12);
}

// On a left-handed grid every face still points out of the cell it is seen from.
BOOST_AUTO_TEST_CASE(leftHandedFaultedSide)
{
    const auto depth = [](int, int j_, int k_) {
        return 2.0*nodeOf(k_) + ((cellOf(j_) >= 1) ? 0.6 : 0.0);
    };
    const auto g = makeGrid({2, 2, 2}, leftHanded, depth);
    Dune::CpGrid grid;
    grid.processEclipseFormat(g.raw(), false);
    refine(grid, g, {box("LGR1", {2,2,2}, {0,1,0}, {2,2,2})});

    for (const auto& element : Dune::elements(grid.leafGridView())) {
        const auto centre = element.geometry().center();
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            auto outward = is.geometry().center();
            outward -= centre;
            BOOST_CHECK_GT(is.centerUnitOuterNormal() * outward, 0.0);
        }
    }
    BOOST_CHECK_SMALL(worstClosure(grid), 1e-12);
}

// Across a fault, a box twice as fine as its neighbour box: each pair of children shares
// one face, and the faces do not overlap.
BOOST_AUTO_TEST_CASE(unequalBoxesAcrossTheFault)
{
    const auto g = makeGrid({4, 2, 2}, vertical, faulted(2, 0.6));
    Dune::CpGrid grid;
    grid.processEclipseFormat(g.raw(), false);
    const double volume = totalVolume(grid);
    refine(grid, g, {box("FINE", {4,4,4}, {0,0,0}, {2,2,2}),
                     box("COARSE", {2,2,2}, {2,0,0}, {4,2,2})});

    BOOST_CHECK_EQUAL(maxFacesPerPair(grid), 1);
    BOOST_CHECK_SMALL(worstClosure(grid), 1e-12);
    BOOST_CHECK_CLOSE(totalVolume(grid), volume, 1e-10);
}

// The same grid written twice: once as is, once with each pillar's end points moved along
// its own line to depths that differ from pillar to pillar. Refined pillars depend on the
// lines only, so a box and the helper grid of its faulted side, or two boxes across a
// fault, compute the same sub-pillars bit for bit either way: the leaf is the same.
BOOST_AUTO_TEST_CASE(faultedSideIndependentOfEndPointDepths)
{
    const auto movedEnds = [](int i, int j) {
        const auto p = fanning(i, j);
        const auto at = [&p](double z) {
            const double t = (z - p[2]) / (p[5] - p[2]);
            return std::array { p[0] + t*(p[3] - p[0]), p[1] + t*(p[4] - p[1]), z };
        };
        const auto top = at(-1.0 - i - 2.0*j);
        const auto bottom = at(100.0 + 3.0*i + j);
        return std::array { top[0], top[1], top[2], bottom[0], bottom[1], bottom[2] };
    };
    const std::vector<std::vector<Opm::Refinement::BlockRefinement>> layouts {
        { box("LGR1", {2,2,2}, {2,0,0}, {4,2,3}) },
        { box("LGR1", {2,2,2}, {0,0,0}, {2,2,3}) },
        { box("LGR1", {2,2,2}, {0,0,0}, {2,2,3}), box("LGR2", {2,2,2}, {2,0,0}, {4,2,3}) },
    };
    for (const auto& requests : layouts) {
        const auto leaf = [&requests](const Pillar& pillar) {
            const auto g = makeGrid({4, 2, 3}, pillar, wedge(2, 0.6));
            Dune::CpGrid grid;
            grid.processEclipseFormat(g.raw(), false);
            refine(grid, g, requests);
            BOOST_CHECK_EQUAL(maxFacesPerPair(grid), 1);
            return std::array { double(grid.size(3)), double(grid.size(0)), totalVolume(grid),
                                worstClosure(grid) };
        };
        const auto written = leaf(fanning);
        const auto moved = leaf(movedEnds);
        BOOST_CHECK_EQUAL(written[0], moved[0]);
        BOOST_CHECK_EQUAL(written[1], moved[1]);
        BOOST_CHECK_CLOSE(written[2], moved[2], 1e-10);
        // Faces on fanning pillars are not planar: the cells close equally well, not exactly.
        BOOST_CHECK_SMALL(written[3] - moved[3], 1e-10);
    }
}
