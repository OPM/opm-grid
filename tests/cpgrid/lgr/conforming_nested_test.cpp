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

#define BOOST_TEST_MODULE ConformingNestedTests
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
#include <stdexcept>
#include <string>
#include <vector>

// A nested box refines cells of another box's level grid. It is assembled with the
// other boxes, so the leaf keeps what they need (faces between boxes, faulted sides)
// and joins the nested box to its parent LGR.

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
int nodeOf(int doubled) { return cellOf(doubled) + sideOf(doubled); }

TestGrdecl makeGrid(const std::array<int,3>& dims, const std::function<double(int,int,int)>& depth)
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

double flat(int, int, int k_) { return 2.0*nodeOf(k_); }

Opm::Refinement::BlockRefinement box(const std::string& name, std::array<int,3> cellsPerDim,
                                     std::array<int,3> start, std::array<int,3> end,
                                     const std::string& parent = "GLOBAL")
{
    Opm::Refinement::BlockRefinement request;
    request.name = name;
    request.parentGridName = parent;
    request.cellsPerDim = cellsPerDim;
    request.startIJK = start;
    request.endIJK = end;
    return request;
}

void build(Dune::CpGrid& grid, const TestGrdecl& g, std::vector<Opm::Refinement::BlockRefinement> requests)
{
    grid.processEclipseFormat(g.raw(), false);
    grid.setLgrBackend(Opm::Refinement::Backend::Conforming);
    grid.setRefinementBuilder(std::make_shared<Opm::Refinement::ConformingBlockBuilder>(
        g.dims, g.coord, g.zcorn, g.actnum));
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

// Faces between leaf cells of the two levels, and the most faces any cell pair has.
struct Contacts
{
    std::map<std::pair<int,int>, int> byLevels;
    int mostPerPair{0};
};

Contacts contacts(const Dune::CpGrid& grid)
{
    Contacts result;
    std::map<std::pair<int,int>, int> pairs;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (!is.neighbor() || (is.inside().index() > is.outside().index())) {
                continue;
            }
            const int a = is.inside().level();
            const int b = is.outside().level();
            ++result.byLevels[{std::min(a, b), std::max(a, b)}];
            result.mostPerPair = std::max(result.mostPerPair,
                                          ++pairs[{is.inside().index(), is.outside().index()}]);
        }
    }
    return result;
}

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(MPIFixture);

// A nested box next to two top-level boxes sharing a face: their faces stay.
BOOST_AUTO_TEST_CASE(nestedBoxKeepsFacesBetweenOtherBoxes)
{
    const auto g = makeGrid({6, 3, 3}, flat);
    const auto A = box("A", {3,3,3}, {1,1,1}, {3,2,2});
    const auto B = box("B", {3,3,3}, {3,1,1}, {5,2,2});

    Dune::CpGrid plain;
    build(plain, g, {A, B});
    Dune::CpGrid withNested;
    build(withNested, g, {A, B, box("NEST", {2,2,2}, {1,1,1}, {5,2,2}, "A")});

    const auto before = contacts(plain);
    const auto after = contacts(withNested);
    BOOST_CHECK_GT(before.byLevels.at({1, 2}), 0);
    BOOST_CHECK_EQUAL(after.byLevels.at({1, 2}), before.byLevels.at({1, 2}));
    BOOST_CHECK_GT(after.byLevels.at({1, 3}), 0);            // the nested box joins A
    BOOST_CHECK_EQUAL(after.mostPerPair, 1);
    BOOST_CHECK_SMALL(worstClosure(withNested), 1e-12);
    BOOST_CHECK_CLOSE(totalVolume(withNested), totalVolume(plain), 1e-10);
}

// A box on a fault keeps its faulted side when it holds a nested box.
BOOST_AUTO_TEST_CASE(nestedBoxInABoxOnAFault)
{
    const auto g = makeGrid({4, 2, 2}, [](int i_, int, int k_) {
        return 2.0*nodeOf(k_) + ((cellOf(i_) >= 2) ? 0.6 : 0.0);
    });
    const auto A = box("A", {3,3,3}, {2,0,0}, {4,2,2});

    Dune::CpGrid plain;
    build(plain, g, {A});
    Dune::CpGrid withNested;
    build(withNested, g, {A, box("NEST", {2,2,2}, {2,1,1}, {4,5,5}, "A")});

    const auto before = contacts(plain);
    const auto after = contacts(withNested);
    BOOST_CHECK_EQUAL(after.byLevels.at({0, 1}), before.byLevels.at({0, 1}));
    BOOST_CHECK_EQUAL(after.mostPerPair, 1);
    BOOST_CHECK_SMALL(worstClosure(withNested), 1e-12);
    BOOST_CHECK_CLOSE(totalVolume(withNested), totalVolume(plain), 1e-10);
}

// Two nested boxes in one LGR: the LGR knows the children of both.
BOOST_AUTO_TEST_CASE(twoNestedBoxesInOneLgr)
{
    const auto g = makeGrid({4, 3, 3}, flat);
    Dune::CpGrid grid;
    build(grid, g, {box("LGR1", {3,3,3}, {1,1,1}, {3,2,2}),
                    box("NEST1", {2,2,2}, {1,1,1}, {3,2,2}, "LGR1"),
                    box("NEST2", {2,2,2}, {4,1,1}, {5,2,2}, "LGR1")});

    int refined = 0;
    int children = 0;
    for (const auto& element : Dune::elements(grid.levelGridView(1))) {
        if (!element.isLeaf()) {
            ++refined;
            for (const auto& child : Dune::descendantElements(element, grid.maxLevel())) {
                children += (child.level() > 1);
            }
        }
    }
    BOOST_CHECK_EQUAL(refined, 3);
    BOOST_CHECK_EQUAL(children, 3*8);
    BOOST_CHECK_SMALL(worstClosure(grid), 1e-12);
}

// Nested boxes sharing a face inside their LGR are joined, with equal or nested subdivisions.
BOOST_AUTO_TEST_CASE(nestedBoxesSharingAFace)
{
    const auto g = makeGrid({4, 3, 3}, flat);
    for (const std::array<int,3> second : { std::array<int,3>{2,2,2}, std::array<int,3>{4,4,4} }) {
        Dune::CpGrid grid;
        build(grid, g, {box("LGR1", {3,3,3}, {1,1,1}, {3,2,2}),
                        box("NEST1", {2,2,2}, {1,1,1}, {3,2,2}, "LGR1"),
                        box("NEST2", second, {3,1,1}, {5,2,2}, "LGR1")});
        const auto c = contacts(grid);
        BOOST_CHECK_GT(c.byLevels.at({2, 3}), 0);
        BOOST_CHECK_EQUAL(c.mostPerPair, 1);
        BOOST_CHECK_SMALL(worstClosure(grid), 1e-12);
    }
}

// A nested box whose side lies on a fault inside its LGR is refused.
BOOST_AUTO_TEST_CASE(nestedSideOnAFaultIsRefused)
{
    const auto g = makeGrid({4, 1, 2}, [](int i_, int, int k_) {
        return 2.0*nodeOf(k_) + ((cellOf(i_) >= 2) ? 0.6 : 0.0);
    });
    Dune::CpGrid grid;
    const auto onTheFault = [](const std::logic_error& e) {
        return std::string(e.what()).find("fault or pinch-out") != std::string::npos;
    };
    BOOST_CHECK_EXCEPTION(build(grid, g, {box("LGR1", {3,3,3}, {0,0,0}, {4,1,2}),
                                          box("NEST", {2,2,2}, {4,1,1}, {6,2,5}, "LGR1")}),
                          std::logic_error, onTheFault);
}

// The nesting rules read only the requests, so they are checked before any level grid
// is built: a nested box on its parent LGR's boundary, and a box nested in a nested box,
// are refused with the grid left unrefined.
BOOST_AUTO_TEST_CASE(nestingRulesCheckedBeforeTheBuild)
{
    const auto g = makeGrid({4, 3, 3}, flat);
    // LGR1 refines cells i = 1..2, j = 1, k = 1 into 3x3x3: its refined index space is 6x3x3.
    const auto lgr1 = box("LGR1", {3,3,3}, {1,1,1}, {3,2,2});
    const auto onTheBoundary = [](const std::logic_error& e) {
        return std::string(e.what()).find("touches the boundary") != std::string::npos;
    };
    const auto twoLevels = [](const std::logic_error& e) {
        return std::string(e.what()).find("only one level of nesting") != std::string::npos;
    };

    // NEST starts at i = 0, on LGR1's boundary.
    Dune::CpGrid onBoundary;
    BOOST_CHECK_EXCEPTION(build(onBoundary, g, {lgr1,
                                                box("NEST", {2,2,2}, {0,1,1}, {2,2,2}, "LGR1")}),
                          std::logic_error, onTheBoundary);
    BOOST_CHECK_EQUAL(onBoundary.maxLevel(), 0);

    // NEST1 (refined index space 12x3x3) keeps one cell inside LGR1, and NEST2 one cell
    // inside NEST1; NEST2 is refused for its depth alone.
    Dune::CpGrid twoDeep;
    BOOST_CHECK_EXCEPTION(build(twoDeep, g, {lgr1,
                                             box("NEST1", {3,3,3}, {1,1,1}, {5,2,2}, "LGR1"),
                                             box("NEST2", {2,2,2}, {1,1,1}, {2,2,2}, "NEST1")}),
                          std::logic_error, twoLevels);
    BOOST_CHECK_EQUAL(twoDeep.maxLevel(), 0);
}

// A nested box over a fault inside its parent LGR that does not cover every LGR1 layer the
// throw joins: across the throw its cells face LGR1 cells outside it, through faces inside
// it. That is refused, naming the cells. A nested box on one side of the fault is built.
BOOST_AUTO_TEST_CASE(nestedThrowJoiningCellsOutsideIsRefused)
{
    // Columns i >= 2 are 0.6 deeper; LGR1 refines every cell into 3x1x3 (12x3x6).
    const auto g = makeGrid({4, 3, 2}, [](int i_, int, int k_) {
        return 2.0*nodeOf(k_) + ((cellOf(i_) >= 2) ? 0.6 : 0.0);
    });
    const auto lgr1 = box("LGR1", {3,1,3}, {0,0,0}, {4,3,2});
    // NEST covers LGR1 cells i = 4..7 across the fault (between i = 5 and 6), k = 1..4 of 0..5.
    const auto outsideTheBox = [](const std::invalid_argument& e) {
        const std::string what = e.what();
        return (what.find("Refinement box 'NEST'") != std::string::npos)
            && (what.find("outside the box (I,J,K from 1 in 'LGR1')") != std::string::npos)
            && (what.find("must cover every layer") != std::string::npos);
    };
    Dune::CpGrid straddling;
    BOOST_CHECK_EXCEPTION(build(straddling, g, {lgr1, box("NEST", {2,2,2}, {4,1,1}, {8,2,5}, "LGR1")}),
                          std::invalid_argument, outsideTheBox);

    // NEST covers LGR1 cells i = 3..4, one cell from the fault.
    Dune::CpGrid oneSide;
    build(oneSide, g, {lgr1, box("NEST", {2,2,2}, {3,1,1}, {5,2,5}, "LGR1")});
    BOOST_CHECK_EQUAL(oneSide.maxLevel(), 2);
    BOOST_CHECK_EQUAL(contacts(oneSide).mostPerPair, 1);
    BOOST_CHECK_SMALL(worstClosure(oneSide), 1e-12);
}

// The nested box's MINPV removes every child: the nested box is an empty LGR.
BOOST_AUTO_TEST_CASE(nestedBoxWithNoChildLeft)
{
    const auto g = makeGrid({4, 3, 3}, flat);
    // LGR1 refines every cell into 3x3x3 (12x9x9); NEST refines LGR1 cells 3..8, 3..5,
    // 3..5 into 2x2x2 (12x6x6).
    auto nest = box("NEST", {2,2,2}, {3,3,3}, {9,6,6}, "LGR1");
    nest.minpvRemoved.assign(12*6*6, 1);
    Dune::CpGrid grid;
    build(grid, g, {box("LGR1", {3,3,3}, {0,0,0}, {4,3,3}), nest});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 2);
    BOOST_CHECK_EQUAL(grid.levelGridView(2).size(0), 0);
    BOOST_CHECK_SMALL(worstClosure(grid), 1e-12);
    // The 6x3x3 LGR1 cells NEST refines are refined into nothing.
    BOOST_CHECK_EQUAL(grid.size(0), 4*3*3*27 - 6*3*3);
}

// The nested box's MINPV removes every child of one LGR1 cell, inside the nested box or on
// its side: that cell is refined into nothing, and the cells next to it keep their faces
// toward it as boundary faces. With one child inside NEST removed, or none, every LGR1
// cell stays.
BOOST_AUTO_TEST_CASE(nestedParentWithNoChildLeft)
{
    const auto g = makeGrid({4, 3, 3}, flat);
    const auto lgr1 = box("LGR1", {3,3,3}, {0,0,0}, {4,3,3});
    // NEST refines LGR1 cells 3..8, 3..5, 3..5 into 2x2x2: 12x6x6 children of volume
    // (2/27)/8, the children of NEST's parent (pi,pj,pk) at 2*pi.., 2*pj.., 2*pk...
    const auto child = [](int i, int j, int k) { return i + 12*(j + 6*k); };
    const auto removeParent = [&child](std::vector<int>& mask, int pi, int pj, int pk) {
        for (int k = 0; k < 2; ++k) {
            for (int j = 0; j < 2; ++j) {
                for (int i = 0; i < 2; ++i) {
                    mask[child(2*pi + i, 2*pj + j, 2*pk + k)] = 1;
                }
            }
        }
    };
    const std::vector<int> none(12*6*6, 0);
    auto interior = none;
    removeParent(interior, 2, 1, 1);          // LGR1 cell (5,4,4)
    auto onTheSide = none;
    removeParent(onTheSide, 0, 0, 0);         // LGR1 cell (3,3,3), on NEST's corner
    auto oneChild = none;
    oneChild[child(4, 2, 2)] = 1;             // inside LGR1 cell (5,4,4)

    for (const auto& [mask, emptied] : std::vector<std::pair<std::vector<int>, std::array<int,3>>>{
             {none, {-1,-1,-1}}, {oneChild, {-1,-1,-1}}, {interior, {5,4,4}}, {onTheSide, {3,3,3}}}) {
        auto nest = box("NEST", {2,2,2}, {3,3,3}, {9,6,6}, "LGR1");
        nest.minpvRemoved = mask;
        Dune::CpGrid grid;
        build(grid, g, {lgr1, nest});

        const auto removed = std::count(mask.begin(), mask.end(), 1);
        BOOST_CHECK_EQUAL(grid.size(0), 4*3*3*27 - 6*3*3 + 12*6*6 - removed);
        BOOST_CHECK_CLOSE(totalVolume(grid), 72.0 - removed*(2.0/27.0)/8.0, 1e-10);
        BOOST_CHECK_SMALL(worstClosure(grid), 1e-12);
        // An LGR1 cell is a leaf unless NEST refines it, into children or into nothing.
        int lgr1Leaves = 0;
        for (const auto& element : Dune::elements(grid.levelGridView(1))) {
            std::array<int,3> ijk;
            grid.currentData()[1]->getIJK(element.index(), ijk);
            lgr1Leaves += element.isLeaf() ? 1 : 0;
            if (ijk == emptied) {
                BOOST_CHECK(!element.isLeaf());
            }
        }
        BOOST_CHECK_EQUAL(lgr1Leaves, 4*3*3*27 - 6*3*3);
    }
}

// stableCellId() cannot tell nested cells apart, and says so.
BOOST_AUTO_TEST_CASE(stableCellIdRefusesNestedGrids)
{
    const auto g = makeGrid({4, 3, 3}, flat);
    Dune::CpGrid grid;
    build(grid, g, {box("LGR1", {3,3,3}, {1,1,1}, {3,2,2}),
                    box("NEST1", {2,2,2}, {1,1,1}, {5,2,2}, "LGR1")});
    BOOST_CHECK_THROW(grid.stableCellId(), std::logic_error);
}
