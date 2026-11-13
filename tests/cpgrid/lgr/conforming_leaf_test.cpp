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

#define BOOST_TEST_MODULE ConformingLeafTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/refinement/RefinementRequest.hpp>

#include <tests/cpgrid/lgr/ConformingTestHelpers.hpp>

#include <array>
#include <cmath>
#include <vector>

// The leaf of the conforming refinement: cells whose children a box's MINPV removes,
// the corner order of a box on a left-handed grid, the level-zero face an LGR
// boundary face comes from, and the leaf index each grid gives its cells.

using namespace Opm::ConformingTest;

BOOST_GLOBAL_FIXTURE(MPIFixture);

// The box's MINPV removes both children of one host: the host is refined into nothing,
// and its coarse neighbour keeps its face toward it as a boundary face.
BOOST_AUTO_TEST_CASE(hostWithNoChildLeft)
{
    const auto g = makeGrid({4, 1, 1}, flat);
    auto request = box("LGR1", {2,1,1}, {1,0,0}, {3,1,1});
    request.minpvRemoved = {1, 1, 0, 0};                 // children of host i = 1
    Dune::CpGrid grid;
    grid.processEclipseFormat(g.raw(), false);
    refine(grid, g, {request});

    BOOST_CHECK_EQUAL(grid.size(0), 1 + 2 + 1);
    BOOST_CHECK_SMALL(worstClosure(grid), 1e-12);
    for (const auto& element : Dune::elements(grid.levelGridView(0))) {
        const int i = grid.currentData().front()->globalCell()[element.index()];
        BOOST_CHECK_EQUAL(element.isLeaf(), (i == 0) || (i == 3));
    }
}

// The box's MINPV removes every child: the box is an empty LGR, its hosts are refined
// into nothing, and their neighbours keep their faces toward them as boundary faces.
BOOST_AUTO_TEST_CASE(boxWithNoChildLeft)
{
    const auto g = makeGrid({4, 3, 3}, flat);
    auto request = box("LGR1", {2,2,2}, {1,1,1}, {3,2,2});  // hosts (1,1,1) and (2,1,1)
    request.minpvRemoved.assign(4*2*2, 1);
    Dune::CpGrid grid;
    grid.processEclipseFormat(g.raw(), false);
    refine(grid, g, {request});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_EQUAL(grid.levelGridView(1).size(0), 0);
    BOOST_CHECK_EQUAL(grid.size(0), 36 - 2);
    BOOST_CHECK_SMALL(worstClosure(grid), 1e-12);
    for (const auto& element : Dune::elements(grid.levelGridView(0))) {
        const int cell = grid.currentData().front()->globalCell()[element.index()];
        BOOST_CHECK_EQUAL(element.isLeaf(), (cell != 1 + 4*1 + 12*1) && (cell != 2 + 4*1 + 12*1));
    }
}

// A box one parent wide in J has a level grid one cell thick in J. Its
// corners must still match its parents' corners for both J directions and both
// ZCORN directions, including when the two winding reversals cancel.
BOOST_AUTO_TEST_CASE(thinBoxCornerOrder)
{
    for (const bool reverseJ : {false, true}) {
        for (const bool upwardZcorn : {false, true}) {
            BOOST_TEST_CONTEXT("reverseJ=" << reverseJ << ", upwardZcorn=" << upwardZcorn) {
                const auto g = makeGrid({3, 2, 2}, reverseJ ? leftHanded : vertical,
                    [upwardZcorn](int, int, int k_) {
                        const double depth = 2.0*nodeOf(k_);
                        return upwardZcorn ? 4.0 - depth : depth;
                    });
                Dune::CpGrid grid;
                grid.processEclipseFormat(g.raw(), false);
                refine(grid, g, {box("LGR1", {2,1,2}, {1,0,0}, {2,1,2})});

                for (const auto& element : Dune::elements(grid.leafGridView())) {
                    const auto& geometry = element.geometry();
                    Dune::FieldVector<double,3> mean(0.0);
                    for (int corner = 0; corner < 8; ++corner) {
                        mean += geometry.corner(corner);
                    }
                    mean /= 8.0;
                    BOOST_CHECK_SMALL((mean - geometry.center()).two_norm(), 1e-12);
                }
                BOOST_CHECK_SMALL(worstClosure(grid), 1e-12);
            }
        }
    }
}

// Each LGR boundary face comes from the level-zero face between the two cells' origins
// that holds it.
BOOST_AUTO_TEST_CASE(parentFaceHoldsTheLgrBoundaryFace)
{
    const auto g = makeGrid({4, 2, 2}, vertical, [](int i_, int, int k_) {
        return 2.0*nodeOf(k_) + ((cellOf(i_) >= 2) ? 0.6 : 0.0);
    });
    Dune::CpGrid grid;
    grid.processEclipseFormat(g.raw(), false);
    refine(grid, g, {box("LGR1", {2,2,2}, {2,0,0}, {4,2,2})});

    int checked = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (!is.neighbor() || (element.level() == is.outside().level())) {
                continue;
            }
            const auto parent = grid.getParentIntersectionFromLgrBoundaryFace(is);
            BOOST_CHECK_EQUAL(parent.inside().index(), element.getOrigin().index());
            BOOST_CHECK_EQUAL(parent.outside().index(), is.outside().getOrigin().index());
            // The leaf face lies on the parent face: its centre is within the parent's area
            // around the parent's centre along the fault plane.
            const auto offset = is.geometry().center() - parent.geometry().center();
            BOOST_CHECK_SMALL(std::abs(offset * parent.centerUnitOuterNormal()), 1e-12);
            ++checked;
        }
    }
    BOOST_CHECK_GT(checked, 0);
}

// The box's MINPV removes children on the box's sides, of a host that keeps others: the
// coarse cells across keep a boundary face where each removed child would be. The two
// children removed together share a corner no other cell has.
BOOST_AUTO_TEST_CASE(removedChildOnTheBoxSide)
{
    const auto g = makeGrid({4, 3, 3}, flat);
    // Hosts (1,1,1) and (2,1,1) into 4x2x2 children of volume 0.25. Child 0 is the
    // corner child on the -I, -J and -K sides; child 8 is above it.
    for (const auto& removed : std::vector<std::vector<int>>{ {0}, {0, 8} }) {
        auto request = box("LGR1", {2,2,2}, {1,1,1}, {3,2,2});
        request.minpvRemoved.assign(16, 0);
        for (const int child : removed) {
            request.minpvRemoved[child] = 1;
        }
        Dune::CpGrid grid;
        grid.processEclipseFormat(g.raw(), false);
        refine(grid, g, {request});

        const int count = static_cast<int>(removed.size());
        BOOST_CHECK_EQUAL(grid.size(0), 36 - 2 + 16 - count);
        BOOST_CHECK_SMALL(worstClosure(grid), 1e-12);
        BOOST_CHECK_CLOSE(totalVolume(grid), 72.0 - 0.25*count, 1e-10);
    }
}

// The box's MINPV removes a child on a side across a throw of two whole layers: the
// side is not faulted, as each host's face is the whole face of the cell across, but
// that cell is two layers above the host. It keeps a boundary face where the child would
// be. The same for a nested box beside such a throw inside its parent LGR.
//
//   Section y = const (depth down), 2x1x4 grid, column i=1 shifted up by 2:
//
//   depth   i=0      i=1
//    -2              +-----+
//                    | k=0 |
//    -1              +-----+
//                    | k=1 |
//     0    +-----+   +-----+
//          | k=0 |===| box |   LGR1 refines (1,0,2) into 1x1x2; its child 0, the
//     1    +-----+   +-----+   upper one, is removed: (0,0,0) keeps the upper half
//          | k=1 |   | k=3 |   of the face '===' as a boundary face.
//     2    +-----+   +-----+
//          | k=2 |
//     3    +-----+
//          | k=3 |
//     4    +-----+
BOOST_AUTO_TEST_CASE(removedChildAcrossAWholeLayerThrow)
{
    {
        const auto g = makeGrid({2, 1, 4}, vertical, [](int i_, int, int k_) {
            return 1.0*nodeOf(k_) - ((cellOf(i_) >= 1) ? 2.0 : 0.0);
        });
        auto request = box("LGR1", {1,1,2}, {1,0,2}, {2,1,3});
        request.minpvRemoved = {1, 0};
        Dune::CpGrid grid;
        grid.processEclipseFormat(g.raw(), false);
        refine(grid, g, {request});

        BOOST_CHECK_EQUAL(grid.size(0), 8 - 1 + 2 - 1);
        BOOST_CHECK_CLOSE(totalVolume(grid), 8.0 - 0.5, 1e-10);
        BOOST_CHECK_SMALL(worstClosure(grid), 1e-12);
        BOOST_CHECK_EQUAL(maxFacesPerPair(grid), 1);
    }
    {
        // LGR1 refines the whole 3x1x6 grid into 2x3x1 cells of volume 1/6; NEST refines
        // its cell (4,1,3), beside the throw between parent columns i=1 and i=2, into
        // 1x1x2 children of volume 1/12, and removes the upper one.
        const auto g = makeGrid({3, 1, 6}, vertical, [](int i_, int, int k_) {
            return 1.0*nodeOf(k_) - ((cellOf(i_) >= 2) ? 2.0 : 0.0);
        });
        auto nest = box("NEST", {1,1,2}, {4,1,3}, {5,2,4}, "LGR1");
        nest.minpvRemoved = {1, 0};
        Dune::CpGrid grid;
        grid.processEclipseFormat(g.raw(), false);
        refine(grid, g, {box("LGR1", {2,3,1}, {0,0,0}, {3,1,6}), nest});

        BOOST_CHECK_EQUAL(grid.size(0), 18*6 - 1 + 2 - 1);
        BOOST_CHECK_CLOSE(totalVolume(grid), 18.0 - 1.0/12.0, 1e-10);
        BOOST_CHECK_SMALL(worstClosure(grid), 1e-12);
        BOOST_CHECK_EQUAL(maxFacesPerPair(grid), 1);
    }
}

// Two boxes share the face x = 2 of a 2x1x1 grid of cells 2 x 4 x 6: LGR1 refines cell
// (0,0,0), LGR2 cell (1,0,0). One box's MINPV removes children on the shared face, or
// every child of its host: the cells of the other box keep a boundary face where each
// removed child would be. With different factors the finer box makes the faces between
// the boxes, so a removal in the finer box leaves the coarser box's child a boundary
// face on part of its side; with equal ones, a face's partner is missing.
//
// The face x = 2 (y right, depth down) as each box divides it, with child 0 marked:
//
//            LGR1 1x2x3       LGR2 1x4x6
//   z=0   +-----+-----+     +--+--+--+--+
//         |  0  |     |     |0 |  |  |  |
//   z=1   |     |     |     +--+--+--+--+
//         |     |     |     |  |  |  |  |
//   z=2   +-----+-----+     +--+--+--+--+
//         |     |     |     |  |  |  |  |
//   z=3   |     |     |     +--+--+--+--+
//         |     |     |     |  |  |  |  |
//   z=4   +-----+-----+     +--+--+--+--+
//         |     |     |     |  |  |  |  |
//   z=5   |     |     |     +--+--+--+--+
//         |     |     |     |  |  |  |  |
//   z=6   +-----+-----+     +--+--+--+--+
//        y=0         y=4   y=0         y=4
//
// LGR1's child 0 (volume 8) faces four children of LGR2 (volume 2 each), over area 4.
BOOST_AUTO_TEST_CASE(removedChildOnASharedFace)
{
    const auto g = makeGrid({2, 1, 1},
        [](int i, int j) { return std::array<double,6>{2.0*i, 4.0*j, 0.0, 2.0*i, 4.0*j, 100.0}; },
        [](int, int, int k_) { return 6.0*nodeOf(k_); });
    // The area of the leaf faces on x = 2 with no cell across.
    const auto boundaryAreaOnTheSharedFace = [](const Dune::CpGrid& grid) {
        double area = 0.0;
        for (const auto& element : Dune::elements(grid.leafGridView())) {
            for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
                if (!is.neighbor() && (std::abs(is.geometry().center()[0] - 2.0) < 1e-12)) {
                    area += is.geometry().volume();
                }
            }
        }
        return area;
    };
    const auto mask = [](int size, const std::vector<int>& removed) {
        std::vector<int> m(size, 0);
        for (const int child : removed) {
            m[child] = 1;
        }
        return m;
    };
    struct Case
    {
        std::array<int,3> factors1;
        std::array<int,3> factors2;
        std::vector<int> removed1;      // LGR1's MINPV mask, empty for none
        std::vector<int> removed2;      // LGR2's MINPV mask, empty for none
        int cells;                      // leaf cells
        double volume;
        double boundaryArea;            // on x = 2
    };
    for (const auto& c : std::vector<Case>{
             // The coarser box removes its child 0: LGR2's 4 sub-faces there face the domain.
             {{1,2,3}, {1,4,6}, mask(6, {0}), {}, 5 + 24, 96.0 - 8.0, 4.0},
             // The finer box removes its child 0, then every child.
             {{1,2,3}, {1,4,6}, {}, mask(24, {0}), 6 + 23, 96.0 - 2.0, 1.0},
             {{1,2,3}, {1,4,6}, {}, std::vector<int>(24, 1), 6, 48.0, 24.0},
             // The finer box first.
             {{1,4,6}, {1,2,3}, mask(24, {0}), {}, 23 + 6, 96.0 - 2.0, 1.0},
             {{1,4,6}, {1,2,3}, {}, mask(6, {0}), 24 + 5, 96.0 - 8.0, 4.0},
             // Equal factors: the child across is the face's partner.
             {{1,2,3}, {1,2,3}, mask(6, {0}), {}, 5 + 6, 96.0 - 8.0, 4.0},
             {{1,2,3}, {1,2,3}, {}, mask(6, {0}), 6 + 5, 96.0 - 8.0, 4.0}}) {
        auto lgr1 = box("LGR1", c.factors1, {0,0,0}, {1,1,1});
        auto lgr2 = box("LGR2", c.factors2, {1,0,0}, {2,1,1});
        lgr1.minpvRemoved = c.removed1;
        lgr2.minpvRemoved = c.removed2;
        Dune::CpGrid grid;
        grid.processEclipseFormat(g.raw(), false);
        refine(grid, g, {lgr1, lgr2});

        BOOST_CHECK_EQUAL(grid.size(0), c.cells);
        BOOST_CHECK_CLOSE(totalVolume(grid), c.volume, 1e-10);
        BOOST_CHECK_SMALL(worstClosure(grid), 1e-12);
        BOOST_CHECK_EQUAL(maxFacesPerPair(grid), 1);
        BOOST_CHECK_CLOSE(boundaryAreaOnTheSharedFace(grid), c.boundaryArea, 1e-10);
    }
}

// Each grid maps its cells to the leaf: a leaf cell is found again through the grid it
// comes from, and the hosts of the box, refined, have no leaf cell.
BOOST_AUTO_TEST_CASE(levelCellsMapToTheLeaf)
{
    const auto g = makeGrid({4, 1, 1}, flat);
    Dune::CpGrid grid;
    grid.processEclipseFormat(g.raw(), false);
    refine(grid, g, {box("LGR1", {2,1,1}, {1,0,0}, {3,1,1})});

    const auto& data = grid.currentData();
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        BOOST_CHECK_EQUAL(data[element.level()]->getLeafIdxFromLevelIdx(element.getLevelElem().index()),
                          element.index());
    }
    for (const int host : {1, 2}) {
        BOOST_CHECK_EQUAL(data.front()->getLeafIdxFromLevelIdx(host), -1);
    }
}
