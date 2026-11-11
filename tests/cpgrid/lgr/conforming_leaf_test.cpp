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
