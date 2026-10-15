/*
  Copyright 2026 Equinor ASA

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
#include "config.h"

#define BOOST_TEST_MODULE PolyhedralGridFaultsTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/utility/platform_dependent/disable_warnings.h>
#include <dune/common/parallel/mpihelper.hh>
#include <opm/grid/polyhedralgrid.hh>
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>

#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>

#include "CornerPointTestGrid.hpp"

#include <array>
#include <cmath>
#include <vector>

using Grid = Dune::PolyhedralGrid<3, 3>;
using Opm::TestGrids::CornerPointGrid;
using Opm::TestGrids::layersWithShift;

namespace
{

struct MpiFixture
{
    MpiFixture()
    {
        int argc = boost::unit_test::framework::master_test_suite().argc;
        char** argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(argc, argv);
    }
};

// Check that the corners of each cell are at the positions of the corner-point grid, in the
// order of the Dune reference cube: the x coordinate varies fastest, then y and then z.
void checkCorners(const Grid& grid, const CornerPointGrid& cp)
{
    BOOST_REQUIRE_EQUAL(grid.size(0), cp.dims[0] * cp.dims[1] * cp.dims[2]);
    for (auto it = grid.leafGridView().begin<0>(); it != grid.leafGridView().end<0>(); ++it) {
        const int cell = grid.leafIndexSet().index(*it);
        BOOST_REQUIRE_EQUAL(it->geometry().corners(), 8);
        for (int corner = 0; corner < 8; ++corner) {
            const auto position = it->geometry().corner(corner);
            const auto& expected = cp.expectedCorners[8 * cell + corner];
            for (int d = 0; d < 3; ++d) {
                BOOST_CHECK_MESSAGE(std::abs(position[d] - expected[d]) < 1e-10,
                                    "cell " << cell << " corner " << corner << " coordinate " << d
                                    << " is " << position[d] << " but should be " << expected[d]);
            }
        }
    }
}

double totalVolume(const Grid& grid)
{
    double volume = 0.0;
    for (auto it = grid.leafGridView().begin<0>(); it != grid.leafGridView().end<0>(); ++it) {
        volume += it->geometry().volume();
    }
    return volume;
}

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(MpiFixture);

BOOST_AUTO_TEST_CASE(RegularGridHasTheCornersInOrder)
{
    const CornerPointGrid cp(2, 2, 2, layersWithShift([](int, int, int, int, int) { return 0.0; }));
    const Opm::EclipseGrid eclipseGrid(cp.dims, cp.coord, cp.zcorn);
    const Grid grid(eclipseGrid);
    checkCorners(grid, cp);
    BOOST_CHECK_CLOSE(totalVolume(grid), 8.0, 1e-8);
}

// The right column is moved down by half a layer, so the faces between the columns are split
// where the layers of the neighbor start. The dividing lines meet the edges of the cells in the
// middle, and no cell corner is on a dividing line.
//
// Seen along y, with the depth z downwards:
//
//        x=0       x=1       x=2
//  z=0   +---------+
//        |         |
//  z=0.5 |  (0,0)  +---------+
//        |         |         |
//  z=1   +---------+  (1,0)  |
//        |         |         |
//  z=1.5 |  (0,1)  +---------+
//        |         |         |
//  z=2   +---------+  (1,1)  |
//                  |         |
//  z=2.5           +---------+
//
// The face at x=1 is split at z=0.5, 1, 1.5 and 2 into three pieces, one for each pair
// of neighbors: (0,0)-(1,0), (0,1)-(1,0) and (0,1)-(1,1), with (i,k) the cell.
BOOST_AUTO_TEST_CASE(FaultWithDividingLinesInTheMiddleOfTheEdges)
{
    const CornerPointGrid cp(2, 1, 2, layersWithShift([](int i, int, int, int, int) { return i == 1 ? 0.5 : 0.0; }));
    const Opm::EclipseGrid eclipseGrid(cp.dims, cp.coord, cp.zcorn);
    const Grid grid(eclipseGrid);
    checkCorners(grid, cp);
    BOOST_CHECK_CLOSE(totalVolume(grid), 4.0, 1e-8);
}

// The right column is moved down by half a layer at one of the pillars and not at the other,
// so the dividing lines are slanted and some of them end in a corner of a cell.
//
// Seen along y, at the pillars y=0 and y=1:
//
//   y=0, right column moved down     y=1, no shift
//
//        x=0       x=1       x=2          x=0       x=1       x=2
//  z=0   +---------+                      +---------+---------+
//        |         |                      |         |         |
//  z=0.5 |         +---------+            |         |         |
//        |         |         |            |         |         |
//  z=1   +---------+         |            +---------+---------+
//        |         |         |            |         |         |
//  z=1.5 |         +---------+            |         |         |
//        |         |         |            |         |         |
//  z=2   +---------+         |            +---------+---------+
//                  |         |
//  z=2.5           +---------+
//
// The dividing lines of the face at x=1 go from the boundary of the left column at z=1 (and 2)
// at y=1, where the layers meet, to the middle of the cell edge at y=0. They end in a corner there.
BOOST_AUTO_TEST_CASE(FaultWithDividingLinesEndingInCorners)
{
    const CornerPointGrid cp(2, 1, 2, layersWithShift([](int i, int, int, int, int jy) { return i == 1 ? 0.5 * (1 - jy) : 0.0; }));
    const Opm::EclipseGrid eclipseGrid(cp.dims, cp.coord, cp.zcorn);
    const Grid grid(eclipseGrid);
    checkCorners(grid, cp);
    BOOST_CHECK_CLOSE(totalVolume(grid), 4.0, 1e-8);
}

// A fault in both horizontal directions, with tilted layers on the other side.
//
// Seen from above, with the shift of each cell. The column i=1 is moved down by 0.35 at the
// pillar x=1 and 0.70 at x=2, and the row j=2 by 0.45 at y=2 and 0.90 at y=3, so the cell (1,2)
// has both shifts added.
//
//  y=3  +---------+---------+---------+
//       |  y-shift|  x+y    |  y-shift|
//  y=2  +---------+---------+---------+
//       |         |  x-shift|         |
//  y=1  +---------+---------+---------+
//       |         |  x-shift|         |
//  y=0  +---------+---------+---------+
//      x=0       x=1       x=2       x=3
BOOST_AUTO_TEST_CASE(FaultsInBothDirections)
{
    const CornerPointGrid cp(3, 3, 2, layersWithShift([](int i, int j, int, int ix, int jy) {
        return (i == 1 ? 0.35 * (1 + ix) : 0.0) + (j == 2 ? 0.45 * (1 + jy) : 0.0);
    }));
    const Opm::EclipseGrid eclipseGrid(cp.dims, cp.coord, cp.zcorn);
    const Grid grid(eclipseGrid);
    checkCorners(grid, cp);
}

// The middle layer of a single column has no thickness at one of the pillars, so two
// corners of that cell are the same node and the cell has six nodes only. The layer
// boundaries are at the depths 0, 1, 2 and 3, except that at the collapsed pillar the
// boundary 2 is at the depth 1, where the boundary 1 is.
//
// Seen along y at y=0, where the pillar x=1 is the collapsed one:
//
//        x=0                 x=1
//  z=0   +-------------------+
//        |        k=0        |
//  z=1   +-------------------+  <- the boundaries 1 and 2 meet here
//        | k=1          ..--''
//        |         ..--''
//        |   ..--''     k=2
//  z=2   +-''                |
//        |                   |
//  z=3   +-------------------+
//
// The cell k=1 is a triangle in this view, and its two corners at x=1 are the same node.
BOOST_AUTO_TEST_CASE(CollapsedEdge)
{
    const CornerPointGrid cp(1, 1, 3, [](int, int, int k, int zs, int ix, int jy) {
        return (ix == 1 && jy == 0 && k + zs == 2) ? 1.0 : double(k + zs);
    });
    const Opm::EclipseGrid eclipseGrid(cp.dims, cp.coord, cp.zcorn);
    const Grid grid(eclipseGrid);
    checkCorners(grid, cp);
}

// The whole side of a cell at x = 1 has collapsed to a line, because the top and the bottom
// are at the same depth at both pillars there. The cell is a wedge. That side has no face, and
// the two nodes of the line are the corners of both the top and the bottom.
//
// Seen along y:
//
//        x=0         x=1
//  z=0   +--.
//        |    ``--.
//  z=0.5 |           ``+  <- top and bottom at x=1 are both at z=0.5
//        |       ..--''
//  z=1   +---''
BOOST_AUTO_TEST_CASE(CollapsedSide)
{
    const CornerPointGrid cp(1, 1, 1, [](int, int, int, int zs, int ix, int) {
        return ix == 1 ? 0.5 : double(zs);
    });
    const Opm::EclipseGrid eclipseGrid(cp.dims, cp.coord, cp.zcorn);
    const Grid grid(eclipseGrid);
    checkCorners(grid, cp);
    // A triangular prism: the area of the triangle is 1/2 * 1 * 1 in the xz plane, and its length 1.
    BOOST_CHECK_CLOSE(totalVolume(grid), 0.5, 1e-8);
}

// The same wedge in a grid, with neighbors in the three directions. Only the corners of
// the cell in the middle are moved, so there are faults between it and the cells around it.
//
// Seen along y, at the row j=1, with the cells (i,k):
//
//        x=0     x=1     x=2     x=3
//  z=0   +-------+-------+-------+
//        |       |       |       |
//  z=1   +-------+-------+-------+
//        |       |``--.  |       |
//  z=1.5 |       | (1,1) `+      |
//        |       |  ..--'|       |
//  z=2   +-------+-------+-------+
//        |       |       |       |
//  z=3   +-------+-------+-------+
BOOST_AUTO_TEST_CASE(CollapsedSideWithNeighbors)
{
    const CornerPointGrid cp(3, 3, 3, [](int i, int j, int k, int zs, int ix, int) {
        // The cell in the middle is collapsed at its side at x = 2, and the cells next to it follow.
        const bool collapsed = (i == 1 && j == 1 && k == 1 && ix == 1);
        return collapsed ? k + 0.5 : double(k + zs);
    });
    const Opm::EclipseGrid eclipseGrid(cp.dims, cp.coord, cp.zcorn);
    const Grid grid(eclipseGrid);
    BOOST_REQUIRE_EQUAL(grid.size(0), 27);
    // Only the corners of the collapsed cell are checked, since the cells around it have faults.
    const int collapsedCell = 1 + 3 * (1 + 3 * 1);
    for (auto it = grid.leafGridView().begin<0>(); it != grid.leafGridView().end<0>(); ++it) {
        if (grid.leafIndexSet().index(*it) != collapsedCell) {
            continue;
        }
        for (int corner = 0; corner < 8; ++corner) {
            const auto position = it->geometry().corner(corner);
            const auto& expected = cp.expectedCorners[8 * collapsedCell + corner];
            for (int d = 0; d < 3; ++d) {
                BOOST_CHECK_CLOSE(position[d] + 1.0, expected[d] + 1.0, 1e-8);
            }
        }
    }
}

// A sliver: the cell has thickness at one pillar only, and its top and bottom are at the same depth
// at the three other pillars. The sides at x = 1 and y = 1 have collapsed, so they have no face, and
// the node at the pillar where they meet is on the top and the bottom only.
//
// Seen from above, with the depths of the top and the bottom at each pillar:
//
//  y=1   o 0.2 ------------ o 0.1
//        |                  |
//        |                  |
//  y=0   o 0 and 0.01 ----- o 0.3
//       x=0                x=1
//
// Only at the pillar x=0, y=0 the cell has thickness. The sides at x=1 and y=1 are triangles
// with no area, so they have no face, and the node at x=1, y=1 is on the top and the bottom only.
BOOST_AUTO_TEST_CASE(CollapsedAtThreePillars)
{
    const CornerPointGrid cp(1, 1, 1, [](int, int, int, int zs, int ix, int jy) {
        if (ix == 0 && jy == 0) {
            return 0.01 * zs;
        }
        return ix == 1 && jy == 1 ? 0.1 : (ix == 1 ? 0.3 : 0.2);
    });
    const Opm::EclipseGrid eclipseGrid(cp.dims, cp.coord, cp.zcorn);
    const Grid grid(eclipseGrid);
    checkCorners(grid, cp);
}
