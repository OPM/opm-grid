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

#include <array>
#include <cmath>
#include <vector>

using Grid = Dune::PolyhedralGrid<3, 3>;

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

// A corner-point grid with vertical pillars at the integer positions x = 0..nx, y = 0..ny.
// The depth of the corner (ix, jy, zs) of the cell (i, j, k), where ix and jy tell which of
// the pillars of the cell and zs whether it is the top or the bottom, is given by depth().
struct CornerPointGrid
{
    std::array<int, 3> dims;
    std::vector<double> coord;
    std::vector<double> zcorn;
    // For each cell its eight corners, in the order of the Dune reference cube: ix + 2 * jy + 4 * zs.
    std::vector<std::array<double, 3>> expectedCorners;

    template <class Depth>
    CornerPointGrid(const int nx, const int ny, const int nz, Depth depth)
        : dims{nx, ny, nz}
    {
        for (int py = 0; py <= ny; ++py) {
            for (int px = 0; px <= nx; ++px) {
                for (const double z : {-10.0, 10.0}) {
                    coord.insert(coord.end(), {double(px), double(py), z});
                }
            }
        }

        // The order of ZCORN: layers from top to bottom, within each the rows of corners in y,
        // and within each row the corners in x. Each cell has two corners in each direction.
        zcorn.resize(8 * nx * ny * nz);
        expectedCorners.resize(8 * nx * ny * nz);
        for (int k = 0; k < nz; ++k) {
            for (int zs = 0; zs < 2; ++zs) {
                for (int j = 0; j < ny; ++j) {
                    for (int jy = 0; jy < 2; ++jy) {
                        for (int i = 0; i < nx; ++i) {
                            for (int ix = 0; ix < 2; ++ix) {
                                const double z = depth(i, j, k, zs, ix, jy);
                                const std::size_t index = ((2 * k + zs) * 2 * ny + 2 * j + jy) * 2 * nx + 2 * i + ix;
                                zcorn[index] = z;
                                const std::size_t cell = i + nx * (j + ny * k);
                                expectedCorners[8 * cell + ix + 2 * jy + 4 * zs] = {double(i + ix), double(j + jy), z};
                            }
                        }
                    }
                }
            }
        }
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

// Layers of unit thickness, with the corners moved down by shift(i, j, k, ix, jy).
template <class Shift>
auto layersWithShift(Shift shift)
{
    return [shift](int i, int j, int k, int zs, int ix, int jy) { return k + zs + shift(i, j, k, ix, jy); };
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
BOOST_AUTO_TEST_CASE(FaultWithDividingLinesEndingInCorners)
{
    const CornerPointGrid cp(2, 1, 2, layersWithShift([](int i, int, int, int, int jy) { return i == 1 ? 0.5 * (1 - jy) : 0.0; }));
    const Opm::EclipseGrid eclipseGrid(cp.dims, cp.coord, cp.zcorn);
    const Grid grid(eclipseGrid);
    checkCorners(grid, cp);
    BOOST_CHECK_CLOSE(totalVolume(grid), 4.0, 1e-8);
}

// A fault in both horizontal directions, with tilted layers on the other side.
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

// The same wedge in a grid, with neighbors in the three directions.
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
