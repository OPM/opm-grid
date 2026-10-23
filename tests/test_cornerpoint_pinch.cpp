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

#define BOOST_TEST_MODULE CornerPointPinchTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/UnstructuredGrid.h>
#include <opm/grid/cornerpoint_grid.h>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <array>
#include <memory>
#include <vector>

namespace
{

using GridPtr = std::unique_ptr<UnstructuredGrid, void (*)(UnstructuredGrid*)>;

// One column of three layers with unit area. The layers are 1 m thick, except for the middle one, which has the thickness given and
// is inactive. The depth of the top of the first layer is zero.
GridPtr column(const double middleThickness, const bool pinchActive)
{
    std::vector<double> coord;
    for (int py = 0; py <= 1; ++py) {
        for (int px = 0; px <= 1; ++px) {
            coord.insert(coord.end(), {double(px), double(py), -10.0, double(px), double(py), 10.0});
        }
    }
    const std::array<std::array<double, 2>, 3> layers{{{0.0, 1.0}, {1.0, 1.0 + middleThickness},
                                                         {1.0 + middleThickness, 2.0 + middleThickness}}};
    std::vector<double> zcorn;
    for (const auto& layer : layers) {
        for (const double z : layer) {
            zcorn.insert(zcorn.end(), 4, z);
        }
    }
    std::vector<int> actnum{1, 0, 1};
    grdecl g{};
    g.dims[0] = g.dims[1] = 1;
    g.dims[2] = 3;
    g.coord = coord.data();
    g.zcorn = zcorn.data();
    g.actnum = actnum.data();
    return GridPtr(create_grid_cornerpoint_pinch(&g, 0.0, 0, pinchActive ? 1 : 0), destroy_grid);
}

// The faces with a cell on both sides.
int connections(const UnstructuredGrid& g)
{
    int count = 0;
    for (int f = 0; f < g.number_of_faces; ++f) {
        count += g.face_cells[2 * f] >= 0 && g.face_cells[2 * f + 1] >= 0;
    }
    return count;
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(CollapsedInactiveLayerIsNotAConnectionWithoutPinch)
{
    const auto grid = column(0.0, false);
    BOOST_REQUIRE(grid);
    BOOST_CHECK_EQUAL(grid->number_of_cells, 2);
    BOOST_CHECK_EQUAL(connections(*grid), 0);
}

BOOST_AUTO_TEST_CASE(CollapsedInactiveLayerIsAConnectionWithPinch)
{
    const auto grid = column(0.0, true);
    BOOST_REQUIRE(grid);
    BOOST_CHECK_EQUAL(grid->number_of_cells, 2);
    BOOST_CHECK_EQUAL(connections(*grid), 1);
    for (int f = 0; f < grid->number_of_faces; ++f) {
        if (grid->face_cells[2 * f] >= 0 && grid->face_cells[2 * f + 1] >= 0) {
            // the face between the two cells is the shared face, with the area of the column
            BOOST_CHECK_CLOSE(grid->face_areas[f], 1.0, 1e-10);
            BOOST_CHECK_EQUAL(grid->face_cells[2 * f], 0);
            BOOST_CHECK_EQUAL(grid->face_cells[2 * f + 1], 1);
        }
    }
    // the volume of the cells is not changed
    BOOST_CHECK_CLOSE(grid->cell_volumes[0], 1.0, 1e-10);
    BOOST_CHECK_CLOSE(grid->cell_volumes[1], 1.0, 1e-10);
}

// An inactive layer with thickness is not a pinch-out, whatever PINCH says.
BOOST_AUTO_TEST_CASE(InactiveLayerWithThicknessIsNotAConnection)
{
    for (const bool pinch : {false, true}) {
        const auto grid = column(0.5, pinch);
        BOOST_REQUIRE(grid);
        BOOST_CHECK_EQUAL(grid->number_of_cells, 2);
        BOOST_CHECK_EQUAL(connections(*grid), 0);
    }
}

// The function that does not take the pinch argument is the one without pinch.
BOOST_AUTO_TEST_CASE(CreateWithoutPinchArgumentHasNoPinchConnections)
{
    std::vector<double> coord;
    for (int py = 0; py <= 1; ++py) {
        for (int px = 0; px <= 1; ++px) {
            coord.insert(coord.end(), {double(px), double(py), -10.0, double(px), double(py), 10.0});
        }
    }
    std::vector<double> zcorn;
    for (const double z : {0.0, 1.0, 1.0, 1.0, 1.0, 2.0}) {
        zcorn.insert(zcorn.end(), 4, z);
    }
    std::vector<int> actnum{1, 0, 1};
    grdecl g{};
    g.dims[0] = g.dims[1] = 1;
    g.dims[2] = 3;
    g.coord = coord.data();
    g.zcorn = zcorn.data();
    g.actnum = actnum.data();
    GridPtr grid(create_grid_cornerpoint(&g, 0.0, 0), destroy_grid);
    BOOST_REQUIRE(grid);
    BOOST_CHECK_EQUAL(connections(*grid), 0);
}
