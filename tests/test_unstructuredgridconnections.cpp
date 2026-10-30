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

#define BOOST_TEST_MODULE UnstructuredGridConnectionsTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/UnstructuredGrid.h>
#include <opm/grid/common/UnstructuredGridConnections.hpp>
#include <opm/grid/cornerpoint_grid.h>
#include <opm/grid/cpgpreprocess/preprocess.h>
#include <opm/grid/polyhedralgrid.hh>

#include <dune/common/parallel/mpihelper.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{

using GridPtr = std::unique_ptr<UnstructuredGrid, void (*)(UnstructuredGrid*)>;

// One column with unit area, with the layers given as the depths of the four top corners and of the four bottom corners of each. A layer
// where the top and the bottom are the same is collapsed and is not a cell.
std::array<double, 8> flat(const double top, const double bottom)
{
    return {top, top, top, top, bottom, bottom, bottom, bottom};
}

GridPtr column(const std::vector<std::array<double, 8>>& layers)
{
    std::vector<double> coord;
    for (int py = 0; py <= 1; ++py) {
        for (int px = 0; px <= 1; ++px) {
            coord.insert(coord.end(), {static_cast<double>(px), static_cast<double>(py), -10.0,
                                       static_cast<double>(px), static_cast<double>(py), 100.0});
        }
    }
    std::vector<double> zcorn;
    for (const auto& layer : layers) {
        zcorn.insert(zcorn.end(), layer.begin(), layer.end());
    }
    std::vector<int> actnum(layers.size(), 1);
    grdecl g{};
    g.dims[0] = g.dims[1] = 1;
    g.dims[2] = static_cast<int>(layers.size());
    g.coord = coord.data();
    g.zcorn = zcorn.data();
    g.actnum = actnum.data();
    return {create_grid_cornerpoint(&g, 0.0, 0), destroy_grid};
}

int connections(const UnstructuredGrid& g)
{
    int count = 0;
    for (int f = 0; f < g.number_of_faces; ++f) {
        count += g.face_cells[2 * f] >= 0 && g.face_cells[2 * f + 1] >= 0;
    }
    return count;
}

int numberOfFaces(const UnstructuredGrid& g, const int cell)
{
    return g.cell_facepos[cell + 1] - g.cell_facepos[cell];
}

// The sum of the normals of the faces of a cell, pointing out of it, is zero if the faces fit.
double normalSum(const UnstructuredGrid& g, const int cell)
{
    std::array<double, 3> sum{};
    for (auto pos = g.cell_facepos[cell]; pos < g.cell_facepos[cell + 1]; ++pos) {
        const int f = g.cell_faces[pos];
        const double sign = g.face_cells[2 * f] == cell ? 1.0 : -1.0;
        for (int d = 0; d < 3; ++d) {
            sum[d] += sign * g.face_normals[3 * f + d];
        }
    }
    return std::abs(sum[2]);
}

struct MpiFixture
{
    MpiFixture()
    {
        auto& suite = boost::unit_test::framework::master_test_suite();
        Dune::MPIHelper::instance(suite.argc, suite.argv);
    }
};

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(MpiFixture);

// The space between the cells has the top face of the first cell as a face of the cell, and the second cell has it as one more.
BOOST_AUTO_TEST_CASE(CellsAcrossSpaceAreConnected)
{
    // A (0, 10), a collapsed layer (10, 10), and C (15, 25): 5 m of space that is not a cell
    const auto grid = column({flat(0.0, 10.0), flat(10.0, 10.0), flat(15.0, 25.0)});
    BOOST_REQUIRE(grid);
    BOOST_REQUIRE_EQUAL(grid->number_of_cells, 2);
    BOOST_CHECK_EQUAL(connections(*grid), 0);
    BOOST_CHECK_EQUAL(numberOfFaces(*grid, 0), 6);
    BOOST_CHECK_EQUAL(numberOfFaces(*grid, 1), 6);
    const int boundaryFaces = grid->number_of_faces;

    BOOST_CHECK_EQUAL(Dune::unstructured::connectCellsAcrossGaps(*grid, {{0, 1}}), 1);
    BOOST_CHECK_EQUAL(connections(*grid), 1);
    BOOST_CHECK_EQUAL(grid->number_of_faces, boundaryFaces);
    BOOST_CHECK_EQUAL(numberOfFaces(*grid, 0), 6);
    BOOST_CHECK_EQUAL(numberOfFaces(*grid, 1), 7);

    // the face of the first cell, from it to the other, with the tag of its top side for the second cell
    int found = 0;
    for (int f = 0; f < grid->number_of_faces; ++f) {
        if (grid->face_cells[2 * f] == 0 && grid->face_cells[2 * f + 1] == 1) {
            ++found;
            BOOST_CHECK_CLOSE(grid->face_areas[f], 1.0, 1e-10);
            BOOST_CHECK_CLOSE(grid->face_centroids[3 * f + 2], 10.0, 1e-10);
        }
    }
    BOOST_CHECK_EQUAL(found, 1);
    const int last = grid->cell_facepos[2] - 1;
    BOOST_CHECK_EQUAL(grid->cell_facetag[last], 4);
    BOOST_CHECK_EQUAL(grid->face_cells[2 * grid->cell_faces[last]], 0);
    // the volumes and the other faces of the cells are as they were
    BOOST_CHECK_CLOSE(grid->cell_volumes[0], 10.0, 1e-10);
    BOOST_CHECK_CLOSE(grid->cell_volumes[1], 10.0, 1e-10);
    BOOST_CHECK_SMALL(normalSum(*grid, 0), 1e-10);
}

BOOST_AUTO_TEST_CASE(ConnectedCellsAndUnknownCellsAreLeftAlone)
{
    const auto grid = column({flat(0.0, 10.0), flat(10.0, 10.0), flat(15.0, 25.0)});
    BOOST_REQUIRE(grid);
    BOOST_CHECK_EQUAL(Dune::unstructured::connectCellsAcrossGaps(*grid, {{0, 1}}), 1);
    // the first cell has no face at the boundary at the bottom any more
    BOOST_CHECK_EQUAL(Dune::unstructured::connectCellsAcrossGaps(*grid, {{0, 1}}), 0);
    BOOST_CHECK_EQUAL(Dune::unstructured::connectCellsAcrossGaps(*grid, {{0, 5}, {-1, 1}, {1, 1}}), 0);
    BOOST_CHECK_EQUAL(connections(*grid), 1);
    BOOST_CHECK_EQUAL(numberOfFaces(*grid, 1), 7);
}

// Cells that are neighbors are connected already.
BOOST_AUTO_TEST_CASE(NeighborsAreNotConnectedAgain)
{
    const auto grid = column({flat(0.0, 10.0), flat(10.0, 20.0)});
    BOOST_REQUIRE(grid);
    BOOST_CHECK_EQUAL(connections(*grid), 1);
    BOOST_CHECK_EQUAL(Dune::unstructured::connectCellsAcrossGaps(*grid, {{0, 1}}), 0);
    BOOST_CHECK_EQUAL(connections(*grid), 1);
}

BOOST_AUTO_TEST_CASE(InvalidGrid)
{
    UnstructuredGrid empty{};
    empty.dimensions = 2;
    BOOST_CHECK_THROW(Dune::unstructured::connectCellsAcrossGaps(empty, {}), std::invalid_argument);
}

// The face to the cell above is not a face of the top of the cell: the cell keeps its corners, with its top at its own depth.
BOOST_AUTO_TEST_CASE(PolyhedralGridKeepsTheCornersOfTheCells)
{
    auto grid = column({flat(0.0, 10.0), flat(10.0, 10.0), flat(15.0, 25.0)});
    BOOST_REQUIRE(grid);
    BOOST_CHECK_EQUAL(Dune::unstructured::connectCellsAcrossGaps(*grid, {{0, 1}}), 1);
    Dune::PolyhedralGrid<3, 3>::UnstructuredGridPtr gridPtr(grid.release());
    Dune::PolyhedralGrid<3, 3> polyhedralGrid(std::move(gridPtr));

    const std::array<std::array<double, 2>, 2> expected{{{0.0, 10.0}, {15.0, 25.0}}};
    int numberOfCells = 0;
    for (const auto& element : elements(polyhedralGrid.leafGridView())) {
        const int cell = polyhedralGrid.leafIndexSet().index(element);
        const auto& geometry = element.geometry();
        BOOST_CHECK_EQUAL(geometry.corners(), 8);
        double lowest = 1e300;
        double highest = -1e300;
        for (int i = 0; i < geometry.corners(); ++i) {
            lowest = std::min(lowest, geometry.corner(i)[2]);
            highest = std::max(highest, geometry.corner(i)[2]);
        }
        BOOST_CHECK_CLOSE(lowest + 100.0, expected[cell][0] + 100.0, 1e-10);
        BOOST_CHECK_CLOSE(highest, expected[cell][1], 1e-10);
        ++numberOfCells;
    }
    BOOST_CHECK_EQUAL(numberOfCells, 2);
}

// The top of the cell below is at the bottom of the cell above in part, so that the faces have nodes in common. The nodes of the face of the cell
// above are not all on the sides of the cell below, and the cell keeps its corners.
BOOST_AUTO_TEST_CASE(PolyhedralGridKeepsTheCornersWithNodesInCommon)
{
    // the top of C is at the bottom of A at three of the corners, and 5 m below it at the fourth
    auto grid = column({flat(0.0, 10.0), flat(10.0, 10.0), {10.0, 10.0, 10.0, 15.0, 25.0, 25.0, 25.0, 25.0}});
    BOOST_REQUIRE(grid);
    BOOST_REQUIRE_EQUAL(grid->number_of_cells, 2);
    BOOST_CHECK_EQUAL(Dune::unstructured::connectCellsAcrossGaps(*grid, {{0, 1}}), 1);
    Dune::PolyhedralGrid<3, 3>::UnstructuredGridPtr gridPtr(grid.release());
    Dune::PolyhedralGrid<3, 3> polyhedralGrid(std::move(gridPtr));
    for (const auto& element : elements(polyhedralGrid.leafGridView())) {
        const auto& geometry = element.geometry();
        BOOST_CHECK_EQUAL(geometry.corners(), 8);
        if (polyhedralGrid.leafIndexSet().index(element) == 1) {
            // the corners of the top are those of the top of the cell
            std::array<double, 4> top{};
            for (int i = 0; i < 4; ++i) {
                top[i] = geometry.corner(i)[2];
            }
            std::sort(top.begin(), top.end());
            BOOST_CHECK_CLOSE(top[0], 10.0, 1e-10);
            BOOST_CHECK_CLOSE(top[1], 10.0, 1e-10);
            BOOST_CHECK_CLOSE(top[2], 10.0, 1e-10);
            BOOST_CHECK_CLOSE(top[3], 15.0, 1e-10);
        }
    }
}
