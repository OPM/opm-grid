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

#define BOOST_TEST_MODULE PolyhedralGridViewsTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/utility/platform_dependent/disable_warnings.h>
#include <dune/common/parallel/mpihelper.hh>
#include <opm/grid/polyhedralgrid.hh>
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>

#include <opm/grid/common/UnstructuredGridPartitioning.hpp>
#include <opm/grid/polyhedralgrid/cartesianindexmapper.hh>

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

// Cells owned by the two leftmost columns of a 4x2x1 grid plus one overlap layer.
std::vector<int> localCells(const Grid& grid)
{
    const auto& g = static_cast<const UnstructuredGrid&>(grid);
    std::vector<int> part(g.number_of_cells);
    for (int c = 0; c < g.number_of_cells; ++c) {
        part[c] = ((c % 4) < 2) ? 0 : 1;
    }
    std::vector<int> cells;
    for (int c = 0; c < g.number_of_cells; ++c) {
        if (part[c] == 0) {
            cells.push_back(c);
        }
    }
    const auto overlap = Dune::unstructured::computeOverlapCells(g, part, 2, 1);
    cells.insert(cells.end(), overlap[0].begin(), overlap[0].end());
    return cells;
}

Grid::UnstructuredGridPtr localGrid(const Grid& grid, const std::vector<int>& cells)
{
    auto local = Dune::unstructured::extractLocalGrid(grid, cells);
    return Grid::UnstructuredGridPtr(local.grid.release());
}

int countCells(const Grid& grid)
{
    int n = 0;
    for (auto it = grid.leafGridView().begin<0>(); it != grid.leafGridView().end<0>(); ++it) {
        ++n;
    }
    return n;
}

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(MpiFixture);

BOOST_AUTO_TEST_CASE(NoDistributedViewByDefault)
{
    Grid grid(std::vector<int>{4, 2, 1}, std::vector<double>{1.0, 1.0, 1.0});
    BOOST_CHECK(!grid.hasDistributedView());
    BOOST_CHECK(!grid.distributedViewActive());
    BOOST_CHECK_THROW(grid.switchToDistributedView(), std::logic_error);
    // Switching to the global view is always fine.
    grid.switchToGlobalView();
    BOOST_CHECK_EQUAL(grid.size(0), 8);
}

BOOST_AUTO_TEST_CASE(SwitchBetweenViews)
{
    Grid grid(std::vector<int>{4, 2, 1}, std::vector<double>{1.0, 1.0, 1.0});
    const auto cells = localCells(grid);   // 4 owned + 2 overlap cells
    BOOST_REQUIRE_EQUAL(cells.size(), 6u);

    // Copy of the global cell centers, since the active grid changes.
    std::vector<std::array<double, 3>> globalCenters(8);
    for (auto it = grid.leafGridView().begin<0>(); it != grid.leafGridView().end<0>(); ++it) {
        const auto c = it->geometry().center();
        globalCenters[grid.leafIndexSet().index(*it)] = {c[0], c[1], c[2]};
    }

    grid.setDistributedView(localGrid(grid, cells));
    BOOST_CHECK(grid.hasDistributedView());
    // Installing a view does not activate it.
    BOOST_CHECK(!grid.distributedViewActive());
    BOOST_CHECK_EQUAL(grid.size(0), 8);
    BOOST_CHECK_EQUAL(countCells(grid), 8);

    grid.switchToDistributedView();
    BOOST_CHECK(grid.distributedViewActive());
    BOOST_CHECK_EQUAL(grid.size(0), 6);
    BOOST_CHECK_EQUAL(countCells(grid), 6);

    // Logical Cartesian size and cell numbering stay global.
    BOOST_CHECK_EQUAL(grid.logicalCartesianSize()[0], 4);
    BOOST_CHECK_EQUAL(grid.logicalCartesianSize()[1], 2);
    Dune::CartesianIndexMapper<Grid> mapper(grid);
    BOOST_CHECK_EQUAL(mapper.cartesianSize(), 8);
    BOOST_CHECK_EQUAL(mapper.compressedSize(), 6);

    double volume = 0.0;
    for (auto it = grid.leafGridView().begin<0>(); it != grid.leafGridView().end<0>(); ++it) {
        const int local = grid.leafIndexSet().index(*it);
        BOOST_CHECK_EQUAL(mapper.cartesianIndex(local), cells[local]);
        BOOST_CHECK_EQUAL(grid.globalIdSet().id(*it), static_cast<std::size_t>(cells[local]));
        const auto c = it->geometry().center();
        for (int d = 0; d < 3; ++d) {
            BOOST_CHECK_CLOSE(c[d] + 1.0, globalCenters[cells[local]][d] + 1.0, 1e-10);
        }
        volume += it->geometry().volume();
    }
    BOOST_CHECK_CLOSE(volume, 6.0, 1e-10);

    // The faces of a cell are available in the local grid.
    for (auto it = grid.leafGridView().begin<0>(); it != grid.leafGridView().end<0>(); ++it) {
        int faces = 0;
        for (auto is = grid.leafGridView().ibegin(*it); is != grid.leafGridView().iend(*it); ++is) {
            ++faces;
        }
        BOOST_CHECK_EQUAL(faces, 6);
    }

    grid.switchToGlobalView();
    BOOST_CHECK(!grid.distributedViewActive());
    BOOST_CHECK_EQUAL(grid.size(0), 8);
    BOOST_CHECK_EQUAL(countCells(grid), 8);
    for (auto it = grid.leafGridView().begin<0>(); it != grid.leafGridView().end<0>(); ++it) {
        const int index = grid.leafIndexSet().index(*it);
        BOOST_CHECK_EQUAL(grid.globalIdSet().id(*it), static_cast<std::size_t>(index));
    }
}

BOOST_AUTO_TEST_CASE(ReplaceDistributedViewWhileActive)
{
    Grid grid(std::vector<int>{4, 2, 1}, std::vector<double>{1.0, 1.0, 1.0});
    auto cells = localCells(grid);
    grid.setDistributedView(localGrid(grid, cells));
    grid.switchToDistributedView();
    BOOST_CHECK_EQUAL(grid.size(0), 6);

    // Replacing the view keeps the distributed view active.
    const std::vector<int> fewer{0, 1};
    grid.setDistributedView(localGrid(grid, fewer));
    BOOST_CHECK(grid.distributedViewActive());
    BOOST_CHECK_EQUAL(grid.size(0), 2);
    BOOST_CHECK_EQUAL(countCells(grid), 2);

    grid.switchToGlobalView();
    BOOST_CHECK_EQUAL(grid.size(0), 8);
}
