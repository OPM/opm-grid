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

#define BOOST_TEST_MODULE PolyhedralGridDistributionTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/utility/platform_dependent/disable_warnings.h>
#include <dune/common/parallel/mpihelper.hh>
#include <opm/grid/polyhedralgrid.hh>
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>

#include <algorithm>
#include <stdexcept>
#include <vector>

// Run with 1, 2 or 4 processes. With one process there is nothing to distribute.

using Grid = Dune::PolyhedralGrid<3, 3>;

namespace
{

constexpr int nx = 8;
constexpr int ny = 3;

struct MpiFixture
{
    MpiFixture()
    {
        int argc = boost::unit_test::framework::master_test_suite().argc;
        char** argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(argc, argv);
    }
};

// Unit sized cells. The processes get contiguous strips in the x direction.
Grid makeGrid()
{
    return Grid(std::vector<int>{nx, ny, 1}, std::vector<double>{1.0, 1.0, 1.0});
}

int stripOwner(const int i, const int numRanks)
{
    return i * numRanks / nx;
}

// The partition, only provided on rank 0.
std::vector<int> makePartition(const Grid& grid)
{
    std::vector<int> parts;
    if (grid.comm().rank() == 0) {
        const int numRanks = grid.comm().size();
        for (int c = 0; c < nx * ny; ++c) {
            parts.push_back(stripOwner(c % nx, numRanks));
        }
    }
    return parts;
}

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(MpiFixture);

BOOST_AUTO_TEST_CASE(DistributeInStrips)
{
    Grid grid = makeGrid();
    const int rank = grid.comm().rank();
    const int numRanks = grid.comm().size();

    const auto parts = makePartition(grid);
    const bool distributed = grid.loadBalance(parts, /* ownersFirst = */ true);
    BOOST_CHECK_EQUAL(distributed, numRanks > 1);
    if (numRanks == 1) {
        return;
    }

    BOOST_CHECK(grid.hasDistributedView());
    BOOST_CHECK(grid.distributedViewActive());

    // Owned cells: my strip; overlap: one column on each side that has a neighbor.
    int expectedOwned = 0;
    for (int i = 0; i < nx; ++i) {
        expectedOwned += (stripOwner(i, numRanks) == rank) ? ny : 0;
    }
    const bool hasLeft = rank > 0;
    const bool hasRight = rank < numRanks - 1;
    const int expectedOverlap = ny * (int(hasLeft) + int(hasRight));
    BOOST_CHECK_EQUAL(grid.size(0), expectedOwned + expectedOverlap);

    const auto& indexSet = grid.getCellIndexSet();
    BOOST_CHECK_EQUAL(static_cast<int>(indexSet.size()), grid.size(0));

    int numOwned = 0;
    int lastOwnerPosition = -1;
    int firstOverlapPosition = grid.size(0);
    for (auto it = grid.leafGridView().begin<0>(); it != grid.leafGridView().end<0>(); ++it) {
        const int local = grid.leafIndexSet().index(*it);
        const int globalCell = grid.globalCell()[local];
        const int i = globalCell % nx;
        const int j = globalCell / nx;

        // The cell is where it should be.
        const auto center = it->geometry().center();
        BOOST_CHECK_CLOSE(center[0], i + 0.5, 1e-10);
        BOOST_CHECK_CLOSE(center[1], j + 0.5, 1e-10);
        BOOST_CHECK_EQUAL(grid.globalIdSet().id(*it), static_cast<std::size_t>(globalCell));

        // Ownership agrees with the partition, in the index set as well.
        const bool owner = (stripOwner(i, numRanks) == rank);
        numOwned += owner;
        if (owner) {
            lastOwnerPosition = std::max(lastOwnerPosition, local);
        }
        else {
            firstOverlapPosition = std::min(firstOverlapPosition, local);
        }
        BOOST_CHECK(std::abs(stripOwner(i, numRanks) - rank) <= 1);
    }
    for (const auto& index : indexSet) {
        const int local = index.local().local();
        const int globalCell = grid.globalCell()[local];
        BOOST_CHECK_EQUAL(static_cast<int>(index.global()), globalCell);
        const bool owner = (stripOwner(globalCell % nx, numRanks) == rank);
        BOOST_CHECK(index.local().attribute() == (owner ? Grid::AttributeSet::owner
                                                         : Grid::AttributeSet::copy));
    }
    BOOST_CHECK_EQUAL(numOwned, expectedOwned);

    // Owners first.
    BOOST_CHECK(lastOwnerPosition < firstOverlapPosition);

    // Every cell is owned exactly once.
    BOOST_CHECK_EQUAL(grid.comm().sum(numOwned), nx * ny);

    // Neighboring strips share their boundary columns: ny owned cells of
    // each of them are overlap cells of the other.
    const auto& remote = grid.getCellRemoteIndices();
    int numNeighbors = 0;
    for (auto it = remote.begin(); it != remote.end(); ++it) {
        ++numNeighbors;
        BOOST_CHECK_EQUAL(std::abs(it->first - rank), 1);
        BOOST_CHECK_EQUAL(static_cast<int>(it->second.first->size()), 2 * ny);
        BOOST_CHECK_EQUAL(static_cast<int>(it->second.second->size()), 2 * ny);
    }
    BOOST_CHECK_EQUAL(numNeighbors, int(hasLeft) + int(hasRight));

    // The logical Cartesian size stays global.
    BOOST_CHECK_EQUAL(grid.logicalCartesianSize()[0], nx);
    BOOST_CHECK_EQUAL(grid.logicalCartesianSize()[1], ny);
}

BOOST_AUTO_TEST_CASE(DistributeInGlobalOrder)
{
    Grid grid = makeGrid();
    if (grid.comm().size() == 1) {
        return;
    }
    BOOST_REQUIRE(grid.loadBalance(makePartition(grid), /* ownersFirst = */ false));

    // Without owners first the cells are ordered by global index.
    int previous = -1;
    for (int local = 0; local < grid.size(0); ++local) {
        BOOST_CHECK(grid.globalCell()[local] > previous);
        previous = grid.globalCell()[local];
    }
}

BOOST_AUTO_TEST_CASE(NoOverlap)
{
    Grid grid = makeGrid();
    if (grid.comm().size() == 1) {
        return;
    }
    BOOST_REQUIRE(grid.loadBalance(makePartition(grid), true, false, /* overlapLayers = */ 0));
    for (const auto& index : grid.getCellIndexSet()) {
        BOOST_CHECK(index.local().attribute() == Grid::AttributeSet::owner);
    }
    BOOST_CHECK_EQUAL(grid.getCellRemoteIndices().neighbours(), 0u);
}

BOOST_AUTO_TEST_CASE(SwitchViewsAfterDistribution)
{
    Grid grid = makeGrid();
    if (grid.comm().size() == 1) {
        return;
    }
    BOOST_REQUIRE(grid.loadBalance(makePartition(grid)));
    const int localCells = grid.size(0);

    // Rank 0 has the whole global grid, and a second load balance is refused.
    grid.switchToGlobalView();
    if (grid.comm().rank() == 0) {
        BOOST_CHECK_EQUAL(grid.size(0), nx * ny);
    }
    BOOST_CHECK(!grid.loadBalance(makePartition(grid)));
    grid.switchToDistributedView();
    BOOST_CHECK_EQUAL(grid.size(0), localCells);
}

BOOST_AUTO_TEST_CASE(InvalidPartitionThrowsOnAllRanks)
{
    Grid grid = makeGrid();
    if (grid.comm().size() == 1) {
        return;
    }
    auto parts = makePartition(grid);
    if (grid.comm().rank() == 0) {
        parts[2] = grid.comm().size();   // no such rank
    }
    BOOST_CHECK_THROW(grid.loadBalance(parts), std::invalid_argument);
    BOOST_CHECK(!grid.hasDistributedView());
}
