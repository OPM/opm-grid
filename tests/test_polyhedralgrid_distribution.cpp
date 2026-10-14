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

#include <dune/grid/common/datahandleif.hh>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <vector>

// Run with 1, 2 or 4 processes. With one process there is nothing to distribute.

using Grid = Dune::PolyhedralGrid<3, 3>;

#if HAVE_MPI

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

namespace
{

template <Dune::PartitionIteratorType pitype>
int countEntities(const Grid& grid)
{
    int count = 0;
    for (auto it = grid.leafGridView().begin<0, pitype>(); it != grid.leafGridView().end<0, pitype>(); ++it) {
        ++count;
    }
    return count;
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(PartitionTypesAndIterators)
{
    Grid grid = makeGrid();
    const int numRanks = grid.comm().size();
    const int rank = grid.comm().rank();
    if (numRanks == 1) {
        // Everything is interior in a serial grid.
        BOOST_CHECK_EQUAL(countEntities<Dune::Interior_Partition>(grid), grid.size(0));
        BOOST_CHECK_EQUAL(countEntities<Dune::Overlap_Partition>(grid), grid.size(0));
        BOOST_CHECK_EQUAL(countEntities<Dune::Ghost_Partition>(grid), 0);
        BOOST_CHECK_EQUAL(grid.overlapSize(0), 0);
        return;
    }
    BOOST_REQUIRE(grid.loadBalance(makePartition(grid), true));
    BOOST_CHECK_EQUAL(grid.overlapSize(0), 1);

    int owned = 0;
    for (int i = 0; i < nx; ++i) {
        owned += (stripOwner(i, numRanks) == rank) ? ny : 0;
    }
    const int all = grid.size(0);

    BOOST_CHECK_EQUAL(countEntities<Dune::Interior_Partition>(grid), owned);
    BOOST_CHECK_EQUAL(countEntities<Dune::InteriorBorder_Partition>(grid), owned);
    BOOST_CHECK_EQUAL(countEntities<Dune::Overlap_Partition>(grid), all);
    BOOST_CHECK_EQUAL(countEntities<Dune::OverlapFront_Partition>(grid), all);
    BOOST_CHECK_EQUAL(countEntities<Dune::All_Partition>(grid), all);
    BOOST_CHECK_EQUAL(countEntities<Dune::Ghost_Partition>(grid), 0);

    // The partition type of each cell follows its ownership, and the
    // partition grid views agree with the iterators.
    int numInterior = 0;
    for (auto it = grid.leafGridView().begin<0>(); it != grid.leafGridView().end<0>(); ++it) {
        const int i = grid.globalCell()[grid.leafIndexSet().index(*it)] % nx;
        const bool isMine = (stripOwner(i, numRanks) == rank);
        BOOST_CHECK(it->partitionType() == (isMine ? Dune::InteriorEntity : Dune::OverlapEntity));
        numInterior += (it->partitionType() == Dune::InteriorEntity);
    }
    BOOST_CHECK_EQUAL(numInterior, owned);

    int viewCount = 0;
    const auto interiorView = grid.leafGridView<Dune::Interior_Partition>();
    for (auto it = interiorView.begin<0>(); it != interiorView.end<0>(); ++it) {
        BOOST_CHECK(it->partitionType() == Dune::InteriorEntity);
        ++viewCount;
    }
    BOOST_CHECK_EQUAL(viewCount, owned);

    // Other codimensions are not partitioned: they are all interior.
    int numFaces = 0;
    for (auto it = grid.leafGridView().begin<3, Dune::Interior_Partition>();
         it != grid.leafGridView().end<3, Dune::Interior_Partition>(); ++it) {
        ++numFaces;
    }
    BOOST_CHECK_EQUAL(numFaces, grid.size(3));
}

BOOST_AUTO_TEST_CASE(IntersectionsAtPartitionBorders)
{
    for (const int layers : {1, 0}) {
        Grid grid = makeGrid();
        const int numRanks = grid.comm().size();
        const int rank = grid.comm().rank();
        if (numRanks == 1) {
            return;
        }
        BOOST_REQUIRE(grid.loadBalance(makePartition(grid), true, false, layers));

        for (auto it = grid.leafGridView().begin<0, Dune::Interior_Partition>();
             it != grid.leafGridView().end<0, Dune::Interior_Partition>(); ++it) {
            const int gc = grid.globalCell()[grid.leafIndexSet().index(*it)];
            const int i = gc % nx;
            const int j = gc / nx;
            // The boundary of the whole domain, including top and bottom.
            int physical = int(i == 0) + int(i == nx - 1) + int(j == 0) + int(j == ny - 1) + 2;
            // Faces to cells that are on another process but not part of the grid.
            int onBorder = 0;
            if (layers == 0) {
                onBorder = int(i > 0 && stripOwner(i - 1, numRanks) != rank)
                         + int(i < nx - 1 && stripOwner(i + 1, numRanks) != rank);
            }

            int boundary = 0;
            int neighbors = 0;
            for (auto is = grid.leafGridView().ibegin(*it); is != grid.leafGridView().iend(*it); ++is) {
                boundary += is->boundary();
                if (is->neighbor()) {
                    ++neighbors;
                    // A cell across the partition border is an overlap cell here.
                    const int ni = grid.globalCell()[grid.leafIndexSet().index(is->outside())] % nx;
                    const bool foreign = (stripOwner(ni, numRanks) != rank);
                    BOOST_CHECK(is->outside().partitionType() == (foreign ? Dune::OverlapEntity
                                                                          : Dune::InteriorEntity));
                }
            }
            BOOST_CHECK_EQUAL(boundary, physical + onBorder);
            BOOST_CHECK_EQUAL(boundary + neighbors, 6);
        }
    }
}

namespace
{

// The Cartesian index of a cell. A serial grid has no global_cell array.
int globalCellOf(const Grid& grid, const int cell)
{
    const int* globalCell = grid.globalCellPtr();
    return globalCell ? globalCell[cell] : cell;
}

// One int per cell. The received value replaces, or is added to, the local one.
class CellValueHandle : public Dune::CommDataHandleIF<CellValueHandle, int>
{
public:
    CellValueHandle(const Grid& grid, std::vector<int>& values, const bool add = false)
        : grid_(grid), values_(values), add_(add)
    {}

    bool contains(int, int codim) const { return codim == 0; }
    bool fixedSize(int, int) const { return true; }

    template <class Entity>
    std::size_t size(const Entity&) const { return 1; }

    template <class Buffer, class Entity>
    void gather(Buffer& buffer, const Entity& entity) const
    {
        buffer.write(values_[grid_.leafIndexSet().index(entity)]);
    }

    template <class Buffer, class Entity>
    void scatter(Buffer& buffer, const Entity& entity, std::size_t)
    {
        int value;
        buffer.read(value);
        auto& target = values_[grid_.leafIndexSet().index(entity)];
        target = add_ ? target + value : value;
    }

private:
    const Grid& grid_;
    std::vector<int>& values_;
    bool add_;
};

// 1 + globalCell % 3 values per cell, which tests communication of variable sizes.
class CellVectorHandle : public Dune::CommDataHandleIF<CellVectorHandle, int>
{
public:
    CellVectorHandle(const Grid& grid, std::vector<std::vector<int>>& values)
        : grid_(grid), values_(values)
    {}

    bool contains(int, int codim) const { return codim == 0; }
    bool fixedSize(int, int) const { return false; }

    template <class Entity>
    std::size_t size(const Entity& entity) const
    {
        return values_[grid_.leafIndexSet().index(entity)].size();
    }

    template <class Buffer, class Entity>
    void gather(Buffer& buffer, const Entity& entity) const
    {
        for (const int v : values_[grid_.leafIndexSet().index(entity)]) {
            buffer.write(v);
        }
    }

    template <class Buffer, class Entity>
    void scatter(Buffer& buffer, const Entity& entity, std::size_t n)
    {
        auto& target = values_[grid_.leafIndexSet().index(entity)];
        target.resize(n);
        for (auto& v : target) {
            buffer.read(v);
        }
    }

private:
    const Grid& grid_;
    std::vector<std::vector<int>>& values_;
};

int globalColumn(const Grid& grid, const int local)
{
    return globalCellOf(grid, local) % nx;
}

// Number of other processes that have the (owned) cell in their overlap, for
// one overlap layer.
int numOverlapCopies(const int i, const int rank, const int numRanks)
{
    return int(i > 0 && stripOwner(i - 1, numRanks) != rank)
         + int(i < nx - 1 && stripOwner(i + 1, numRanks) != rank);
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(CommunicateOwnerToAll)
{
    Grid grid = makeGrid();
    const int rank = grid.comm().rank();
    const int numRanks = grid.comm().size();
    if (numRanks > 1) {
        BOOST_REQUIRE(grid.loadBalance(makePartition(grid), true));
    }

    std::vector<int> values(grid.size(0), -1);
    for (int c = 0; c < grid.size(0); ++c) {
        if (grid.entityPartitionType<0>(c) == Dune::InteriorEntity) {
            values[c] = globalCellOf(grid, c);
        }
    }

    // Through the grid view, as simulators do.
    CellValueHandle handle(grid, values);
    grid.leafGridView().communicate(handle, Dune::InteriorBorder_All_Interface,
                                    Dune::ForwardCommunication);

    for (int c = 0; c < grid.size(0); ++c) {
        BOOST_CHECK_EQUAL(values[c], globalCellOf(grid, c));
    }
    static_cast<void>(rank);
}

BOOST_AUTO_TEST_CASE(CommunicateAllToOwnerAdds)
{
    Grid grid = makeGrid();
    const int rank = grid.comm().rank();
    const int numRanks = grid.comm().size();
    if (numRanks == 1) {
        return;
    }
    BOOST_REQUIRE(grid.loadBalance(makePartition(grid), true));

    std::vector<int> values(grid.size(0), 1);
    CellValueHandle handle(grid, values, /* add = */ true);
    grid.communicate(handle, Dune::InteriorBorder_All_Interface, Dune::BackwardCommunication);

    for (int c = 0; c < grid.size(0); ++c) {
        if (grid.entityPartitionType<0>(c) == Dune::InteriorEntity) {
            // The cell itself plus its copies on the other processes.
            BOOST_CHECK_EQUAL(values[c], 1 + numOverlapCopies(globalColumn(grid, c), rank, numRanks));
        }
        else {
            BOOST_CHECK_EQUAL(values[c], 1);
        }
    }
}

BOOST_AUTO_TEST_CASE(CommunicateAllToAll)
{
    Grid grid = makeGrid();
    const int rank = grid.comm().rank();
    const int numRanks = grid.comm().size();
    if (numRanks == 1) {
        return;
    }
    BOOST_REQUIRE(grid.loadBalance(makePartition(grid), true));

    // Every cell that is on two processes gets the value of the other one.
    std::vector<int> values(grid.size(0), rank);
    CellValueHandle handle(grid, values);
    grid.communicate(handle, Dune::All_All_Interface, Dune::ForwardCommunication);

    for (int c = 0; c < grid.size(0); ++c) {
        const int i = globalColumn(grid, c);
        int expected = rank;
        if (grid.entityPartitionType<0>(c) == Dune::OverlapEntity) {
            expected = stripOwner(i, numRanks);
        }
        else if (numOverlapCopies(i, rank, numRanks) == 1) {
            expected = (i > 0 && stripOwner(i - 1, numRanks) != rank) ? stripOwner(i - 1, numRanks)
                                                                     : stripOwner(i + 1, numRanks);
        }
        BOOST_CHECK_EQUAL(values[c], expected);
    }
}

BOOST_AUTO_TEST_CASE(CommunicateBetweenOverlapCells)
{
    Grid grid = makeGrid();
    const int rank = grid.comm().rank();
    if (grid.comm().size() != 4) {
        return;
    }
    // With two layers, the cells in the columns 2 and 3 are overlap cells on
    // rank 0 and 2, and the ones in the columns 4 and 5 on rank 1 and 3.
    BOOST_REQUIRE(grid.loadBalance(makePartition(grid), true, false, 2));

    std::vector<int> values(grid.size(0), rank);
    CellValueHandle handle(grid, values);
    grid.communicate(handle, Dune::Overlap_OverlapFront_Interface, Dune::ForwardCommunication);

    for (int c = 0; c < grid.size(0); ++c) {
        const int i = globalColumn(grid, c);
        int expected = rank;
        if (grid.entityPartitionType<0>(c) == Dune::OverlapEntity) {
            if ((i == 2 || i == 3) && (rank == 0 || rank == 2)) {
                expected = 2 - rank;
            }
            if ((i == 4 || i == 5) && (rank == 1 || rank == 3)) {
                expected = 4 - rank;
            }
        }
        BOOST_CHECK_EQUAL(values[c], expected);
    }
}

BOOST_AUTO_TEST_CASE(CommunicateVariableSize)
{
    Grid grid = makeGrid();
    if (grid.comm().size() > 1) {
        BOOST_REQUIRE(grid.loadBalance(makePartition(grid), true));
    }

    auto expected = [&](const int c) {
        std::vector<int> v(1 + globalCellOf(grid, c) % 3);
        for (std::size_t k = 0; k < v.size(); ++k) {
            v[k] = 10 * globalCellOf(grid, c) + static_cast<int>(k);
        }
        return v;
    };

    std::vector<std::vector<int>> values(grid.size(0));
    for (int c = 0; c < grid.size(0); ++c) {
        if (grid.entityPartitionType<0>(c) == Dune::InteriorEntity) {
            values[c] = expected(c);
        }
    }
    CellVectorHandle handle(grid, values);
    grid.communicate(handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
    for (int c = 0; c < grid.size(0); ++c) {
        BOOST_CHECK(values[c] == expected(c));
    }
}

BOOST_AUTO_TEST_CASE(CommunicateInGlobalViewDoesNothing)
{
    Grid grid = makeGrid();
    if (grid.comm().size() > 1) {
        BOOST_REQUIRE(grid.loadBalance(makePartition(grid), true));
        grid.switchToGlobalView();
    }
    std::vector<int> values(grid.size(0), 7);
    CellValueHandle handle(grid, values, true);
    BOOST_CHECK_NO_THROW(grid.communicate(handle, Dune::All_All_Interface, Dune::ForwardCommunication));
    BOOST_CHECK(std::all_of(values.begin(), values.end(), [](int v) { return v == 7; }));
}

namespace
{

// As in the simulators, only rank 0 holds the global grid. The other ranks
// have a grid without cells.
std::unique_ptr<Grid> makeGridOnRoot()
{
    if (Dune::MPIHelper::getCommunication().rank() == 0) {
        return std::make_unique<Grid>(std::vector<int>{nx, ny, 1}, std::vector<double>{1.0, 1.0, 1.0});
    }
    Grid::UnstructuredGridPtr empty(allocate_grid(3, 0, 0, 0, 0, 0));
    empty->cartdims[0] = nx;
    empty->cartdims[1] = ny;
    empty->cartdims[2] = 1;
    return std::make_unique<Grid>(std::move(empty));
}

// Moves variable sized data per cell: gather() reads the values of the cell in
// `from`, and scatter() writes the received values to the cell in `to`. Which
// grid view the cells are of is up to the caller.
class FromToHandle : public Dune::CommDataHandleIF<FromToHandle, int>
{
public:
    FromToHandle(const Grid& grid,
                 std::vector<std::vector<int>>& from,
                 std::vector<std::vector<int>>& to)
        : grid_(grid), from_(from), to_(to)
    {}

    bool contains(int, int codim) const { return codim == 0; }
    bool fixedSize(int, int) const { return false; }

    template <class Entity>
    std::size_t size(const Entity& entity) const
    {
        return from_[grid_.leafIndexSet().index(entity)].size();
    }

    template <class Buffer, class Entity>
    void gather(Buffer& buffer, const Entity& entity) const
    {
        for (const int v : from_[grid_.leafIndexSet().index(entity)]) {
            buffer.write(v);
        }
    }

    template <class Buffer, class Entity>
    void scatter(Buffer& buffer, const Entity& entity, std::size_t n)
    {
        auto& target = to_[grid_.leafIndexSet().index(entity)];
        target.resize(n);
        for (auto& v : target) {
            buffer.read(v);
        }
    }

private:
    const Grid& grid_;
    std::vector<std::vector<int>>& from_;
    std::vector<std::vector<int>>& to_;
};

// The data of a global cell: 1 + c % 3 values, which depend on the cell and a tag.
std::vector<int> cellValues(const int globalCell, const int tag)
{
    std::vector<int> v(1 + globalCell % 3);
    for (std::size_t k = 0; k < v.size(); ++k) {
        v[k] = tag + 10 * globalCell + static_cast<int>(k);
    }
    return v;
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(ScatterDataFromRoot)
{
    auto grid = makeGridOnRoot();
    const int rank = grid->comm().rank();
    if (grid->comm().size() == 1) {
        return;
    }
    BOOST_REQUIRE(grid->loadBalance(makePartition(*grid), true));

    // Rank 0 has the data of the global grid; the other ranks have none.
    grid->switchToGlobalView();
    std::vector<std::vector<int>> globalValues(rank == 0 ? nx * ny : 0);
    for (int c = 0; c < static_cast<int>(globalValues.size()); ++c) {
        globalValues[c] = cellValues(c, 100000);
    }
    std::vector<std::vector<int>> localValues;

    grid->switchToDistributedView();
    localValues.resize(grid->size(0));
    grid->switchToGlobalView();
    FromToHandle handle(*grid, globalValues, localValues);
    grid->scatterData(handle);

    // The view is as it was before.
    BOOST_CHECK(!grid->distributedViewActive());
    grid->switchToDistributedView();
    BOOST_REQUIRE_EQUAL(static_cast<int>(localValues.size()), grid->size(0));
    for (int c = 0; c < grid->size(0); ++c) {
        BOOST_CHECK(localValues[c] == cellValues(globalCellOf(*grid, c), 100000));
    }
}

BOOST_AUTO_TEST_CASE(GatherDataToRoot)
{
    auto grid = makeGridOnRoot();
    const int rank = grid->comm().rank();
    if (grid->comm().size() == 1) {
        return;
    }
    BOOST_REQUIRE(grid->loadBalance(makePartition(*grid), true));

    // Owned cells have data, overlap cells do not.
    std::vector<std::vector<int>> localValues(grid->size(0));
    for (int c = 0; c < grid->size(0); ++c) {
        if (grid->entityPartitionType<0>(c) == Dune::InteriorEntity) {
            localValues[c] = cellValues(globalCellOf(*grid, c), 200000);
        }
    }
    grid->switchToGlobalView();
    std::vector<std::vector<int>> globalValues(rank == 0 ? nx * ny : 0);
    grid->switchToDistributedView();

    FromToHandle handle(*grid, localValues, globalValues);
    grid->gatherData(handle);

    BOOST_CHECK(grid->distributedViewActive());
    if (rank == 0) {
        for (int c = 0; c < nx * ny; ++c) {
            BOOST_CHECK(globalValues[c] == cellValues(c, 200000));
        }
    }
}

BOOST_AUTO_TEST_CASE(LoadBalanceWithData)
{
    auto grid = makeGridOnRoot();
    const int rank = grid->comm().rank();
    if (grid->comm().size() == 1) {
        return;
    }

    std::vector<std::vector<int>> globalValues(rank == 0 ? nx * ny : 0);
    for (int c = 0; c < static_cast<int>(globalValues.size()); ++c) {
        globalValues[c] = cellValues(c, 300000);
    }
    std::vector<std::vector<int>> localValues;

    // The distributed grid does not exist yet, so there is room for all cells.
    localValues.resize(nx * ny);
    FromToHandle handle(*grid, globalValues, localValues);
    BOOST_REQUIRE(grid->loadBalance(handle, makePartition(*grid), true));
    BOOST_CHECK(grid->distributedViewActive());

    for (int c = 0; c < grid->size(0); ++c) {
        BOOST_CHECK(localValues[c] == cellValues(globalCellOf(*grid, c), 300000));
    }
}

#else

// The distribution needs MPI, so there is nothing to test without it.
BOOST_AUTO_TEST_CASE(RequiresMpi)
{
}

#endif // HAVE_MPI
