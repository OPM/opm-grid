/*
  Copyright 2026 NORCE Norwegian Research Centre AS.

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

// Distributing the grid of a dual-continuum run. A matrix cell and its fracture
// twin have to be on the same process, both owned or both overlap.

#include <config.h>

#define BOOST_TEST_MODULE DualContinuumDistributionTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>
#include <opm/input/eclipse/EclipseState/Aquifer/NumericalAquifer/NumericalAquifers.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Python/Python.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <algorithm>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct Fixture
{
    Fixture()
    {
        int argc = boost::unit_test::framework::master_test_suite().argc;
        char** argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(argc, argv);
    }
};

BOOST_GLOBAL_FIXTURE(Fixture);

constexpr int nx = 8;
constexpr int ny = 6;
constexpr int nz = 4; // two matrix layers, two fracture layers

// The matrix has the first two layers. With 'holes', some cells are inactive in one continuum
// only, and with 'sigma' the continua are coupled.
Opm::Deck makeDeck(const bool holes, const double sigma)
{
    const int cells = nx * ny * nz;
    std::ostringstream actnum;
    for (int c = 0; c < cells; ++c) {
        const bool hole = holes && (c == 10 || c == 45 || c == cells / 2 + 20 || c == cells / 2 + 77);
        actnum << (hole ? 0 : 1) << ' ';
    }
    std::ostringstream deck;
    deck << "RUNSPEC\nOIL\nWATER\nDIMENS\n " << nx << ' ' << ny << ' ' << nz << " /\nDUALPORO\nGRID\n"
         << "DX\n " << cells << "*100 /\nDY\n " << cells << "*100 /\nDZ\n " << cells << "*10 /\n"
         << "TOPS\n " << nx * ny << "*2000 /\nACTNUM\n" << actnum.str() << "/\n"
         << "PORO\n " << cells << "*0.2 /\nPERMX\n " << cells << "*10 /\nSIGMA\n " << sigma << " /\n";
    return Opm::Parser{}.parseString(deck.str());
}

// Cartesian index and whether the cell is owned, of every cell of the distributed grid.
std::map<int, bool> localCells(const Dune::CpGrid& grid)
{
    std::map<int, bool> cells;
    const auto& globalCell = grid.globalCell();
    for (const auto& element : elements(grid.leafGridView())) {
        const auto index = grid.leafGridView().indexSet().index(element);
        cells[globalCell[index]] = element.partitionType() == Dune::InteriorEntity;
    }
    return cells;
}

void checkTwinsTogether(const bool holes, const double sigma)
{
    // Distribution needs more than one process, see the _parallel ctest.
    if (Dune::CpGrid().comm().size() == 1) {
        return;
    }
    const auto deck = makeDeck(holes, sigma);
    Opm::EclipseState eclState(deck);

    Dune::CpGrid grid;
    const bool root = grid.comm().rank() == 0;
    grid.processEclipseFormat(root ? &eclState.getInputGrid() : nullptr, &eclState, false, false, false);

    const int numCells = root ? grid.numCells() : 0;
    const int total = grid.comm().sum(numCells);
    BOOST_REQUIRE(total > 0);

    BOOST_REQUIRE(grid.loadBalance(/* overlapLayers = */ 1, Dune::PartitionMethod::zoltan));

    const auto cells = localCells(grid);
    const int half = nx * ny * nz / 2;

    // The distributed grid knows the offset between a cell and its twin, as the global grid does.
    BOOST_CHECK_EQUAL(grid.currentData().front()->continuumCartOffset(), static_cast<std::size_t>(half));

    int owned = 0;
    for (const auto& [cartesian, isOwned] : cells) {
        owned += isOwned;
        const int twin = cartesian < half ? cartesian + half : cartesian - half;
        // A twin that is active is here, and with the same ownership.
        const auto it = cells.find(twin);
        if (it == cells.end()) {
            // The twin of an owned cell is here, unless it is inactive in the deck. The twin of an
            // overlap cell is missing if the neighbours of the twin are inactive.
            BOOST_CHECK_MESSAGE(!isOwned || !eclState.getInputGrid().cellActive(twin),
                                "the active twin " << twin << " of the owned cell " << cartesian << " is missing");
            continue;
        }
        BOOST_CHECK_MESSAGE(it->second == isOwned,
                            "cell " << cartesian << " and its twin " << twin << " differ in ownership");
    }

    const int allOwned = grid.comm().sum(owned);
    BOOST_CHECK_EQUAL(allOwned, total);
    // Every process has cells, and the cells are spread over the processes.
    BOOST_CHECK(owned > 0);
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(TwinsAreOnTheSameProcess)
{
    checkTwinsTogether(/* holes = */ false, 0.1);
}

BOOST_AUTO_TEST_CASE(TwinsAreOnTheSameProcessWithInactiveCells)
{
    checkTwinsTogether(/* holes = */ true, 0.1);
}

BOOST_AUTO_TEST_CASE(TwinsAreOnTheSameProcessWithoutCoupling)
{
    // No coupling, so there is no connection between the twins that the partitioner sees.
    checkTwinsTogether(/* holes = */ false, 0.0);
}

BOOST_AUTO_TEST_CASE(ExtendedGridCanBeDistributed)
{
    if (Dune::CpGrid().comm().size() == 1) {
        return;
    }
    // Two Cartesian grids joined into one, with the cells of both in each column on the same process.
    Dune::CpGrid grid;
    grid.createCartesian({4, 3, 2}, {1.0, 1.0, 1.0});
    Dune::CpGrid other;
    other.createCartesian({4, 3, 2}, {1.0, 1.0, 1.0});
    std::vector<int> parts;
    if (grid.comm().rank() == 0) {
        grid.appendDisjointGrid(other);
        for (int c = 0; c < grid.numCells(); ++c) {
            parts.push_back((c % 12) % grid.comm().size());
        }
    }
    BOOST_REQUIRE(grid.loadBalance(parts, /* ownersFirst = */ false, /* addCornerCells = */ false, 1));
    BOOST_CHECK(grid.numCells() > 0);
    BOOST_CHECK_EQUAL(grid.comm().sum(grid.numCells()) >= 48, true);
}

namespace {

// Cartesian index of the cell with the one-based indices (i, j, k).
int cartesian(const int i, const int j, const int k)
{
    return (i - 1) + (j - 1) * nx + (k - 1) * nx * ny;
}

// Whether one process owns all the cells.
bool ownedOnOneProcess(const Dune::CpGrid& grid, const std::map<int, bool>& cells, const std::set<int>& wanted)
{
    int owned = 0;
    for (const auto c : wanted) {
        const auto it = cells.find(c);
        owned += (it != cells.end() && it->second);
    }
    return grid.comm().max(owned) == static_cast<int>(wanted.size());
}

} // anonymous namespace

static void checkWellsAndTwins(const int partitionMethod)
{
    if (Dune::CpGrid().comm().size() == 1) {
        return;
    }

    // W1 is in the matrix, W3 in the fracture continuum. W2 has a cell in the fracture twin
    // of a cell of W1, as a well can perforate both continua, so that the partition that
    // keeps W1 together and the one that keeps W2 together have to meet.
    const std::string deckText =
        "RUNSPEC\nOIL\nWATER\nDIMENS\n 8 6 4 /\nDUALPORO\nSTART\n 1 'JAN' 2015 /\n"
        "WELLDIMS\n 3 8 1 3 /\n"
        "GRID\nDX\n 192*100 /\nDY\n 192*100 /\nDZ\n 192*10 /\nTOPS\n 48*2000 /\n"
        "PORO\n 192*0.2 /\nPERMX\n 192*10 /\nPERMY\n 192*10 /\nPERMZ\n 192*10 /\nSIGMA\n 0.1 /\n"
        "SCHEDULE\n"
        "WELSPECS\n 'W1' 'G' 1 1 2005 'OIL' /\n 'W2' 'G' 1 1 2025 'OIL' /\n 'W3' 'G' 4 3 2035 'OIL' /\n/\n"
        "COMPDAT\n"
        " 'W1' 1 1 1 1 'OPEN' 1* 1* 0.5 /\n 'W1' 8 1 1 1 'OPEN' 1* 1* 0.5 /\n"
        " 'W2' 1 1 3 3 'OPEN' 1* 1* 0.5 /\n 'W2' 8 6 3 3 'OPEN' 1* 1* 0.5 /\n"
        " 'W3' 4 3 4 4 'OPEN' 1* 1* 0.5 /\n 'W3' 5 3 4 4 'OPEN' 1* 1* 0.5 /\n/\n"
        "END\n";
    const auto deck = Opm::Parser{}.parseString(deckText);
    Opm::EclipseState eclState(deck);
    Opm::EclipseGrid eclGrid = eclState.getInputGrid();
    const Opm::TableManager table(deck);
    const Opm::FieldPropsManager fp(deck, Opm::Phases{true, true, false}, eclGrid, table);
    const Opm::Runspec runspec(deck);
    const Opm::Schedule schedule{deck, eclGrid, fp, Opm::NumericalAquifers{}, runspec,
                                 std::make_shared<Opm::Python>()};
    const auto wells = schedule.getWellsatEnd();
    BOOST_REQUIRE_EQUAL(wells.size(), 3U);

    // Cells of a well that may get a connection later.
    const std::unordered_map<std::string, std::set<int>> future{
        {"W3", {cartesian(2, 5, 2), cartesian(7, 2, 4)}}};

    Dune::CpGrid grid;
    const bool root = grid.comm().rank() == 0;
    grid.processEclipseFormat(root ? &eclGrid : nullptr, &eclState, false, false, false);
    BOOST_REQUIRE(grid.loadBalance(&wells, future, nullptr, 1, partitionMethod).first);

    const auto cells = localCells(grid);
    const int half = nx * ny * nz / 2;

    // The cells of a well, and their twins, are owned by one process.
    std::map<std::string, std::set<int>> wellCells;
    for (const auto& well : wells) {
        for (const auto& connection : well.getConnections()) {
            const auto c = static_cast<int>(connection.global_index());
            const auto twin = c < half ? c + half : c - half;
            wellCells[well.name()].insert({c, twin});
        }
    }
    for (const auto& [name, c] : future) {
        for (const auto cell : c) {
            wellCells[name].insert({cell, cell < half ? cell + half : cell - half});
        }
    }
    for (const auto& [name, wanted] : wellCells) {
        BOOST_CHECK_MESSAGE(ownedOnOneProcess(grid, cells, wanted),
                            "the cells of well " << name << " and their twins are not on one process");
    }

    // W1 and W2 meet in a twin pair, so they are on the same process.
    auto both = wellCells["W1"];
    both.insert(wellCells["W2"].begin(), wellCells["W2"].end());
    BOOST_CHECK_MESSAGE(ownedOnOneProcess(grid, cells, both),
                        "the wells W1 and W2 share a twin pair, and are not on one process");
}

BOOST_AUTO_TEST_CASE(WellsAndTheirTwinsAreOnOneProcess)
{
    checkWellsAndTwins(Dune::PartitionMethod::zoltan);
}

BOOST_AUTO_TEST_CASE(WellsAndTheirTwinsAreOnOneProcessWithoutZoltan)
{
    // The partitioner does not know the twins, so they are put together afterwards,
    // together with the cells of the wells that they join.
    checkWellsAndTwins(Dune::PartitionMethod::simple);
}
