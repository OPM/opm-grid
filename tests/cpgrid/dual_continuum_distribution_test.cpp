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
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <algorithm>
#include <map>
#include <sstream>
#include <string>
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
