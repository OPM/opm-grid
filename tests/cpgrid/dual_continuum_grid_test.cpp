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

// Building the CpGrid of a dual-continuum run from the EclipseGrid of the
// deck: one CpGrid per continuum, from EclipseGrid::continuumGrid(), joined
// with CpGrid::extendGrid().

#include <config.h>

#define BOOST_TEST_MODULE DualContinuumGridTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <sstream>
#include <stdexcept>
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

constexpr int nx = 3;
constexpr int ny = 2;
constexpr int nz = 4; // two matrix layers, two stacked fracture layers

// Corner-point deck whose fracture layers are stacked below the matrix
// layers, with non-parallel pillars. Without DUALPORO it is an ordinary
// single-porosity grid with four layers.
Opm::Deck makeDeck(const bool dual = true, const std::string& extraGridKeywords = "")
{
    constexpr double top = 2000.0;
    constexpr double dz = 10.0;

    std::ostringstream coord;
    int pillar = 0;
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i, ++pillar) {
            const double slope = 0.03 * pillar;
            coord << 100.0 * i << ' ' << 100.0 * j << ' ' << top << ' '
                  << 100.0 * i + slope * nz * dz << ' ' << 100.0 * j << ' ' << top + nz * dz << '\n';
        }
    }

    std::ostringstream zcorn;
    for (int k = 0; k < nz; ++k) {
        for (int face = 0; face < 2; ++face) {
            for (int j = 0; j < ny; ++j) {
                for (int jj = 0; jj < 2; ++jj) {
                    for (int i = 0; i < 2 * nx; ++i) {
                        zcorn << top + (k + face) * dz << ' ';
                    }
                    zcorn << '\n';
                }
            }
        }
    }

    const int cells = nx * ny * nz;
    const std::string deckData =
        "RUNSPEC\nOIL\nWATER\nDIMENS\n " + std::to_string(nx) + ' ' + std::to_string(ny) + ' '
        + std::to_string(nz) + " /\n" + (dual ? "DUALPORO\n" : "")
        + "GRID\n" + extraGridKeywords + "COORD\n" + coord.str() + "/\nZCORN\n" + zcorn.str() + "/\n"
        "PORO\n " + std::to_string(cells) + "*0.2 /\nPERMX\n " + std::to_string(cells) + "*1.0 /\n";
    return Opm::Parser{}.parseString(deckData);
}

// Matrix cells 1 and 8 and fracture cells 12+2 and 12+9 inactive: the two
// continua have different activity.
std::vector<int> actnum()
{
    std::vector<int> a(nx * ny * nz, 1);
    for (const int c : {1, 8, 14, 21}) {
        a[c] = 0;
    }
    return a;
}

std::unique_ptr<Dune::CpGrid> makeCpGrid(const Opm::EclipseGrid& eclGrid)
{
    auto grid = std::make_unique<Dune::CpGrid>();
    grid->processEclipseFormat(&eclGrid, nullptr, false, false, false);
    return grid;
}

// processEclipseFormat() builds the dual-continuum grid from a
// dual-continuum EclipseGrid.
std::unique_ptr<Dune::CpGrid> makeDualCpGrid(const Opm::EclipseGrid& eclGrid)
{
    return makeCpGrid(eclGrid);
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(CellsAndCartesianIndicesMatchEclipseGrid)
{
    const auto deck = makeDeck();
    const auto active = actnum();
    const Opm::EclipseGrid eclGrid(deck, active.data());
    const auto grid = makeDualCpGrid(eclGrid);

    const std::array<int, 3> size{nx, ny, nz};
    BOOST_CHECK(grid->logicalCartesianSize() == size);
    BOOST_REQUIRE_EQUAL(grid->numCells(), static_cast<int>(eclGrid.getNumActive()));

    // Cell c of the CpGrid is active cell c of the EclipseGrid, so everything
    // indexed by active cell, like the field properties, lines up.
    const auto& activeMap = eclGrid.getActiveMap();
    BOOST_REQUIRE_EQUAL(activeMap.size(), grid->globalCell().size());
    for (std::size_t c = 0; c < activeMap.size(); ++c) {
        BOOST_CHECK_EQUAL(grid->globalCell()[c], static_cast<int>(activeMap[c]));
    }
}

BOOST_AUTO_TEST_CASE(TwinsAreColocatedWithoutAnyDepthOverride)
{
    const auto deck = makeDeck();
    const auto active = actnum();
    const Opm::EclipseGrid eclGrid(deck, active.data());
    const auto grid = makeDualCpGrid(eclGrid);

    std::unordered_map<int, int> cellOfGlobal;
    for (int c = 0; c < grid->numCells(); ++c) {
        cellOfGlobal[grid->globalCell()[c]] = c;
    }

    int numTwins = 0;
    for (std::size_t g = 0; g < eclGrid.getCartesianSize() / 2; ++g) {
        const auto twin = eclGrid.fractureTwin(g);
        if (!eclGrid.cellActive(g) || !eclGrid.cellActive(twin)) {
            continue;
        }
        ++numTwins;

        const int matrix = cellOfGlobal.at(static_cast<int>(g));
        const int fracture = cellOfGlobal.at(static_cast<int>(twin));
        BOOST_CHECK(eclGrid.isTwinPair(g, twin));
        BOOST_CHECK_CLOSE_FRACTION(grid->cellVolume(fracture), grid->cellVolume(matrix), 1e-12);
        for (int d = 0; d < 3; ++d) {
            BOOST_CHECK_CLOSE_FRACTION(grid->cellCentroid(fracture)[d] + 1.0,
                                       grid->cellCentroid(matrix)[d] + 1.0, 1e-12);
        }
    }
    // 12 pairs, minus the 4 with an inactive cell in either half.
    BOOST_CHECK_EQUAL(numTwins, 12 - 4);
}

BOOST_AUTO_TEST_CASE(MatrixAndFractureCellsAreNeverNeighbours)
{
    const auto deck = makeDeck();
    const Opm::EclipseGrid eclGrid(deck);

    // Stacked in a single grid, the last matrix layer and the first fracture
    // layer share faces; this is what the simulator had to zero out.
    const Opm::EclipseGrid stackedEclGrid(makeDeck(/*dual=*/false));
    const auto stacked = makeCpGrid(stackedEclGrid);
    const int half = static_cast<int>(eclGrid.getCartesianSize() / 2);
    const auto crossing = [&](const Dune::CpGrid& g)
    {
        int count = 0;
        for (int f = 0; f < g.numFaces(); ++f) {
            const int a = g.faceCell(f, 0);
            const int b = g.faceCell(f, 1);
            if (a >= 0 && b >= 0 && ((g.globalCell()[a] < half) != (g.globalCell()[b] < half))) {
                ++count;
            }
        }
        return count;
    };
    BOOST_CHECK(crossing(*stacked) > 0);

    // With one grid per continuum no face connects the continua.
    const auto dual = makeDualCpGrid(eclGrid);
    BOOST_CHECK_EQUAL(crossing(*dual), 0);

    // The same cells, just not connected across the continua.
    BOOST_CHECK_EQUAL(dual->numCells(), stacked->numCells());
}

BOOST_AUTO_TEST_CASE(MatchesJoiningTheContinuumGridsByHand)
{
    const auto deck = makeDeck();
    const auto active = actnum();
    const Opm::EclipseGrid eclGrid(deck, active.data());

    const auto automatic = makeCpGrid(eclGrid);

    auto byHand = makeCpGrid(eclGrid.continuumGrid(Opm::EclipseGrid::Continuum::Matrix));
    byHand->extendGrid(*makeCpGrid(eclGrid.continuumGrid(Opm::EclipseGrid::Continuum::Fracture)));

    BOOST_REQUIRE_EQUAL(automatic->numCells(), byHand->numCells());
    BOOST_CHECK_EQUAL(automatic->numFaces(), byHand->numFaces());
    BOOST_CHECK_EQUAL(automatic->numVertices(), byHand->numVertices());
    BOOST_CHECK(automatic->globalCell() == byHand->globalCell());
    BOOST_CHECK(automatic->logicalCartesianSize() == byHand->logicalCartesianSize());
}

BOOST_AUTO_TEST_CASE(SinglePorosityGridIsNotSplit)
{
    const Opm::EclipseGrid eclGrid(makeDeck(/*dual=*/false));
    const auto grid = makeCpGrid(eclGrid);
    BOOST_CHECK_EQUAL(grid->numCells(), nx * ny * nz);
    const std::array<int, 3> size{nx, ny, nz};
    BOOST_CHECK(grid->logicalCartesianSize() == size);
}

BOOST_AUTO_TEST_CASE(MinpvIsNotYetSupported)
{
    const Opm::EclipseGrid eclGrid(makeDeck(/*dual=*/true, "MINPV\n 1 /\n"));
    Dune::CpGrid grid;
    BOOST_CHECK_THROW(grid.processEclipseFormat(&eclGrid, nullptr, false, false, false),
                      std::invalid_argument);
}
