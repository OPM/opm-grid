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
#include <opm/grid/cpgrid/GridHelpers.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/Grid/NNC.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
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
Opm::Deck makeDeck(const bool dual = true, const std::string& extraGridKeywords = "",
                   const std::map<int, double>& poroOverride = {})
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
    std::ostringstream poro;
    for (int c = 0; c < cells; ++c) {
        const auto it = poroOverride.find(c);
        poro << (it != poroOverride.end() ? it->second : 0.2) << ' ';
    }
    const std::string deckData =
        "RUNSPEC\nOIL\nWATER\nDIMENS\n " + std::to_string(nx) + ' ' + std::to_string(ny) + ' '
        + std::to_string(nz) + " /\n" + (dual ? "DUALPORO\n" : "")
        + "GRID\n" + extraGridKeywords + "COORD\n" + coord.str() + "/\nZCORN\n" + zcorn.str() + "/\n"
        "PORO\n " + poro.str() + "/\nPERMX\n " + std::to_string(cells) + "*1.0 /\n"
        + (dual ? "SIGMA\n 0.12 /\n" : "");
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

BOOST_AUTO_TEST_CASE(FacesAndCellsReferToEachOther)
{
    for (const bool coupled : {false, true}) {
        const auto deck = makeDeck(true, "");
        Opm::EclipseState eclState(deck);
        Dune::CpGrid grid;
        grid.processEclipseFormat(&eclState.getInputGrid(), coupled ? &eclState : nullptr, false, false, false);

        // Every face of a cell has the cell among its cells, and every cell of a face has the
        // face among its faces.
        for (int c = 0; c < grid.numCells(); ++c) {
            const auto faces = grid.cellFaceRow(c);
            for (int i = 0; i < faces.size(); ++i) {
                const int f = faces[i].index();
                BOOST_CHECK_MESSAGE(grid.faceCell(f, 0) == c || grid.faceCell(f, 1) == c,
                                    "cell " << c << " has face " << f << " that does not have the cell");
            }
        }
        for (int f = 0; f < grid.numFaces(); ++f) {
            for (const int side : {0, 1}) {
                const int c = grid.faceCell(f, side);
                if (c < 0) {
                    continue;
                }
                const auto faces = grid.cellFaceRow(c);
                bool found = false;
                for (int i = 0; i < faces.size(); ++i) {
                    found = found || faces[i].index() == f;
                }
                BOOST_CHECK_MESSAGE(found, "face " << f << " has cell " << c << " that does not have the face");
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(SinglePorosityGridIsNotSplit)
{
    const Opm::EclipseGrid eclGrid(makeDeck(/*dual=*/false));
    const auto grid = makeCpGrid(eclGrid);
    BOOST_CHECK_EQUAL(grid->numCells(), nx * ny * nz);
    const std::array<int, 3> size{nx, ny, nz};
    BOOST_CHECK(grid->logicalCartesianSize() == size);
}

namespace {

// Number of faces, of any kind, shared by two active cells given by Cartesian index.
int sharedFaces(const Dune::CpGrid& grid, const int cartesian1, const int cartesian2)
{
    int c1 = -1;
    int c2 = -1;
    for (int c = 0; c < grid.numCells(); ++c) {
        if (grid.globalCell()[c] == cartesian1) { c1 = c; }
        if (grid.globalCell()[c] == cartesian2) { c2 = c; }
    }
    if (c1 < 0 || c2 < 0) {
        return -1;
    }
    int count = 0;
    const auto faces = grid.cellFaceRow(c1);
    for (int i = 0; i < faces.size(); ++i) {
        const int face = faces[i].index();
        if (grid.faceCell(face, 0) == c2 || grid.faceCell(face, 1) == c2) {
            ++count;
        }
    }
    return count;
}

// Two columns, three matrix layers on top of three fracture layers. The middle
// layer of the first column is very thin, in the matrix geometry that both
// continua share.
Opm::Deck makeThinLayerDeck(const std::string& pinch, const std::string& actnum,
                            const std::map<int, double>& poroOverride = {})
{
    constexpr int nxp = 2;
    constexpr int nzp = 6;
    std::ostringstream coord;
    for (int j = 0; j <= 1; ++j) {
        for (int i = 0; i <= nxp; ++i) {
            coord << 100.0 * i << ' ' << 100.0 * j << " 2000 " << 100.0 * i << ' ' << 100.0 * j << " 2200\n";
        }
    }
    // Depths of the layer boundaries of the matrix, per column. The middle layer of the first
    // column is thinner than the PINCH threshold. The fracture layers are the matrix layers 60 m below.
    const std::array<std::array<double, 4>, 2> boundary{{{2000.0, 2010.0, 2010.2, 2020.0},
                                                         {2000.0, 2010.0, 2020.0, 2030.0}}};
    const auto top = [&boundary](const int k, const int column) {
        return boundary[column][k % 3] + 60.0 * (k / 3);
    };
    const auto bottom = [&boundary](const int k, const int column) {
        return boundary[column][k % 3 + 1] + 60.0 * (k / 3);
    };
    std::ostringstream zcorn;
    for (int k = 0; k < nzp; ++k) {
        for (int face = 0; face < 2; ++face) {
            for (int jj = 0; jj < 2; ++jj) {
                for (int i = 0; i < nxp; ++i) {
                    for (int ii = 0; ii < 2; ++ii) {
                        zcorn << (face == 0 ? top(k, i) : bottom(k, i)) << ' ';
                    }
                }
                zcorn << '\n';
            }
        }
    }
    const int cells = nxp * 1 * nzp;
    std::ostringstream poro;
    for (int c = 0; c < cells; ++c) {
        const auto it = poroOverride.find(c);
        poro << (it != poroOverride.end() ? it->second : 0.2) << ' ';
    }
    return Opm::Parser{}.parseString(
        "RUNSPEC\nOIL\nWATER\nDIMENS\n 2 1 6 /\nDUALPORO\nGRID\n" + pinch + "COORD\n" + coord.str()
        + "/\nZCORN\n" + zcorn.str() + "/\nACTNUM\n" + actnum + "/\nPORO\n " + poro.str()
        + "/\nPERMX\n " + std::to_string(cells) + "*1.0 /\nSIGMA\n 0.12 /\n");
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(MinpvRemovesCellsOfEachContinuum)
{
    // Matrix cell 3 and fracture cell 12+4 have too little pore volume.
    const auto deck = makeDeck(true, "MINPV\n 1e4 /\n", {{3, 0.001}, {16, 0.001}});
    Opm::EclipseState eclState(deck);
    const auto& eclGrid = eclState.getInputGrid();

    Dune::CpGrid grid;
    auto removed = grid.processEclipseFormat(&eclGrid, &eclState, false, false, false);
    std::sort(removed.begin(), removed.end());

    // The cells are reported with the Cartesian index of the dual-continuum grid.
    BOOST_CHECK((removed == std::vector<std::size_t>{3, 16}));
    BOOST_CHECK_EQUAL(grid.numCells(), nx * ny * nz - 2);
    for (const int cell : {3, 16}) {
        BOOST_CHECK(std::find(grid.globalCell().begin(), grid.globalCell().end(), cell)
                    == grid.globalCell().end());
    }

    // The twins of the removed cells stay. Their coupling is dropped with the cell.
    for (const int cell : {3 + 12, 16 - 12}) {
        BOOST_CHECK(std::find(grid.globalCell().begin(), grid.globalCell().end(), cell)
                    != grid.globalCell().end());
    }
    BOOST_CHECK_EQUAL(sharedFaces(grid, 0, 12), 1);
    BOOST_CHECK_EQUAL(sharedFaces(grid, 3, 15), -1);
    BOOST_CHECK_EQUAL(sharedFaces(grid, 4, 16), -1);
}

BOOST_AUTO_TEST_CASE(MinpvOfEachContinuumFollowsItsOwnPoreVolume)
{
    // Only the pore volume of the continuum's own cells counts: with a small
    // value in the matrix cell, the fracture twin is kept.
    const auto deck = makeDeck(true, "MINPV\n 1e4 /\n", {{3, 0.001}});
    Opm::EclipseState eclState(deck);
    Dune::CpGrid grid;
    auto removed = grid.processEclipseFormat(&eclState.getInputGrid(), &eclState, false, false, false);
    BOOST_CHECK((removed == std::vector<std::size_t>{3}));
    BOOST_CHECK_EQUAL(grid.numCells(), nx * ny * nz - 1);
}

BOOST_AUTO_TEST_CASE(PinchedOutLayerIsBridgedInEachContinuum)
{
    // The zero-thickness layer is inactive in both continua. With the default
    // PINCH options the cells above and below it are connected by a face.
    const std::string actnum = " 1 1  0 1  1 1   1 1  0 1  1 1 ";
    const auto deck = makeThinLayerDeck("PINCH\n 0.5 /\n", actnum);
    Opm::EclipseState eclState(deck);

    Dune::CpGrid grid;
    grid.processEclipseFormat(&eclState.getInputGrid(), &eclState, false, false, false);

    BOOST_CHECK_EQUAL(grid.numCells(), 12 - 2);
    // Column 0: layers 0 and 2 of the matrix (cells 0 and 4) and of the fracture (6 and 10).
    BOOST_CHECK_EQUAL(sharedFaces(grid, 0, 4), 1);
    BOOST_CHECK_EQUAL(sharedFaces(grid, 6, 10), 1);
    // The continua are not connected through the pinched layers.
    BOOST_CHECK_EQUAL(sharedFaces(grid, 4, 6), 0);
    BOOST_CHECK(eclState.getPinchNNC().empty());
}

BOOST_AUTO_TEST_CASE(PinchOutConnectionsOfBothContinuaUseTheDualGridNumbering)
{
    // Item 4 of PINCH set to ALL: the transmissibility across the pinched layer is
    // given to the simulator as a connection, one per continuum.
    const std::string actnum = " 1 1  0 1  1 1   1 1  0 1  1 1 ";
    const auto deck = makeThinLayerDeck("PINCH\n 0.5 GAP 1* ALL /\n", actnum);
    Opm::EclipseState eclState(deck);

    Dune::CpGrid grid;
    grid.processEclipseFormat(&eclState.getInputGrid(), &eclState, false, false, false);

    const auto& pinch = eclState.getPinchNNC();
    BOOST_REQUIRE_EQUAL(pinch.size(), std::size_t{2});
    BOOST_CHECK_EQUAL(pinch[0].cell1, std::size_t{0});
    BOOST_CHECK_EQUAL(pinch[0].cell2, std::size_t{4});
    BOOST_CHECK_EQUAL(pinch[1].cell1, std::size_t{6});
    BOOST_CHECK_EQUAL(pinch[1].cell2, std::size_t{10});
    BOOST_CHECK(pinch[0].trans > 0.0);
    // Same permeability and geometry in both continua.
    BOOST_CHECK_CLOSE(pinch[0].trans, pinch[1].trans, 1e-10);
}

BOOST_AUTO_TEST_CASE(ContinuumWithoutActiveCellsIsRejected)
{
    const auto deck = makeDeck();
    for (const bool matrixEmpty : {true, false}) {
        std::vector<int> active(nx * ny * nz, 1);
        const int half = nx * ny * nz / 2;
        for (int c = 0; c < half; ++c) {
            active[c + (matrixEmpty ? 0 : half)] = 0;
        }
        const Opm::EclipseGrid eclGrid(deck, active.data());
        Dune::CpGrid grid;
        BOOST_CHECK_THROW(grid.processEclipseFormat(&eclGrid, nullptr, false, false, false),
                          std::invalid_argument);
    }
}

BOOST_AUTO_TEST_CASE(EclipseGridOfProcessedGridKeepsTheGeometryOfRemovedMatrixCells)
{
    // The middle matrix cell of the second column (3) has too little pore volume, and is
    // removed and bridged. The processed grid then has its corners collapsed. The fracture
    // twin (9) is not removed, and has to keep the geometry of the input.
    const std::string actnum = " 1 1 1 1 1 1  1 1 1 1 1 1 ";
    const auto deck = makeThinLayerDeck("PINCH\n 0.5 /\nMINPV\n 1e4 /\n", actnum, {{3, 0.001}});
    Opm::EclipseState eclState(deck);
    const auto& inputGrid = eclState.getInputGrid();

    Dune::CpGrid grid;
    auto removed = grid.processEclipseFormat(&inputGrid, &eclState, false, false, false);
    std::sort(removed.begin(), removed.end());
    // The thin layer of the first column has too little pore volume as well.
    BOOST_CHECK((removed == std::vector<std::size_t>{2, 3, 8}));
    BOOST_CHECK(!grid.zcornData().empty());

    const auto output = Opm::UgGridHelpers::createEclipseGrid(grid, inputGrid);
    BOOST_CHECK(!output.cellActive(3));
    BOOST_REQUIRE(output.cellActive(9));
    BOOST_CHECK_CLOSE(output.getCellThickness(std::size_t{9}), 10.0, 1e-8);
    BOOST_CHECK_CLOSE(output.getCellDepth(std::size_t{9}), inputGrid.getCellDepth(std::size_t{3}), 1e-8);
}

BOOST_AUTO_TEST_CASE(NumericalAquifersAreRejected)
{
    const auto deck = Opm::Parser{}.parseString(
        "RUNSPEC\nOIL\nWATER\nDIMENS\n 2 1 2 /\nAQUDIMS\n 1 1 1* 1* 1* 1* 1* 1* /\nDUALPORO\nGRID\n"
        "DX\n 4*100 /\nDY\n 4*100 /\nDZ\n 4*10 /\nTOPS\n 2*2000 /\nPORO\n 4*0.2 /\nPERMX\n 4*1 /\n"
        "SIGMA\n 0.1 /\n"
        "AQUNUM\n 1 1 1 1 1.0e3 100 0.3 100 1* 1* 1* 1* /\n/\n"
        "AQUCON\n 1 2 2 1 1 1 1 'I-' 1* 1* 'NO' /\n/\n");
    Opm::EclipseState eclState(deck);
    Dune::CpGrid grid;
    BOOST_CHECK_THROW(grid.processEclipseFormat(&eclState.getInputGrid(), &eclState, false, false, false),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(EclipseStateNncGiveTheMatrixFractureCouplingFaces)
{
    // From an EclipseState, whose input NNCs hold the coupling of each twin pair.
    const auto deck = makeDeck();
    Opm::EclipseState eclState(deck);
    const auto& eclGrid = eclState.getInputGrid();
    const auto numNnc = eclState.getInputNNC().input().size();
    BOOST_REQUIRE_EQUAL(numNnc, eclGrid.getCartesianSize() / 2);

    Dune::CpGrid grid;
    grid.processEclipseFormat(&eclGrid, &eclState, false, false, false);
    const auto withoutNnc = makeCpGrid(eclGrid);

    BOOST_REQUIRE_EQUAL(grid.numCells(), static_cast<int>(eclGrid.getNumActive()));
    BOOST_CHECK_EQUAL(grid.numFaces(), withoutNnc->numFaces() + static_cast<int>(numNnc));

    // Every twin pair is connected by one NNC face, and these come first.
    std::unordered_map<int, int> cellOfGlobal;
    for (int c = 0; c < grid.numCells(); ++c) {
        cellOfGlobal[grid.globalCell()[c]] = c;
    }
    const int half = static_cast<int>(eclGrid.getCartesianSize() / 2);
    for (std::size_t f = 0; f < numNnc; ++f) {
        const int a = grid.faceCell(static_cast<int>(f), 0);
        const int b = grid.faceCell(static_cast<int>(f), 1);
        BOOST_REQUIRE(a >= 0 && b >= 0);
        BOOST_CHECK_EQUAL(grid.globalCell()[b], grid.globalCell()[a] + half);
        BOOST_CHECK_EQUAL(grid.numFaceVertices(static_cast<int>(f)), 0);
    }
    for (int g = 0; g < half; ++g) {
        BOOST_CHECK(cellOfGlobal.count(g) == 1 && cellOfGlobal.count(g + half) == 1);
    }
}
