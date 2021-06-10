/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2021 Equinor ASA.

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

#include <config.h>

#define NVERBOSE

#define BOOST_TEST_MODULE CpGridPINCH

#include <boost/test/unit_test.hpp>

#include <dune/common/version.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/grid/CpGrid.hpp>
#include <vector>
#include <utility>

struct Fixture
{
    Fixture()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(m_argc, m_argv);
        Opm::OpmLog::setupSimpleDefaultLogging();
    }

    static int rank()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        return Dune::MPIHelper::instance(m_argc, m_argv).rank();
    }
};

BOOST_GLOBAL_FIXTURE(Fixture);

BOOST_AUTO_TEST_CASE(NoPinchWithThickness)
{
    BOOST_TEST_MESSAGE("Testing that pinch thickness does not create NNCs");
    // Specify a stack of three cells where the middle
    // cells is less thick than the PINCH thickness
    // Should still result in a grid with 3 cells.
    const char *deckString =
        "RUNSPEC\n"
        "DIMENS\n"
        "1  1  3 /"
        "GRID\n"
        "COORD\n"
        "0 0 0  0 0 3\n"
        "1 0 0  1 0 3\n"
        "0 1 0  0 1 3\n"
        "1 1 0  1 1 3\n"
        "/\n"
        "ZCORN\n"
        "4*0\n"
        "8*1\n"
        "8*1.1\n"
        "4*2.1\n"
        "/\n"
        "PINCH\n"
        "0.02   NOGAP   1*   1*\n"
        "/\n";

    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);

    Dune::CpGrid grid;
    Opm::EclipseGrid ecl_grid(deck);

    grid.processEclipseFormat(&ecl_grid, nullptr, false, false, false);

    if(Fixture::rank() == 0)
    {
        BOOST_CHECK(grid.size(0)==3);
    }
}


BOOST_AUTO_TEST_CASE(PinchZeroVolumeNoBarrier)
{
    BOOST_TEST_MESSAGE("Testing that cells with zero volume present no barriers"
                       " with PINCH.");
    // Specify a stack of five cells where the middle
    // where two will present a barrier because of zero
    // value
    // cells is less thick than the PINCH thickness
    // Should still result in a grid with 3 cells.
    const char *deckString =
        "RUNSPEC\n"
        "DIMENS\n"
        "1  1  5 /"
        "GRID\n"
        "COORD\n"
        "0 0 0  0 0 3\n"
        "1 0 0  1 0 3\n"
        "0 1 0  0 1 3\n"
        "1 1 0  1 1 3\n"
        "/\n"
        "ZCORN\n"
        "4*0\n"
        "8*1\n"
        "8*1\n"
        "8*1.1\n"
        "8*1.1\n"
        "4*2.1\n"
        "/\n"
        "PINCH\n"
        "0.01   NOGAP   1*   1*\n"
        "/\n";

    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);

    Dune::CpGrid grid;
    Opm::EclipseGrid ecl_grid(deck);

    grid.processEclipseFormat(&ecl_grid, nullptr, false, false, false);

    if(Fixture::rank() == 0)
    {
        BOOST_CHECK(grid.size(0)==3);
    }

    // Check that cells 2 and 4 present no barriers
    // i.d. each cell has neighbors
    const auto& gridView = grid.leafGridView();
    for(const auto& element : elements(gridView))
    {
        std::size_t neighbors = 0;
        for(auto it = gridView.ibegin(element), endIt = gridView.iend(element);
            it != endIt; ++it)
        {
            neighbors += it.neighbor();
        }
        BOOST_CHECK(neighbors > 0);
    }
}


BOOST_AUTO_TEST_CASE(NoPinchZeroVolumeBarrier)
{
    BOOST_TEST_MESSAGE("Testing that cells with zero volume present barriers"
                       " without PINCH.");
    // Specify a stack of five cells where the middle
    // where two will present a barrier because of zero
    // value
    // cells is less thick than the PINCH thickness
    // Should still result in a grid with 3 cells.
    const char *deckString =
        "RUNSPEC\n"
        "DIMENS\n"
        "1  1  5 /"
        "GRID\n"
        "COORD\n"
        "0 0 0  0 0 3\n"
        "1 0 0  1 0 3\n"
        "0 1 0  0 1 3\n"
        "1 1 0  1 1 3\n"
        "/\n"
        "ZCORN\n"
        "4*0\n"
        "8*1\n"
        "8*1\n"
        "8*1.1\n"
        "8*1.1\n"
        "4*2.1\n"
        "/\n";

    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);

    Dune::CpGrid grid;
    Opm::EclipseGrid ecl_grid(deck);

    grid.processEclipseFormat(&ecl_grid, nullptr, false, false, false);

    if(Fixture::rank() == 0)
    {
        BOOST_CHECK(grid.size(0)==3);
    }

    // Check that cells 2 and 4 present no barriers
    // i.d. each cell has neighbors
    const auto& gridView = grid.leafGridView();
    for(const auto& element : elements(gridView))
    {
        for(auto it = gridView.ibegin(element), endIt = gridView.iend(element);
            it != endIt; ++it)
        {
            BOOST_CHECK(!it.neighbor());
        }
    }
}
