//===========================================================================
//
// File: inactiveCells_lgr_test.cpp
//
// Created:  July 04 2024 11:06:00
//
// Author(s): Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2024 Equinor ASA.

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

#define BOOST_TEST_MODULE InactiveCellsLgr

#include <boost/test/unit_test.hpp>

#include <dune/common/version.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/grid/CpGrid.hpp>

#include <iostream>
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

void testInactiveCellsLgrs(const std::string& deckString,
                           const std::vector<std::array<int,3>>& cells_per_dim_vec,
                           const std::vector<std::array<int,3>>& startIJK_vec,
                           const std::vector<std::array<int,3>>& endIJK_vec,
                           const std::vector<std::string>& lgr_name_vec)
{
    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);
    Opm::EclipseState es(deck);

    Dune::CpGrid grid;
    grid.processEclipseFormat(&es.getInputGrid(), &es, false, false, false);

    try {
        grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
    }
    catch (const std::exception& e) {
        std::cout << e.what() << "'\n";
    }
}

BOOST_GLOBAL_FIXTURE(Fixture);


BOOST_AUTO_TEST_CASE(inactiveCells_in_at_least_one_lgr)
{
    const std::string deckString =
        R"( RUNSPEC
        DIMENS
        1  1  5 /
        GRID
        COORD
        0 0 0
        0 0 1
        1 0 0
        1 0 1
        0 1 0
        0 1 1
        1 1 0
        1 1 1
        /
        ZCORN
        4*0
        8*1
        8*2
        8*3
        8*4
        4*5
        /
        ACTNUM
        0
        1
        1
        1
        1
        /
        PORO
        5*0.15
        /)";

    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}};
    // LGR1 cell indices: {0,1}, 0 INACTIVE CELL, 1 ACTIVE CELL
    // LGR2 cell indices: {4} ACTIVE CELL
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,4}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,2}, {1,1,5}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    testInactiveCellsLgrs(deckString, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
}

BOOST_AUTO_TEST_CASE(inactiveCells_outside_lgr)
{
    const std::string deckString =
        R"( RUNSPEC
        DIMENS
        1  1  5 /
        GRID
        COORD
        0 0 0
        0 0 1
        1 0 0
        1 0 1
        0 1 0
        0 1 1
        1 1 0
        1 1 1
        /
        ZCORN
        4*0
        8*1
        8*2
        8*3
        8*4
        4*5
        /
        ACTNUM
        1
        1
        1
        0
        1
        /
        PORO
        5*0.15
        /)";

    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}, {4,4,4}};
    // LGR1 cell indices: {0}
    // LGR2 cell indices: {2}
    // LGR3 cell indices: {4}
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,2}, {0,0,4}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,1}, {1,1,3}, {1,1,5}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3"};
    testInactiveCellsLgrs(deckString, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
}

BOOST_AUTO_TEST_CASE(inactiveCells_in_lgrs)
{

    const std::string deckString =
        R"( RUNSPEC
        DIMENS
        1  1  5 /
        GRID
        COORD
        0 0 0
        0 0 1
        1 0 0
        1 0 1
        0 1 0
        0 1 1
        1 1 0
        1 1 1
        /
        ZCORN
        4*0
        8*1
        8*2
        8*3
        8*4
        4*5
        /
        ACTNUM
        0
        1
        1
        1
        0
        /
        PORO
        5*0.15
        /)";

    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,3}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,2}, {1,1,5}};
    // LGR1 cell indices = {0,1}, LGR2 cell indices = {3,4}.
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    testInactiveCellsLgrs(deckString, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
}



BOOST_AUTO_TEST_CASE(inactiveCells_in_and_outside_lgrs)
{

    const std::string deckString =
        R"( RUNSPEC
        DIMENS
        1  1  5 /
        GRID
        COORD
        0 0 0
        0 0 1
        1 0 0
        1 0 1
        0 1 0
        0 1 1
        1 1 0
        1 1 1
        /
        ZCORN
        4*0
        8*1
        8*2
        8*3
        8*4
        4*5
        /
        ACTNUM
        1
        0
        0
        1
        0
        /
        PORO
        5*0.15
        /)";

    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,3}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,2}, {1,1,5}};
    // LGR1 cell indices = {0,1}, LGR2 cell indices = {3,4}.
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    testInactiveCellsLgrs(deckString, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
}

BOOST_AUTO_TEST_CASE(inactiveCells_only_outside_lgrs)
{

    const std::string deckString =
        R"( RUNSPEC
        DIMENS
        1  1  5 /
        GRID
        COORD
        0 0 0
        0 0 1
        1 0 0
        1 0 1
        0 1 0
        0 1 1
        1 1 0
        1 1 1
        /
        ZCORN
        4*0
        8*1
        8*2
        8*3
        8*4
        4*5
        /
        ACTNUM
        1
        1
        0
        0
        1
        /
        PORO
        5*0.15
        /)";

    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,4}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,2}, {1,1,5}};
    // LGR1 cell indices = {0,1}, LGR2 cell indices = {4}.
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    testInactiveCellsLgrs(deckString, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
}


