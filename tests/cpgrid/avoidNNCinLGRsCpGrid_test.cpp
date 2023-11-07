//===========================================================================
//
// File: avoidNNCinLGRsCpGrid_test.cpp
//
// Created:  Nov 09 2023 15:25:00
//
// Author(s): Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2023 Equinor ASA.

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

#define BOOST_TEST_MODULE CpGridWithNNC

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

void testCase(const std::string& deckString,
              const Opm::NNC& nnc,
              const std::vector<std::array<int,3>>& cells_per_dim_vec,
              const std::vector<std::array<int,3>>& startIJK_vec,
              const std::vector<std::array<int,3>>& endIJK_vec,
              const std::vector<std::string>& lgr_name_vec,
              bool hasNNC)
{

    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);
    Opm::EclipseState es(deck);
    es.appendInputNNC(nnc.input());

    Dune::CpGrid grid;
    grid.processEclipseFormat(&es.getInputGrid(), &es, false, false, false);
    
    if(hasNNC){
        BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec), std::logic_error);
    }
    else{
        grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
    }
}

BOOST_GLOBAL_FIXTURE(Fixture);

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
        5*1
        /
        PORO
        5*0.15
        /)";

BOOST_AUTO_TEST_CASE(NNCatAnLgr)
{
    Opm::NNC nnc;
    nnc.addNNC(2, 4, 1.0); // connect cell 2 (does not belong to any LGR) and cell 4 (belongs to LGR2)
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,4}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,1,1}, {1,1,5}};
    // LGR1 cell indices = {0,1}, LGR2 cell indices = {4}.
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    testCase(deckString, nnc,  cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec, true);
}

BOOST_AUTO_TEST_CASE(NNCAtSeveralLgrs)
{
    Opm::NNC nnc;
    nnc.addNNC(0, 1, 1.0); // connect cell 0 and cell 1 (both belong to LGR1)
    nnc.addNNC(2, 4, 1.0); // connect cell 2 (belongs to LGR2) and cell 4 (belongs to LGR3)
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}, {4,4,4}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,2}, {0,0,4}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,1}, {1,1,3}, {1,1,5}};
    // LGR1 cell indices = {0}, LGR2 cell indices = {2}, LGR3 cell indices = {4}.
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3"};
    testCase(deckString, nnc,  cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec, true);
}

BOOST_AUTO_TEST_CASE(NNCoutsideLgrs)
{
    Opm::NNC nnc;
    nnc.addNNC(1, 3, 1.0); // connect cell 1 and cell 3 (both do NOT belong to any LGR)
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}, {4,4,4}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,2}, {0,0,4}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,1}, {1,1,3}, {1,1,5}};
    // LGR1 cell indices = {0}, LGR2 cell indices = {2}, LGR3 cell indices = {4}.
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3"};
    testCase(deckString, nnc,  cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec, false);
}

