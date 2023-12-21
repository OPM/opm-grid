//===========================================================================
//
// File: cuboidShape_test.cpp
//
// Created: Tue 05.12.2023 11:00:00
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
#include "config.h"

#define BOOST_TEST_MODULE CuboidShapeTest
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif

#include <dune/common/version.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/LookUpCellCentroid.hh>

#include <sstream>
#include <iostream>

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

void createEclGridCpGrid_and_checkCuboidShape(const std::string& deckString)
{
    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);

    Dune::CpGrid grid;
    Opm::EclipseGrid eclGrid(deck);

    grid.processEclipseFormat(&eclGrid, nullptr, false, false, false);

    const auto& leafGridView = grid.chooseData()[0]; // Leaf Grid coincides with level 0, in this case. 
    BOOST_CHECK_THROW(leafGridView->checkCuboidShape({0}), std::logic_error);
}

BOOST_AUTO_TEST_CASE(nonCuboidCell)
{/*
   Cell corners:                     COORD
   0 {0, 0, 0}                       line 1: corners 0 and 4
   1 {1, 0, 0}                       line 2: corners 1 and 5
   2 {0, 1, 0}                       line 3: corners 2 and 6
   3 {1, 1, 0}                       line 4: corners 3 and 7
   4 {0, 0, 1}
   5 {1, 0, 0}  coincides with 1
   6 {0, 1, 1}
   7 {1, 1, 0}  coincides with 3 */
    const std::string deckString =
        R"(RUNSPEC
        DIMENS
        1  1  1 /
        GRID
        COORD
        0 0 0  0 0 1
        1 0 0  1 0 0
        0 1 0  0 1 1
        1 1 0  1 1 0
        /
        ZCORN
        6*0
        2*1
        /
        ACTNUM
        1*1
        /
        )";

    createEclGridCpGrid_and_checkCuboidShape(deckString);
}
