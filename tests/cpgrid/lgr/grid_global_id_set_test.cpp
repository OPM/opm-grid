//===========================================================================
//
// File: grid_global_id_set_test.cpp
//
// Created: Friday 28.02.2025 11:09:00
//
// Author(s): Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================
/*
  Copyright 2025 Equinor ASA.

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

#define BOOST_TEST_MODULE GridGlobalIdSetTests
#include <boost/test/unit_test.hpp>

#include <dune/common/version.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/grid/CpGrid.hpp>

#include <array>
#include <stdexcept>
#include <string>
#include <vector>

struct Fixture {
    Fixture()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(m_argc, m_argv);
        Opm::OpmLog::setupSimpleDefaultLogging();
    }
};

BOOST_GLOBAL_FIXTURE(Fixture);

// Create a test grid
Dune::CpGrid createTestGrid()
{
    Dune::CpGrid grid;
    // Create the starting grid (before adding LGRs)
    Opm::Parser parser;
    const std::string deck_string = R"(
RUNSPEC
DIMENS
  3 3 1 /
GRID
CARFIN
-- NAME I1-I2 J1-J2 K1-K2 NX NY NZ
'LGR1'  2  3  2  3  1  1  6  6  3/
ENDFIN
DX
  9*1000 /
DY
	9*1000 /
DZ
	9*20 /
TOPS
	9*8325 /
 ACTNUM
        1 1 1
        1 1 1
        1 1 1
        /
PORO
  9*0.15 /
PERMX
  9*1 /
COPY
  PERMX PERMZ /
  PERMX PERMY /
/
EDIT
OIL
GAS
TITLE
The title
START
16 JUN 1988 /
PROPS
REGIONS
SOLUTION
SCHEDULE
)";
    const auto deck = parser.parseString(deck_string);
    Opm::EclipseState ecl_state(deck);
    Opm::EclipseGrid eclipse_grid = ecl_state.getInputGrid();
    grid.processEclipseFormat(&eclipse_grid, &ecl_state, false, false, false);
    return grid;
}


BOOST_AUTO_TEST_CASE(allLevelCellsHaveIdOnGridGlobalIdSet)
{
    auto grid = createTestGrid();
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{3, 3, 3}},
                               /* startIJK_vec = */ {{1, 1, 0}},
                               /* endIJK_vec = */ {{3, 3, 1}},
                               /* lgr_name_vec = */ {"LGR1"});
    const auto& gridGlobalIdSet = grid.globalIdSet();

    const int maxLevel = grid.maxLevel();
    for (int level = 0; level <= maxLevel; ++level) {
        const auto& levelGlobalIdSet = grid.currentData()[level]->globalIdSet();
        const auto& levelElements = Dune::elements(grid.levelGridView(level));
        for (const auto& element : levelElements) {
            BOOST_CHECK_EQUAL( gridGlobalIdSet.id(element), levelGlobalIdSet.id(element));
        }
    }
}

BOOST_AUTO_TEST_CASE(leafGlobalIdCellAndGridGlobalIdsCoincide)
{
    auto grid = createTestGrid();
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{3, 3, 3}},
                               /* startIJK_vec = */ {{1, 1, 0}},
                               /* endIJK_vec = */ {{3, 3, 1}},
                               /* lgr_name_vec = */ {"LGR1"});
    const auto& gridGlobalIdSet = grid.globalIdSet();

    const auto& leafGlobalIdSet = grid.currentData().back()->globalIdSet();
    const auto& leafElements = Dune::elements(grid.leafGridView());
    for (const auto& element : leafElements) {
        BOOST_CHECK_EQUAL( gridGlobalIdSet.id(element), leafGlobalIdSet.id(element));
    }
}

BOOST_AUTO_TEST_CASE(allLevelCellsHaveIdOnGridGlobalIdSetWithDistributedLevelZero)
{
    auto grid = createTestGrid();
    if (grid.comm().size()>1) {
        grid.loadBalance();
        grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{3, 3, 3}},
                                   /* startIJK_vec = */ {{1, 1, 0}},
                                   /* endIJK_vec = */ {{3, 3, 1}},
                                   /* lgr_name_vec = */ {"LGR1"});
        const auto& gridGlobalIdSet = grid.globalIdSet();

        const int maxLevel = grid.maxLevel();
        for (int level = 0; level <= maxLevel; ++level) {
            const auto& levelGlobalIdSet = grid.currentData()[level]->globalIdSet();
            const auto& levelElements = Dune::elements(grid.levelGridView(level));
            for (const auto& element : levelElements) {
                BOOST_CHECK_EQUAL( gridGlobalIdSet.id(element), levelGlobalIdSet.id(element));
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(leafGlobalIdCellAndGridGlobalIdsCoincideWithDistributedLevelZero)
{
    auto grid = createTestGrid();
    if (grid.comm().size()>1) {
        grid.loadBalance();
        grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{3, 3, 3}},
                                   /* startIJK_vec = */ {{1, 1, 0}},
                                   /* endIJK_vec = */ {{3, 3, 1}},
                                   /* lgr_name_vec = */ {"LGR1"});
        const auto& gridGlobalIdSet = grid.globalIdSet();

        const auto& leafGlobalIdSet = grid.currentData().back()->globalIdSet();
        const auto& leafElements = Dune::elements(grid.leafGridView());
        for (const auto& element : leafElements) {
            BOOST_CHECK_EQUAL( gridGlobalIdSet.id(element), leafGlobalIdSet.id(element));
        }
    }
}

