//===========================================================================
//
// File: distribute_level_zero_from_grid_with_lgrs_test.cpp
//
// Created: Thursday 17.02.2025 10:51:00
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

#define BOOST_TEST_MODULE DistributeLevelZeroFromGridWithLgrsTests
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif


#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/CommunicationUtils.hpp>

#include <dune/common/version.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

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

// Create a grid and add one test LGR with dimension 6x6x3.
Dune::CpGrid
createGridAndAddTestLgr(const std::string& deck_string)
{
    Dune::CpGrid grid;
    // Create the starting grid (before adding LGRs)
    Opm::Parser parser;
    const auto deck = parser.parseString(deck_string);
    Opm::EclipseState ecl_state(deck);
    Opm::EclipseGrid eclipse_grid = ecl_state.getInputGrid();

    grid.processEclipseFormat(&eclipse_grid, &ecl_state, false, false, false);


    // Add LGR1 and update grid view
    const std::vector<std::array<int, 3>> cells_per_dim_vec
        = {{3, 3, 3}}; // 3x3x3 child cells in x-,y-, and z-direction per ACTIVE parent cell
    const std::vector<std::array<int, 3>> startIJK_vec
        = {{1, 1, 0}}; // starts at (1,1,0) in coarse grid - equivalent to (I1-1, J1-1, K1-1) from its CARFIN block
    const std::vector<std::array<int, 3>> endIJK_vec
        = {{3, 3, 1}}; // ends at (3,3,1) in coarse grid - equivalent to (I2, J2, K2) from its CARFIN block
    const std::vector<std::string> lgr_name_vec = {"LGR1"};

    grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
    return grid;
}

BOOST_AUTO_TEST_CASE(fullActiveParentCellsBlock)
{
    const std::string deck_string = R"(
RUNSPEC
DIMENS
  4 4 1 /
GRID
CARFIN
-- NAME I1-I2 J1-J2 K1-K2 NX NY NZ
'LGR1'  2  3  2  3  1  1  6  6  3/
ENDFIN
DX
  16*1000 /
DY
	16*1000 /
DZ
	16*20 /
TOPS
	16*8325 /
ACTNUM
        16*1 /
PORO
  16*0.15 /
PERMX
  16*1 /
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

    auto grid = createGridAndAddTestLgr(deck_string);

    // loadBalance( overlapSize, partitionMethod, imbalanceTolerance, level grid to be distributed )
    grid.loadBalance(1, Dune::PartitionMethod::zoltanGoG, 1.1, 0);
    // grid.loadBalance();
}
