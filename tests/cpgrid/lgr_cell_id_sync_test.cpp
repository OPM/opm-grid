//===========================================================================
//
// File: lgr_cell_id_sync_test.cpp
//
// Created: Wednesday 05.03.2025 08:47:00
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

#define BOOST_TEST_MODULE LgrCellIdSyncTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/ParentToChildrenCellGlobalIdHandle.hpp>
#include <opm/grid/common/CommunicationUtils.hpp>
#include <tests/cpgrid/LgrChecks.hpp>

#include <dune/common/version.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <array>
#include <stdexcept>
#include <string>
#include <map>
#include <vector>
#include <tuple>

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


void createTestGrid(Dune::CpGrid& grid)
{
    Opm::Parser parser;
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

    const auto deck = parser.parseString(deck_string);
    Opm::EclipseState ecl_state(deck);
    Opm::EclipseGrid eclipse_grid = ecl_state.getInputGrid();
    grid.processEclipseFormat(&eclipse_grid, &ecl_state, false, false, false);
}

void checkChildGlobalIdsTest(const Dune::CpGrid& grid)
{
    // first_child_ids [ parent cell global id ] =  first child cell global id from undistributed view.
    const std::unordered_map<int,int> first_child_ids = { {5, 66}, {6, 93}, {9, 120}, {10, 147} };
    // parent cell with global id equal to 5 has children cell ids 66, ..., 92 (= 66 + 26).
    // parent cell with global id equal to 6 has children cell ids 93, ..., 119 (= 119 + 26).
    // parent cell with global id equal to 9 has children cell ids 120, ..., 146 (= 120 + 26).
    // parent cell with global id equal to 10 has children cell ids 147, ..., 173 (= 147 + 26).

    for (const auto& element : elements(grid.leafGridView())) {
        if (!element.hasFather()) {
            continue;
        }
        const auto& parent_globalId = grid.globalIdSet().id(element.father());
        const auto& idx_in_parent = element.getIdxInParentCell();

        const auto& expected_elem_globalId = first_child_ids.at(parent_globalId) + idx_in_parent;
        const auto& actual_elem_globalId = grid.globalIdSet().id(element);

        BOOST_CHECK_EQUAL( expected_elem_globalId, actual_elem_globalId);
    }
}

BOOST_AUTO_TEST_CASE(cellIdSyncWhenLgrsAddedFirstInUndistributedViewThenInDistributedView)
{
    // Create the grid and add the LGRs in the global-view (non-distributed view)
    Dune::CpGrid grid;
    createTestGrid(grid);
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{3, 3, 3}},
                               /* startIJK_vec = */ {{1, 1, 0}},
                               /* endIJK_vec = */  {{3, 3, 1}},
                               /* lgr_name_vec = */  {"LGR1"});

    // Per parent cell, communicating the first child global id on the list is enough.
    // For the other siblings, first_child_global_idx + idx, with idx=1,..., children_list.size() -1.
    Opm::checkConsecutiveChildGlobalIdsPerParent(grid);

    // Note: parent cell global id from undistirbuted and distributed view coincide.
    //       CpGrid::syncCellIds() rewrite the cell ids of refined cells to coincide with the ones from
    //       the undistributed view.
    checkChildGlobalIdsTest(grid);

    const int maxLevelBeforeLoadBalance = grid.maxLevel();

    if (grid.comm().size()>1) {

        // Load balance the grid and add the LGRs in the distributed view.
        grid.loadBalance(/*overlapLayers*/ 1,
                         /*partitionMethod*/ Dune::PartitionMethod::zoltanGoG,
                         /*imbalanceTol*/ 1.1,
                         /*level*/ 0);
        grid.addLgrsUpdateLeafView(/*cells_per_dim_vec*/ {{3, 3, 3}},
                                   /*startIJK_vec*/ {{1, 1, 0}},
                                   /*endIJK_vec*/ {{3, 3, 1}},
                                   /*lgr_name_vec*/ {"LGR1"});


        const int maxLevelAfterLoadBalance = grid.maxLevel();
        BOOST_CHECK_EQUAL( maxLevelBeforeLoadBalance, maxLevelAfterLoadBalance);

        Opm::checkConsecutiveChildGlobalIdsPerParent(grid);

        grid.syncDistributedGlobalCellIds();
        const auto& data = grid.currentData();
        Opm::checkCellGlobalIdUniquenessForInteriorCells(grid, data);

        checkChildGlobalIdsTest(grid);
    }
}

BOOST_AUTO_TEST_CASE(cellIdSyncWhenLgrsAddedFirstInDistributedViewThenInUndistributedView)
{
    // Create the grid and add the LGRs in the global-view (non-distributed view)
    Dune::CpGrid grid;
    createTestGrid(grid);

    if (grid.comm().size()>1) {

        // Load balance the grid and add the LGRs in the distributed view.
        grid.loadBalance(/*overlapLayers*/ 1,
                         /*partitionMethod*/ Dune::PartitionMethod::zoltanGoG,
                         /*imbalanceTol*/ 1.1,
                         /*level*/ 0);
        grid.addLgrsUpdateLeafView(/*cells_per_dim_vec*/ {{3, 3, 3}},
                                   /*startIJK_vec*/ {{1, 1, 0}},
                                   /*endIJK_vec*/ {{3, 3, 1}},
                                   /*lgr_name_vec*/ {"LGR1"});
    }

    grid.switchToGlobalView();
    grid.addLgrsUpdateLeafView(/*cells_per_dim_vec*/ {{3, 3, 3}},
                               /*startIJK_vec*/ {{1, 1, 0}},
                               /*endIJK_vec*/ {{3, 3, 1}},
                               /*lgr_name_vec*/ {"LGR1"});
    Opm::checkConsecutiveChildGlobalIdsPerParent(grid);
    checkChildGlobalIdsTest(grid);

    if (grid.comm().size()>1) {
        grid.switchToDistributedView();

        Opm::checkConsecutiveChildGlobalIdsPerParent(grid);

        grid.syncDistributedGlobalCellIds();
        const auto& data = grid.currentData();
        Opm::checkCellGlobalIdUniquenessForInteriorCells(grid, data);

        checkChildGlobalIdsTest(grid);
    }
}
