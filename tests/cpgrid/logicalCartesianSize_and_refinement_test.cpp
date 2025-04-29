//===========================================================================
//
// File: logicalCartesianSize_and_refinement_test.cpp
//
// Created: Tuesday 29.04.2025 09:50:00
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

#define BOOST_TEST_MODULE LogicalCartesianSizeTests
#include <boost/test/unit_test.hpp>

#include <tests/cpgrid/LgrChecks.hpp>


#include <array>
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


void areEqual(const std::array<int,3>& expected_logicalCartesianSize,
              const std::array<int,3>& actual_logicalCartesianSize)
{
    BOOST_CHECK_EQUAL(expected_logicalCartesianSize[0], actual_logicalCartesianSize[0]);
    BOOST_CHECK_EQUAL(expected_logicalCartesianSize[1], actual_logicalCartesianSize[1]);
    BOOST_CHECK_EQUAL(expected_logicalCartesianSize[2], actual_logicalCartesianSize[2]);
}

BOOST_AUTO_TEST_CASE(lgrLogicalCartesianSize_afterAddLgrsUpdateLeafView_makesSense)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    grid.addLgrsUpdateLeafView( /* cells_per_dim = */ {{3,3,3}, {3,3,3}},
                                /* startIJK_vec = */ {{0,0,1}, {2,2,2}},
                                /* endIJK_vec = */ {{3,2,2}, {4,3,3}},
                                /* lgr_name_vec = */ {"LGR1", "LGR2"});

    // Block shaped parent cells of LGR1 dimensions (3-0)x(2-0)x(2-1). Number of subdivisions per cell, per direction {3,3,3}.
    areEqual( /* expected_logicalCartisianSize = */  {9,6,3},  // LGR1 dimensions {(3-0)*3, (2-0)*3, (2-1)*3}.
              /* LGR1 logicalCartesianSize = */ grid.currentData()[1]->logicalCartesianSize());

    // Block shaped parent cells of LGR2 dimensions (4-2)x(3-2)x(3-2). Number of subdivisions per cell, per direction {3,3,3}.
    areEqual( /* expected_logicalCartisianSize = */ {6,3,3}, // LGR2 dimensions {(4-2)*3, (3-2)*3, (3-2)*3}.
              /* LGR2 logicalCartesianSize = */ grid.currentData()[2]->logicalCartesianSize());

}


BOOST_AUTO_TEST_CASE(gridLogicalCartesianSize_afterStrictLocalRefinementWith_addLgrsUpdateLeafView_isACopyOfLevelZeroLogicalCartesianSize)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    grid.addLgrsUpdateLeafView(/* cells_per_dim = */ {{3,3,3}, {3,3,3}},
                               /* startIJK_vec = */ {{0,0,1}, {2,2,2}},
                               /* endIJK_vec = */ {{3,2,2}, {4,3,3}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2"});

    areEqual(/* grid dimensions before refinement = */ {4,3,3},
             /* level 0 logicalCartesianSize = */ grid.currentData().front()->logicalCartesianSize());

    areEqual(/* level 0 logicalCartesianSize = */ grid.currentData().front()->logicalCartesianSize(),
             grid.logicalCartesianSize());
}

BOOST_AUTO_TEST_CASE(gridLogicalCartesianSize_afterHiddenGlobalRefinementWith_addLgrsUpdateLeafView_makesSense)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    grid.addLgrsUpdateLeafView(/* cells_per_dim = */ {{3,3,3}},
                               /* startIJK_vec = */ {{0,0,0}},
                               /* endIJK_vec = */ {{4,3,3}},
                               /* lgr_name_vec = */ {"LGR1"});

    areEqual(/* grid dimensions before refinement = */ {4,3,3},
             /* level 0 logicalCartesianSize = */ grid.currentData().front()->logicalCartesianSize());

    // Block shaped parent cells of LGR1 is the entire level zero grid, dimensions (4-0)x(3-0)x(3-1).
    // Number of subdivisions per cell, per direction {3,3,3}.
    areEqual(/* expected logicalCartesianSize = */ {12, 9, 9},  // LGR1 dimensions {4*3, 3*3, 3*3}.
             /* LGR1 logicalCartesianSize = */ grid.currentData()[1]->logicalCartesianSize());

    areEqual(/* expected logicalCartesianSize = */ {12, 9, 9},  // LGR1 dimensions {4*3, 3*3, 3*3}.
             grid.logicalCartesianSize());
}


BOOST_AUTO_TEST_CASE(lgrAndGridLogicalCartesianSize_afterStrictLocalRefinementWith_adapt_areACopyOfLevelZeroLogicalCartesianSize)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    std::vector<int> markedCells = {17,18,21,22};
    Opm::adaptGrid(grid, markedCells); // Default subdivisions per cell 2x2x2 in x-,y-, and z-direction.

    areEqual(/* grid dimensions before refinement = */ {4,3,3},
             /* level 0 logicalCartesianSize = */ grid.currentData().front()->logicalCartesianSize());

    // Even though the marked cells form a 2x2x1 block, the logicalCartesianSize of LGR1 is NOT {4,4,2}.
    areEqual(/* level 0 logicalCartesianSize = */ grid.currentData().front()->logicalCartesianSize(), // {4,3,3}
             /* LGR1 logicalCartesianSize = */  grid.currentData()[1]->logicalCartesianSize());

    areEqual(/* level 0 logicalCartesianSize = */ grid.currentData().front()->logicalCartesianSize(),
             grid.logicalCartesianSize());
}

BOOST_AUTO_TEST_CASE(lgrAndGridLogicalCartesianSize_afterHiddenGlobalRefinementWith_adapt_makeSense)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    std::vector<int> markedCells(36); // 36 = 4x3x3
    std::iota(markedCells.begin(), markedCells.end(), 0);
    Opm::adaptGrid(grid, markedCells); // Default subdivisions per cell 2x2x2 in x-,y-, and z-direction.

    areEqual(/* grid dimensions before refinement = */ {4,3,3},
             /* level 0 logicalCartesianSize = */ grid.currentData().front()->logicalCartesianSize());

    areEqual(/* expected logicalCartesianSize = */ {4*2, 3*2, 3*2},
             /* LGR1 logicalCartesianSize = */  grid.currentData()[1]->logicalCartesianSize());

    areEqual(/* expected logicalCartesianSize = */ {4*2, 3*2, 3*2},
             grid.logicalCartesianSize());
}

BOOST_AUTO_TEST_CASE(lgrAndGridLogicalCartesianSize_after_globalRefine_makeSense)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    grid.globalRefine(1); // Default subdivisions per cell 2x2x2 in x-,y-, and z-direction.

    areEqual(/* expected logicalCartesianSize = */ {4*2, 3*2, 3*2}, grid.logicalCartesianSize());
    // The refined level grid is a "copy" of the leaf grid view, if globalRefine has been invoked.
    // TODO: remove the refined level grid in this case.
    areEqual(/* expected logicalCartesianSize = */ {4*2, 3*2, 3*2}, grid.currentData()[1]->logicalCartesianSize());
}
