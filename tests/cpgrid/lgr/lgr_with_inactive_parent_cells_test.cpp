/*
  Copyright 2024, 2025 Equinor ASA.

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

#define BOOST_TEST_MODULE LgrWithInactiveParentCells

#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <tests/cpgrid/lgr/LgrChecks.hpp>

#include <array>
#include <string>
#include <vector>

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

void createTestGridWithLgrsSerial(Dune::CpGrid& grid,
                                  const std::string& deckString)
{
    Opm::createGridAndAddLgrs(grid,
                              deckString,
                              /* cells_per_dim_vec = */  {{2,2,2}, {3,3,3}},
                              /* startIJK_vec = */   {{0,0,0}, {2,3,1}},
                              /* endIJK_vec = */{{3,3,1}, {4,5,2}},
                              /* lgr_name_vec = */ {"LGR1", "LGR2"});
}

// Distribute level-zero grid and add LGRs to its distributed view,
// of the test CpGrid with existing LGRs in the global view.
void createTestGridWithLgrsParallel(Dune::CpGrid& grid)
{
    grid.loadBalance(/*overlapLayers*/ 1,
                     /*partitionMethod*/ Dune::PartitionMethod::zoltanGoG,
                     /*imbalanceTol*/ 1.1,
                     /*level*/ 0);

    grid.addLgrsUpdateLeafView( /* cells_per_dim_vec = */  {{2,2,2}, {3,3,3}},
                                /* startIJK_vec = */   {{0,0,0}, {2,3,1}},
                                /* endIJK_vec = */{{3,3,1}, {4,5,2}},
                                /* lgr_name_vec = */ {"LGR1", "LGR2"});
}

void checkGridInactiveCellsCount(const Dune::CpGrid& grid,
                                 int expected_maxLevel,
                                 const std::vector<int>& expected_global_cells)
{
    BOOST_CHECK_EQUAL( grid.maxLevel(), expected_maxLevel);
    Opm::checkGlobalActiveCellsCountInGridWithLgrs(grid, expected_global_cells);

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{2,2,2}, {3,3,3}},
                           /* lgr_name_vec = */  {"LGR1", "LGR2"});
}

// This test reuses the same deck, grid, and LGRs across cases.
// Only the ACTNUM block (active/inactive parent cells) differs.
// We avoid repetition and test serial/parallel behavior by distributing
// the level-zero grid and adding LGRs to its distributed view,
// even when LGRs exist in the global CpGrid view.
BOOST_AUTO_TEST_CASE(refinementDoesNotOccurIfAllParentCellsAreInactive)
{

    const std::string deckString =
        R"( RUNSPEC
  DIMENS
 -- NX NY NZ cells per x-,y-, and z-direction
     4 5 2 /
  GRID
  COORD -- grid corner coordinates (bounding box), 6*(NX +1)*(NY +1) values
  0 0 0   0 0 1 -- bottom of the pillar 0 | top of the pillar 0
  1 0 0   1 0 1 -- bottom of the pillar 1 | top of the pillar 1 (...)
  2 0 0   2 0 1
  3 0 0   3 0 1
  4 0 0   4 0 1
  0 1 0   0 1 1
  1 1 0   1 1 1
  2 1 0   2 1 1
  3 1 0   3 1 1
  4 1 0   4 1 1
  0 2 0   0 2 1
  1 2 0   1 2 1
  2 2 0   2 2 1
  3 2 0   3 2 1
  4 2 0   4 2 1
  0 3 0   0 3 1
  1 3 0   1 3 1
  2 3 0   2 3 1
  3 3 0   3 3 1
  4 3 0   4 3 1
  0 4 0   0 4 1
  1 4 0   1 4 1
  2 4 0   2 4 1
  3 4 0   3 4 1
  4 4 0   4 4 1
  0 5 0   0 5 1
  1 5 0   1 5 1
  2 5 0   2 5 1
  3 5 0   3 5 1
  4 5 0   4 5 1 -- bottom of the pillar 29 | top of the pillar 29
  /
  ZCORN
-- to easily deduce total active cells, no pinch-outs (consistent values to avoid collapsed cells, flat layers)
  80*0  -- top layer    k = 0
  80*1  -- bottom layer k = 0
  80*1  -- top layer    k = 1
  80*2  -- bottom layer k = 1
  /
  ACTNUM
-- i = 0 1 2 3
       0 0 0 1 -- layer k = 0     j = 0
       0 0 0 1 --                 j = 1
       0 0 0 1 --                 j = 2
       1 1 1 1 --                 j = 3
       1 1 1 0 --                 j = 4
       1 1 1 1 -- layer k = 1     j = 0
       0 1 1 1 --                 j = 1
       1 1 1 1 --                 j = 2
       1 1 0 0 --                 j = 3
       0 0 0 0 --                 j = 4
  /
  PORO
  40*0.15
  /)";

    Dune::CpGrid grid;
    createTestGridWithLgrsSerial(grid, deckString);
    // LGR1: 9 inactive parent cells.      |   LGR2: 4 inactive parent cells.
    // i=0 i=1 i=2          layer k = 0    |   i=2 i=3             layer k = 1
    //  0   0   0    j = 0                 |    0   0     j = 3
    //  0   0   0    j = 1                 |    0   0     j = 4
    //  0   0   0    j = 2

    // Refinement does not occur (since all parent cells are inactive, for all LGRs).
    BOOST_CHECK_EQUAL( grid.maxLevel(), 0);
    BOOST_CHECK_EQUAL( grid.levelGridView(0).size(0), grid.size(0));

    bool isSerial = grid.comm().size() == 1;
    // In serial, total active coarse cells (from level zero grid): 40 - 17(0's in ACTNUM block) = 23.
    if(isSerial) {
        checkGridInactiveCellsCount(grid, /* expected_maxLevel = */ 0, /* expected_active_cells = */ {23});
    }
    else { // Distribute level zero grid and add the same LGRs to the distributed view of the grid.
        createTestGridWithLgrsParallel(grid);
        checkGridInactiveCellsCount(grid, /* expected_maxLevel = */ 0, /* expected_global_cells = */ {23});
    }
}

BOOST_AUTO_TEST_CASE(refinementOccursIfAtLeastOneLgrHasAtLeastOneActiveParentCell)
{
    const std::string deckString =
        R"( RUNSPEC
  DIMENS
 -- NX NY NZ cells per x-,y-, and z-direction
     4 5 2 /
  GRID
  COORD -- grid corner coordinates (bounding box), 6*(NX +1)*(NY +1) values
  0 0 0   0 0 1 -- bottom of the pillar 0 | top of the pillar 0
  1 0 0   1 0 1 -- bottom of the pillar 1 | top of the pillar 1 (...)
  2 0 0   2 0 1
  3 0 0   3 0 1
  4 0 0   4 0 1
  0 1 0   0 1 1
  1 1 0   1 1 1
  2 1 0   2 1 1
  3 1 0   3 1 1
  4 1 0   4 1 1
  0 2 0   0 2 1
  1 2 0   1 2 1
  2 2 0   2 2 1
  3 2 0   3 2 1
  4 2 0   4 2 1
  0 3 0   0 3 1
  1 3 0   1 3 1
  2 3 0   2 3 1
  3 3 0   3 3 1
  4 3 0   4 3 1
  0 4 0   0 4 1
  1 4 0   1 4 1
  2 4 0   2 4 1
  3 4 0   3 4 1
  4 4 0   4 4 1
  0 5 0   0 5 1
  1 5 0   1 5 1
  2 5 0   2 5 1
  3 5 0   3 5 1
  4 5 0   4 5 1 -- bottom of the pillar 29 | top of the pillar 29
  /
  ZCORN
-- to easily deduce total active cells, no pinch-outs (consistent values to avoid collapsed cells, flat layers)
  80*0  -- top layer    k = 0
  80*1  -- bottom layer k = 0
  80*1  -- top layer    k = 1
  80*2  -- bottom layer k = 1
  /
  ACTNUM
-- i = 0 1 2 3
       0 0 0 1 -- layer k = 0     j = 0
       0 0 0 1 --                 j = 1
       0 0 1 1 --                 j = 2
       1 1 1 1 --                 j = 3
       1 1 1 0 --                 j = 4
       1 1 1 1 -- layer k = 1     j = 0
       0 1 1 1 --                 j = 1
       1 1 1 1 --                 j = 2
       1 1 0 0 --                 j = 3
       0 0 0 0 --                 j = 4
  /
  PORO
  40*0.15
  /)";

    Dune::CpGrid grid;
    createTestGridWithLgrsSerial(grid, deckString);
    // LGR1: 1 active, 8 inactive parent cells    |   LGR2: 4 inactive parent cells.
    // i=0 i=1 i=2          layer k = 0           |   i=2 i=3             layer k = 1
    //  0   0   0    j = 0                        |    0   0     j = 3
    //  0   0   0    j = 1                        |    0   0     j = 4
    //  0   0   1    j = 2

    bool isSerial = grid.comm().size() == 1;
    // In serial, total active coarse cells (from level zero grid): 40 - 14(0's in ACTNUM block) = 24.
    //            total active refined cells in LGR1: 1 active parent cells, number subd 2x2x2 -> 1*(2*2*2) = 8.
    //            total active refined cells in LGR2: 0 active parent cells, number subd 3x3x3 -> 0*(3*3*3) = 0.
    //            total active leaf cells: 24 in level zero - 1 parent cell + 8 new refined cells = 31.
    if (isSerial) {
        checkGridInactiveCellsCount(grid, /* expected_maxLevel = */ 2, /* expected_active_cells = */ {24,8,0,31});
    }
    else { // Distribute level zero grid and add the same LGRs to the distributed view of the grid.
        createTestGridWithLgrsParallel(grid);
        checkGridInactiveCellsCount(grid, /* expected_maxLevel = */ 2, /* expected_global_cells = */ {24,8,0,31});
    }
}

BOOST_AUTO_TEST_CASE(refineBlocksWithAtLeastOneActiveCell)
{

    const std::string deckString =
        R"( RUNSPEC
  DIMENS
 -- NX NY NZ cells per x-,y-, and z-direction
     4 5 2 /
  GRID
  COORD -- grid corner coordinates (bounding box), 6*(NX +1)*(NY +1) values.
  0 0 0   0 0 1 -- bottom of the pillar 0 | top of the pillar 0
  1 0 0   1 0 1 -- bottom of the pillar 1 | top of the pillar 1 (...)
  2 0 0   2 0 1
  3 0 0   3 0 1
  4 0 0   4 0 1
  0 1 0   0 1 1
  1 1 0   1 1 1
  2 1 0   2 1 1
  3 1 0   3 1 1
  4 1 0   4 1 1
  0 2 0   0 2 1
  1 2 0   1 2 1
  2 2 0   2 2 1
  3 2 0   3 2 1
  4 2 0   4 2 1
  0 3 0   0 3 1
  1 3 0   1 3 1
  2 3 0   2 3 1
  3 3 0   3 3 1
  4 3 0   4 3 1
  0 4 0   0 4 1
  1 4 0   1 4 1
  2 4 0   2 4 1
  3 4 0   3 4 1
  4 4 0   4 4 1
  0 5 0   0 5 1
  1 5 0   1 5 1
  2 5 0   2 5 1
  3 5 0   3 5 1
  4 5 0   4 5 1 -- bottom of the pillar 29 | top of the pillar 29
  /
  ZCORN
-- to easily deduce total active cells, no pinch-outs (consistent values to avoid collapsed cells, flat layers)
  80*0  -- top layer    k = 0
  80*1  -- bottom layer k = 0
  80*1  -- top layer    k = 1
  80*2  -- bottom layer k = 1
  /
  ACTNUM
-- i = 0 1 2 3
       0 0 1 1 -- layer k = 0     j = 0
       0 0 1 1 --                 j = 1
       1 1 1 1 --                 j = 2
       1 1 1 1 --                 j = 3
       1 1 1 0 --                 j = 4
       1 1 1 1 -- layer k = 1     j = 0
       0 1 1 1 --                 j = 1
       1 1 1 1 --                 j = 2
       1 1 1 1 --                 j = 3
       0 0 0 0 --                 j = 4
  /
  PORO
  40*0.15
  /)";

    Dune::CpGrid grid;
    createTestGridWithLgrsSerial(grid, deckString);
    // LGR1: 5 active, 4 inactive parent cells    |   LGR2: 2 inactive, 2 active parent cells.
    // i=0 i=1 i=2          layer k = 0           |   i=2 i=3             layer k = 1
    //  0   0   1    j = 0                        |    1   1     j = 3
    //  0   0   1    j = 1                        |    0   0     j = 4
    //  1   1   1    j = 2

    bool isSerial = grid.comm().size() == 1;
    // In serial, total active coarse cells (from level zero grid): 40 - 10(0's in ACTNUM block) = 30.
    //            total active refined cells in LGR1: 5 active parent cells, number subd 2x2x2 -> 5*(2*2*2) = 40.
    //            total active refined cells in LGR2: 2 active parent cells, number subd 3x3x3 -> 2*(3*3*3) = 54.
    //            total active leaf cells: 30 in level zero - 7 parent cells + (40+54) new refined cells = 117.
    if (isSerial) {
        checkGridInactiveCellsCount(grid, /* expected_maxLevel = */ 2, /* expected_active_cells = */  {30, 40, 54, 117});
    }
    else { // Distribute level zero grid and add the same LGRs to the distributed view of the grid.
        createTestGridWithLgrsParallel(grid);
        checkGridInactiveCellsCount(grid, /* expected_maxLevel = */ 2, /* expected_global_cells = */ {30, 40, 54, 117});
    }
}
