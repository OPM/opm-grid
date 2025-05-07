//===========================================================================
//
// File: lgr_with_inactive_parent_cells_test.cpp
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

#define BOOST_TEST_MODULE LgrWithInactiveParentCells

#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <tests/cpgrid/LgrChecks.hpp>

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
  ZCORN -- no pinch-outs: consistent ZCORN values for each cell corner to avoid collapsed cells (flat layers).
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
    Opm::createGridAndAddLgrs(grid,
                              deckString,
                              /* cells_per_dim_vec = */  {{2,2,2}, {3,3,3}},
                              /* startIJK_vec = */   {{0,0,0}, {2,3,1}},
                              /* endIJK_vec = */{{3,3,1}, {4,5,2}},
                              /* lgr_name_vec = */ {"LGR1", "LGR2"});

    // LGR1: 9 inactive parent cells.
    // i=0 i=1 i=2          layer k = 0
    //  0   0   0    j = 0
    //  0   0   0    j = 1
    //  0   0   0    j = 2

    // LGR2: 4 inactive parent cells.
    // i=2 i=3             layer k = 1
    //  0   0     j = 3
    //  0   0     j = 4

    if (grid.comm().size() == 1) {
        // In serial, total active coarse cells (from level zero grid): 40 - 17(0's in ACTNUM block) = 23.
        BOOST_CHECK_EQUAL( grid.levelGridView(0).size(0), 23);
        // Refinement does not occur (since all parent cells are inactive, for all LGRs).
        BOOST_CHECK_EQUAL( grid.maxLevel(), 0);
        BOOST_CHECK_EQUAL( grid.levelGridView(0).size(0), grid.size(0));

        Opm::checkGridWithLgrs(grid, /* cells_per_dim_vec = */ {{2,2,2}, {3,3,3}},   /* lgr_name_vec = */  {"LGR1", "LGR2"});
    }

    if (grid.comm().size()>1)
    {
        grid.loadBalance(/*overlapLayers*/ 1,
                         /*partitionMethod*/ Dune::PartitionMethod::zoltanGoG,
                         /*imbalanceTol*/ 1.1,
                         /*level*/ 0);

        grid.addLgrsUpdateLeafView( /* cells_per_dim_vec = */  {{2,2,2}, {3,3,3}},
                                    /* startIJK_vec = */   {{0,0,0}, {2,3,1}},
                                    /* endIJK_vec = */{{3,3,1}, {4,5,2}},
                                    /* lgr_name_vec = */ {"LGR1", "LGR2"});

        // In serial, total active coarse cells (from level zero grid): 40 - 17(0's in ACTNUM block) = 23.
        // auto global_active_coarse_cells = grid.comm().sum(grid.levelGridView(0).size(0)); count only interior!
        //  BOOST_CHECK_EQUAL( grid.levelGridView(0).size(0), 23);

        // Refinement does not occur (since all parent cells are inactive, for all LGRs).
        BOOST_CHECK_EQUAL( grid.maxLevel(), 0);
        BOOST_CHECK_EQUAL( grid.levelGridView(0).size(0), grid.size(0));

        Opm::checkGridWithLgrs(grid, /* cells_per_dim_vec = */ {{2,2,2}, {3,3,3}},   /* lgr_name_vec = */  {"LGR1", "LGR2"});
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
  ZCORN -- no pinch-outs: consistent ZCORN values for each cell corner to avoid collapsed cells (flat layers).
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
    Opm::createGridAndAddLgrs(grid,
                              deckString,
                              /* cells_per_dim_vec = */  {{2,2,2}, {3,3,3}},
                              /* startIJK_vec = */   {{0,0,0}, {2,3,1}},
                              /* endIJK_vec = */{{3,3,1}, {4,5,2}},
                              /* lgr_name_vec = */ {"LGR1", "LGR2"});

    // LGR1: 1 active, 8 inactive parent cells.
    // i=0 i=1 i=2          layer k = 0
    //  0   0   0    j = 0
    //  0   0   0    j = 1
    //  0   0   1    j = 2

    // LGR2: 4 inactive parent cells.
    // i=2 i=3             layer k = 1
    //  0   0     j = 3
    //  0   0     j = 4

    if (grid.comm().size() == 1) {
        // In serial, total active coarse cells (from level zero grid): 40 - 14(0's in ACTNUM block) = 24.
        BOOST_CHECK_EQUAL( grid.levelGridView(0).size(0), 24);
        BOOST_CHECK_EQUAL( grid.levelGridView(1).size(0), 1*(2*2*2)); // LGR1 1 active parent cells, number subd 2x2x2
        // 'Refined' level grid with all inactive parent cells is empty.
        BOOST_CHECK_EQUAL( grid.levelGridView(2).size(0), 0*(3*3*3)); // LGR2 all inactive parent cells, number subd 3x3x3
        // 24 active cells in level zero - 1 parent cells + 1*(2*2*2) new refined cells = 31
        BOOST_CHECK_EQUAL( grid.size(0), 31);

        Opm::checkGridWithLgrs(grid, /* cells_per_dim_vec = */ {{2,2,2}, {3,3,3}},   /* lgr_name_vec = */  {"LGR1", "LGR2"});
    }


    if (grid.comm().size()>1)
    {
        grid.loadBalance(/*overlapLayers*/ 1,
                         /*partitionMethod*/ Dune::PartitionMethod::zoltanGoG,
                         /*imbalanceTol*/ 1.1,
                         /*level*/ 0);

        grid.addLgrsUpdateLeafView( /* cells_per_dim_vec = */  {{2,2,2}, {3,3,3}},
                                    /* startIJK_vec = */   {{0,0,0}, {2,3,1}},
                                    /* endIJK_vec = */{{3,3,1}, {4,5,2}},
                                    /* lgr_name_vec = */ {"LGR1", "LGR2"});

        // In serial, total active coarse cells (from level zero grid): 40 - 17(0's in ACTNUM block) = 23.
        // auto global_active_coarse_cells = grid.comm().sum(grid.levelGridView(0).size(0)); count only interior!
        //  BOOST_CHECK_EQUAL( grid.levelGridView(0).size(0), 23);

        // Refinement does not occur (since all parent cells are inactive, for all LGRs).
        // BOOST_CHECK_EQUAL( grid.maxLevel(), 0);
        // BOOST_CHECK_EQUAL( grid.levelGridView(0).size(0), grid.size(0));

        Opm::checkGridWithLgrs(grid, /* cells_per_dim_vec = */ {{2,2,2}, {3,3,3}},   /* lgr_name_vec = */  {"LGR1", "LGR2"});
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
  ZCORN -- no pinch-outs: consistent ZCORN values for each cell corner to avoid collapsed cells (flat layers).
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
    Opm::createGridAndAddLgrs(grid,
                              deckString,
                              /* cells_per_dim_vec = */  {{2,2,2}, {3,3,3}},
                              /* startIJK_vec = */   {{0,0,0}, {2,3,1}},
                              /* endIJK_vec = */{{3,3,1}, {4,5,2}},
                              /* lgr_name_vec = */ {"LGR1", "LGR2"});

    // LGR1: 4 inactive, 5 inactive parent cells.
    // i=0 i=1 i=2          layer k = 0
    //  0   0   1    j = 0
    //  0   0   1    j = 1
    //  1   1   1    j = 2

    // LGR2: 2 inactive, 2 active parent cells.
    // i=2 i=3             layer k = 1
    //  1   1     j = 3
    //  0   0     j = 4

    if (grid.comm().size() == 1) {
        // In serial, total active coarse cells (from level zero grid): 40 - 10(0's in ACTNUM block) = 30.
        BOOST_CHECK_EQUAL( grid.levelGridView(0).size(0), 30);
        BOOST_CHECK_EQUAL( grid.maxLevel(), 2);
        BOOST_CHECK_EQUAL( grid.levelGridView(1).size(0), 5*(2*2*2)); // 40. LGR1 5 active parent cells, number subd per cell 2x2x2
        BOOST_CHECK_EQUAL( grid.levelGridView(2).size(0), 2*(3*3*3)); // 54. LGR2 2 active parent cells, number subd per cell 3x3x3
        // 30 active coarse cells - 7 active parent cells + new refined cells 5*(2*2*2) + 2*(3*3*3) = 30-7+40+54 = 117
        BOOST_CHECK_EQUAL( grid.size(0), 117);

        Opm::checkGridWithLgrs(grid, /* cells_per_dim_vec = */ {{2,2,2}, {3,3,3}},   /* lgr_name_vec = */  {"LGR1", "LGR2"});
    }

    if (grid.comm().size()>1)
    {

        grid.loadBalance(/*overlapLayers*/ 1,
                         /*partitionMethod*/ Dune::PartitionMethod::zoltanGoG,
                         /*imbalanceTol*/ 1.1,
                         /*level*/ 0);

        grid.addLgrsUpdateLeafView( /* cells_per_dim_vec = */  {{2,2,2}, {3,3,3}},
                                    /* startIJK_vec = */   {{0,0,0}, {2,3,1}},
                                    /* endIJK_vec = */{{3,3,1}, {4,5,2}},
                                    /* lgr_name_vec = */ {"LGR1", "LGR2"});

        //  BOOST_CHECK_EQUAL( grid.levelGridView(0).size(0), grid.size(0));
        BOOST_CHECK_EQUAL( grid.maxLevel(), 2);

        Opm::checkGridWithLgrs(grid, /* cells_per_dim_vec = */ {{2,2,2}, {3,3,3}},   /* lgr_name_vec = */  {"LGR1", "LGR2"});
    }
}
