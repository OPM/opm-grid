/*
  Copyright 2026 Equinor ASA.

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

#define BOOST_TEST_MODULE AquiferCellsAndConnectionsNotRefinedTests

#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <tests/cpgrid/lgr/LgrChecks.hpp>

#include <string>

struct Fixture
{
    Fixture()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(m_argc, m_argv);
        Opm::OpmLog::setupSimpleDefaultLogging();
    }
};

BOOST_GLOBAL_FIXTURE(Fixture);

// Local Grid Refinement (LGR) in CpGrid does not support refinement of
// aquifer cells or cells that have neighboring aquifer connections.
//
// The following methods will throw an exception if aquifer cells or
// aquifer connections occur inside an LGR region:
//
//  - CpGrid::addLgrsUpdateLeafView(...)  (used by the CARFIN keyword)
//  - CpGrid::autoRefine(...)             (used by the AUTOREF keyword)
//
// In contrast, DUNE Grid refinement methods:
//
//  - CpGrid::globalRefine(...)
//  - CpGrid::adapt()
//
// ignore aquifer cells and aquifer connections during
// refinement and perform the refinement only on the remaining cells.

// To create a CpGrid containing aquifer cells and connections,
// the input deck must include the following keywords:
//
//   - AQUDIMS  : defines aquifer dimensions
//   - AQUNUM   : defines the number of aquifers
//   - AQUCON   : defines aquifer connections

// ** AQUDIMS **
// MXNAQN MXNAQC NIFTBL NRIFTB NANAQ NCAMAX MXNALI MXNAQL
// - MXNAQN maximum number of aquifers
// - MXNAQC maximum number of aquifer connections
// - NIFTBL maximum number of Carter-Tracy aquifer tables associated with AQUTAB
// - NRIFTB maximum number of rows in the Carter-Tracy aquifer tables associated with AQUTAB
// - NANAQ  maximum number of analytical aquifers defined by AQUFETP, AQUFLUX and AQUCT
// - NCAMAX maximum number of cells connected to an analytical aquifer
// - MXNALI maximum number of aquifer lists
// - MXNAQL maximum number of analytical aquifers in any single aquifer list

// ** AQUNUM **  Numerical aquifer description
// - AQUID   id number >= 1 and <= maximum number of aquifers (MXNAQN in AQUDIMS)
// - I, J, K location 1<= I, J, K <= NX, NY, NZ respectively.
// - AREA    cross-sectional area of the aquifer used in calculating the aquifer connection transmissibility
// - LENGTH  length of the numerical aquifer
// - PORO
// - PERM
// - DATUM   reference datum depth of the numerical aquifer
// - PRESS   numerial aquifer pressure at DATUM
// - PVTNUM
// - SATNUM

// ** AQUCON **   Numerial aquifer connections
// - id number
// - box I1 I2 J1 J2 K1 K2
// - connect face
// - trans mult
// - trans optn
// - adjoin cells

const std::string deckString =
    R"( RUNSPEC
  DIMENS
 -- NX NY NZ cells per x-,y-, and z-direction
     4 5 2 /
AQUDIMS
--  MXNAQN MXNAQC NIFTBL NRIFTB NANAQ NCAMAX MXNALI MXNAQL
    1      1      1*      1*     1*     1*      1*      1* /
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
       1 0 1 1 --                 j = 1
       1 1 1 1 --                 j = 2
       1 1 1 1 --                 j = 3
       1 1 1 0 --                 j = 4
       1 1 1 1 -- layer k = 1     j = 0
       1 1 1 1 --                 j = 1
       1 1 1 1 --                 j = 2
       1 1 1 1 --                 j = 3
       0 0 0 0 --                 j = 4
  /
AQUNUM
-- ID  I J K AREA LENGTH PORO PERM DEPTH PRESS PVTNUM SATNUM
   1   1 3 1 0.5  0.5    0.15 0.3  1*    1*    1*     1*     /
/
AQUCON
-- ID  I1 I2 J1 J2 K1 K2  ConnectFace  TransMult TransOptn AdjoinCells
   1   1  1  1  2  1  1   'J+'         1*        1*        'NO' /
/
  PORO
  40*0.15
  /
PERMX
40*100
/
PERMY
40*100
/
PERMZ
40*100
/)";

// The grid has an aquifer cell at ijk = {0,2,0} in CpGrid (IJK_deck = {1,3,1}).
//
// The grid has an aquifer connection:
// connected grid cells with ijk_1 = {0,0,0} and ijk_2 = {0,1,0}
//                     (IJK_deck_1 = {1,1,1} and IJK_deck_2 = {1,2,1})
// with J+
//
//         Cartesian level zero indices
// k = 0    16  17  18  19
//          12  13  14  15
//          |8|  9  10  11       Cell with index 8 is aquifer cell (ijk = {0,2,0}).
//          |4|  5   6   7       Cells 0 (ijk = {0,0,0}) and 4 (ijk= {0,1,0}) connect
//          |0|  1   2   3       with aquifer cell in J+.
// -----------------------

void aquiferCellsAndConnsIgnoredInRefinement(const Dune::CpGrid& grid)
{
    const auto& levelZeroData = *grid.currentData().front();
    const auto& levelZeroAquCells = levelZeroData.sortedNumAquiferCells();

    for (const auto& aquCellIdx : levelZeroAquCells) {

        const auto& aquElem = Dune::cpgrid::Entity<0>(levelZeroData, aquCellIdx, true);
        BOOST_CHECK( aquElem.isLeaf()); // aquifer cells are not involved in refinement

        const auto& leafAquCellIdx = levelZeroData.getLeafIdxFromLevelIdx(aquCellIdx);
        const auto& leafAquElem = Dune::cpgrid::Entity<0>(grid.currentLeafData(), leafAquCellIdx, true);

        for (const auto& intersection : Dune::intersections(grid.leafGridView(), leafAquElem)){
            if (intersection.neighbor()) {
                BOOST_CHECK(intersection.outside().isLeaf()); // connections to aquifers not involved in refinement
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(addLgrsUpdateLeafViewThrowsIfLgrHasAquiferCellOrConnention)
{
    Dune::CpGrid grid;
    Opm::createGridFromDeckString(grid, deckString);

    // ** LGR with Aquifer Cell**
    //
    // LGR1: 6 active parent cells                |   Cartesian indices of parent cells in LGR1
    // i=0 i=1 i=2          layer k = 0           |   12  13  14
    //  1   1   1    j = 2                        |   |8|  9  10      Cell 8 is aquifer cell.
    //  1   1   1    j = 3                        |
    // LGR1 parent cell block contains an         |
    // aquifer cell ijk = {0,2,0} (in deckString  |
    // IJK_deck = {1,3,1}).                       |
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView(  /* cells_per_dim_vec = */  {{2,2,2}},
                                                   /* startIJK_vec = */       {{0,2,0}},
                                                   /* endIJK_vec = */         {{3,4,1}},
                                                   /* lgr_name_vec = */       {"LGR1"}), std::logic_error);


    // ** LGR with Aquifer Connection **
    //
    // LGR1: 2 active parent cells                |   Cartesian indices of parent cells in LGR1
    // i=0 i=1 i=2          layer k = 0           |   |4|  5  6     Cell 4 is connected to the aquifer cell |8|.
    //  1   0   1    j = 1                        |
    // LGR1 parent cell block contains an         |
    // aquifer connection ijk = {0,1,0}           |
    // (IJK_deck = {1,2,1}).                      |
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView(  /* cells_per_dim_vec = */  {{2,2,2}},
                                                   /* startIJK_vec = */       {{0,1,0}},
                                                   /* endIJK_vec = */         {{3,2,1}},
                                                   /* lgr_name_vec = */       {"LGR1"}), std::logic_error);
}

BOOST_AUTO_TEST_CASE(refinementViaAddLgrsUpdateLeafViewOccursIfLgrHasNoAquiferCellOrConnention)
{
    Dune::CpGrid grid;
    Opm::createGridFromDeckString(grid, deckString);

    // LGR1: 3 active parent cells                |   Cartesian indices of parent cells in LGR1
    // i=0 i=1 i=2          layer k = 1           |   24  25  26
    //  1   1   1    j = 3                        |   20  21  22      LGR does not contain aquifer data.
    //  0   0   0    j = 4                        |
    grid.addLgrsUpdateLeafView(  /* cells_per_dim_vec = */  {{2,2,2}},
                                 /* startIJK_vec = */       {{0,3,1}},
                                 /* endIJK_vec = */         {{3,4,2}},
                                 /* lgr_name_vec = */       {"LGR1"});

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{2,2,2}},
                           /* lgr_name_vec = */ {"LGR1"});

    aquiferCellsAndConnsIgnoredInRefinement(grid); // grid has aquifer data but outside LGRs
}

BOOST_AUTO_TEST_CASE(autoRefineThrowsIfGridHasAquiferCellOrConnections)
{
    Dune::CpGrid grid;
    Opm::createGridFromDeckString(grid, deckString);

    BOOST_CHECK_THROW(grid.autoRefine(/* nxnynz = */ std::array{3,3,5}), std::logic_error);
}

BOOST_AUTO_TEST_CASE(refinemenViaGlobalRefineOccursIgnoringAquiferCellsAndConnections)
{
    Dune::CpGrid grid;
    Opm::createGridFromDeckString(grid, deckString);

    grid.globalRefine(2);

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{2,2,2}, {2,2,2}},
                           /* lgr_name_vec = */ {"GR1", "GR2"}, // aquifer cells and connections are ignored
                           /* preRefinedMaxLevel = */ 0,
                           /* isNested = */ true);

    aquiferCellsAndConnsIgnoredInRefinement(grid);
}


BOOST_AUTO_TEST_CASE(refinementViaAdaptOccursIgnoringAquiferCellsAndConnections)
{
    Dune::CpGrid grid;
    Opm::createGridFromDeckString(grid, deckString);

    for (const auto& element : Dune::elements(grid.leafGridView())) {
        grid.mark(/* refine mark = */ 1, element);
    }

    grid.preAdapt();
    grid.adapt();
    grid.postAdapt();

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{2,2,2}},
                           /* lgr_name_vec = */ {"GR1"}); // aquifer cells and connections are ignored

    aquiferCellsAndConnsIgnoredInRefinement(grid);
}
