//===========================================================================
//
// File: adapt_cpgrid_test.cpp
//
// Created: Apr 22 2024 16:45
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
#include "config.h"

#define BOOST_TEST_MODULE AdaptTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <tests/cpgrid/LgrChecks.hpp>

#include <array>
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
};

BOOST_GLOBAL_FIXTURE(Fixture);

void adaptGridWithParams(Dune::CpGrid& grid,
                         const std::array<int,3>& cells_per_dim,
                         const std::vector<int>& markedCells)
{
    const int startingGridIdx = grid.currentData().size() -1; // size before calling adapt
    std::vector<int> assignRefinedLevel(grid.currentData()[startingGridIdx]->size(0));

    for (const auto& elemIdx : markedCells)
    {
        const auto& elem =  Dune::cpgrid::Entity<0>(*(grid.currentData()[startingGridIdx]), elemIdx, true);
        grid.mark(1, elem);
        assignRefinedLevel[elemIdx] = grid.maxLevel() + 1;
        BOOST_CHECK( grid.getMark(elem) == 1);
        BOOST_CHECK( elem.mightVanish() == true);
    }
    grid.preAdapt();
    if( grid.preAdapt() ) { // markedCells can be empty
        grid.adapt({cells_per_dim}, assignRefinedLevel, {"LGR"+std::to_string(grid.maxLevel() +1)});
        grid.postAdapt();
    }
}

void adaptGrid(Dune::CpGrid& grid,
               const std::vector<int>& markedCells)
{
    const auto& leafGridView = grid.currentData().back();
    for (const auto& elemIdx : markedCells)
    {
        const auto& elem =  Dune::cpgrid::Entity<0>(*leafGridView, elemIdx, true);
        grid.mark(1, elem);
    }
    grid.preAdapt();
    grid.adapt();
    grid.postAdapt();
}

void compareGrids(const Dune::CpGrid& grid,
                  const Dune::CpGrid& other_grid,
                  bool isBlockShape,
                  bool hasBeenRefinedAtLeastOnce,
                  bool isGlobalRefinement)
{
    if(isBlockShape) { // For a mixed grid that gets refined a second time, isBlockShape == false, even though the marked elements form a block.
        Opm::checkCellBlockRefinements(grid, other_grid);

        if (isGlobalRefinement) {
            Opm::checkLeafGridGeometryEquality(grid, other_grid);
        }
    }
}

void checkAdaptedGrid(Dune::CpGrid& grid,
                      const std::array<int,3>& cells_per_dim,
                      bool isBlockShape,
                      bool hasBeenRefinedAtLeastOnce,
                      bool isGlobalRefinement)
{
    // BOOST_CHECK( grid.preAdapt() );
    const auto& data = grid.currentData();
    BOOST_CHECK(static_cast<int>(data.size()) == grid.maxLevel() +2);

    Opm::checkVertexAndFaceIndexAreNonNegative(grid);
    Opm::checkGridBasicHiearchyInfo(grid, {cells_per_dim});
    Opm::checkGridLocalAndGlobalIdConsistency(grid, data);

    /* for (int level = 1; level < coarse_grid.maxLevel(); ++level)
       {
       auto itMinLevel = std::min_element((data[level] -> global_cell_).begin(),  (data[level] -> global_cell_).end());
       auto itMaxLevel = std::max_element((data[level] -> global_cell_).begin(),  (data[level] -> global_cell_).end());
       BOOST_CHECK_EQUAL( *itMinLevel, 0);
       const auto& maxCartesianIdxLevel = data[level]->logical_cartesian_size_[0]*data[level]->logical_cartesian_size_[1]* data[level]->logical_cartesian_size_[2] -1;
       BOOST_CHECK_EQUAL( *itMaxLevel, maxCartesianIdxLevel);
       }*/

    //Opm::checkGlobalCellBounds(coarse_grid, data);

    /* auto itMin = std::min_element((data.back() -> global_cell_).begin(),  (data.back()-> global_cell_).end());
       auto itMax = std::max_element((data.back() -> global_cell_).begin(),  (data.back() -> global_cell_).end());
       BOOST_CHECK_EQUAL( *itMin, 0);
       const auto& maxCartesianIdx = coarse_grid.logicalCartesianSize()[0]*coarse_grid.logicalCartesianSize()[1]*coarse_grid.logicalCartesianSize()[2] -1;
       BOOST_CHECK_EQUAL( *itMax, maxCartesianIdx);*/
}

BOOST_AUTO_TEST_CASE(emptyMarkedElemForRefinementSetDoesNothingToTheGrid)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    adaptGridWithParams(grid, /* cells_per_dim = */ {2,2,2}, /* markedCells = */ {});

    // Create other grid for comparison
    Dune::CpGrid other_grid;
    other_grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    // isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement
    //compareGrids(grid, other_grid, true, false, false);
}


BOOST_AUTO_TEST_CASE(markNoElemForRefinementIsEquivalentToCallGlobalRefineWithZero)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    adaptGridWithParams(grid, /* cells_per_dim = */ {2,2,2}, /* markedCells = */ {});

    // Create other grid for comparison
    Dune::CpGrid other_grid;
    other_grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    other_grid.globalRefine(0);

    // isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement
    // compareGrids(grid, other_grid, true, false, true);
}


BOOST_AUTO_TEST_CASE(markAllElementsForRefinementIsEquivalentToCallGlobalRefinementWithOne)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    std::vector<int> markedCells(36); // 36 = 4x3x3
    std::iota(markedCells.begin(), markedCells.end(), 0);
    adaptGridWithParams(grid, /* cells_per_dim = */ {2,2,2}, markedCells);

    // We set isBlockShape as false, even though global-refinement implies refinement of a block of cells.
    //  isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    checkAdaptedGrid(grid, /* cells_per_dim = */ {2,2,2}, false, false, true);


    // Create other grid for comparison
    Dune::CpGrid other_grid;
    other_grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    other_grid.globalRefine(1);

    // We set isBlockShape as false, even though global-refinement implies refinement of a block of cells.
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    // isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement
    compareGrids(grid, other_grid, false, false, true);
}

BOOST_AUTO_TEST_CASE(markCellBlockForRefinementIsEquivalentToCallAddLgrsUpdateLeafView)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    std::vector<int> markedCells = {17,18};
    adaptGridWithParams(grid, /* cells_per_dim = */ {2,2,2}, markedCells);

    checkAdaptedGrid(grid, /* cells_per_dim = */ {2,2,2}, true, false, false);

    // Create other grid for comparison
    Dune::CpGrid other_grid;
    other_grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    other_grid.addLgrsUpdateLeafView( /* cells_per_dim = */ {{2,2,2}},
                                      /* startIJK = */ {{1,1,1}},
                                      /* endIJK = */  {{3,2,2}}, // block cell indices = {17, 18}
                                      /* lgr_name = */  {"LGR1"});


    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    compareGrids(grid, other_grid, true, false, false);
}

BOOST_AUTO_TEST_CASE(markCellsNotFormingABlock)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    std::vector<int> markedCells = {0,4,5,17,18,30,31,35};

    // k = 0   8  9  10  11 | k = 1  20   21   22  23 | k = 2  32 33  34  (35) |
    //       (4) (5)  6   7 |        16  (17) (18) 19 |        28 29 (30) (31) |
    //       (0)  1   2   3 |        12   13   14  15 |        24 25   26   27 |

    adaptGridWithParams(grid, /* cells_per_dim = */ {2,3,4}, markedCells);
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    checkAdaptedGrid(grid, /* cells_per_dim = */ {2,3,4}, false, false, false);
}

/*BOOST_AUTO_TEST_CASE(markNonBlockCells_compareAdapt)
  {
  // Create a grid
  Dune::CpGrid coarse_grid;
  const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
  const std::array<int, 3> grid_dim = {4,3,3};
  const std::array<int, 3> cells_per_dim = {3,3,3}; // WHY THIS DOES NOT FAIL??
  std::vector<int> markedCells = {1,4,6,9,17,22,28,32,33};
  coarse_grid.createCartesian(grid_dim, cell_sizes);

  // Create a grid
  Dune::CpGrid other_grid;
  other_grid.createCartesian(grid_dim, cell_sizes);
  for (const auto& elemIdx : markedCells)
  {
  const auto& elem =  Dune::cpgrid::Entity<0>(*(other_grid.currentData()[0]), elemIdx, true);
  other_grid.mark(1, elem);
  }
  other_grid.preAdapt();
  other_grid.adapt();
  other_grid.postAdapt();

  // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
  markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, other_grid, false, false, false);
  }*/

BOOST_AUTO_TEST_CASE(callAdaptMultipleTimes)
{

    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    std::vector<int> markedCells = {1,4,6,17,22,28,32};
    adaptGridWithParams(grid, /* cells_per_dim = */ {2,2,2}, markedCells);

    // k = 0  8   9  10  11 | k = 1  20  21  (22) 23 | k = 2  (32) 33 34 35 |
    //       (4)  5  (6)  7 |        16 (17)  18  19 |        (28) 29 30 31 |
    //        0  (1)  2   3 |        12  13   14  15 |          24 25 26 27 |


    std::vector<int> markedCells1 = {1,4,6};
    std::vector<int> markedCells2 = {38, 43}; // Equivalent cells to level 0 cells with indices {17,22};
    std::vector<int> markedCells3 = {63, 67}; // Equivalent cells to level 0 cells with indices {28,32};

    // Create a grid
    Dune::CpGrid other_grid;
    other_grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    adaptGrid(other_grid, markedCells1);
    adaptGrid(other_grid, markedCells2);
    adaptGrid(other_grid, markedCells3);

    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    // INFO ABOUR GRID (NOT OTHER_GRID)
    compareGrids(grid, other_grid, false, false, false);
}

/*BOOST_AUTO_TEST_CASE(markCoarseCellsOnLgrBoundaryOfGridWithLgrsThrows)
  {
  // Create a grid
  Dune::CpGrid coarse_grid;
  const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
  const std::array<int, 3> grid_dim = {4,3,3};
  const std::array<int, 3> cells_per_dim = {3,3,3};
  coarse_grid.createCartesian(grid_dim, cell_sizes);

  // LGR1 marked elements with elemIdx = 17 and 18, refined into 27 children cells each,
  // with leaf indices 17+0,...,17+26 = 43 (children {level 0, cell index 17}),44,...,70 (children {level 0, cell index 18}).
  const std::array<int, 3> startIJK = {1,1,1};
  const std::array<int, 3> endIJK = {3,2,2};
  const std::string lgr_name = {"LGR1"};
  coarse_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

  // - Coarse cells touching the LGR1 on its boundary.
  // Cell 5, 14, 16, 71, 73, and 82, touching the bottom, front, left, right, back, and the top of LGR1, respectively.
  std::vector<int> markedCells = {5,14,16,71,73,82};
  // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
  BOOST_CHECK_THROW(markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, true, true, false), std::logic_error);
  }*/


// Nested cases -

/*BOOST_AUTO_TEST_CASE(refineCoarseCells_in_mixedGrid_notTouchingLgrBoundary) {
// Create a grid
Dune::CpGrid coarse_grid;
const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
const std::array<int, 3> grid_dim = {4,3,3};
const std::array<int, 3> cells_per_dim = {2,2,2};
coarse_grid.createCartesian(grid_dim, cell_sizes);

// LGR1 marked element with elemIdx = 3, refined into 8 children cells with leaf indices 3,...,10.
const std::array<int, 3> startIJK = {3,0,0};
const std::array<int, 3> endIJK = {4,1,1};
const std::string lgr_name = {"LGR1"};
coarse_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

std::vector<int> markedCells = {0,1,11,15}; // coarse cells (in level 0 grid, this cell has index 8)
// The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, false, true, false);
}

BOOST_AUTO_TEST_CASE(refineInteriorRefinedCells_in_mixedGrid) {
// Create a grid
Dune::CpGrid coarse_grid;
const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
const std::array<int, 3> grid_dim = {4,3,3};
const std::array<int, 3> cells_per_dim = {3,3,3};
coarse_grid.createCartesian(grid_dim, cell_sizes);

// LGR1 marked elements with elemIdx = 17 and 18, refined into 27 children cells each,
// with leaf indices 17+0,...,17+26 = 43 (children {level 0, cell index 17}),44,...,70 (children {level 0, cell index 18}).
const std::array<int, 3> startIJK = {1,1,1};
const std::array<int, 3> endIJK = {3,2,2};
const std::string lgr_name = {"LGR1"};
coarse_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

// Cells 30, 31, 56, and 57 are refined cells, located in the interior of the refined-level-grid-1 (lgr 1 / level 1).
// Therefore, cell_to_face_ for all of them has size 6. (Their faces have all 2 refined neigboring cells - (not one coarse cell, and one refined)).
std::vector<int> markedCells = {30,31, 56,57};
// The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, false, true, false);
}


BOOST_AUTO_TEST_CASE(refineMixedCells_in_mixedGrid) {
// Create a grid
Dune::CpGrid coarse_grid;
const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
const std::array<int, 3> grid_dim = {4,3,3};
const std::array<int, 3> cells_per_dim = {3,3,3};
coarse_grid.createCartesian(grid_dim, cell_sizes);

// LGR1 marked elements with elemIdx = 17 and 18, refined into 27 children cells each,
// with leaf indices 17+0,...,17+26 = 43 (children {level 0, cell index 17}),44,...,70 (children {level 0, cell index 18}).
const std::array<int, 3> startIJK = {1,1,1};
const std::array<int, 3> endIJK = {3,2,2};
const std::string lgr_name = {"LGR1"};
coarse_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

// - Cells 30, 31, 56, and 57 are refined cells, located in the interior of the refined-level-grid-1 (lgr 1 / level 1).
// Therefore, cell_to_face_ for all of them has size 6. (Their faces have all 2 refined neigboring cells - (not one coarse cell, and one refined)).
// - Cells 0,1,2,12, and 15 are coarse cells, not touching the boundary of the LGR1 (cells 12 and 15 do share corners with LGR1 but do not share
// any face. Therefore, the faces of cells 0,1,2,12,and 15 have all 1 or 2 neighboring coarse cells).
std::vector<int> markedCells = {0,1,2,12,15,30,31,56,57};
// The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, false, true, false);
}


BOOST_AUTO_TEST_CASE(refineMixedCells_in_multiLevelGrid) {
// Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {3,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);

    // LGR1: element with elemIdx = 17, refined into 27 children cells with leaf indices 17+0,...,17+26 = 43 (children {level 0, cell index 17}).
    // LGR2: element with elemIdx = 18, refined into 27 children cells with leaf indices 44,...,70 (children {level 0, cell index 18}).
    const std::vector<std::array<int, 3>> startIJK_vec = {{1,1,1}, {2,1,1}};
    const std::vector<std::array<int, 3>> endIJK_vec = {{2,2,2}, {3,2,2}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    coarse_grid.addLgrsUpdateLeafView({cells_per_dim, cells_per_dim}, startIJK_vec, endIJK_vec, lgr_name_vec);

    // - Cells 30, 31, 56, and 57 are refined cells, located in the interior of the refined-level-grid-1 (lgr 1 / level 1).
    // Therefore, cell_to_face_ for all of them has size 6. (Their faces have all 2 refined neigboring cells - (not one coarse cell, and one refined)).
    // - Cells 0,1,2,12, and 15 are coarse cells, not touching the boundary of the LGR1 (cells 12 and 15 do share corners with LGR1 but do not share
    // any face. Therefore, the faces of cells 0,1,2,12,and 15 have all 1 or 2 neighboring coarse cells).
    std::vector<int> markedCells = {0,1,2,12,15,30,31,56,57};
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, false, true, false);
}


BOOST_AUTO_TEST_CASE(refineMixedCells_in_mixedGrid_II)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {3,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);

    // LGR1 marked elements with elemIdx = 17 and 18, refined into 27 children cells each,
    // with leaf indices 17+0,...,17+26 = 43 (children {level 0, cell index 17}),44,...,70 (children {level 0, cell index 18}).
    const std::array<int, 3> startIJK = {1,1,1};
    const std::array<int, 3> endIJK = {3,2,2};
    const std::string lgr_name = {"LGR1"};
    coarse_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

    // - Cells 72 (in level 0 with cell index 20) and 84 (in level 0 with cell index 32) are coarse cells,
    // sharing one K_FACE, that do not share faces with LGR1 (they do share corners).
    // - Cells 25,34,43 are refined cells, children of {level 0, cell index 17}, forming a collum.
    // - Cells 50,59,68 are refined cells, children of {level 0, cell index 18}, forming a collum.
    // Cells 25 and 50, 34 and 59, 43 and 68, share a face (the collums are next to each other).
    std::vector<int> markedCells = {25,34,43,50,59,68, 72, 84};
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, true, true, false);
}
*/





