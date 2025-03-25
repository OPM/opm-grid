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

#include <algorithm>
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
    grid.adapt({cells_per_dim}, assignRefinedLevel, {"LGR"+std::to_string(grid.maxLevel() +1)});
    grid.postAdapt();
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
                  const Dune::CpGrid& equivalent_grid,
                  bool gridLgrsHaveBlockShape,
                  bool gridHasBeenGlobalRefined)
{
    if(gridLgrsHaveBlockShape) { // For a mixed grid that gets refined a second time, isBlockShape == false, even though the marked elements form a block.
        Opm::checkCellBlockRefinements(grid, equivalent_grid);

        if (gridHasBeenGlobalRefined) {
            Opm::checkLeafGridGeometryEquality(grid, equivalent_grid);
        }
    }
}

void checkAdaptedGrid(Dune::CpGrid& grid,
                      const std::array<int,3>& cells_per_dim,
                      bool lgrsHaveBlockShape,
                      bool gridHasBeenGlobalRefined,
                      int preAdaptMaxLevel)
{
    const auto& data = grid.currentData();
    BOOST_CHECK(static_cast<int>(data.size()) == grid.maxLevel() +2);

    Opm::checkVertexAndFaceIndexAreNonNegative(grid);
    Opm::checkGridBasicHiearchyInfo(grid, {cells_per_dim}, preAdaptMaxLevel);
    Opm::checkGridLocalAndGlobalIdConsistency(grid, data);

    if (lgrsHaveBlockShape) {
        Opm::checkGlobalCellBounds(grid, data);
    }
}

BOOST_AUTO_TEST_CASE(emptyMarkedElemForRefinementSetDoesNothingToTheGrid)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    // If no elements are marked for refinement, calling adapt(/*args*/) or adapt()
    // will not modify the grid.
    adaptGridWithParams(grid, /* cells_per_dim = */ {2,2,2}, /* markedCells = */ {});
    adaptGrid(grid, /* markedCells = */ {});

    BOOST_CHECK_EQUAL( grid.maxLevel(), 0);
    BOOST_CHECK_EQUAL( grid.preAdapt(), false);
    BOOST_CHECK_EQUAL( grid.adapt(), false);
    BOOST_CHECK_EQUAL( grid.leafGridView().size(0), grid.levelGridView(0).size(0));
}

BOOST_AUTO_TEST_CASE(markAllElementsForRefinementIsEquivalentToCallGlobalRefinementWithOne)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    std::vector<int> markedCells(36); // 36 = 4x3x3
    std::iota(markedCells.begin(), markedCells.end(), 0);
    adaptGridWithParams(grid, /* cells_per_dim = */ {2,2,2}, markedCells);

    /** Why? check if we can change what's said below */
    // We set isBlockShape as false, even though global-refinement implies refinement of a block of cells..
    checkAdaptedGrid(grid,
                     /* cells_per_dim = */ {2,2,2},
                     /* lgrsHaveBlockShape = */ false,
                     /* gridHasBeenGlobalRefined = */ true,
                     /* preAdaptMaxLevel = */ 0);


    // Create other grid for comparison
    Dune::CpGrid equivalent_grid;
    equivalent_grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    equivalent_grid.globalRefine(1);

    // We set isBlockShape as false, even though global-refinement implies refinement of a block of cells.
    compareGrids(grid,
                 equivalent_grid,
                 /* lgrsHaveBlockShape = */ false,
                 /* gridHasBeenGlobalRefined = */ true);
}

BOOST_AUTO_TEST_CASE(markCellBlockForRefinementIsEquivalentToCallAddLgrsUpdateLeafView)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    std::vector<int> markedCells = {17,18};
    adaptGridWithParams(grid, /* cells_per_dim = */ {2,2,2}, markedCells);

    checkAdaptedGrid(grid,
                     /* cells_per_dim = */ {2,2,2},
                     /* lgrsHaveBlockShape = */ true,
                     /* gridHasBeenGlobalRefined = */ false,
                     /* preAdaptMaxLevel = */ 0);

    // Create other grid for comparison
    Dune::CpGrid equivalent_grid;
    equivalent_grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    equivalent_grid.addLgrsUpdateLeafView( /* cells_per_dim = */ {{2,2,2}},
                                      /* startIJK = */ {{1,1,1}},
                                      /* endIJK = */  {{3,2,2}}, // block cell indices = {17, 18}
                                      /* lgr_name = */  {"LGR1"});
    
    compareGrids(grid, equivalent_grid, /* lgrsHaveBlockShape = */ true, /* gridHasBeenGlobalRefined = */ false);
}

BOOST_AUTO_TEST_CASE(refinementOfCellsNotFormingABlockIsSupported)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    std::vector<int> markedCells = {0,4,5,17,18,30,31,35};

    // level zero grid view - cells marked for refinement are denoted with (/*index*/).
    // k = 0   8  9  10  11 | k = 1  20   21   22  23 | k = 2  32 33  34  (35) | 
    //       (4) (5)  6   7 |        16  (17) (18) 19 |        28 29 (30) (31) |
    //       (0)  1   2   3 |        12   13   14  15 |        24 25   26   27 |

    adaptGridWithParams(grid, /* cells_per_dim = */ {2,3,4}, markedCells);
  
    checkAdaptedGrid(grid,
                     /* cells_per_dim = */ {2,3,4},
                     /* lgrsHaveBlockShape = */ false,
                     /* gridHasBeenGlobalRefined = */ false,
                     /* preAdaptMaxLevel = */ 0);
}

BOOST_AUTO_TEST_CASE(callAdaptMultipleTimesAsLongAsCoarseMarkedElementsAreNotAtLgrBoundaries)
{

    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    std::vector<int> markedCells = {1,4,6,17,22,28,32};
    adaptGridWithParams(grid, /* cells_per_dim = */ {2,2,2}, markedCells);

    // level zero grid view - cells marked for refinement are denoted with (/*index*/).
    // k = 0  8   9  10  11 | k = 1  20  21  (22) 23 | k = 2  (32) 33 34 35 |
    //       (4)  5  (6)  7 |        16 (17)  18  19 |        (28) 29 30 31 |
    //        0  (1)  2   3 |        12  13   14  15 |          24 25 26 27 |


    std::vector<int> markedCells1 = {1,4,6};
    std::vector<int> markedCells2 = {38, 43}; // Equivalent cells to level 0 cells with indices {17,22};
    std::vector<int> markedCells3 = {63, 67}; // Equivalent cells to level 0 cells with indices {28,32};

    // Create a grid for comparison
    Dune::CpGrid equivalent_grid;
    equivalent_grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    adaptGrid(equivalent_grid, markedCells1);
    adaptGrid(equivalent_grid, markedCells2);
    adaptGrid(equivalent_grid, markedCells3);
    
    compareGrids(grid, equivalent_grid, /* lgrsHaveBlockShape = */ false, /* gridHasBeenGlobalRefined = */ false);
}

BOOST_AUTO_TEST_CASE(refinedCoarseOrRefinedCellsOnLgrBoundaryInAGridWithLgrsThrows)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    // level zero grid view - cells marked for refinement are denoted with (/*index*/).
    // k = 0  8   9  10  11 | k = 1  20  21   22  23 | k = 2  32 33 34 35 | 
    //        4   5   6   7 |        16 (17) (18) 19 |        28 29 30 31 |
    //        0   1   2   3 |        12  13   14  15 |        24 25 26 27 |
    grid.addLgrsUpdateLeafView( /* cells_per_dim = */ {{3,3,3}},
                                /* startIJK_vec = */ {{1,1,1}},
                                /* endIJK_vec = */ {{3,2,2}},
                                /* lgr_name_vec = */ {"LGR1"});

    // level 1 grid view, dimension 6x6x3 (leaf cell indices).
    //
    // LGR1 marked elements with elemIdx = 17 and 18, refined into 27 children cells each,
    // with leaf indices 17+0,...,17+26 = 43 (children {level 0, cell index 17}),
    //                   44,...,70 (children {level 0, cell index 18}).
    //
    // k = 2  41  42  43  | 68  69  70 |
    //        38  39  40  | 65  66  67 |
    //        35  36  37  | 62  63  64 |
    // ---------------------------------
    // k = 1  32  33  34  | 59  60  61 |
    //        29  30  31  | 56  57  58 |
    //        26  27  28  | 53  54  55 |
    // ---------------------------------
    // k = 0  23  24  25  | 50  51  52 |
    //        20  21  22  | 47  48  49 |
    //        17  18  19  | 44  45  46 |

    // leaf grid view
    // k = 0  8   9  10  11 | k = 1  72  73   74  75 | k = 2  84 85 86 87 | 
    //        4   5   6   7 |        16  ()* ()** 71 |        80 81 82 83 |
    //        0   1   2   3 |        12  13   14  15 |        76 77 78 79 |
    // *  indices 17, 18, ..., 43 (children of parent cell with index 17 in level zero grid)
    // ** indices 44, 45, ..., 70 (children of parent cell with index 18 in level zero grid)

    // - Examples of coarse cells touching the LGR1 on its boundary.
    // Cell 5, 14, 16, 71, 73, and 82, touching the bottom, front, left, right, back, and the top of LGR1, respectively.
    // - Example of refined cells with neighboring cell on a different level (here, coarser neighboring cell in level zero).
    // Cell 21 (coarse neighbor cell index 5) bottom lgr boundary,
    //      45 (coarse neighbor cell index 14) front lgr boundary,
    //      29 (coarse neighbor cel index 16) left lgr boundary,
    //      67 (coarse neighbor cell index 71) right lgr boundary,
    //      33 (coarse neighbor cell index 73) back lgr boundary,
    //      66 (coarse neighbor cell index 82) top lgr boundary.
    std::vector<int> throwingCells = {5, 14, 16, 71, 73, 82, 21, 45, 29, 67, 33, 66};
    std::for_each( throwingCells.begin(), throwingCells.end(), [&]( int elemIdx) {
        BOOST_CHECK_THROW( adaptGrid(grid, {elemIdx}) , std::logic_error);
        BOOST_CHECK_THROW( adaptGridWithParams(grid, /* cells_per_dim = */ {3,3,3}, {elemIdx}); , std::logic_error);
    });
}

BOOST_AUTO_TEST_CASE(refineCoarseCellsNotTouchingLgrBoundaryInAGridWithLgrsIsSupported)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    // level zero grid - cells marked for refinement are denoted with (/*index*/).
    // k = 0  8   9  10  11 | k = 1  20  21   22  23 | k = 2  32 33 34 35 | 
    //        4   5   6   7 |        16  17   18  19 |        28 29 30 31 |
    //        0   1   2  (3)|        12  13   14  15 |        24 25 26 27 |
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ { {2,2,2}},
                               /* startIJK_vec = */ {{3,0,0}},
                               /* endIJK_vec = */ {{4,1,1}},
                               /* lgr_name_vec = */ {"LGR1"});

    // leaf grid view (after adding LGR1 and before calling adapt with parameters).
    // k = 0 [15]  16  17  18 | k = 1  27  28   29  30 | k = 2  39 40 41 42 | 
    //       [11]  12  13  14 |        23  24   25  26 |        35 36 37 38 |
    //        [0]  [1]  2  ()*|        19  20   21  22 |        31 32 33 34 |
    //* indices 3-10 of refined cells (children of parent cell with index 3 in level zero grid).
    // [/*index*/] coarse marked cells for refinement, (not forming a block) not touching the lgr boundary.
    
    std::vector<int> markedCells = {0,1,11,15}; 
    adaptGridWithParams(grid, /* cells_per_dim = */ {2,3,4}, markedCells);
    
    checkAdaptedGrid(grid,
                     /* cells_per_dim = */ {2,3,4},
                     /* lgrsHaveBlockShape = */ false,
                     /* gridHasBeenGlobalRefined = */ false,
                     /* preAdaptMaxLevel = */ 1);
}

BOOST_AUTO_TEST_CASE(refineRefinedCellsInTheInteriorOfAnLgrIsSupported)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    // level zero grid view - cells marked for refinement are denoted with (/*index*/).
    // k = 0  8   9  10  11 | k = 1  20  21   22  23 | k = 2  32 33 34 35 | 
    //        4   5   6   7 |        16 (17) (18) 19 |        28 29 30 31 |
    //        0   1   2   3 |        12  13   14  15 |        24 25 26 27 |

    // LGR1 marked elements with elemIdx = 17 and 18, refined into 27 children cells each,
    // with leaf indices 17+0,...,17+26 = 43 (children {level 0, cell index 17}),44,...,70 (children {level 0, cell index 18}).
    grid.addLgrsUpdateLeafView( /* cells_per_dim = */ {{3,3,3}},
                                /* startIJK_vec = */ {{1,1,1}},
                                /* endIJK_vec = */ {{3,2,2}},
                                /* lgr_name_vec = */ {"LGR1"});


    // level 1 grid view, dimension 6x6x3 (leaf cell indices).
    // k = 2  41  42  43  | 68  69  70 |
    //        38  39  40  | 65  66  67 |
    //        35  36  37  | 62  63  64 |
    // ---------------------------------
    // k = 1  32  33  34  | 59  60  61 |
    //        29 [30][31] |[56][57] 58 |
    //        26  27  28  | 53  54  55 |
    // ---------------------------------
    // k = 0  23  24  25  | 50  51  52 |
    //        20  21  22  | 47  48  49 |
    //        17  18  19  | 44  45  46 |

    // Cells 30, 31, 56, and 57 are refined cells, located in the interior of the refined-level-grid-1 (lgr 1 / level 1).
    // Therefore, cell_to_face_ for all of them has size 6. (Their faces have all 2 refined neigboring cells - (not one coarse cell, and one refined)).
    std::vector<int> markedCells = {30,31, 56,57};
    adaptGridWithParams(grid, /* cells_per_dim = */ {2,3,4}, markedCells);

    checkAdaptedGrid(grid,
                     /* cells_per_dim = */ {2,3,4},
                     /* lgrsHaveBlockShape = */ false,
                     /* gridHasBeenGlobalRefined = */ false,
                     /* preAdaptMaxLevel = */ 1);
}

BOOST_AUTO_TEST_CASE(refineCoarseAndRefinedCellsAwayFromLgrBondaryIsSupported) {

    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    // level zero grid view - cells marked for refinement are denoted with (/*index*/).
    // k = 0  8   9  10  11 | k = 1  20  21   22  23 | k = 2  32 33 34 35 | 
    //        4   5   6   7 |        16 (17) (18) 19 |        28 29 30 31 |
    //        0   1   2   3 |        12  13   14  15 |        24 25 26 27 |
  
   // LGR1 marked elements with elemIdx = 17 and 18, refined into 27 children cells each,
    // with leaf indices 17+0,...,17+26 = 43 (children {level 0, cell index 17}),44,...,70 (children {level 0, cell index 18}).
    grid.addLgrsUpdateLeafView( /* cells_per_dim = */ {{3,3,3}},
                                /* startIJK_vec = */ {{1,1,1}},
                                /* endIJK_vec = */ {{3,2,2}},
                                /* lgr_name_vec = */ {"LGR1"});

    // level 1 grid view, dimension 6x6x3 (leaf cell indices).
    // k = 2  41  42  43  | 68  69  70 |
    //        38  39  40  | 65  66  67 |
    //        35  36  37  | 62  63  64 |
    // ---------------------------------
    // k = 1  32  33  34  | 59  60  61 |
    //        29 [30][31] |[56][57] 58 |
    //        26  27  28  | 53  54  55 |
    // ---------------------------------
    // k = 0  23  24  25  | 50  51  52 |
    //        20  21  22  | 47  48  49 |
    //        17  18  19  | 44  45  46 |
    //
    // Cells 30, 31, 56, and 57 are refined cells, located in the interior of the refined-level-grid-1 (lgr 1 / level 1).
    // Therefore, cell_to_face_ for all of them has size 6. (Their faces have all 2 refined neigboring cells - (not one coarse cell, and one refined)).


    // leaf grid view
    // k = 0  8  9 10 11 | k = 1  72  73   74  75 | k = 2  84 85 [86][87]| 
    //        4  5  6  7 |        16  ()* ()** 71 |        80 81  82  83 |
    //       [0][1] 2  3 |        12  13   14  15 |        76 77  78  79 |
    //
    // *  indices 17, 18, ..., 43 (children of parent cell with index 17 in level zero grid)
    // ** indices 44, 45, ..., 70 (children of parent cell with index 18 in level zero grid)
   
    // - Cells 30, 31, 56, and 57 are refined cells, located in the interior of the refined-level-grid-1 (lgr 1 / level 1).
    // Therefore, cell_to_face_ for all of them has size 6. (Their faces have all 2 refined neigboring cells - (not one coarse cell, and one refined)).
    // - Cells 0,1,2,86, and 87 are coarse cells, not touching the boundary of the LGR1 (cell 86 shares corners with LGR1 but do not share
    // any face. Therefore, the faces of cells 0,1,2,86, and 87 have all 1 or 2 neighboring coarse cells).
    std::vector<int> markedCells = {0,1,2,86,87,30,31,56,57};
    adaptGridWithParams(grid, /* cells_per_dim = */ {2,3,4}, markedCells);

    checkAdaptedGrid(grid,
                     /* cells_per_dim = */ {2,3,4},
                     /* lgrsHaveBlockShape = */ false,
                     /* gridHasBeenGlobalRefined = */ false,
                     /* preAdaptMaxLevel = */ 1);
}


/*BOOST_AUTO_TEST_CASE(refineMixedCells_in_multiLevelGrid) {
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





