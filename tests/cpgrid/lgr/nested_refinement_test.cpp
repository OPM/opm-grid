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

#define BOOST_TEST_MODULE NestedRefinementTests
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
};

BOOST_GLOBAL_FIXTURE(Fixture);


BOOST_AUTO_TEST_CASE(ifNonParentGridNameProvidedDefaultIsAllChildGridsFromGlobal) {

    const std::vector<std::array<int,3>> cells_per_dim_vec = {{3,3,1}, {3,3,1}};
    const std::vector<std::array<int,3>> startIJK_vec = {{1,1,0}, {2,2,0}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,2,1}, {3,3,1}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    const std::vector<std::string> lgr_parent_grid_name_vec = {"GLOBAL", "GLOBAL"};

    //   GLOBAL            The grids are stored: GLOBAL,
    //   |    |                                  LGR1, LGR2 (child grids from GLOBAL),
    // LGR1  LGR2                                leaf grid view (without name).
    
    Dune::CpGrid grid, equiv_grid;
    grid.createCartesian(/* grid_dim = */ {3,3,1}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    equiv_grid.createCartesian(/* grid_dim = */ {3,3,1}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    grid.addLgrsUpdateLeafView(cells_per_dim_vec,
                               startIJK_vec,
                               endIJK_vec,
                               lgr_name_vec,
                               lgr_parent_grid_name_vec);
    
    equiv_grid.addLgrsUpdateLeafView(cells_per_dim_vec,
                                     startIJK_vec,
                                     endIJK_vec,
                                     lgr_name_vec);

    Opm::checkLeafGridGeometryEquality(grid,
                                       equiv_grid);

    BOOST_CHECK_EQUAL(grid.maxLevel(), 2);
    BOOST_CHECK_EQUAL(grid.levelGridView(1).size(0), 9); // 1 parent cell from GLOBAL, refined into 3x3x1 children
    BOOST_CHECK_EQUAL(grid.levelGridView(2).size(0), 9); // 1 parent cell from LGR1, refined into 3x3x1 children
    BOOST_CHECK_EQUAL(grid.leafGridView().size(0), 25); // 3x3x1 level zero - 1 parent cell + 9 LGR1 - 1 LGR1-parent + 9 LGR2 = 25

    for (int level = 1; level <=2; ++level) { // unique parent grid = "GLOBAL"->level zero
        for (const auto& element : Dune::elements(grid.levelGridView(1))) {
            BOOST_CHECK_EQUAL( element.father().level(), 0);
        }
    }

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{3,3,1}, {3,3,1}},
                           /* lgr_name_vec = */ {"LGR1", "LGR2"},
                           /* gridHasBeenGlobalRefined = */ false,
                           /* preRefineMaxLevel = */ 0,
                           /* isNested = */ false);
}

BOOST_AUTO_TEST_CASE(nestedRefinementOnly) {

    const std::vector<std::array<int,3>> cells_per_dim_vec = {{3,3,1}, {3,3,1}, {3,3,1}, {3,3,1}};
    const std::vector<std::array<int,3>> startIJK_vec = {{1,1,0},{1,1,0}, {1,1,0}, {1,1,0}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,2,1},{2,2,1}, {2,2,1}, {2,2,1}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3", "LGR4"};
    const std::vector<std::string> lgr_parent_grid_name_vec = {"GLOBAL","LGR1","LGR2","LGR3"};

    // GLOBAL            The grids are stored: GLOBAL,
    //  |                                      LGR1 (child grid from GLOBAL),
    // LGR1                                    LGR2 (child grid from LGR1),
    //  |                                      LGR3 (child grid from LGR2),
    // LGR2                                    LGR4 (child grid from LGR3),
    //  |                                      leaf grid view (without name).
    // LGR3                                    
    //  |
    // LGR4
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {3,3,1}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    grid.addLgrsUpdateLeafView(cells_per_dim_vec,
                               startIJK_vec,
                               endIJK_vec,
                               lgr_name_vec,
                               lgr_parent_grid_name_vec);

    BOOST_CHECK_EQUAL(grid.maxLevel(), 4);
    BOOST_CHECK_EQUAL(grid.levelGridView(1).size(0), 9); // 1 parent cell from GLOBAL, refined into 3x3x1 children
    BOOST_CHECK_EQUAL(grid.levelGridView(2).size(0), 9); // 1 parent cell from LGR1, refined into 3x3x1 children
    BOOST_CHECK_EQUAL(grid.levelGridView(3).size(0), 9); // 1 parent cell from GLOBAL, refined into 3x3x1 children
    BOOST_CHECK_EQUAL(grid.levelGridView(4).size(0), 9); // 1 parent cell from LGR1, refined into 3x3x1 children
    BOOST_CHECK_EQUAL(grid.leafGridView().size(0), 41);
    // 3x3x1 level zero - 1 parent cell + 9 LGR1 - 1 LGR1-parent + 9 LGR2 - 1 l0-parent-cell + 9 LGR3 - 1 LGR-parent + 9 LGR4 = 41

    for (const auto& element : Dune::elements(grid.levelGridView(1))) {
        BOOST_CHECK_EQUAL( element.father().level(), 0);
    }

    for (const auto& element : Dune::elements(grid.levelGridView(2))) { 
        BOOST_CHECK_EQUAL( element.father().level(), 1);
    }

    for (const auto& element : Dune::elements(grid.levelGridView(3))) {
        BOOST_CHECK_EQUAL( element.father().level(), 2);
    }

    for (const auto& element : Dune::elements(grid.levelGridView(4))) {
        BOOST_CHECK_EQUAL( element.father().level(), 3);
    }

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{3,3,1}, {3,3,1}, {3,3,1}, {3,3,1}},
                           /* lgr_name_vec = */ {"LGR1", "LGR2", "LGR3", "LGR4"},
                           /* gridHasBeenGlobalRefined = */ false,
                           /* preRefineMaxLevel = */ 0,
                           /* isNested = */ true);
}

BOOST_AUTO_TEST_CASE(mixNameOrderAndNestedRefinement){

    const std::vector<std::array<int,3>> cells_per_dim_vec = {{3,3,1}, {3,3,1}, {3,3,1}, {3,3,1}};
    const std::vector<std::array<int,3>> startIJK_vec = {{1,1,0},{0,0,0}, {0,0,0}, {1,1,0}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,2,1},{1,1,1}, {1,1,1}, {2,2,1}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3", "LGR4"};
    const std::vector<std::string> lgr_parent_grid_name_vec = {"GLOBAL","LGR1","GLOBAL","LGR3"};

    //   GLOBAL            The grids are stored: GLOBAL,
    //   |    |                                  LGR1, LGR3 (child grid from GLOBAL),
    // LGR1  LGR3                                LGR2 (child grid from LGR1),
    //   |    |                                  LGR4 (child grid from LGR3),
    // LGR2  LGR4                                leaf grid view (without name).
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {3,3,1}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    grid.addLgrsUpdateLeafView(cells_per_dim_vec,
                               startIJK_vec,
                               endIJK_vec,
                               lgr_name_vec,
                               lgr_parent_grid_name_vec);

    BOOST_CHECK_EQUAL(grid.maxLevel(), 4);
    BOOST_CHECK_EQUAL(grid.levelGridView(0).size(0), 9);
    BOOST_CHECK_EQUAL(grid.levelGridView(1).size(0), 9); // 1 parent cell from GLOBAL, refined into 3x3x1 children
    BOOST_CHECK_EQUAL(grid.levelGridView(2).size(0), 9); // 2 parent cell from LGR1, refined into 3x3x1 children
    BOOST_CHECK_EQUAL(grid.levelGridView(3).size(0), 9); // 1 parent cell from LGR2, refined into 3x3x1 children
    BOOST_CHECK_EQUAL(grid.levelGridView(4).size(0), 9); // 1 parent cell from LGR2, refined into 3x3x1 children
    BOOST_CHECK_EQUAL(grid.leafGridView().size(0), 41);
    // 3x3x1 level zero - 2 parent cells + 9 LGR1 - 1 LGR1-parent + 9 LGR2 - 1 l0-parent-cell + 9 LGR3 - 1 LGR-parent + 9 LGR4 = 41
    
    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{3,3,1}, {3,3,1}, {3,3,1}, {3,3,1}},
                           /* lgr_name_vec = */ {"LGR1", "LGR3", "LGR2", "LGR4"},
                           /* gridHasBeenGlobalRefined = */ false,
                           /* preRefineMaxLevel = */ 0,
                           /* isNested = */ true);

    for (const auto& element : Dune::elements(grid.levelGridView(1))) { // LGR1
        BOOST_CHECK_EQUAL( element.father().level(), 0); // LGR1 parent grid is GLOBAL
    }

    for (const auto& element : Dune::elements(grid.levelGridView(2))) { // LGR3 
        BOOST_CHECK_EQUAL( element.father().level(), 0); // LGR3 parent grid is GLOBAL
    }

    for (const auto& element : Dune::elements(grid.levelGridView(3))) { // LGR2
        BOOST_CHECK_EQUAL( element.father().level(), 1); // LGR2 parent grid is LGR1
    }

    for (const auto& element : Dune::elements(grid.levelGridView(4))) { // LGR4
        BOOST_CHECK_EQUAL( element.father().level(), 2); // LGR4 parent grid is LGR3-> grid index 2
    }
}

BOOST_AUTO_TEST_CASE(throwIfParentGridNameDoesNotExitBeforeItsLgrs){

    const std::vector<std::array<int,3>> cells_per_dim_vec = {{3,3,1}, {2,2,2}, {3,2,1}, {3,4,1}};
    const std::vector<std::array<int,3>> startIJK_vec = {{1,1,0},{0,0,0}, {0,0,0}, {1,1,0}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,2,1},{1,1,1}, {1,1,1}, {2,2,1}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3", "LGR4"};
    const std::vector<std::string> lgr_parent_grid_name_vec = {"GLOBAL","LGR3","GLOBAL","LGR1"};

    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {3,3,1}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    BOOST_CHECK_THROW( grid.addLgrsUpdateLeafView(cells_per_dim_vec,
                                                  startIJK_vec,
                                                  endIJK_vec,
                                                  lgr_name_vec,
                                                  lgr_parent_grid_name_vec), std::invalid_argument);
}

