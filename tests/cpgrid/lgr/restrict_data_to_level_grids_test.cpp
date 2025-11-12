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

#define BOOST_TEST_MODULE RestrictDataToLevelGridsTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CartesianIndexMapper.hpp>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>
#include <opm/grid/cpgrid/LgrOutputHelpers.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>
#include <opm/output/data/Cells.hpp>
#include <opm/output/data/Solution.hpp>

#include <array>
#include <string>
#include <utility>  // for std::move
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

void restrictFakeLeafDataToLevelGrids(const Dune::CpGrid& grid,
                                      const std::vector<std::vector<double>>& expected_data_levels)
{
    // Create a fake leaf solution
    Opm::data::Solution leafSolution{};
    BOOST_CHECK(!leafSolution.has("NOTHING")); // dummy check

    const auto& leafView = grid.leafGridView();

    const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapper(grid);
    
    std::vector<double> leafVec{};
    leafVec.resize(leafView.size(0));
    for (const auto& element : Dune::elements(leafView)) {
        leafVec[element.index()] = cartMapper.cartesianIndex(element.index())*0.1;
    }

    leafSolution.insert("FAKEPROP",
                        Opm::UnitSystem::measure::liquid_surface_volume, // just one possible meassure
                        std::move(leafVec),
                        Opm::data::TargetType::RESTART_OPM_EXTENDED);    // just one possible TargetType

    BOOST_CHECK(leafSolution.has("FAKEPROP"));

    const auto& leafFakePropData = leafSolution.data<double>("FAKEPROP");
    BOOST_CHECK_EQUAL(leafFakePropData.size(), leafView.size(0));

    std::vector<Opm::data::Solution> levelSolutions{};
    Opm::Lgr::extractSolutionLevelGrids(grid,
                                        leafSolution,
                                        levelSolutions);

    // By now, all level cells have data assigned.
    // Notice that the cells that vanished (i.e. do not appear on the leaf grid view,
    // they were parent cells) got the value "rubbish = -1"
    for (int level = 0; level <= grid.maxLevel(); ++level) {
        for (const auto& element : Dune::elements(grid.levelGridView(level))) {
            if (element.isLeaf())
                continue;
            BOOST_CHECK_EQUAL(levelSolutions[level].data<double>("FAKEPROP")[element.index()], -1);
        }
    }
     
    const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
    for (int level = 0; level <= grid.maxLevel(); ++level) {

        BOOST_CHECK(levelSolutions[level].has("FAKEPROP"));

        const auto& levelFakePropData = levelSolutions[level].data<double>("FAKEPROP");
        BOOST_CHECK_EQUAL(levelFakePropData.size(), grid.levelGridView(level).size(0));
        BOOST_CHECK_EQUAL_COLLECTIONS(levelFakePropData.begin(), levelFakePropData.end(),
                                      expected_data_levels[level].begin(), expected_data_levels[level].end());
    }

    Opm::data::Wells dummyWells{};
    Opm::data::GroupAndNetworkValues dummyGroupAndNetworkValues{};
    Opm::data::Aquifers dummyAquifer{};
    Opm::RestartValue leafRestartValue(leafSolution, dummyWells, dummyGroupAndNetworkValues, dummyAquifer);
    std::vector<Opm::RestartValue> restartValue_levels{};
    Opm::Lgr::extractRestartValueLevelGrids<Dune::CpGrid>(grid, leafRestartValue, restartValue_levels);
}

BOOST_AUTO_TEST_CASE(restrictDataForNonNestedLgrsSharingEdges)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.,1.,1.});
    //                          LGR1 parent cells          LGR2 parent cells
    // k = 2   32 33 34 35 |
    //         28 29 30 31 |          29 30
    //         24 25 26 27 |          25 26
    // --------------------|    ------------------       ------------------
    // k = 1   20 21 22 23 |
    //         16 17 18 19 |          17 18
    //         12 13 14 15 |          13 14
    // --------------------|    ------------------       ------------------
    // k = 0    8  9 10 11 |                                   9 10
    //          4  5  6  7 |
    //          0  1  2  3 |
    //---------------------|
    // Refine a few cells, each into 2x2x2 children (2 subdivisions per x-,y-,z- direction).
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}, {2,3,2}},
                               /* startIJK_vec = */ {{1,0,1}, {1,2,0}},
                               /* endIJK_vec = */ {{3,2,3}, {3,3,1}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2"});

    std::vector<std::vector<double>> expected_data_levels{};
    expected_data_levels.resize(grid.maxLevel()+1);

    // Expected data is cartesianIndex/10, expect for parent cells.
    expected_data_levels[0] = std::vector<double>(grid.levelGridView(0).size(0));
    std::iota(expected_data_levels[0].begin(), expected_data_levels[0].end(), 0);

    expected_data_levels[1] = {13,13,14,14,13,13,14,14,  // layer 0
                               17,17,18,18,17,17,18,18,
                               13,13,14,14,13,13,14,14,
                               17,17,18,18,17,17,18,18,
                               25,25,26,26,25,25,26,26,  // layer 1
                               29,29,30,30,29,29,30,30,
                               25,25,26,26,25,25,26,26,
                               29,29,30,30,29,29,30,30};
    expected_data_levels[2] = {  9, 9, 10, 10, 9, 9, 10, 10, 9, 9, 10, 10,  // layer 0
                                 9, 9, 10, 10, 9, 9, 10, 10, 9, 9, 10, 10}; // layer 1
     for (int level = 0; level <= grid.maxLevel(); ++level) {
        for (auto& v : expected_data_levels[level]){
            v *= 0.1;
        }
    }
    // Parent cells that get "rubbish = -1"
    for (const auto& idx : std::vector<int>{9,10,13,14,17,18,25,26,29,30}) {
        expected_data_levels[0][idx] = -1;
    }

    restrictFakeLeafDataToLevelGrids(grid, expected_data_levels);
}

BOOST_AUTO_TEST_CASE(restrictDataForNestedRefinementOnly) {

    const std::vector<std::array<int,3>>  cells_per_dim_vec = {{3,2,2}, {2,2,2}, {2,2,1}, {3,4,2}};
    const std::vector<std::array<int,3>>       startIJK_vec = {{1,1,0}, {1,1,0}, {1,1,0}, {1,1,0}};
    const std::vector<std::array<int,3>>         endIJK_vec = {{2,2,1}, {2,2,1}, {2,2,1}, {2,2,1}};
    const std::vector<std::string>             lgr_name_vec = { "LGR1",  "LGR2",  "LGR3",  "LGR4"};
    const std::vector<std::string> lgr_parent_grid_name_vec = {"GLOBAL", "LGR1",  "LGR2",  "LGR3"};

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

    std::vector<std::vector<double>> expected_data_levels{};
    expected_data_levels.resize(grid.maxLevel()+1);

    // Expected data is cartesianIndex/10, expect for parent cells.
    expected_data_levels[0] = std::vector<double>(grid.levelGridView(0).size(0));
    std::iota(expected_data_levels[0].begin(), expected_data_levels[0].end(), 0);
    for (int level = 0; level <= grid.maxLevel(); ++level) {
        if (level) {
            // cartesian index parent cell/origin cell  = 4
            expected_data_levels[level] =  std::vector<double>(grid.levelGridView(level).size(0), 4);
        }
        for (auto& v : expected_data_levels[level]){
            v *= 0.1;
        }
    }
    // Parent cells get "rubbish = -1"
    expected_data_levels[0][4] = -1; // parent index LGR1 = {4}
    expected_data_levels[1][4] = -1; // {level 1, level element index 4} parent cell (children in LGR2)
    expected_data_levels[2][3] = -1; // {level 2, level element index 3} parent cell (children in LGR3)
    expected_data_levels[3][3] = -1; // {level 3, level element index 3} parent cell (children in LGR4)

    restrictFakeLeafDataToLevelGrids(grid, expected_data_levels);
}

BOOST_AUTO_TEST_CASE(restrictDataForMixNameOrderAndNestedRefinement){

    const std::vector<std::array<int,3>>  cells_per_dim_vec = { {3,2,2}, {2,2,2}, {2,2,1}, {3,4,2}};
    const std::vector<std::array<int,3>>       startIJK_vec = { {1,1,0}, {0,0,0}, {0,0,0}, {1,1,0}};
    const std::vector<std::array<int,3>>         endIJK_vec = { {2,2,1}, {1,1,1}, {1,1,1}, {2,2,1}};
    const std::vector<std::string>             lgr_name_vec = {  "LGR1", "LGR2",   "LGR3", "LGR4"};
    const std::vector<std::string> lgr_parent_grid_name_vec = {"GLOBAL", "LGR1", "GLOBAL", "LGR3"};

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

    std::vector<std::vector<double>> expected_data_levels{};
    expected_data_levels.resize(grid.maxLevel()+1);

    /** level 1 contains LGR1, whose parent grid is GLOBAL
        level 2 contains LGR3, whose parent grid is GLOBAL
        level 3 contains LGR2, whose parent grid is LGR1
        level 4 contains LGR4, whose parent grid is LGR3   */

    // Expected data is cartesianIndex/10, expect for parent cells.
    expected_data_levels[0] = std::vector<double>(grid.levelGridView(0).size(0));
    std::iota(expected_data_levels[0].begin(), expected_data_levels[0].end(), 0);

    // cartesian index parent cell/origin cell  = 4
    expected_data_levels[1] =  std::vector<double>(grid.levelGridView(1).size(0), 4);    // LGR1
    // cartesian index parent cell/origin cell  = 0
    expected_data_levels[2] =  std::vector<double>(grid.levelGridView(2).size(0));       // LGR3
    // cartesian index parent cell/origin cell  = 4
    expected_data_levels[3] =  std::vector<double>(grid.levelGridView(3).size(0), 4);  // LGR2
    // cartesian index parent cell/origin cell  = 0
    expected_data_levels[4] =  std::vector<double>(grid.levelGridView(4).size(0)); // LGR4

    for (int level = 0; level <= grid.maxLevel(); ++level) {
        for (auto& v : expected_data_levels[level]){
            v *= 0.1;
        }
    }
    // Parent cells that get "rubbish = -1"
    expected_data_levels[0][4] = -1; // parent index LGR1 = {4}
    expected_data_levels[0][0] = -1; // parent index LGR3 = {0}
    expected_data_levels[1][0] = -1; // {level 1, level element index 0} parent cell (children in LGR2)
    expected_data_levels[2][3] = -1; // {level 2, level element index 3} parent cell (children in LGR4)

    restrictFakeLeafDataToLevelGrids(grid, expected_data_levels);
}
