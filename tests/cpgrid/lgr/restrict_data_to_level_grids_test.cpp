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
#include <opm/grid/cpgrid/CpGridData.hpp>

#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>
#include <opm/grid/cpgrid/LgrHelpers.hpp>

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>
#include <opm/output/data/Solution.hpp>
#include <opm/output/data/Cells.hpp>

#include <dune/grid/common/mcmgmapper.hh>

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

template <typename Scalar>
std::vector<Opm::data::Solution> extractDataForLevelGrids(const Dune::CpGrid& grid,
                                                          const Opm::data::Solution& leafSolution)
{
    int maxLevel = grid.maxLevel();

    // To restrict/create the level cell data, based on the leaf cells and the hierarchy
    std::vector<Opm::data::Solution> levelSolutions{};
    levelSolutions.resize(maxLevel+1);

    for (const auto& [name, leafCellData] : leafSolution)
    {
        const auto& leafVector = leafCellData.template data<Scalar>();
        const auto& messure =  leafCellData.dim;  // Opm::UnitSystem::measure;
        const auto& target = leafCellData.target; // Opm::data::TargetType>;

        if (leafVector.empty()) {
            continue;
        }

        std::vector<std::vector<Scalar>> levelVectors{};
        levelVectors.resize(maxLevel+1);

        const auto rubbish = -1;
        for (int level = 0; level <= maxLevel; ++level) {
            levelVectors[level].resize(grid.levelGridView(level).size(0), rubbish);
        }

        // For level cells that appear in the leaf, extract the data value from leafVector
        // and assign it the the equivalent level cell.
        // Notice that cells that vanished (parent cells) get the rubbish value.
        for (const auto& element : Dune::elements(grid.leafGridView())) {
            levelVectors[element.level()][element.getLevelElem().index()] = leafVector[element.index()];
        }

        // Reorder the containers in the order expected by outout files (increasing level Cartesian indices)
        const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
        for (int level = 0; level <= maxLevel; ++level) {

            const auto toOutput = Opm::Lgr::mapLevelIndicesToCartesianOutputOrder(grid, levelCartMapp, level);
            auto outputContainer = Opm::Lgr::reorderForOutput( levelVectors[level],
                                                               toOutput);

            levelSolutions[level].insert(name,
                                         messure,
                                         std::move(outputContainer),
                                         target);
        }
    }
    return levelSolutions;
}

void restrictFakeLeafDataToLevelGrids(const Dune::CpGrid& grid,
                                      const std::vector<std::vector<double>>& expected_data_levels)
{
    // Create a fake leaf solution
    Opm::data::Solution leafSolution{};
    BOOST_CHECK(!leafSolution.has("NOTHING")); // dummy check

    const auto& leafView = grid.leafGridView();

    using LeafMapper = Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LeafGridView>;
    const LeafMapper leafMapper(leafView, Dune::mcmgElementLayout());

    std::vector<double> leafVec{};
    leafVec.resize(leafView.size(0));
    for (const auto& element : Dune::elements(leafView)) {
        leafVec[leafMapper.index(element)] = element.hasFather()? 1 : 0;
    }
    // leafVec is a vector of 0's and 1's.
    // Set 1 for those cells that have a father, 0 otherwise.

    leafSolution.insert("FAKEPROP",
                        Opm::UnitSystem::measure::liquid_surface_volume, // just one possible meassure
                        std::move(leafVec),
                        Opm::data::TargetType::RESTART_OPM_EXTENDED);    // just one possible TargetType

    BOOST_CHECK(leafSolution.has("FAKEPROP"));

    const auto& leafFakePropData = leafSolution.data<double>("FAKEPROP");
    BOOST_CHECK_EQUAL(leafFakePropData.size(), leafView.size(0));

    const auto levelSolutions = extractDataForLevelGrids<double>(grid, leafSolution);

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

    // Expected data for level grids is 1 for parent and children cells, 0 otherwise.

    // level 0 data vector size 36, with all entries 0 except for parent cells that get "rubbish = -1"
    expected_data_levels[0] = std::vector<double>(36,0);
    for (const auto& idx : std::vector<int>{9,10,13,14,17,18,25,26,29,30}) {
        expected_data_levels[0][idx] = -1;
    }

    // level 1 data vector of 1's, size 64 (cells_per_dim = {2,2,2}, block parent cells dim = {2,2,2})
    expected_data_levels[1] = std::vector<double>(64,1);

    // level 2 data vector of 1's, size 24 (cells_per_dim = {2,3,2}, block parent cells dim = {2,1,1})
    expected_data_levels[2] = std::vector<double>(24,1);

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

    // Expected data for level grids is 1 for parent and children cells, 0 otherwise.

    // level 0 data vector size 36, with all entries 0 except for parent cells, that get "rubbish = -1"
    expected_data_levels[0] = std::vector<double>(9,0);
    expected_data_levels[0][4] = -1; // parent index LGR1 = {4}

    // level 1 data vector size 12 (cells_per_dim = {3,2,2}, block parent cells dim = {1,1,1}).
    // All cells have a parent and only cell with level 1 index 4 is as well a parent cell, refined
    // in LGR2. Therefore, cell with level 1 index 4 gets "rubbish = -1"
    expected_data_levels[1] = std::vector<double>(12,1);
    expected_data_levels[1][4] = -1;

    // level 2 data vector size 8 (cells_per_dim = {2,2,2}, block parent cells dim = {1,1,1})
    // All cells have a parent and only cell with level 2 index 3 is as well a parent cell, refined
    // in LGR3. Therefore, cell with level 2 index 3 gets "rubbish = -1"
    expected_data_levels[2] = std::vector<double>(8,1);
    expected_data_levels[2][3] = -1;

    // level 3 data vector size 4 (cells_per_dim = {2,2,1}, block parent cells dim = {1,1,1})
    // All cells have a parent and only cell with level 3 index 3 is as well a parent cell, refined
    // in LGR4. Therefore, cell with level 3 index 3 gets "rubbish = -1"
    expected_data_levels[3] = std::vector<double>(4,1);
    expected_data_levels[3][3] = -1;

    // level 4 data vector of 1's, size 24 (cells_per_dim = {3,4,2}, block parent cells dim = {1,1,1})
    expected_data_levels[4] = std::vector<double>(24,1);

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

    // Expected data for level grids is 1 for parent and children cells, 0 otherwise.

    // level 0 data vector size 36, with all entries 0 except for parent cells that get "rubbish = -1"
    expected_data_levels[0] = std::vector<double>(9,0);
    expected_data_levels[0][4] = -1; // parent index LGR1 = {4}
    expected_data_levels[0][0] = -1; // parent index LGR3 = {0}
    /** level 1 contains LGR1, whose parent grid is GLOBAL
        level 2 contains LGR3, whose parent grid is GLOBAL
        level 3 contains LGR2, whose parent grid is LGR1
        level 4 contains LGR4, whose parent grid is LGR3   */

    // level 1 (LGR1) data vector size 12 (cells_per_dim = {3,2,2}, block parent cells dim = {1,1,1})
    // All cells have a parent and only cell with level 1 index 0 is as well a parent cell, refined
    // in LGR2. Therefore, cell with level 1 index 0 gets "rubbish = -1"
    expected_data_levels[1] = std::vector<double>(12,1);
    expected_data_levels[1][0] = -1;

    // level 2 (LGR3) data vector of size 4 (cells_per_dim = {2,2,1}, block parent cells dim = {1,1,1})
    // All cells have a parent and only cell with level 2 index 3 is as well a parent cell, refined
    // in LGR4. Therefore, cell with level 2 index 3 gets "rubbish = -1"
    expected_data_levels[2] = std::vector<double>(4,1);
    expected_data_levels[2][3] = -1;

    // level 3 (LGR2) data vector of 1's, size 8 (cells_per_dim = {2,2,2}, block parent cells dim = {1,1,1})
    expected_data_levels[3] = std::vector<double>(8,1);

    // level 4 (LGR4) data vector of 1's, size 24 (cells_per_dim = {3,4,2}, block parent cells dim = {1,1,1})
    expected_data_levels[4] = std::vector<double>(24,1);

    restrictFakeLeafDataToLevelGrids(grid, expected_data_levels);
}
