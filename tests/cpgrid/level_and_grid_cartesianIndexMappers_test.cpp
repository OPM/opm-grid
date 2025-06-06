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

#define BOOST_TEST_MODULE LevelAndGridCartesianIndexMappersTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/cpgrid/CartesianIndexMapper.hpp>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>
#include <tests/cpgrid/LgrChecks.hpp>


#include <array>
#include <unordered_set>
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

// This test reuses in each case the same grid and LGRs, to check
// serial and parallel bahavior. The difference is how refinement
// gets trigered, namemly, by calling addLgrsUpdateLeafView(...),
// adapt(), or globalRefine(..).
BOOST_AUTO_TEST_CASE(levalAndGridCartesianIndexMappers_afterAddLgrsUpdateLeafView)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    bool isParallel = grid.comm().size() > 1;
    if (isParallel) {
        grid.loadBalance();
    }
    grid.addLgrsUpdateLeafView( /* cells_per_dim = */ {{3,3,3}, {3,3,3}},
                                /* startIJK_vec = */ {{0,0,1}, {2,2,2}},
                                /* endIJK_vec = */ {{3,2,2}, {4,3,3}},
                                /* lgr_name_vec = */ {"LGR1", "LGR2"});

    // Take all the level grids and the leaf grid view - Remove? Does it help readability?
    const auto& level_0_grid_data = grid.currentData().front();
    const auto& level_1_grid_data = grid.currentData()[1];
    const auto& level_2_grid_data = grid.currentData()[2];
    const auto& leaf_grid_data = grid.currentData().back();
    const int unexisting_level = 3;

    const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapp(grid);
    const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
    // I would like to create Opm::LeafCartesianIndexMapper<Dune::CpGrid> leafCartMapp(grid);

    // ---------------------------------------------------------- cartesianDimensions()
    areEqual(/* levelZeroCartesianDimensions = */ {4,3,3}, level_0_grid_data->logicalCartesianSize());
    areEqual(level_0_grid_data->logicalCartesianSize(), cartMapp.cartesianDimensions());
    areEqual(level_0_grid_data->logicalCartesianSize(), levelCartMapp.cartesianDimensions(0));
    areEqual(cartMapp.cartesianDimensions(), levelCartMapp.cartesianDimensions(0));
    areEqual(grid.logicalCartesianSize(), cartMapp.cartesianDimensions());
    
    // Block shaped parent cells of LGR1 dimensions (3-0)x(2-0)x(2-1). Number of subdivisions per cell, per direction {3,3,3}.
    areEqual(/* LGR1 dimensions {(3-0)*3, (2-0)*3, (2-1)*3} */ {9,6,3}, level_1_grid_data->logicalCartesianSize());
    areEqual(level_1_grid_data->logicalCartesianSize(), levelCartMapp.cartesianDimensions(1));

    // Block shaped parent cells of LGR2 dimensions (4-2)x(3-2)x(3-2). Number of subdivisions per cell, per direction {3,3,3}.
    areEqual(/* LGR2 dimensions {(4-2)*3, (3-2)*3, (3-2)*3} */ {6,3,3}, level_2_grid_data->logicalCartesianSize());
    areEqual(level_2_grid_data->logicalCartesianSize(), levelCartMapp.cartesianDimensions(2));

    BOOST_CHECK_THROW(levelCartMapp.cartesianDimensions(unexisting_level), std::logic_error);

    // ---------------------------------------------------------- cartesianSize()
    BOOST_CHECK_EQUAL(cartMapp.cartesianSize(), /* level zero Cartesian size */ 4*3*3);
    BOOST_CHECK_EQUAL(cartMapp.cartesianSize(), level_0_grid_data->size(0)); // Unfortunately, not grid.size(0)
    
    BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(0), /* level zero Cartesian size */ 4*3*3);
    BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(0), level_0_grid_data->size(0)); // BECAUSE ALL ACTIVE CELLS
     
    BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(1), /* level 1 Cartesian size */ 9*6*3);
    BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(1), level_1_grid_data->size(0)); // BECAUSE ALL ACTIVE CELLS
    
    BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(2), /* level 2 Cartesian size */ 6*3*3);
    BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(2), level_2_grid_data->size(0)); // BECAUSE ALL ACTIVE CELLS
    
    BOOST_CHECK_THROW(levelCartMapp.cartesianSize(unexisting_level), std::logic_error);

    // ---------------------------------------------------------- compressedSize()
    BOOST_CHECK_EQUAL(cartMapp.compressedSize(), 244); // 244 = 4*3*3 - (parent-cells (3*2*1) + (2*1*1)) + (9*6*3) + (6*3*3)
    BOOST_CHECK_EQUAL(levelCartMapp.compressedSize(0), 4*3*3);
    BOOST_CHECK_EQUAL(levelCartMapp.compressedSize(1), 9*6*3);
    BOOST_CHECK_EQUAL(levelCartMapp.compressedSize(2), 6*3*3);
    BOOST_CHECK_THROW(levelCartMapp.compressedSize(unexisting_level), std::logic_error);

    // ---------------------------------------------------------- cartesianIndex()
    const int dummy_compressedElementIndex = 0;
    BOOST_CHECK_EQUAL(cartMapp.cartesianIndex(dummy_compressedElementIndex), 0); //To be improved
    BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(dummy_compressedElementIndex, 0), 0);
    BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(dummy_compressedElementIndex, 1), 0);
    BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(dummy_compressedElementIndex, 2), 0);
    //BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(dummy_compressedElementIndex, 0), 35);
    //BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(dummy_compressedElementIndex, 1), 161);
    // BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(dummy_compressedElementIndex, 2), 53);
    BOOST_CHECK_THROW(levelCartMapp.cartesianIndex(dummy_compressedElementIndex, unexisting_level), std::logic_error);

    const int dummy_compressedElementIndexOnLevel = 0;
    std::array<int,3> dummy_coordsOnLevel{};
    BOOST_CHECK_THROW(levelCartMapp.cartesianCoordinate(dummy_compressedElementIndexOnLevel,
                                                        dummy_coordsOnLevel,
                                                        unexisting_level), std::logic_error);
}

BOOST_AUTO_TEST_CASE(gridLogCartSize_afterStrictLocalRefinementWith_addLgrsUpdateLeafView_isACopyOfLevelZeroLogCartSize)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    bool isParallel = grid.comm().size() > 1;
    if (isParallel) {
        grid.loadBalance();
    }
    grid.addLgrsUpdateLeafView(/* cells_per_dim = */ {{3,3,3}, {3,3,3}},
                               /* startIJK_vec = */ {{0,0,1}, {2,2,2}},
                               /* endIJK_vec = */ {{3,2,2}, {4,3,3}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2"});

    areEqual(/* grid dimensions before refinement = */ {4,3,3},
             /* level 0 logicalCartesianSize = */ grid.currentData().front()->logicalCartesianSize());

    areEqual(/* level 0 logicalCartesianSize = */ grid.currentData().front()->logicalCartesianSize(),
             grid.logicalCartesianSize());
}

BOOST_AUTO_TEST_CASE(gridLogCartSize_afterHiddenGlobalRefinementWith_addLgrsUpdateLeafView_makesSense)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    bool isParallel = grid.comm().size() > 1;
    if (isParallel) {
        grid.loadBalance();
    }
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

BOOST_AUTO_TEST_CASE(lgrAndGridLogCartSize_afterStrictLocalRefinementWith_adapt_areACopyOfLevelZeroLogCartSize)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    bool isParallel = grid.comm().size() > 1;
    if (isParallel){
        grid.loadBalance();
    }

    std::unordered_set<int> markedCells = {17,18,21,22}; // parent cell global ids
    // Mark selected elements for refinement
    for (const auto& element : elements(grid.leafGridView())) {
        const auto& id = grid.globalIdSet().id(element);
        if (markedCells.count(id) > 0) {
            grid.mark(1, element);
        }
    }

    grid.preAdapt();
    grid.adapt(); // Default subdivisions per cell 2x2x2 in x-,y-, and z-direction.
    grid.postAdapt();

    areEqual(/* grid dimensions before refinement = */ {4,3,3},
             /* level 0 logicalCartesianSize = */ grid.currentData().front()->logicalCartesianSize());

    // Even though the marked cells form a 2x2x1 block, the logicalCartesianSize of LGR1 is NOT {4,4,2}.
    areEqual(/* level 0 logicalCartesianSize = */ grid.currentData().front()->logicalCartesianSize(), // {4,3,3}
             /* LGR1 logicalCartesianSize = */  grid.currentData()[1]->logicalCartesianSize());

    areEqual(/* level 0 logicalCartesianSize = */ grid.currentData().front()->logicalCartesianSize(),
             grid.logicalCartesianSize());
}

BOOST_AUTO_TEST_CASE(lgrAndGridLogCartSize_afterHiddenGlobalRefinementWith_adapt_makeSense)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    areEqual(/* grid dimensions before refinement = */ {4,3,3},
             /* level 0 logicalCartesianSize = */ grid.currentData().front()->logicalCartesianSize());

    bool isParallel = grid.comm().size() > 1;
    if (isParallel) {
        grid.loadBalance();
    }
    // Mark all elements -> 'indirect' globalRefine
    for (const auto& element : elements(grid.leafGridView())) {
        grid.mark(1, element);
    }
    grid.preAdapt();
    grid.adapt(); // Default subdivisions per cell 2x2x2 in x-,y-, and z-direction.
    grid.postAdapt();

    areEqual(/* expected logicalCartesianSize = */ {4*2, 3*2, 3*2},
             /* LGR1 logicalCartesianSize = */ grid.currentData()[1]->logicalCartesianSize());

    areEqual(/* expected logicalCartesianSize = */ {4*2, 3*2, 3*2},
             grid.logicalCartesianSize());
}

BOOST_AUTO_TEST_CASE(lgrAndGridLogCartSize_after_globalRefine_makeSense)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    bool isParallel = grid.comm().size() > 1;
    if (isParallel) {
        grid.loadBalance();
    }
    grid.globalRefine(1); // Default subdivisions per cell 2x2x2 in x-,y-, and z-direction.

    areEqual(/* expected logicalCartesianSize = */ {4*2, 3*2, 3*2},
             grid.logicalCartesianSize());
    // The refined level grid is a "copy" of the leaf grid view, if globalRefine has been invoked.
    // TODO: remove the refined level grid in this case.
    areEqual(/* expected logicalCartesianSize = */ {4*2, 3*2, 3*2},
             grid.currentData()[1]->logicalCartesianSize());
}
