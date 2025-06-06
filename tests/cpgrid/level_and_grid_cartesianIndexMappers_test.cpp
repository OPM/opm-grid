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

void cartesianDimensionsCoincideWithLevelZeroOnes(const Dune::CpGrid& grid,
                                                  const Dune::CartesianIndexMapper<Dune::CpGrid>& cartMapp,
                                                  const Opm::LevelCartesianIndexMapper<Dune::CpGrid>& levelCartMapp,
                                                  const std::array<int,3>& levelZero_cartDims)
{
    Opm::areEqual(cartMapp.cartesianDimensions(), levelZero_cartDims);
    Opm::areEqual(cartMapp.cartesianDimensions(), levelCartMapp.cartesianDimensions(0));
    Opm::areEqual(cartMapp.cartesianDimensions(), grid.currentData().front()->logicalCartesianSize());
}

void cartesianSizeCoincidesWithLevelZeroOne(const Dune::CpGrid& grid,
                                            const Dune::CartesianIndexMapper<Dune::CpGrid>& cartMapp,
                                            const Opm::LevelCartesianIndexMapper<Dune::CpGrid>& levelCartMapp,
                                            int levelZero_cartSize)
{
    BOOST_CHECK_EQUAL(cartMapp.cartesianSize(), levelZero_cartSize);
    BOOST_CHECK_EQUAL(cartMapp.cartesianSize(), levelCartMapp.cartesianSize(0));
    BOOST_CHECK_EQUAL(cartMapp.cartesianSize(), grid.currentData().front()->size(0)); // only for ALL active & serial?
}


void compressedSizeCoincidesWithLeafGridCellCount(const Dune::CpGrid& grid,
                                                  const Dune::CartesianIndexMapper<Dune::CpGrid>& cartMapp,
                                                  int leafGridCellCount)
{
    BOOST_CHECK_EQUAL(cartMapp.compressedSize(), leafGridCellCount);
    BOOST_CHECK_EQUAL(cartMapp.compressedSize(), grid.size(0));
    BOOST_CHECK_EQUAL(cartMapp.compressedSize(), grid.currentData().back()->size(0)); // only for ALL active & serial?
}

void checkLevels(const Dune::CpGrid& grid,
                 const Opm::LevelCartesianIndexMapper<Dune::CpGrid>& levelCartMapp,
                 const std::vector<std::array<int,3>>& expected_level_cartDims,
                 const std::vector<int>& expected_level_cartSizes,
                 const std::vector<int>& expected_level_compressedSizes)
{
    for (int level = 0; level <= grid.maxLevel(); ++level) {
        Opm::areEqual( levelCartMapp.cartesianDimensions(level), expected_level_cartDims[level] );
        BOOST_CHECK_EQUAL( levelCartMapp.cartesianSize(level), expected_level_cartSizes[level] );
        BOOST_CHECK_EQUAL( levelCartMapp.compressedSize(level), expected_level_compressedSizes[level] );
    }
}

void allThrow(const Opm::LevelCartesianIndexMapper<Dune::CpGrid>& levelCartMapp,
              int unexisting_level)
{
    BOOST_CHECK_THROW(levelCartMapp.cartesianDimensions(unexisting_level), std::logic_error);
    BOOST_CHECK_THROW(levelCartMapp.cartesianSize(unexisting_level), std::logic_error);
    BOOST_CHECK_THROW(levelCartMapp.compressedSize(unexisting_level), std::logic_error);

    const int dummy_compressedElementIndex = 0;
    BOOST_CHECK_THROW(levelCartMapp.cartesianIndex(dummy_compressedElementIndex, unexisting_level), std::logic_error);

    const int dummy_compressedElementIndexOnLevel = 0;
    std::array<int,3> dummy_coordsOnLevel{};
    BOOST_CHECK_THROW(levelCartMapp.cartesianCoordinate(dummy_compressedElementIndexOnLevel,
                                                        dummy_coordsOnLevel,
                                                        unexisting_level), std::logic_error);
}

// This test reuses in each case the same grid and LGRs, to check
// serial and parallel bahavior. The difference is how refinement
// gets trigered, namemly, by calling addLgrsUpdateLeafView(...),
// adapt(), or globalRefine(..).
BOOST_AUTO_TEST_CASE(level_and_grid_cartesianIndexMapper_afterStrictLocalRefinementWith_addLgrsUpdateLeafView)
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

    // Block shaped parent cells of LGR1 dimensions (3-0)x(2-0)x(2-1).
    // Number of subdivisions per cell, per direction {3,3,3}.       -> LGR1 dims = {9,6,3}
    //
    // Block shaped parent cells of LGR2 dimensions (4-2)x(3-2)x(3-2).
    // Number of subdivisions per cell, per direction {3,3,3}.       -> LGR2 dims = {6,3,3}

    const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapp(grid);
    const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
    // I would like to create Opm::LeafCartesianIndexMapper<Dune::CpGrid> leafCartMapp(grid);

    // Check level Cartesian dimensions and size, and compressed size, for level grids.
    checkLevels(grid,
                levelCartMapp,
                {{4,3,3}, {9,6,3}, {6,3,3}}, // expected Cartesian dimensions per level grid
                {4*3*3, 9*6*3, 6*3*3},       // expected Cartesian sizes per level grid
                {4*3*3, 9*6*3, 6*3*3});      // expected compressed sizes per level grid

    // In this case, all parent cells are active, therefore level Cartesian size
    // and compressed size coincide.
    for (int level = 0; level <= grid.maxLevel(); ++level) {
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(level), levelCartMapp.compressedSize(level));
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(level), grid.currentData()[level]->size(0));
    }

    /*// ---------------------------------------------------------- cartesianIndex() TO BE IMPROVED
      const int dummy_compressedElementIndex = 0;
      const int dummy_compressedElementIndexOnLevel = 0;
      std::array<int,3> dummy_coordsOnLevel{};

      BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(dummy_compressedElementIndex, 0), 0);
      BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(dummy_compressedElementIndex, 1), 0);
      BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(dummy_compressedElementIndex, 2), 0);
      //BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(dummy_compressedElementIndex, 0), 35);
      //BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(dummy_compressedElementIndex, 1), 161);
      // BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(dummy_compressedElementIndex, 2), 53);
      levelCartMapp.cartesianCoordinate(dummy_compressedElementIndexOnLevel, dummy_coordsOnLevel, 0);
      Opm::areEqual(dummy_coordsOnLevel, {0,0,0});
      levelCartMapp.cartesianCoordinate(dummy_compressedElementIndexOnLevel, dummy_coordsOnLevel, 1);
      Opm::areEqual(dummy_coordsOnLevel, {0,0,0});
      levelCartMapp.cartesianCoordinate(dummy_compressedElementIndexOnLevel, dummy_coordsOnLevel, 2);
      Opm::areEqual(dummy_coordsOnLevel, {0,0,0});*/

    allThrow(levelCartMapp, /* unexisting_level = */ grid.maxLevel()+1);

    // CartesianIndexMapper checks - how it relates to LevelCartesianIndexMapper
    cartesianDimensionsCoincideWithLevelZeroOnes(grid,
                                                 cartMapp,
                                                 levelCartMapp,
                                                 /* level zero Cartesian dims = */ {4,3,3});
    // For strict-locally-refined grid, CartesianIndexMapper::cartesianDimensions and
    // grid.logicalCartesianSize() coincide (with level zero Cartesian dimensions/logical Cartesian size).
    Opm::areEqual(cartMapp.cartesianDimensions(), grid.logicalCartesianSize());

    cartesianSizeCoincidesWithLevelZeroOne(grid,
                                           cartMapp,
                                           levelCartMapp,
                                           /* level zero Cartesian size = */ 4*3*3);

    int leafGridCellCount = 244; // 4*3*3 levelZeroCells - [(3*2*1) + (2*1*1) parentCells] + (9*6*3 LGR1Cells) + (6*3*3 LGR2Cells)
    compressedSizeCoincidesWithLeafGridCellCount(grid, cartMapp, leafGridCellCount);

    // CartesianIndexMapper -cartesianIndex(...) TO BE IMPROVED
    //BOOST_CHECK_EQUAL(cartMapp.cartesianIndex(dummy_compressedElementIndex), 0); //To be improved
    // CartesianIndexMapper - cartesianCoordinate(...)
    //  cartMapp.cartesianCoordinate(dummy_compressedElementIndexOnLevel, dummy_coordsOnLevel);
    //Opm::areEqual(dummy_coordsOnLevel, {0,0,0});
}

BOOST_AUTO_TEST_CASE(level_and_grid_cartesianIndexMapper_afterHiddenGlobalRefinementWith_addLgrsUpdateLeafView)
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

    // Block shaped parent cells of LGR1 is the entire level zero grid, dimensions (4-0)x(3-0)x(3-1).
    // Number of subdivisions per cell, per direction {3,3,3}.  -> LGR1 dims = {12,9,9}

    const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapp(grid);
    const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
    // I would like to create Opm::LeafCartesianIndexMapper<Dune::CpGrid> leafCartMapp(grid);

    // Check level Cartesian dimensions and size, and compressed size, for level grids.
    checkLevels(grid,
                levelCartMapp,
                {{4,3,3}, {12,9,9}}, // expected Cartesian dimensions per level grid
                {4*3*3, 12*9*9},     // expected Cartesian sizes per level grid
                {4*3*3, 12*9*9});    // expected compressed sizes per level grid

    // In this case, all parent cells are active, therefore level Cartesian size
    // and compressed size coincide.
    for (int level = 0; level <= grid.maxLevel(); ++level) {
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(level), levelCartMapp.compressedSize(level));
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(level), grid.currentData()[level]->size(0));
    }

    allThrow(levelCartMapp, /* unexisting_level = */ grid.maxLevel()+1);

    cartesianDimensionsCoincideWithLevelZeroOnes(grid,
                                                 cartMapp,
                                                 levelCartMapp,
                                                 /* level zero Cartesian dims = */ {4,3,3});
    // For hidden-globally-refined grid via addLgrsUpdateLeafView, CartesianIndexMapper::cartesianDimensions
    // and grid.logicalCartesianSize() DO NOT coincide
    Opm::areEqual(cartMapp.cartesianDimensions(), {4,3,3} /* level zero grid Cartesian dimensions */);
    Opm::areEqual(grid.logicalCartesianSize(), {12,9,9} /* LEAF grid view Cartesian dimensions */);

    cartesianSizeCoincidesWithLevelZeroOne(grid,
                                           cartMapp,
                                           levelCartMapp,
                                           /* level zero Cartesian size = */ 4*3*3);

    int leafGridCellCount = 972; // 4*3*3 levelZeroCells - [4*3*3 parentCells] + (12*9*9 LGR1Cells)
    compressedSizeCoincidesWithLeafGridCellCount(grid, cartMapp, leafGridCellCount);
}

BOOST_AUTO_TEST_CASE(level_and_grid_cartesianIndexMapper_afterStrictLocalRefinementWith_adapt)
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

    const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapp(grid);
    const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
    // I would like to create Opm::LeafCartesianIndexMapper<Dune::CpGrid> leafCartMapp(grid);

    // Level Cartesian dimensions and size, and compressed size, for refined level grids
    // take the values from level zero grid.
    // Even though the marked cells form a 2x2x1 block, the Cartesian dimensions and
    // Cartesian size of LGR1 is NOT {4,4,2} and 4*4*2 = 32, but {4,3,3} and 4*3*3 respectively.
    checkLevels(grid,
                levelCartMapp,
                {{4,3,3}, {4,3,3}}, // expected Cartesian dimensions per level grid
                {4*3*3, 4*3*3},     // expected Cartesian sizes per level grid
                {4*3*3, 32});    // expected compressed sizes per level grid

    // In this case, all parent cells are active, however, level Cartesian size
    // and compressed size do not coincide for the refined level grid.
    BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(1), /* levelZeroGridCellCount = */ 4*3*3);
    BOOST_CHECK_EQUAL(levelCartMapp.compressedSize(1), /* LGR1CellCount = */ grid.currentData()[1]->size(0));

    allThrow(levelCartMapp, /* unexisting_level = */ grid.maxLevel()+1);

    cartesianDimensionsCoincideWithLevelZeroOnes(grid,
                                                 cartMapp,
                                                 levelCartMapp,
                                                 /* level zero Cartesian dims = */ {4,3,3});
    // For strict-local-refined grid via adapt(), CartesianIndexMapper::cartesianDimensions
    // and grid.logicalCartesianSize() coincide.
    Opm::areEqual(cartMapp.cartesianDimensions(), grid.logicalCartesianSize());

    cartesianSizeCoincidesWithLevelZeroOne(grid,
                                           cartMapp,
                                           levelCartMapp,
                                           /* level zero Cartesian size = */ 4*3*3);

    int leafGridCellCount = 64; // 4*3*3 levelZeroCells - [2*2*1 parentCells] + (4*4*2 LGR1Cells)
    compressedSizeCoincidesWithLeafGridCellCount(grid, cartMapp, leafGridCellCount);
}

BOOST_AUTO_TEST_CASE(level_and_grid_cartesianIndexMapper_afterHiddenGlobalRefinementWith_adapt)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    Opm::areEqual(/* grid dimensions before refinement = */ {4,3,3},
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

    const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapp(grid);
    const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
    // I would like to create Opm::LeafCartesianIndexMapper<Dune::CpGrid> leafCartMapp(grid);

    // Check level Cartesian dimensions and size, and compressed size, for level grids.
    checkLevels(grid,
                levelCartMapp,
                {{4,3,3}, {8,6,6}}, // expected Cartesian dimensions per level grid
                {4*3*3, 8*6*6},     // expected Cartesian sizes per level grid
                {4*3*3, 8*6*6});    // expected compressed sizes per level grid

    // In this case, all parent cells are active, therefore level Cartesian size
    // and compressed size coincide.
    for (int level = 0; level <= grid.maxLevel(); ++level) {
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(level), levelCartMapp.compressedSize(level));
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(level), grid.currentData()[level]->size(0));
    }

    allThrow(levelCartMapp, /* unexisting_level = */ grid.maxLevel()+1);

    cartesianDimensionsCoincideWithLevelZeroOnes(grid,
                                                 cartMapp,
                                                 levelCartMapp,
                                                 /* level zero Cartesian dims = */ {4,3,3});
    // For hidden-globally-refined grid via adapt(), CartesianIndexMapper::cartesianDimensions
    // and grid.logicalCartesianSize() DO NOT coincide
    Opm::areEqual(cartMapp.cartesianDimensions(), {4,3,3} /* level zero grid Cartesian dimensions */);
    Opm::areEqual(grid.logicalCartesianSize(), {8,6,6} /* LEAF grid view Cartesian dimensions */);

    cartesianSizeCoincidesWithLevelZeroOne(grid,
                                           cartMapp,
                                           levelCartMapp,
                                           /* level zero Cartesian size = */ 4*3*3);

    int leafGridCellCount = 288; // 4*3*3 levelZeroCells - [4*3*3 parentCells] + (8*6*6 LGR1Cells)
    compressedSizeCoincidesWithLeafGridCellCount(grid, cartMapp, leafGridCellCount);
}

BOOST_AUTO_TEST_CASE(level_and_grid_cartesianIndexMapper_after_globalRefine)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    bool isParallel = grid.comm().size() > 1;
    if (isParallel) {
        grid.loadBalance();
    }
    grid.globalRefine(1); // Default subdivisions per cell 2x2x2 in x-,y-, and z-direction.

    // The refined level grid is a "copy" of the leaf grid view, if globalRefine has been invoked.
    // TODO: remove the refined level grid in this case.

    const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapp(grid);
    const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
    // I would like to create Opm::LeafCartesianIndexMapper<Dune::CpGrid> leafCartMapp(grid);

    // Check level Cartesian dimensions and size, and compressed size, for level grids.
    checkLevels(grid,
                levelCartMapp,
                {{4,3,3}, {8,6,6}}, // expected Cartesian dimensions per level grid
                {4*3*3, 8*6*6},     // expected Cartesian sizes per level grid
                {4*3*3, 8*6*6});    // expected compressed sizes per level grid

    // In this case, all parent cells are active, therefore level Cartesian size
    // and compressed size coincide.
    for (int level = 0; level <= grid.maxLevel(); ++level) {
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(level), levelCartMapp.compressedSize(level));
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(level), grid.currentData()[level]->size(0));
    }

    allThrow(levelCartMapp, /* unexisting_level = */ grid.maxLevel()+1);

    cartesianDimensionsCoincideWithLevelZeroOnes(grid,
                                                 cartMapp,
                                                 levelCartMapp,
                                                 /* level zero Cartesian dims = */ {4,3,3});
    // For hidden-globally-refined grid via adapt(), CartesianIndexMapper::cartesianDimensions
    // and grid.logicalCartesianSize() DO NOT coincide
    Opm::areEqual(cartMapp.cartesianDimensions(), {4,3,3} /* level zero grid Cartesian dimensions */);
    Opm::areEqual(grid.logicalCartesianSize(), {8,6,6} /* LEAF grid view Cartesian dimensions */);

    cartesianSizeCoincidesWithLevelZeroOne(grid,
                                           cartMapp,
                                           levelCartMapp,
                                           /* level zero Cartesian size = */ 4*3*3);

    int leafGridCellCount = 288; // 4*3*3 levelZeroCells - [4*3*3 parentCells] + (8*6*6 LGR1Cells)
    compressedSizeCoincidesWithLeafGridCellCount(grid, cartMapp, leafGridCellCount);
}
