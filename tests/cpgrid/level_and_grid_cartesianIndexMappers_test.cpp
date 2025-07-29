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

#include <algorithm>
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
    if (grid.comm().size() == 1) { // only for grid with ALL active cells and serial runs
        BOOST_CHECK_EQUAL(cartMapp.cartesianSize(), grid.currentData().front()->size(0));
    }
}

void compressedSizeCoincidesWithLeafGridCellCount(const Dune::CpGrid& grid,
                                                  const Dune::CartesianIndexMapper<Dune::CpGrid>& cartMapp,
                                                  int serialLeafGridCellCount)
{
    if (grid.comm().size() == 1) { // only for grid with ALL active cells and serial runs
        BOOST_CHECK_EQUAL(cartMapp.compressedSize(), serialLeafGridCellCount);
    }
    BOOST_CHECK_EQUAL(cartMapp.compressedSize(), grid.size(0));
    BOOST_CHECK_EQUAL(cartMapp.compressedSize(), grid.currentData().back()->size(0));
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
        if (grid.comm().size() == 1) { // serial run
            BOOST_CHECK_EQUAL( levelCartMapp.compressedSize(level), expected_level_compressedSizes[level] );
        }
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

void levelCartSizeEqualsCompressedSizeIfAllActiveAndSerial(const Dune::CpGrid& grid,
                                                           const Opm::LevelCartesianIndexMapper<Dune::CpGrid>& levelCartMapp,
                                                           bool isParallel)
{
    if (!isParallel) {
        for (int level = 0; level <= grid.maxLevel(); ++level) {
            BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(level), levelCartMapp.compressedSize(level));
            BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(level), grid.currentData()[level]->size(0));
        }
    }
}

void checkDimsAndSizes(const Dune::CpGrid& grid,
                       const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapp,
                       const Opm::LevelCartesianIndexMapper<Dune::CpGrid>& levelCartMapp,
                       const std::vector<std::array<int,3>>& expectedLevelCartDims,
                       const std::vector<int>& expectedLevelCartSizes,
                       const std::vector<int>&  expectedSerialLevelCompressedSizes,
                       int serialLeafGridCellCount)
{
    checkLevels(grid, levelCartMapp, expectedLevelCartDims, expectedLevelCartSizes, expectedSerialLevelCompressedSizes);
    allThrow(levelCartMapp, /* unexisting_level = */ grid.maxLevel()+1);
    cartesianDimensionsCoincideWithLevelZeroOnes(grid, cartMapp, levelCartMapp, expectedLevelCartDims[0]);
    cartesianSizeCoincidesWithLevelZeroOne(grid, cartMapp, levelCartMapp, expectedLevelCartSizes[0]);
    compressedSizeCoincidesWithLeafGridCellCount(grid, cartMapp, serialLeafGridCellCount);
}

void checkSerialMinMaxCompressedIdx(int compressedIndex,
                                    int serialMinCompressedIndex,
                                    int serialMaxCompressedIndex,
                                    bool isParallel)
{
    if (!isParallel) {
        BOOST_CHECK(compressedIndex >= serialMinCompressedIndex);
        BOOST_CHECK(compressedIndex <= serialMaxCompressedIndex);
    }
}

void computeAndCompareCartCoords(const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapp,
                                 int compressedIdx,
                                 const std::array<int,3>& expectedCoords)
{
    std::array<int,3> coords{};
    cartMapp.cartesianCoordinate(compressedIdx, coords);
    Opm::areEqual(coords, expectedCoords);
}

void checkCartMappIdxAndCoords(const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapp,
                               int compressedIndex,
                               int expectedCartesianIndex,
                               const std::array<int,3>& expectedCoords)
{
    BOOST_CHECK_EQUAL(cartMapp.cartesianIndex(compressedIndex), expectedCartesianIndex);
    computeAndCompareCartCoords(cartMapp, compressedIndex, expectedCoords);
}

void computeAndCompareLevelCartCoords(const Opm::LevelCartesianIndexMapper<Dune::CpGrid>& levelCartMapp,
                                      int level,
                                      int levelCompressedIdx,
                                      const std::array<int,3>& expectedLevelCoords)
{
    std::array<int,3> coords{};
    levelCartMapp.cartesianCoordinate(levelCompressedIdx, coords, level);
    Opm::areEqual(coords, expectedLevelCoords);
}

void checkLevelCartMappIdxAndCoords(const Opm::LevelCartesianIndexMapper<Dune::CpGrid>& levelCartMapp,
                                    int levelCompressedIndex,
                                    int expectedLevelCartesianIndex,
                                    const std::array<int,3>& expectedLevelCoords,
                                    int level)
{
    BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(levelCompressedIndex, level), expectedLevelCartesianIndex);
    computeAndCompareLevelCartCoords(levelCartMapp, level, levelCompressedIndex, expectedLevelCoords);
}

void checkLeafElement(const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapp,
                      const Dune::cpgrid::Entity<0>& element,
                      int serialMinCompressedIndex,
                      int serialMaxCompressedIndex,
                      int expectedCartesianIndex,
                      const std::array<int,3>& expectedCoords,
                      bool isParallel)
{
    checkSerialMinMaxCompressedIdx(element.index(), serialMinCompressedIndex, serialMaxCompressedIndex, isParallel);
    checkCartMappIdxAndCoords(cartMapp, element.index(), expectedCartesianIndex, expectedCoords);
}

void checkLevelElement(const Opm::LevelCartesianIndexMapper<Dune::CpGrid>& levelCartMapp,
                       int levelCompressedIndex,
                       int serialLevelCompressedIndex,
                       int expectedCartesianIndex,
                       const std::array<int,3>& expectedCoords,
                       bool isParallel,
                       int level)
{
    checkSerialMinMaxCompressedIdx(levelCompressedIndex, serialLevelCompressedIndex, serialLevelCompressedIndex, isParallel);
    checkLevelCartMappIdxAndCoords(levelCartMapp, levelCompressedIndex, expectedCartesianIndex, expectedCoords, level);
}

void checkAFewCartIdxAndCoordsForGloballyRefinedTestGrids(const Dune::CpGrid& grid,
                                                          const Opm::LevelCartesianIndexMapper<Dune::CpGrid>& levelCartMapp,
                                                          const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapp,
                                                          bool isParallel)
{
    // Level zero grid for all test cases has dimension 4x3x3
    // --------------------------------------------------------------------
    // k = 0  8   9  10  11 |  k = 1 20  21  22  23 | k = 2 32  33  34  35 |
    //        4   5   6   7 |        16  17  18  19 |       28  29  30  31 |
    //        0   1   2   3 |        12  13  14  15 |       24  25  26  27 |

    // Leaf grid local indices (after hidden-global-refinement via addLgrsUpdateLeafView(...),
    // or adapt() with all elements marked, or explicitly calling globalRefine().)
    // "k = 2"  [256-263] [264-271] [272-279] [280-287] |
    //          [224-231] [232-239] [240-247] [248-255] |
    //          [192-199] [200-207] [208-215] [216-223] |
    // ------------------------------------------------
    // "k = 1"  [160-167] [168-175] [176-183] [184-191] |
    //          [128-135] [136-143] [144-151] [152-159] |
    //          [ 96-103] [104-111] [112-119] [120-127] |
    // ------------------------------------------------
    // "k = 0"  [ 64-71 ] [ 72-79 ] [ 80-87 ] [ 88-95 ] |
    //          [ 32-39 ] [ 40-47 ] [ 48-55 ] [ 56-63 ] |
    //          [  0-7  ] [  8-15 ] [ 16-23 ] [ 24-31 ] |
    // --------------------------------------------------

    // In order to test serial and parallel reusing the same grid.
    bool foundId0 = false;
    bool foundId18 = false;
    bool foundId33 = false;
    for (const auto& element : Dune::elements(grid.leafGridView())) {

        const auto& originId = grid.globalIdSet().id(element.getOrigin());

        // Serial:   Leaf refined cells (born in LGR1/GR) with compressedIndex = 0,1,...,7
        //           have parent cell in level zero with Cartesian index 0 and ijk = {0, 0, 0}.
        // Parallel: seach for the element whose origin cell (parent cell in level zero)
        //           has global id equal to 0.
        bool originHasId0 = (originId == 0);

        // Serial:   Leaf refined cells (born in LGR1/GR) with compressedIndex = 144,145,....,151
        //           have parent cell in level zero with Cartesian index 18 and ijk = {2, 1, 1}.
        // Parallel: seach for the element whose origin cell (parent cell in level zero)
        //           has global id equal to 18.
        bool originHasId18 = (originId == 18);

        // Serial:   Leaf refined cells (born in LGR1/GR) with compressedIndex = 264,265,...,271
        //           have parent cell with Cartesian index 33 and ijk = {1, 2, 2}.
        // Parallel: seach for the element whose origin cell (parent cell in level zero)
        //           has global id equal to 33.
        bool originHasId33 = (originId == 33);

        if (originHasId0) {
            checkLeafElement(cartMapp, element, /* serialMinCompressedIndex = */ 0, /* serialMaxCompressedIndex = */ 7,
                             /* expectedCartesianIndex = */ 0, /* expectedCoords = */ {0, 0, 0}, isParallel);
            foundId0 = true;
        }
        else if (originHasId18) {
            checkLeafElement(cartMapp, element, /* serialMinCompressedIndex = */ 144, /* serialMaxCompressedIndex = */ 151,
                             /* expectedCartesianIndex = */ 18, /* expectedCoords = */ {2, 1, 1}, isParallel);
            foundId18 = true;
        }
        else if (originHasId33) {
            checkLeafElement(cartMapp, element, /* serialMinCompressedIndex = */ 264, /* serialMaxCompressedIndex = */ 271,
                             /* expectedCartesianIndex = */ 33, /* expectedCoords = */ {1, 2, 2}, isParallel);
            foundId33 = true;
        }
    }
    BOOST_CHECK(grid.comm().max(foundId0));
    BOOST_CHECK(grid.comm().max(foundId18));
    BOOST_CHECK(grid.comm().max(foundId33));
    
    // For simplicity, illustration only of layers k = 0 and k = 1 on the LGR1 grid,
    // which corresponds to refined cells with parent cells in coarse (level zero) grid
    // parent cells k = 0   |  8 |  9 | 10 | 11 |   (level zero global ids layer k = 0)
    //                      |  4 |  5 |  6 |  7 |
    //                      |  0 |  1 |  2 |  3 |
    // -------------------------------------------------------------------------------------
    // LGR1 local indices                                 LGR1 (level) Cartesian indices
    // k = 1   70  71 | 78  79 | 86  87 | 94  95 |        88  89 | 90  91 | 92  93 | 94  95 |
    //         68  69 | 76  77 | 84  85 | 92  93 |        80  81 | 82  83 | 84  85 | 86  87 |
    //         -----------------------------------        -----------------------------------
    //         38  39 | 46  47 | 54  55 | 62  63 |        72  73 | 74  75 | 76  77 | 78  79 |
    //         36  37 | 44  45 | 52  53 | 60  61 |        64  65 | 66  67 | 68  69 | 70  71 |
    //         -----------------------------------        -----------------------------------
    //          6   7 | 14  15 | 22  23 | 30  31 |        56  57 | 58  59 | 60  61 | 62  63 |
    //          4   5 | 12  13 | 20  21 | 28  29 |        48  49 | 50  51 | 52  53 | 54  55 |
    // -------------------------------------------------------------------------------------
    // k = 0   66  67 | 74  75 | 82  83 | 90  91 |        40  41 | 42  43 | 44  45 | 46  47 |
    //         64  65 | 72  73 | 80  81 | 88  89 |        32  33 | 34  35 | 36  37 | 38  39 |
    //         -----------------------------------        -----------------------------------
    //         34  35 | 42  43 | 50  51 | 58  59 |        24  25 | 26  27 | 28  29 | 30  31 |
    //         32  33 | 40  41 | 48  49 | 56  57 |        16  17 | 18  19 | 20  21 | 22  23 |
    //         -----------------------------------        -----------------------------------
    //          2   3 | 10  11 | 18  19 | 26  27 |         8   9 | 10  11 | 12  13 | 14  15 |
    //          0   1 |  8   9 | 16  17 | 24  25 |         0   1 |  2   3 |  4   5 |  6   7 |

    // In order to test serial and parallel reusing the same grid.
    
    // Serial: {0,1,2,3, 4,5,6,7} local indices -> {0,1,8,9, 48,49,56,57} level Cartesian indices
    std::vector<int> levelCartIndicesOriginId0 = {0,1,8,9, 48,49,56,57};
    bool found_9 = false; // serial local index 3 -> level Cartesian index 9 
    bool found_57 = false; // serial local index 7 -> level Cartesian index 57
    
    // Serial: {48,49,50,51, 52,53,54,55} local indices -> {20,21,28,29, 68,69,76,77} level Cartesian indices
    std::vector<int> levelCartIndicesOriginId6 = {20,21,28,29, 68,69,76,77};
    bool found_28 = false; // serial local index 50 -> level Cartesian index 28 
    bool found_69 = false; // serial local index 53 -> level Cartesian index 69 
    
    // serial: {88,89,90,91, 92,93,94,95} local indices -> {38,39,40,41, 86,87,94,95} level Cartesian indices
    std::vector<int> levelCartIndicesOriginId11 = {38,39,46,47, 86,87,94,95};
    bool found_39 = false; // serial local index 89 -> level Cartesian index 39 
    bool found_86 = false; // serial local index 92 -> level Cartesian index 86 
    
    for (const auto& element : Dune::elements(grid.levelGridView(1))) {

        const auto& originId = grid.globalIdSet().id(element.getOrigin());

        if (originId == 0){ 
            // Check level Cartesian index is contained in {0,1,8,9, 48,49,56,57}
            const auto it = std::find(levelCartIndicesOriginId0.begin(), levelCartIndicesOriginId0.end(),
                                      levelCartMapp.cartesianIndex(element.index(), 1) );
            BOOST_CHECK( it !=  levelCartIndicesOriginId0.end() );

            // Check a few level Cartesian coordinates (e.g., one in layer k =0, and one in layer k = 1).
            // Note: take into account the "serial illustration".
            // 
            // Serial: LGR1 refined cell with lgr1CompressedIndex = 3 has level Cartesian index 9 and ijk = {1, 1, 0}.
            if (*it == 9) {
                checkLevelElement(levelCartMapp, /* levelCompressedIndex = */ element.index(), /* serialLevelCompressedIndex = */ 3,
                                  /* expectedCartesianIndex = */ 9, /* expectedCoords = */ {1, 1, 0}, isParallel, /*level =*/ 1);
                found_9 = true;
            }
            // Serial: LGR1 refined cell with lgr1CompressedIndex = 7 has level Cartesian index 57 and ijk = {1, 1, 1}.
            if (*it == 57) {
                checkLevelElement(levelCartMapp, /* levelCompressedIndex = */ element.index(), /* serialLevelCompressedIndex = */ 7,
                                  /* expectedCartesianIndex = */ 57, /* expectedCoords = */ {1, 1, 1}, isParallel, /*level =*/ 1);
                found_57 = true;
            }
        }
        else if (originId == 6) {
            // Check level Cartesian index is contained in {20,21,28,29, 68,69,76,77}
            const auto it = std::find(levelCartIndicesOriginId6.begin(), levelCartIndicesOriginId6.end(),
                                      levelCartMapp.cartesianIndex(element.index(), 1) );
            BOOST_CHECK( it !=  levelCartIndicesOriginId6.end() );

            // Check a few level Cartesian coordinates (e.g., one in layer k =0, and one in layer k = 1).
            // Note: take into account the "serial illustration".
            //
            // Serial: LGR1 refined cell with lgr1CompressedIndex = 50 has level Cartesian index 28 and ijk = {4, 3, 0}.
            if (*it == 28) {
                checkLevelElement(levelCartMapp, /* levelCompressedIndex = */ element.index(), /* serialLevelCompressedIndex = */ 50,
                                  /* expectedCartesianIndex = */ 28, /* expectedCoords = */ {4, 3, 0}, isParallel, /*level =*/ 1);
                found_28 = true;
            }
            // Serial: LGR1 refined cell with lgr1CompressedIndex = 53 has level Cartesian index 69 and ijk = {5, 2, 1}.
            if (*it == 69) {
                checkLevelElement(levelCartMapp, /* levelCompressedIndex = */ element.index(), /* serialLevelCompressedIndex = */ 53,
                                  /* expectedCartesianIndex = */ 69, /* expectedCoords = */ {5, 2, 1}, isParallel, /*level =*/ 1);
                found_69 = true;
            }
        }
        else if (originId == 11) {
            // Check level Cartesian index is contained in {38,39,46,47, 86,87,94,95}
            const auto it = std::find(levelCartIndicesOriginId11.begin(), levelCartIndicesOriginId11.end(),
                                      levelCartMapp.cartesianIndex(element.index(), 1) );
            BOOST_CHECK( it !=  levelCartIndicesOriginId11.end() );

            // Check a few level Cartesian coordinates (e.g., one in layer k =0, and one in layer k = 1).
            // Note: take into account the "serial illustration".
            //
            // Serial: LGR1 refined cell with lgr1CompressedIndex = 89 has level Cartesian index 39 and ijk = {7, 4, 0}.
            if (*it == 39) {
                checkLevelElement(levelCartMapp, /* levelCompressedIndex = */ element.index(), /* serialLevelCompressedIndex = */ 89,
                                  /* expectedCartesianIndex = */ 39, /* expectedCoords = */ {7, 4, 0}, isParallel, /*level =*/ 1);
                found_39 = true;
            }
            // Serial: LGR1 refined cell with lgr1CompressedIndex = 92 has level Cartesian index 69 and ijk = {6, 4, 1}.
            if (*it == 86) {
                checkLevelElement(levelCartMapp, /* levelCompressedIndex = */ element.index(), /* serialLevelCompressedIndex = */ 92,
                                  /* expectedCartesianIndex = */ 86, /* expectedCoords = */ {6, 4, 1}, isParallel, /*level =*/ 1);
                found_86 = true;
            }
        }
    }
    BOOST_CHECK(grid.comm().max(found_9));
    BOOST_CHECK(grid.comm().max(found_57));
    BOOST_CHECK(grid.comm().max(found_28));
    BOOST_CHECK(grid.comm().max(found_69));
    BOOST_CHECK(grid.comm().max(found_39));
    BOOST_CHECK(grid.comm().max(found_86));
}

void checkGloballyRefinedTestGrids(const Dune::CpGrid& grid,
                                   const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapp,
                                   const Opm::LevelCartesianIndexMapper<Dune::CpGrid>& levelCartMapp,
                                   bool isParallel)
{
    int serialLeafGridCellCount = 288; // 4*3*3 levelZeroCells - [4*3*3 parentCells] + (8*6*6 LGR1Cells)
    checkDimsAndSizes(grid, cartMapp,  levelCartMapp,
                      {{4,3,3}, {8,6,6}},    // expected Cartesian dimensions per level grid
                      {4*3*3, 8*6*6},        // expected Cartesian sizes per level grid
                      {4*3*3, 8*6*6},        // expected compressed sizes per level grid
                      serialLeafGridCellCount);

    levelCartSizeEqualsCompressedSizeIfAllActiveAndSerial(grid, levelCartMapp, isParallel);

    // For hidden-globally-refined grid via addLgrsUpdateLeafView(...), or adapt(), or explicit
    // globalRefine(...),
    // CartesianIndexMapper::cartesianDimensions and grid.logicalCartesianSize() DO NOT coincide
    Opm::areEqual(cartMapp.cartesianDimensions(), {4,3,3} /* level zero grid Cartesian dimensions */);
    Opm::areEqual(grid.logicalCartesianSize(), {8,6,6} /* LEAF grid view Cartesian dimensions */);

    checkAFewCartIdxAndCoordsForGloballyRefinedTestGrids(grid, levelCartMapp, cartMapp, isParallel);
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

    // Block shaped parent cells of LGR1 dimensions (3-0)x(2-0)x(2-1) = 3x2x1.
    // Number of subdivisions per cell, per direction {3,3,3}.       -> LGR1 dims = {9,6,3}
    //
    // Block shaped parent cells of LGR2 dimensions (4-2)x(3-2)x(3-2) = 2x1x1.
    // Number of subdivisions per cell, per direction {3,3,3}.       -> LGR2 dims = {6,3,3}

    const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapp(grid);
    const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
    // I would like to create Opm::LeafCartesianIndexMapper<Dune::CpGrid> leafCartMapp(grid);

    int serialLeafGridCellCount = 244; // 4*3*3 levelZeroCells - [(3*2*1) + (2*1*1) parentCells] + (9*6*3 LGR1Cells) + (6*3*3 LGR2Cells)
    checkDimsAndSizes(grid, cartMapp,  levelCartMapp,
                      {{4,3,3}, {9,6,3}, {6,3,3}},  // expected Cartesian dimensions per level grid
                      {4*3*3, 9*6*3, 6*3*3},        // expected Cartesian sizes per level grid
                      {4*3*3, 9*6*3, 6*3*3},        // expected compressed sizes per level grid
                      serialLeafGridCellCount);

    // For strict-locally-refined grid, CartesianIndexMapper::cartesianDimensions and
    // grid.logicalCartesianSize() coincide (with level zero Cartesian dimensions/logical Cartesian size).
    Opm::areEqual(cartMapp.cartesianDimensions(), grid.logicalCartesianSize());

    levelCartSizeEqualsCompressedSizeIfAllActiveAndSerial(grid, levelCartMapp, isParallel);

    // Level zero grid
    // --------------------
    // k = 0  8   9  10  11 |  k = 1 20  21  22  23 | k = 2 32  33  34  35 |
    //        4   5   6   7 |        16  17  18  19 |       28  29  30  31 |
    //        0   1   2   3 |        12  13  14  15 |       24  25  26  27 |
    // Notice that all cells are active and the grid is not distributed.

    // In order to test serial and parallel reusing the same grid.
    bool l0_foundId0 = false;
    bool l0_foundId17 = false;
    bool l0_foundId35 = false;
    for (const auto& element : Dune::elements(grid.levelGridView(0))) {

        const auto& elemId = grid.globalIdSet().id(element);

        // Serial:   level zero element compressedIdx = 0 -> level zero Cartesian index = 0, level ijk = {0, 0, 0}
        // Parallel: seach for the element whose global id equal to 0.
        bool hasId0 = (elemId == 0);

        // Serial:   level zero element compressedIdx = 17 -> level zero Cartesian index = 17, level ijk = {1, 1, 1}
        // Parallel: seach for the element whose global id equal to 17.
        bool hasId17 = (elemId == 17);

        // Serial: level zero element compressedIdx = 35 -> level zero Cartesian index = 35, level ijk = {3, 2, 2}
        // Parallel: seach for the element whose global id equal to 35.
        bool hasId35 = (elemId == 35);

        if (hasId0) {
            checkLevelElement(levelCartMapp, element.index(), /* serialLevelCompressedIndex = */ 0,
                              /* expectedCartesianIndex = */ 0, /* expectedCoords = */ {0, 0, 0}, isParallel, /* level = */ 0);
            l0_foundId0 = true;
        }
        else if (hasId17) {
            checkLevelElement(levelCartMapp, element.index(), /* serialLevelCompressedIndex = */ 17,
                              /* expectedCartesianIndex = */ 17, /* expectedCoords = */ {1, 1, 1}, isParallel, /* level = */ 0);
            l0_foundId17 = true;
        }
        else if (hasId35) {
            checkLevelElement(levelCartMapp, element.index(), /* serialLevelCompressedIndex = */ 35,
                              /* expectedCartesianIndex = */ 35, /* expectedCoords = */ {3, 2, 2}, isParallel, /* level = */ 0);
            l0_foundId35 = true;
        }
    }
    BOOST_CHECK(grid.comm().max(l0_foundId0));
    BOOST_CHECK(grid.comm().max(l0_foundId17));
    BOOST_CHECK(grid.comm().max(l0_foundId35));

    // LGR1 parent cell global ids   | 16 17 18 |
    //                               | 12 13 14 |
    // LGR1 local/compressed indices                |     LGR1 level Cartesian indices
    // k = 2  105 106 107 |132 133 134 |159 160 161 |     153 154 155 |156 157 158 |159 160 161 |
    //        102 103 104 |129 130 131 |156 157 158 |     144 145 146 |147 148 149 |150 151 152 |
    //         99 100 101 |126 127 128 |153 154 155 |     135 136 137 |138 139 140 |141 142 143 |
    //         --------------------------------------     ---------------------------------------
    //         24  25  26 | 51  52  53 | 78  79  80 |     126 127 128 |129 130 131 |132 133 134 |
    //         21  22  23 | 48  49  50 | 75  76  77 |     117 118 119 |120 121 122 |123 124 125 |
    //         18  19  20 | 45  46  47 | 72  73  74 |     108 109 110 |111 112 113 |114 115 116 |
    // ---------------------------------------------      ---------------------------------------
    // k = 1   96  97  98 |123 124 125 |150 151 152 |      99 100 101 |102 103 104 |105 106 107 |
    //         93  94  95 |120 121 122 |147 148 149 |      90  91  92 | 93  94  95 | 96  97  98 |
    //         90  91  92 |117 118 119 |144 145 146 |      81  82  83 | 84  85  86 | 87  88  89 |
    //         -------------------------------------      ---------------------------------------
    //         15  16  17 | 42  43  44 | 69  70  71 |      72  73  74 | 75  76  77 | 78  79  80 |
    //         12  13  14 | 39  40  41 | 66  67  68 |      63  64  65 | 66  67  68 | 69  70  71 |
    //          9  10  11 | 36  37  38 | 63  64  65 |      54  55  56 | 57  58  59 | 60  61  62 |
    // ----------------------------------------------     ---------------------------------------
    // k = 0   87  88  89 |114 115 116 |141 142 143 |      45  46  47 | 48  49  50 | 51  52  53 |
    //         84  85  86 |111 112 113 |138 139 140 |      36  37  38 | 39  40  41 | 42  43  44 |
    //         81  82  83 |108 109 110 |135 136 137 |      27  28  29 | 30  31  32 | 33  34  35 |
    //         --------------------------------------     ---------------------------------------
    //          6   7   8 | 33  34  35 | 60  61  62 |      18  19  20 | 21  22  23 | 24  25  26 |
    //          3   4   5 | 30  31  32 | 57  58  59 |       9  10  11 | 12  13  14 | 15  16  17 |
    //          0   1   2 | 27  28  29 | 54  55  56 |       0   1   2 |  3   4   5 |  6   7   8 |
    // ------------------------------------------------------------------------------------------

    // Serial Local indices of refined    |   Level Cartesian indices of refined
    // cells children from parent cell    |   cells children from parent cell 
    // with global id 17                  |   with global id 17
    // k = 2   ...|132 133 134 |...       |    ...|156 157 158 |...
    //         ...|129 130 131 |...       |    ...|147 148 149 |...
    //         ...|126 127 128 |...       |    ...|138 139 140 |...
    // ------------------------------------------------------------
    // k = 1   ...|123 124 125 |...       |    ...|102 103 104 |...
    //         ...|120 121 122 |...       |    ...| 93  94  95 |...
    //         ...|117 118 119 |...       |    ...| 84  85  86 |...
    // ------------------------------------------------------------
    // k = 0   ...|114 115 116 |...       |    ...| 48  49  50 |... 
    //         ...|111 112 113 |...       |    ...| 39  40  41 |... 
    //         ...|108 109 110 |...       |    ...| 30  31  32 |...
    std::vector<int> levelCartIndicesOriginId17 = {30,31,32,39,40,41,48,49,50,           // k = 0
                                                   84,85,86,93,94,95,102,103,104,        // k = 1
                                                   138,139,140,147,148,149,156,157,158}; // k = 2
    bool found_39 = false; // serial local index 111 -> level Cartesian index 39, level ijk = {3,4,0}
    bool found_94 = false; // serial local index 121 -> level Cartesian index 94, level ijk = {4,4,1}
    bool found_158 = false; // serial local index 134 -> level Cartesian index 158, level ijk = {5,5,2}
    
    for (const auto& element : Dune::elements(grid.levelGridView(1))) {

        const auto& originId = grid.globalIdSet().id(element.getOrigin());

        if (originId == 17){ 
            // Check level Cartesian index is contained in levelCartIndicesOriginId17
            const auto it = std::find(levelCartIndicesOriginId17.begin(), levelCartIndicesOriginId17.end(),
                                      levelCartMapp.cartesianIndex(element.index(), 1) );
            BOOST_CHECK( it !=  levelCartIndicesOriginId17.end() );

            // Check a few level Cartesian coordinates (e.g., one per k-layer).
            // Note: take into account the "serial illustration".
            if (*it == 39) {
                // Serial: LGR1 compressedIdx = 111 -> LGR1 Cartesian index = 39, level ijk = {3,4,0}
                checkLevelElement(levelCartMapp, /* levelCompressedIndex = */ element.index(), /* serialLevelCompressedIndex = */ 111,
                                  /* expectedCartesianIndex = */ 39, /* expectedCoords = */ {3,4,0}, isParallel, /*level =*/ 1);
                found_39 = true;
            }
             if (*it == 94) {
                // Serial: LGR1 compressedIdx = 121 -> LGR1 Cartesian index = 94, level ijk = {4,4,1}
                checkLevelElement(levelCartMapp, /* levelCompressedIndex = */ element.index(), /* serialLevelCompressedIndex = */ 121,
                                  /* expectedCartesianIndex = */ 94, /* expectedCoords = */ {4,4,1}, isParallel, /*level =*/ 1);
                found_94 = true;
            }
              if (*it == 158) {
                // Serial: LGR1 compressedIdx = 134 -> LGR1 Cartesian index = 158, level ijk = {5,5,2}
                checkLevelElement(levelCartMapp, /* levelCompressedIndex = */ element.index(), /* serialLevelCompressedIndex = */ 134,
                                  /* expectedCartesianIndex = */ 158, /* expectedCoords = */ {5,5,2}, isParallel, /*level =*/ 1);
                found_158 = true;
              }
        }
    }
    BOOST_CHECK(grid.comm().max(found_39));
    BOOST_CHECK(grid.comm().max(found_94));
    BOOST_CHECK(grid.comm().max(found_158));

    // LGR2 parent cell global ids   | 34 35 |
    
    // LGR2 local/compressed indices   |     LGR2 level Cartesian indices
    // k = 2   24  25  26 | 51  52  53 |     48  49  50 | 51  52  53 |
    //         21  22  23 | 48  49  50 |     42  43  44 | 45  46  47 |
    //         18  19  20 | 45  46  47 |     36  37  38 | 39  40  41 |
    // ---------------------------------     -------------------------
    // k = 1   15  16  17 | 42  43  44 |     30  31  32 | 33  34  35 |
    //         12  13  14 | 39  40  41 |     24  25  26 | 27  28  29 |
    //          9  10  11 | 36  37  38 |     18  19  20 | 21  22  23 |
    // ---------------------------------     -------------------------
    // k = 0    6   7   8 | 33  34  35 |     12  13  14 | 15  16  17 |
    //          3   4   5 | 30  31  32 |      6   7   8 |  9  10  11 |
    //          0   1   2 | 27  28  29 |      0   1   2 |  3   4   5 |
    // ---------------------------------     -------------------------
    std::vector<int> levelCartIndicesOriginId35 = {3,4,5,9,10,11,15,16,17,      // k = 0
                                                   21,22,23,27,28,29,33,34,35,  // k = 1
                                                   39,40,41,45,46,47,51,52,53}; // k = 2
    bool found_3 = false; // serial local index 27 -> level Cartesian index 3, level ijk = {3,0,0}
    bool found_28 = false; // serial local index 40 -> level Cartesian index 28, level ijk = {4,1,1}
    bool found_45 = false; // serial local index 48 -> level Cartesian index 45, level ijk = {3,1,2}

    for (const auto& element : Dune::elements(grid.levelGridView(2))) {

        const auto& originId = grid.globalIdSet().id(element.getOrigin());

        if (originId == 35){
            // Check level Cartesian index is contained in levelCartIndicesOriginId17
            const auto it = std::find(levelCartIndicesOriginId35.begin(), levelCartIndicesOriginId35.end(),
                                      levelCartMapp.cartesianIndex(element.index(), 2) );
            BOOST_CHECK( it !=  levelCartIndicesOriginId35.end() );

            // Check a few level Cartesian coordinates (e.g., one per k-layer).
            // Note: take into account the "serial illustration".
            if (*it == 3) {
                // LGR2 compressedIdx = 27 -> LGR2 Cartesian index = 3,  level ijk = {3, 0, 0}
                checkLevelElement(levelCartMapp, /* levelCompressedIndex = */ element.index(), /* serialLevelCompressedIndex = */ 27,
                                  /* expectedCartesianIndex = */ 3, /* expectedCoords = */ {3, 0, 0}, isParallel, /*level =*/ 2);
                found_3 = true;
            }
            else if (*it == 28) {
                // LGR2 compressedIdx = 40 -> LGR2 Cartesian index = 28,  level ijk = {4,1,1}
                checkLevelElement(levelCartMapp, /* levelCompressedIndex = */ element.index(), /* serialLevelCompressedIndex = */ 40,
                                  /* expectedCartesianIndex = */ 28, /* expectedCoords = */ {4, 1, 1}, isParallel, /*level =*/ 2);
                found_28 = true;
            }
            else if (*it == 45) {
                // LGR2 compressedIdx = 48 -> LGR2 Cartesian index = 45,  level ijk = {3, 1, 2}
                checkLevelElement(levelCartMapp, /* levelCompressedIndex = */ element.index(), /* serialLevelCompressedIndex = */ 48,
                                  /* expectedCartesianIndex = */ 45, /* expectedCoords = */ {3, 1, 2}, isParallel, /*level =*/ 2);
                found_45 = true;
            }
        }
    }
    BOOST_CHECK(grid.comm().max(found_3));
    BOOST_CHECK(grid.comm().max(found_28));
    BOOST_CHECK(grid.comm().max(found_45));

    // Level zero grid
    // --------------------
    // k = 0  8   9  10  11 |  k = 1 20  21  22  23 | k = 2 32  33  34  35 |
    //        4   5   6   7 |        16  17  18  19 |       28  29  30  31 |
    //        0   1   2   3 |        12  13  14  15 |       24  25  26  27 |


    // Leaf grid local indices ( ** and ++ represents refined cells from LGR1 and LGR2 respectively)
    // -----------------------------------------------------------------------
    // k = 0  8   9  10  11 |  k = 1 176 177 178 179 | k = 2 188 189  ++  ++ |
    //        4   5   6   7 |         **  **  ** 175 |       184 185 186 187 |
    //        0   1   2   3 |         **  **  **  93 |       180 181 182 183 |
    // -----------------------------------------------------------------------
    // [] and () represent leaf local indices of refined cells born in LGR1 and LGR2 respectively.
    // -----------------------------------------------------------------------
    // k = 0  8   9  10  11 |  k = 1    176       177       178     179 | k = 2 188 189 (190-216) (217-243) |
    //        4   5   6   7 |        [ 94-120] [121-147] [148-174]  175 |       184 185    186      187     |
    //        0   1   2   3 |        [ 12-38 ] [ 39-65 ] [ 66-92 ]   93 |       180 181    182      183     |
    // ------------------------------------------------------------------------------------------------------
    // In order to test serial and parallel reusing the same grid.
    bool foundId5 = false;
    bool foundId15 = false;
    bool foundId33 = false;
    bool foundId13= false;
    bool foundId16 = false;
    bool foundId34 = false;
    bool foundId35 = false;
    for (const auto& element : Dune::elements(grid.leafGridView())) {

        const auto& originId = grid.globalIdSet().id(element.getOrigin());
        // Serial:   leaf coarse cells in the layer k = 0 have same compressedIndex and CartesianIndex
        //           Leaf coarse cell with compressedIndex = 5 is equivalent to its origin cell from
        //           level zero, with level zero Cartesian Index = 5 and level zero ijk = {1, 1, 0}.
        // Parallel: seach for the element whose origin cell (equivalent cell in level zero)
        //           has global id equal to 5.
        bool originHasId5 = (originId == 5);

        // Serial:   Leaf coarse cells on the boundary of LGR1
        //           Leaf coarse cell with compressedIndex = 93 is equivalent to its origin cell from
        //           level zero, with level zero Cartesian Index = 15 and level zero ijk = {3, 0, 1}.
        // Parallel: seach for the element whose origin cell (equivalent cell in level zero)
        //           has global id equal to 15.
        bool originHasId15 = (originId == 15);

        // Serial:  Leaf coarse cells on the boundary of LGR2
        //          Leaf coarse cell with compressedIndex = 189 is equivalent to its origin cell from
        //          level zero, with level zero Cartesian Index = 33 and level zero ijk = {1, 2, 2}.
        // Parallel: seach for the element whose origin cell (equivalent cell in level zero)
        //           has global id equal to 33.
        bool originHasId33 = (originId == 33);

        // Serial:  Leaf refined cells inherit Cartesian index and coordinates of their parent cell.
        //          Leaf refined cell born in LGR1 with compressedIndex = 39,40,..., 65 have
        //          parent cell in level zero with Cartesian index 13 and ijk = {1, 0, 1}.
        // Parallel: seach for the element whose origin cell (father cell in level zero)
        //           has global id equal to 13.
        bool originHasId13 = (originId == 13);

        // Serial:   Leaf refined cells born in LGR1 with compressedIndex = 94, ..., 120 have
        //           parent cell in level zero with Cartesian index 16 and ijk = {0, 1, 1}.
        // Parallel: seach for the element whose origin cell (father cell in level zero)
        //           has global id equal to 16.
        bool originHasId16 = (originId == 16);

        // Serial:  Leaf refined cell born in LGR2 with compressedIndex = 190,191,..., 216 have
        //          parent cell in level zero with Cartesian index 34 and ijk = {2, 2, 2}.
        // Parallel: seach for the element whose origin cell (father cell in level zero)
        //           has global id equal to 34.
        bool originHasId34 = (originId == 34);

        // Serial: Leaf refined cells born in LGR2 with compressedIndex = 217,218,..., 243 have
        //         parent cell in level zero with Cartesian index 35 and ijk = {3, 2, 2}.
        // Parallel: seach for the element whose origin cell (father cell in level zero)
        //           has global id equal to 35.
        bool originHasId35 = (originId == 35);

        if (originHasId5) {
            checkLeafElement(cartMapp, element, /* serialMinCompressedIndex = */ 5, /* serialMaxCompressedIndex = */ 5,
                             /* expectedCartesianIndex = */ 5, /* expectedCoords = */ {1, 1, 0}, isParallel);
            foundId5 = true;
        }
        else if (originHasId15) {
            checkLeafElement(cartMapp, element, /* serialMinCompressedIndex = */ 93, /* serialMaxCompressedIndex = */ 93,
                             /* expectedCartesianIndex = */ 15, /* expectedCoords = */ {3, 0, 1}, isParallel);
            foundId15 = true;
        }
        else if (originHasId33) {
            checkLeafElement(cartMapp, element, /* serialMinCompressedIndex = */ 189, /* serialMaxCompressedIndex = */ 189,
                             /* expectedCartesianIndex = */ 33, /* expectedCoords = */ {1, 2, 2}, isParallel);
            foundId33 = true;
        }
        else if (originHasId13) {
            checkLeafElement(cartMapp, element, /* serialMinCompressedIndex = */ 39, /* serialMaxCompressedIndex = */ 65,
                             /* expectedCartesianIndex = */ 13, /* expectedCoords = */ {1, 0, 1}, isParallel);
            foundId13 = true;
        }
        else if (originHasId16) {
            checkLeafElement(cartMapp, element, /* serialMinCompressedIndex = */ 94, /* serialMaxCompressedIndex = */ 120,
                             /* expectedCartesianIndex = */ 16, /* expectedCoords = */ {0, 1, 1}, isParallel);
            foundId16 = true;
        }
        else if (originHasId34) {
            checkLeafElement(cartMapp, element, /* serialMinCompressedIndex = */ 190, /* serialMaxCompressedIndex = */ 216,
                             /* expectedCartesianIndex = */ 34, /* expectedCoords = */ {2, 2, 2}, isParallel);
            foundId34 = true;
        }
        else if (originHasId35) {
            checkLeafElement(cartMapp, element, /* serialMinCompressedIndex = */ 217, /* serialMaxCompressedIndex = */ 243,
                             /* expectedCartesianIndex = */ 35, /* expectedCoords = */ {3, 2, 2}, isParallel);
            foundId35 = true;
        }
    }
    BOOST_CHECK(grid.comm().max(foundId5));
    BOOST_CHECK(grid.comm().max(foundId15));
    BOOST_CHECK(grid.comm().max(foundId33));
    BOOST_CHECK(grid.comm().max(foundId13));
    BOOST_CHECK(grid.comm().max(foundId16));
    BOOST_CHECK(grid.comm().max(foundId34));
    BOOST_CHECK(grid.comm().max(foundId35));
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

    int serialLeafGridCellCount = 64; // 4*3*3 levelZeroCells - [2*2*1 parentCells] + (4*4*2 LGR1Cells)
    // Level Cartesian dimensions and size, and compressed size, for refined level grids
    // take the values from level zero grid.
    // Even though the marked cells form a 2x2x1 block, the Cartesian dimensions and
    // Cartesian size of LGR1 is NOT {4,4,2} and 4*4*2 = 32, but {4,3,3} and 4*3*3 respectively.
    checkDimsAndSizes(grid, cartMapp,  levelCartMapp,
                      {{4,3,3}, {4,3,3}},    // expected Cartesian dimensions per level grid
                      {4*3*3, 4*3*3},        // expected Cartesian sizes per level grid
                      {4*3*3, 4*4*2},        // expected compressed sizes per level grid
                      serialLeafGridCellCount);

    // In this case, all parent cells are active, however, level Cartesian size
    // and compressed size do not coincide for the refined level grid.
    BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(1), /* levelZeroGridCellCount = */ 4*3*3);
    BOOST_CHECK_EQUAL(levelCartMapp.compressedSize(1), /* LGR1CellCount = */ grid.currentData()[1]->size(0));

    // For strict-local-refined grid via adapt(), CartesianIndexMapper::cartesianDimensions
    // and grid.logicalCartesianSize() coincide.
    Opm::areEqual(cartMapp.cartesianDimensions(), grid.logicalCartesianSize());

    // Level zero grid
    // --------------------
    // k = 0  8   9  10  11 |  k = 1 20  21  22  23 | k = 2 32  33  34  35 |
    //        4   5   6   7 |        16  17  18  19 |       28  29  30  31 |
    //        0   1   2   3 |        12  13  14  15 |       24  25  26  27 |

    // Leaf grid local indices ([] represent leaf local indices of refined cells born in LGR1)
    // -----------------------------------------------------------------------
    // k = 0  8   9  10  11 |  k = 1 34  [35-42] [43-50]  51 | k = 2 60  61  62  63 |
    //        4   5   6   7 |        16  [17-24] [25-32]  33 |       56  57  58  59 |
    //        0   1   2   3 |        12     13      14    15 |       52  53  54  55 |
    // ------------------------------------------------------------------------------
    // In order to test serial and parallel reusing the same grid.
    bool foundId5 = false;
    bool foundId19 = false;
    bool foundId18 = false;
    bool foundId21 = false;
    for (const auto& element : Dune::elements(grid.leafGridView())) {

        const auto& originId = grid.globalIdSet().id(element.getOrigin());

        // Serial:   Leaf coarse cells in the layer k = 0 have same compressedIndex and CartesianIndex
        //           Leaf coarse cell with compressedIndex = 5 is equivalent to its origin cell from
        //           level zero, with level zero Cartesian Index = 5 and level zero ijk = {1, 1, 0}.
        // Parallel: seach for the element whose origin cell (equivalent cell in level zero)
        //           has global id equal to 5.
        bool originHasId5 = (originId == 5);

        // Serial:   Leaf coarse cells on the boundary of LGR1
        //           Leaf coarse cell with compressedIndex = 33 is equivalent to its origin cell from
        //           level zero, with level zero Cartesian Index = 19 and level zero ijk = {3, 1, 1}.
        // Parallel: seach for the element whose origin cell (equivalent cell in level zero)
        //           has global id equal to 19.
        bool originHasId19 = (originId == 19);

        // Serial:   Leaf refined cells inherit Cartesian index and coordinates of their parent cell.
        //           Leaf refined cells born in LGR1 with compressedIndex = 25, 26, ..., 32 have
        //           parent cell in level zero with Cartesian index 18 and ijk = {2, 1, 1}.
        // Parallel: seach for the element whose origin cell (father cell in level zero)
        //           has global id equal to 18.
        bool originHasId18 = (originId == 18);

        // Serial:   Leaf refined cell born in LGR1 with compressedIndex = 35, 36, ..., 42 have
        //           parent cell in level zero with Cartesian index 21 and ijk = {1, 2, 1}.
        // Parallel: seach for the element whose origin cell (father cell in level zero)
        //           has global id equal to 21.
        bool originHasId21 = (originId == 21);

        if (originHasId5) {
            checkLeafElement(cartMapp, element, /* serialMinCompressedIndex = */ 5, /* serialMaxCompressedIndex = */ 5,
                             /* expectedCartesianIndex = */ 5, /* expectedCoords = */ {1, 1, 0}, isParallel);
            foundId5 = true;
        }
        else if (originHasId19) {
            checkLeafElement(cartMapp, element, /* serialMinCompressedIndex = */ 33, /* serialMaxCompressedIndex = */ 33,
                             /* expectedCartesianIndex = */ 19, /* expectedCoords = */ {3, 1, 1}, isParallel);
            foundId19 = true;
        }
        else if (originHasId18){
            checkLeafElement(cartMapp, element, /* serialMinCompressedIndex = */ 25, /* serialMaxCompressedIndex = */ 32,
                             /* expectedCartesianIndex = */ 18, /* expectedCoords = */ {2, 1, 1}, isParallel);
            foundId18 = true;
        }
        else if (originHasId21) {
            checkLeafElement(cartMapp, element, /* serialMinCompressedIndex = */ 35, /* serialMaxCompressedIndex = */ 42,
                             /* expectedCartesianIndex = */ 21, /* expectedCoords = */ {1, 2, 1}, isParallel);
            foundId21 = true;
        }
    }
    BOOST_CHECK(grid.comm().max(foundId5));
    BOOST_CHECK(grid.comm().max(foundId19));
    BOOST_CHECK(grid.comm().max(foundId18));
    BOOST_CHECK(grid.comm().max(foundId21));

    if (!isParallel) {
        // For simplicity, illustration only of layer k = 0 on the LGR1 grid
        // LGR1 local indices                LGR1 PARENT CELLS level zero Cartesian indices
        // k = 0   18  19 | 26  27 |           21  |  22  |
        //         16  17 | 24  25 |               |      |
        //         -----------------        ---------------
        //          2   3 | 10  11 |           17  |  18  |
        //          0   1 |  8   9 |               |      |
        // LGR1 refined cell with compressedIndex = 3 has parent cell in level zero
        // with Cartesian index 17 and ijk = {1, 1, 1}.
        checkLevelElement(levelCartMapp, /* levelCompressedIndex = */ 3, /* serialLevelCompressedIndex = */ 3,
                          /* expectedCartesianIndex = */ 17, /* expectedCoords = */ {1, 1, 1}, isParallel, /*level = */ 1);
        // LGR1 refined cell with compressedIndex = 26 has parent cell in level zero
        // with Cartesian index 22 and ijk = {2, 2, 1}.
        checkLevelElement(levelCartMapp, /* levelCompressedIndex = */ 26, /* serialLevelCompressedIndex = */ 26,
                          /* expectedCartesianIndex = */ 22, /* expectedCoords = */ {2, 2, 1}, isParallel, /*level =*/ 1);
    }
}

BOOST_AUTO_TEST_CASE(level_and_grid_cartesianIndexMapper_afterHiddenGlobalRefinementWith_addLgrsUpdateLeafView)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    bool isParallel = grid.comm().size() > 1;
    if (isParallel) {
        grid.loadBalance();
    }
    grid.addLgrsUpdateLeafView(/* cells_per_dim = */ {{2,2,2}},
                               /* startIJK_vec = */ {{0,0,0}},
                               /* endIJK_vec = */ {{4,3,3}},
                               /* lgr_name_vec = */ {"LGR1"});

    // Block shaped parent cells of LGR1 is the entire level zero grid, dimensions (4-0)x(3-0)x(3-1).
    // Number of subdivisions per cell, per direction {2,2,2}.  -> LGR1 dims = {8,6,6}

    const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapp(grid);
    const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
    // I would like to create Opm::LeafCartesianIndexMapper<Dune::CpGrid> leafCartMapp(grid);

    checkGloballyRefinedTestGrids(grid, cartMapp, levelCartMapp, isParallel);
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

    checkGloballyRefinedTestGrids(grid, cartMapp, levelCartMapp, isParallel);
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

    checkGloballyRefinedTestGrids(grid, cartMapp, levelCartMapp, isParallel);
}
/** TODO: Add case with inactive cells */
/** TODO: Define class LeafCartesianIndexMapper and include it in the test */
