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

void checkAFewCartIdxAndCoordsForGloballyRefinedGrid(const Dune::CpGrid& grid,
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

        // Serial:   Leaf refined cells (born in LGR1/GR) with compressedIndex = 0,1,...,7
        //           have parent cell in level zero with Cartesian index 0 and ijk = {0, 0, 0}.
        // Parallel: seach for the element whose origin cell (parent cell in level zero)
        //           has global id equal to 0.
        bool originHasId0 = (grid.globalIdSet().id(element.getOrigin()) == 0);

        // Serial:   Leaf refined cells (born in LGR1/GR) with compressedIndex = 144,145,....,151
        //           have parent cell in level zero with Cartesian index 18 and ijk = {2, 1, 1}.
        // Parallel: seach for the element whose origin cell (parent cell in level zero)
        //           has global id equal to 18.
        bool originHasId18 = (grid.globalIdSet().id(element.getOrigin()) == 18);

        // Serial:   Leaf refined cells (born in LGR1/GR) with compressedIndex = 264,265,...,271
        //           have parent cell with Cartesian index 33 and ijk = {1, 2, 2}.
        // Parallel: seach for the element whose origin cell (parent cell in level zero)
        //           has global id equal to 33.
        bool originHasId33 = (grid.globalIdSet().id(element.getOrigin()) == 33);

        if (originHasId0) {
            if (!isParallel) { // ALL CHILDREN of globalId cell 0 = {0,1, ..., 7}
                BOOST_CHECK( element.index() >= 0);
                BOOST_CHECK( element.index() <= 7);
            }
            BOOST_CHECK_EQUAL(cartMapp.cartesianIndex(/* compressedElementIndex = */ element.index()), 0);

            std::array<int,3> coords{};
            cartMapp.cartesianCoordinate(/* compressedElementIndex = */ element.index(), coords);
            Opm::areEqual( coords, {0, 0, 0});

            foundId0 = true;
        }
        else if (originHasId18) {
            if (!isParallel) { // ALL CHILDREN of globalId cell 18 = {144,145,...,151}
                BOOST_CHECK( element.index() >= 144);
                BOOST_CHECK( element.index() <= 151);
            }

            BOOST_CHECK_EQUAL(cartMapp.cartesianIndex(/* compressedElementIndex = */ element.index()), 18);

            std::array<int,3> coords{};
            cartMapp.cartesianCoordinate(/* compressedElementIndex = */ element.index(), coords);
            Opm::areEqual( coords, {2, 1, 1});

            foundId18 = true;
        }
        else if (originHasId33) {
            if (!isParallel) { // ALL CHILDREN of globalId cell 33 = {264,265,...,271}
                BOOST_CHECK( element.index() >= 264);
                BOOST_CHECK( element.index() <= 271);
            }

            BOOST_CHECK_EQUAL(cartMapp.cartesianIndex(/* compressedElementIndex = */ element.index()), 33);

            std::array<int,3> coords{};
            cartMapp.cartesianCoordinate(/* compressedElementIndex = */ element.index(), coords);
            Opm::areEqual( coords, {1, 2, 2});

            foundId33 = true;
        }
    }
    BOOST_CHECK(grid.comm().max(foundId0));
    BOOST_CHECK(grid.comm().max(foundId18));
    BOOST_CHECK(grid.comm().max(foundId33));

    // For simplicity, illustration only of layer k = 0 on the LGR1 grid
    // LGR1 local indices                                 LGR1 (level) Cartesian indices
    // k = 0   66  67 | 74  75 | 82  83 | 90  91 |        40  41 | 42  43 | 44  45 | 46  47 |
    //         64  65 | 72  73 | 80  81 | 88  89 |        32  33 | 34  35 | 36  37 | 38  39 |
    //         -----------------------------------        -----------------------------------
    //         34  35 | 42  43 | 50  51 | 58  59 |        24  25 | 26  27 | 28  29 | 30  31 |
    //         32  33 | 40  41 | 48  49 | 56  57 |        16  17 | 18  19 | 20  21 | 22  23 |
    //         -----------------------------------        -----------------------------------
    //          2   3 | 10  11 | 18  19 | 26  27 |         8   9 | 10  11 | 12  13 | 14  15 |
    //          0   1 |  8   9 | 16  17 | 24  25 |         0   1 |  2   3 |  4   5 |  6   7 |

    /** TODO: Test in parallel levelCartMapp cartesianIndex and cartesianCoordinate */
    if (!isParallel) {
        // LGR1 refined cell with lgr1CompressedIndex = 3 has level Cartesian index 9 and ijk = {1, 1, 0}.
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(/* lgr1CompressedElementIndex = */ 3,  /* level = */ 1), 9);
        std::array<int,3> lgr1Coords3{};
        levelCartMapp.cartesianCoordinate(/* lgr1CompressedElementIndex = */ 3, lgr1Coords3, /* level = */ 1);
        Opm::areEqual( lgr1Coords3, {1, 1, 0});

        // LGR1 refined cell with lgr1CompressedIndex = 50 has level Cartesian index 28 and ijk = {4, 3, 0}.
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(/* compressedElementIndex = */ 50,  /* level = */ 1), 28);
        std::array<int,3> lgr1Coords50{};
        levelCartMapp.cartesianCoordinate(/* compressedElementIndex = */ 50, lgr1Coords50, /* level = */ 1);
        Opm::areEqual( lgr1Coords50, {4, 3, 0});

        // LGR1 refined cell with lgr1CompressedIndex = 90 has level Cartesian index 46 and ijk = {6, 5, 0}.
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(/* compressedElementIndex = */ 90,  /* level = */ 1), 46);
        std::array<int,3> lgr1Coords90{};
        levelCartMapp.cartesianCoordinate(/* compressedElementIndex = */ 90, lgr1Coords90,  /* level = */ 1);
        Opm::areEqual( lgr1Coords90, {6, 5, 0});
    }
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

    // Check level Cartesian dimensions and size, and compressed size, for level grids.
    checkLevels(grid,
                levelCartMapp,
                {{4,3,3}, {9,6,3}, {6,3,3}}, // expected Cartesian dimensions per level grid
                {4*3*3, 9*6*3, 6*3*3},       // expected Cartesian sizes per level grid
                {4*3*3, 9*6*3, 6*3*3});      // expected compressed sizes per level grid

    // In this case, all parent cells are active, therefore level Cartesian size
    // and compressed size coincide.
    if (!isParallel) {
        for (int level = 0; level <= grid.maxLevel(); ++level) {
            BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(level), levelCartMapp.compressedSize(level));
            BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(level), grid.currentData()[level]->size(0));
        }
    }

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

        // Serial:   level zero element compressedIdx = 0 -> level zero Cartesian index = 0, level ijk = {0, 0, 0}
        // Parallel: seach for the element whose global id equal to 0.
        bool hasId0 = (grid.globalIdSet().id(element) == 0);

        // Serial:   level zero element compressedIdx = 17 -> level zero Cartesian index = 17, level ijk = {1, 1, 1}
        // Parallel: seach for the element whose global id equal to 17.
        bool hasId17 = (grid.globalIdSet().id(element) == 17);

        // Serial: level zero element compressedIdx = 35 -> level zero Cartesian index = 35, level ijk = {3, 2, 2}
        // Parallel: seach for the element whose global id equal to 35.
        bool hasId35 = (grid.globalIdSet().id(element) == 35);

        if (hasId0) {
            if (!isParallel) {
                BOOST_CHECK_EQUAL( element.index(), 0);
            }

            BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(/* levelZeroCompressedElementIndex  = */ element.index(), /* level = */ 0), 0);

            std::array<int,3> level0Coords0{};
            levelCartMapp.cartesianCoordinate(/*levelZeroMinCompressedElementIndex = */ element.index(), level0Coords0, /* level = */ 0);
            Opm::areEqual(level0Coords0, {0, 0, 0});

            l0_foundId0 = true;
        }
        else if (hasId17) {
            if (!isParallel) {
                BOOST_CHECK_EQUAL( element.index(), 17);
            }

            BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(/* levelZeroCompressedElementIndex  = */ element.index(), /* level = */ 0), 17);

            std::array<int,3> level0Coords17{};
            levelCartMapp.cartesianCoordinate(/*levelZeroCompressedElementIndex = */ element.index(), level0Coords17, /* level = */ 0);
            Opm::areEqual(level0Coords17, {1, 1, 1});

            l0_foundId17 = true;
        }
        else if (hasId35) {
            if (!isParallel) {
                BOOST_CHECK_EQUAL( element.index(), 35);
            }

            BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(/* levelZeroCompressedElementIndex  = */ element.index(), /* level = */ 0), 35);

            std::array<int,3> level0Coords35{};
            levelCartMapp.cartesianCoordinate(/*levelZeroMaxCompressedElementIndex = */ element.index(), level0Coords35, /* level = */ 0);
            Opm::areEqual(level0Coords35, {3, 2, 2});

            l0_foundId35 = true;
        }
    }
    BOOST_CHECK(grid.comm().max(l0_foundId0));
    BOOST_CHECK(grid.comm().max(l0_foundId17));
    BOOST_CHECK(grid.comm().max(l0_foundId35));


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
    /** TODO: Test in parallel levelCartMapp cartesianIndex and cartesianCoordinate */
    if (!isParallel) {
        // Serial:   LGR1 compressedIdx = 0 -> LGR1 Cartesian index = 0, level ijk = {0, 0, 0}
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(/* lgr1MinCompressedElementIndex  = */ 0,  /* level = */ 1), 0);
        std::array<int,3> lgr1Coords0{};
        levelCartMapp.cartesianCoordinate(/*lgr1MinCompressedElementIndex = */ 0, lgr1Coords0, /* level = */ 1);
        Opm::areEqual(lgr1Coords0, {0, 0, 0});

        // LGR1 compressedIdx = 161 -> LGR1 Cartesian index = 161, level ijk = {8, 5, 2}
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(/* lgr1MaxCompressedElementIndex  = */ 161, /* level = */ 1), 161);
        std::array<int,3> lgr1Coords161{};
        levelCartMapp.cartesianCoordinate(/*lgr1MaxCompressedElementIndex = */ 161, lgr1Coords161, /* level = */ 1);
        Opm::areEqual(lgr1Coords161, {8, 5, 2});

        // LGR1 compressedIdx = 121 -> LGR1 Cartesian index = 94, level ijk = {4, 4, 1}
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(/* lgr1CompressedElementIndex  = */ 121, /* level = */ 1), 94);
        std::array<int,3> lgr1Coords121{};
        levelCartMapp.cartesianCoordinate(/*lgr1CompressedElementIndex = */ 121, lgr1Coords121, /* level = */ 1);
        Opm::areEqual(lgr1Coords121, {4, 4, 1});

        // LGR1 compressedIdx = 104 -> LGR1 Cartesian index = 146, level ijk = {2, 4, 2}
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(/* lgr1CompressedElementIndex  = */ 104, /* level = */ 1), 146);
        std::array<int,3> lgr1Coords104{};
        levelCartMapp.cartesianCoordinate(/*lgr1CompressedElementIndex = */ 104, lgr1Coords104, /* level = */ 1);
        Opm::areEqual(lgr1Coords104, {2, 4, 2});

        // LGR1 compressedIdx = 71 -> LGR1 Cartesian index = 80
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(/* lgr1CompressedElementIndex  = */ 71, /* level = */ 1), 80);
        std::array<int,3> lgr1Coords71{};
        levelCartMapp.cartesianCoordinate(/*lgr1CompressedElementIndex = */ 71, lgr1Coords71, /* level = */ 1);
        Opm::areEqual(lgr1Coords71, {8, 2, 1});


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
        // LGR2 compressedIdx = 0 -> LGR2 Cartesian index = 0, level ijk = {0, 0, 0}
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(/* lgr2MinCompressedElementIndex  = */ 0, /* level = */ 2), 0);
        std::array<int,3> lgr2Coords0{};
        levelCartMapp.cartesianCoordinate(/*lgr2MinCompressedElementIndex = */ 0, lgr2Coords0, /* level = */ 2);
        Opm::areEqual(lgr2Coords0, {0, 0, 0});

        // LGR2 compressedIdx = 53 -> LGR2 Cartesian index = 53, level ijk = {5, 2, 2}
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(/* lgr2MaxCompressedElementIndex  = */ 53, /* level = */ 2), 53);
        std::array<int,3> lgr2Coords53{};
        levelCartMapp.cartesianCoordinate(/*lgr2CompressedElementIndex = */ 53, lgr2Coords53, /* level = */ 2);
        Opm::areEqual(lgr2Coords53, {5, 2, 2});

        // LGR2 compressedIdx = 27 -> LGR2 Cartesian index = 3,  level ijk = {3, 0, 0}
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(/* lgr2CompressedElementIndex  = */ 27, /* level = */ 2), 3);
        std::array<int,3> lgr2Coords27{};
        levelCartMapp.cartesianCoordinate(/*lgr2CompressedElementIndex = */ 27, lgr2Coords27, /* level = */ 2);
        Opm::areEqual(lgr2Coords27, {3, 0, 0});

        // LGR2 compressedIdx = 13 -> LGR2 Cartesian index = 25, level ijk = {1, 1, 1}
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(/* lgr2CompressedElementIndex  = */ 13, /* level = */ 2), 25);
        std::array<int,3> lgr2Coords13{};
        levelCartMapp.cartesianCoordinate(/*lgr2CompressedElementIndex = */ 13, lgr2Coords13, /* level = */ 2);
        Opm::areEqual(lgr2Coords13, {1, 1, 1});

        // LGR2 compressedIdx = 48 -> LGR2 Cartesian index = 45,  level ijk = {3, 1, 2}
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(/* lgr2CompressedElementIndex  = */ 48, /* level = */ 2), 45);
        std::array<int,3> lgr2Coords48{};
        levelCartMapp.cartesianCoordinate(/*lgr2CompressedElementIndex = */ 48, lgr2Coords48, /* level = */ 2);
        Opm::areEqual(lgr2Coords48, {3, 1, 2});
    }


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

    int serialLeafGridCellCount = 244; // 4*3*3 levelZeroCells - [(3*2*1) + (2*1*1) parentCells] + (9*6*3 LGR1Cells) + (6*3*3 LGR2Cells)
    compressedSizeCoincidesWithLeafGridCellCount(grid, cartMapp, serialLeafGridCellCount);

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

        // Serial:   leaf coarse cells in the layer k = 0 have same compressedIndex and CartesianIndex
        //           Leaf coarse cell with compressedIndex = 5 is equivalent to its origin cell from
        //           level zero, with level zero Cartesian Index = 5 and level zero ijk = {1, 1, 0}.
        // Parallel: seach for the element whose origin cell (equivalent cell in level zero)
        //           has global id equal to 5.
        bool originHasId5 = (grid.globalIdSet().id(element.getOrigin()) == 5);

        // Serial:   Leaf coarse cells on the boundary of LGR1
        //           Leaf coarse cell with compressedIndex = 93 is equivalent to its origin cell from
        //           level zero, with level zero Cartesian Index = 15 and level zero ijk = {3, 0, 1}.
        // Parallel: seach for the element whose origin cell (equivalent cell in level zero)
        //           has global id equal to 15.
        bool originHasId15 = (grid.globalIdSet().id(element.getOrigin()) == 15);

        // Serial:  Leaf coarse cells on the boundary of LGR2
        //          Leaf coarse cell with compressedIndex = 189 is equivalent to its origin cell from
        //          level zero, with level zero Cartesian Index = 33 and level zero ijk = {1, 2, 2}.
        // Parallel: seach for the element whose origin cell (equivalent cell in level zero)
        //           has global id equal to 33.
        bool originHasId33 = (grid.globalIdSet().id(element.getOrigin()) == 33);

        // Serial:  Leaf refined cells inherit Cartesian index and coordinates of their parent cell.
        //          Leaf refined cell born in LGR1 with compressedIndex = 39,40,..., 65 have
        //          parent cell in level zero with Cartesian index 13 and ijk = {1, 0, 1}.
        // Parallel: seach for the element whose origin cell (father cell in level zero)
        //           has global id equal to 13.
        bool originHasId13 = (grid.globalIdSet().id(element.getOrigin()) == 13);

        // Serial:   Leaf refined cells born in LGR1 with compressedIndex = 94, ..., 120 have
        //           parent cell in level zero with Cartesian index 16 and ijk = {0, 1, 1}.
        // Parallel: seach for the element whose origin cell (father cell in level zero)
        //           has global id equal to 16.
        bool originHasId16 = (grid.globalIdSet().id(element.getOrigin()) == 16);

        // Serial:  Leaf refined cell born in LGR2 with compressedIndex = 190,191,..., 216 have
        //          parent cell in level zero with Cartesian index 34 and ijk = {2, 2, 2}.
        // Parallel: seach for the element whose origin cell (father cell in level zero)
        //           has global id equal to 34.
        bool originHasId34 = (grid.globalIdSet().id(element.getOrigin()) == 34);

        // Serial: Leaf refined cells born in LGR2 with compressedIndex = 217,218,..., 243 have
        //         parent cell in level zero with Cartesian index 35 and ijk = {3, 2, 2}.
        // Parallel: seach for the element whose origin cell (father cell in level zero)
        //           has global id equal to 35.
        bool originHasId35 = (grid.globalIdSet().id(element.getOrigin()) == 35);

        if (originHasId5) {
            if (!isParallel) {
                BOOST_CHECK_EQUAL( element.index(), 5);
            }
            BOOST_CHECK_EQUAL(cartMapp.cartesianIndex(/* compressedElementIndex = */ element.index()), 5);

            std::array<int,3> coords{};
            cartMapp.cartesianCoordinate(/* compressedElementIndex = */ element.index(), coords);
            Opm::areEqual( coords, {1, 1, 0});

            foundId5 = true;
        }
        else if (originHasId15) {
            if (!isParallel) {
                BOOST_CHECK_EQUAL( element.index(), 93);
            }
            BOOST_CHECK_EQUAL(cartMapp.cartesianIndex(/* compressedElementIndex = */ element.index()), 15);
            std::array<int,3> coords{};
            cartMapp.cartesianCoordinate(/* compressedElementIndex = */ element.index(), coords);
            Opm::areEqual( coords, {3, 0, 1});

            foundId15 = true;
        }
        else if (originHasId33) {
            if (!isParallel) {
                BOOST_CHECK_EQUAL( element.index(), 189);
            }
            BOOST_CHECK_EQUAL(cartMapp.cartesianIndex(/* compressedElementIndex = */ element.index()), 33);

            std::array<int,3> coords{};
            cartMapp.cartesianCoordinate(/* compressedElementIndex = */ element.index(), coords);
            Opm::areEqual( coords, {1, 2, 2});

            foundId33 = true;
        }
        else if (originHasId13) {
            if (!isParallel) { // ALL CHILDREN of globalId cell 13 = {39,40,..., 65}
                BOOST_CHECK( element.index() >= 39);
                BOOST_CHECK( element.index() <= 65);
            }
            BOOST_CHECK_EQUAL(cartMapp.cartesianIndex(/* compressedElementIndex = */ element.index()), 13);

            std::array<int,3> coords{};
            cartMapp.cartesianCoordinate(/* compressedElementIndex = */ element.index(), coords);
            Opm::areEqual( coords, {1, 0, 1});

            foundId13 = true;
        }
        else if (originHasId16) {
            if (!isParallel) { // ALL CHILDREN of globalId cell 16 = {94,95,..., 120}
                BOOST_CHECK( element.index() >= 94);
                BOOST_CHECK( element.index() <= 120);
            }
            BOOST_CHECK_EQUAL(cartMapp.cartesianIndex(/* compressedElementIndex = */ element.index()), 16);

            std::array<int,3> coords{};
            cartMapp.cartesianCoordinate(/* compressedElementIndex = */ element.index(), coords);
            Opm::areEqual( coords, {0, 1, 1});

            foundId16 = true;
        }
        else if (originHasId34) {
            if (!isParallel) { // ALL CHILDREN of globalId cell 34 = {190,191,..., 216}
                BOOST_CHECK( element.index() >= 190);
                BOOST_CHECK( element.index() <= 216);
            }
            BOOST_CHECK_EQUAL(cartMapp.cartesianIndex(/* compressedElementIndex = */ element.index()), 34);

            std::array<int,3> coords{};
            cartMapp.cartesianCoordinate(/* compressedElementIndex = */ element.index(), coords);
            Opm::areEqual( coords, {2, 2, 2});

            foundId34 = true;
        }
        else if (originHasId35) {
            if (!isParallel) { // ALL CHILDREN of globalId cell 35 = {217,218,..., 243}
                BOOST_CHECK( element.index() >= 217);
                BOOST_CHECK( element.index() <= 243);
            }
            BOOST_CHECK_EQUAL(cartMapp.cartesianIndex(/* compressedElementIndex = */ element.index()), 35);

            std::array<int,3> coords{};
            cartMapp.cartesianCoordinate(/* compressedElementIndex = */ element.index(), coords);
            Opm::areEqual( coords, {3, 2, 2});

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

    int serialLeafGridCellCount = 64; // 4*3*3 levelZeroCells - [2*2*1 parentCells] + (4*4*2 LGR1Cells)
    compressedSizeCoincidesWithLeafGridCellCount(grid, cartMapp, serialLeafGridCellCount);

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
        // Serial:   Leaf coarse cells in the layer k = 0 have same compressedIndex and CartesianIndex
        //           Leaf coarse cell with compressedIndex = 5 is equivalent to its origin cell from
        //           level zero, with level zero Cartesian Index = 5 and level zero ijk = {1, 1, 0}.
        // Parallel: seach for the element whose origin cell (equivalent cell in level zero)
        //           has global id equal to 5.
        bool originHasId5 = (grid.globalIdSet().id(element.getOrigin()) == 5);

        // Serial:   Leaf coarse cells on the boundary of LGR1
        //           Leaf coarse cell with compressedIndex = 33 is equivalent to its origin cell from
        //           level zero, with level zero Cartesian Index = 19 and level zero ijk = {3, 1, 1}.
        // Parallel: seach for the element whose origin cell (equivalent cell in level zero)
        //           has global id equal to 19.
        bool originHasId19 = (grid.globalIdSet().id(element.getOrigin()) == 19);

        // Serial:   Leaf refined cells inherit Cartesian index and coordinates of their parent cell.
        //           Leaf refined cells born in LGR1 with compressedIndex = 25, 26, ..., 32 have
        //           parent cell in level zero with Cartesian index 18 and ijk = {2, 1, 1}.
        // Parallel: seach for the element whose origin cell (father cell in level zero)
        //           has global id equal to 18.
        bool originHasId18 = (grid.globalIdSet().id(element.getOrigin()) == 18);

        // Serial:   Leaf refined cell born in LGR1 with compressedIndex = 35, 36, ..., 42 have
        //           parent cell in level zero with Cartesian index 21 and ijk = {1, 2, 1}.
        // Parallel: seach for the element whose origin cell (father cell in level zero)
        //           has global id equal to 21.
        bool originHasId21 = (grid.globalIdSet().id(element.getOrigin()) == 21);


        if (originHasId5) {
            if (!isParallel) {
                BOOST_CHECK_EQUAL( element.index(), 5);
            }
            BOOST_CHECK_EQUAL(cartMapp.cartesianIndex(/* compressedElementIndex = */ element.index()), 5);

            std::array<int,3> coords{};
            cartMapp.cartesianCoordinate(/* compressedElementIndex = */ element.index(), coords);
            Opm::areEqual( coords, {1, 1, 0});

            foundId5 = true;
        }
        else if (originHasId19) {
            if (!isParallel) {
                BOOST_CHECK_EQUAL( element.index(), 33);
            }

            BOOST_CHECK_EQUAL(cartMapp.cartesianIndex(/* compressedElementIndex = */ element.index()), 19);

            std::array<int,3> coords{};
            cartMapp.cartesianCoordinate(/* compressedElementIndex = */ element.index(), coords);
            Opm::areEqual( coords, {3, 1, 1});

            foundId19 = true;
        }
        else if (originHasId18){
            if (!isParallel) {  // ALL CHILDREN of globalId cell 18 = {25,26, ..., 32}
                BOOST_CHECK( element.index() >= 25);
                BOOST_CHECK( element.index() <= 32);
            }

            BOOST_CHECK_EQUAL(cartMapp.cartesianIndex(/* compressedElementIndex = */ element.index()), 18);

            std::array<int,3> coords{};
            cartMapp.cartesianCoordinate(/* compressedElementIndex = */ element.index(), coords);
            Opm::areEqual( coords, {2, 1, 1});

            foundId18 = true;
        }
        else if (originHasId21) {
            if (!isParallel) { // ALL CHILDREN of globalId cell 21 = {35,36, ..., 42}
                BOOST_CHECK( element.index() >= 35);
                BOOST_CHECK( element.index() <= 42);
            }

            BOOST_CHECK_EQUAL(cartMapp.cartesianIndex(/* compressedElementIndex = */ element.index()), 21);

            std::array<int,3> coords{};
            cartMapp.cartesianCoordinate(/* compressedElementIndex = */ element.index(), coords);
            Opm::areEqual( coords, {1, 2, 1});

            foundId21 = true;
        }
    }
    BOOST_CHECK(grid.comm().max(foundId5));
    BOOST_CHECK(grid.comm().max(foundId19));
    BOOST_CHECK(grid.comm().max(foundId18));
    BOOST_CHECK(grid.comm().max(foundId21));

    // For simplicity, illustration only of layer k = 0 on the LGR1 grid
    // LGR1 local indices                LGR1 PARENT CELLS level zero Cartesian indices
    // k = 0   18  19 | 26  27 |           21  |  22  |
    //         16  17 | 24  25 |               |      |
    //         -----------------        ---------------
    //          2   3 | 10  11 |           17  |  18  |
    //          0   1 |  8   9 |               |      |

    if (!isParallel) {

        // LGR1 refined cell with compressedIndex = 3 has parent cell in level zero
        // with Cartesian index 17 and ijk = {1, 1, 1}.
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(/* lgr1CompressedElementIndex = */ 3,  /* level = */ 1), 17);
        std::array<int,3> lgr1Coords3{};
        levelCartMapp.cartesianCoordinate(/* lgr1CompressedElementIndex = */ 3, lgr1Coords3, /* level = */ 1);
        Opm::areEqual( lgr1Coords3, {1, 1, 1});

        // LGR1 refined cell with compressedIndex = 26 has parent cell in level zero
        // with Cartesian index 22 and ijk = {2, 2, 1}.
        BOOST_CHECK_EQUAL(levelCartMapp.cartesianIndex(/* compressedElementIndex = */ 26,  /* level = */ 1), 22);
        std::array<int,3> lgr1Coords26{};
        levelCartMapp.cartesianCoordinate(/* compressedElementIndex = */ 26, lgr1Coords26, /* level = */ 1);
        Opm::areEqual( lgr1Coords26, {2, 2, 1});

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

    // Check level Cartesian dimensions and size, and compressed size, for level grids.
    checkLevels(grid,
                levelCartMapp,
                {{4,3,3}, {8,6,6}}, // expected Cartesian dimensions per level grid
                {4*3*3, 8*6*6},     // expected Cartesian sizes per level grid
                {4*3*3, 8*6*6});    // expected compressed sizes per level grid

    // In this case, all parent cells are active, therefore level Cartesian size
    // and compressed size coincide.
    if (!isParallel) {
        for (int level = 0; level <= grid.maxLevel(); ++level) {
            BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(level), levelCartMapp.compressedSize(level));
            BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(level), grid.currentData()[level]->size(0));
        }
    }

    allThrow(levelCartMapp, /* unexisting_level = */ grid.maxLevel()+1);

    cartesianDimensionsCoincideWithLevelZeroOnes(grid,
                                                 cartMapp,
                                                 levelCartMapp,
                                                 /* level zero Cartesian dims = */ {4,3,3});
    // For hidden-globally-refined grid via addLgrsUpdateLeafView, CartesianIndexMapper::cartesianDimensions
    // and grid.logicalCartesianSize() DO NOT coincide
    Opm::areEqual(cartMapp.cartesianDimensions(), {4,3,3} /* level zero grid Cartesian dimensions */);
    Opm::areEqual(grid.logicalCartesianSize(), {8,6,6} /* LEAF grid view Cartesian dimensions */);

    cartesianSizeCoincidesWithLevelZeroOne(grid,
                                           cartMapp,
                                           levelCartMapp,
                                           /* level zero Cartesian size = */ 4*3*3);

    int serialLeafGridCellCount = 288; // 4*3*3 levelZeroCells - [4*3*3 parentCells] + (8*6*6 LGR1Cells)
    compressedSizeCoincidesWithLeafGridCellCount(grid, cartMapp, serialLeafGridCellCount);

    checkAFewCartIdxAndCoordsForGloballyRefinedGrid(grid, levelCartMapp, cartMapp, isParallel);
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
    if (!isParallel) {
        for (int level = 0; level <= grid.maxLevel(); ++level) {
            BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(level), levelCartMapp.compressedSize(level));
            BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(level), grid.currentData()[level]->size(0));
        }
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

    int serialLeafGridCellCount = 288; // 4*3*3 levelZeroCells - [4*3*3 parentCells] + (8*6*6 LGR1Cells)
    compressedSizeCoincidesWithLeafGridCellCount(grid, cartMapp, serialLeafGridCellCount);

    checkAFewCartIdxAndCoordsForGloballyRefinedGrid(grid, levelCartMapp, cartMapp, isParallel);
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
    if (!isParallel) {
        for (int level = 0; level <= grid.maxLevel(); ++level) {
            BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(level), levelCartMapp.compressedSize(level));
            BOOST_CHECK_EQUAL(levelCartMapp.cartesianSize(level), grid.currentData()[level]->size(0));
        }
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

    int serialLeafGridCellCount = 288; // 4*3*3 levelZeroCells - [4*3*3 parentCells] + (8*6*6 LGR1Cells)
    compressedSizeCoincidesWithLeafGridCellCount(grid, cartMapp, serialLeafGridCellCount);

    checkAFewCartIdxAndCoordsForGloballyRefinedGrid(grid, levelCartMapp, cartMapp, isParallel);
}
