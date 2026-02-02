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

#define BOOST_TEST_MODULE LevelCartToLevelCompressedTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>
#include <opm/grid/cpgrid/LgrOutputHelpers.hpp>
#include <tests/cpgrid/lgr/LgrChecks.hpp>

#include <string>
#include <unordered_map>

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

BOOST_AUTO_TEST_CASE(simpleGridWith1InactiveParentCell_serial)
{
    const std::string deckString = R"( RUNSPEC
  DIMENS
 -- NX NY NZ cells per x-,y-, and z-direction
     3 1 1 /
  GRID
  COORD -- grid corner coordinates (bounding box), 6*(NX +1)*(NY +1) values.
  0 0 0   0 0 2 -- bottom/top of pillar 0
  1 0 0   1 0 2
  2 0 0   2 0 2
  3 0 0   3 0 2
  0 1 0   0 1 2
  1 1 0   1 1 2
  2 1 0   2 1 2
  3 1 0   3 1 2 -- bottom/top of pillar 7
  /
  ZCORN
-- flat layers, no pinch-outs; NZ=1 top and bottom per layer, 8*NX*NY*NZ values
  0 0 0 0  0 0 0 0  0 0 0 0  -- top layer k=0
  2 2 2 2  2 2 2 2  2 2 2 2  -- bottom layer k=0
  /
  ACTNUM
-- i=0,1,2, j=0
  1 0 1  -- layer k=0, j=0
  /
  PORO
  3*0.15
  /)";

    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckString,
                              /* cells_per_dim = */ {{2,2,1}},
                              /* startIJK_vec = */ {{0,0,0}},
                              /* endIJK_vec = */ {{3,1,1}},
                              /* lgr_name_vec = */ {"LGR1"});


    // LGR1 parent cell global ids = {0, 1}, inactive cell in the middle "is not seen"
    // Level zero grid 3x1x1: k = 0 |  0   -   1
    const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
    const auto levelCartToLevelCompressed = Opm::Lgr::levelCartesianToLevelCompressedMaps(grid, levelCartMapp);

    // Level zero grid           Level zero Cartesian indices
    // --------------------
    // k = 0  | 0   -   1 |      | 0   -   2 |
    // Check the size of the map (2 active cells in a level zero grid of dimensions 3x1x1)
    BOOST_CHECK_EQUAL(levelCartToLevelCompressed[/* level */ 0].size(), /* level active cells count */ 2);

    // Check the level Cartesian indices match the level active compressed indices
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 0].at(/* levelCartIdx */ 0), /* levelCompressedIdx */ 0);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 0].at(/* levelCartIdx */ 2), /* levelCompressedIdx */ 1);

    // Check non-existing key level Cartesian values due to inactive cell
    BOOST_CHECK( levelCartToLevelCompressed[/* level */ 0].find(1) ==  levelCartToLevelCompressed[/* level */ 0].end());


    //-----------------------------------------------------------------------
    // LGR1 serial compressed idxs  | LGR1 Cartesian indices, []: inactive  |
    // ----------------------------------------------------------------------
    // k = 0     2  3 | INAC | 6  7 |       k = 0     6  7 | [8] [9] |10 11 |
    //           0  1 |      | 4  5 |                 0  1 | [2] [3] | 4  5 |

    // Check the size of the map (8 active cells in a lgr with dimensions 6x2x2)
    BOOST_CHECK_EQUAL(levelCartToLevelCompressed[/* level */ 1].size(), /* level active cells count */ 8);

    // Check the level Cartesian indices match the level active compressed indices
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 1].at(/* levelCartIdx */ 0), /* levelCompressedIdx */ 0);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 1].at(/* levelCartIdx */ 1), /* levelCompressedIdx */ 1);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 1].at(/* levelCartIdx */ 4), /* levelCompressedIdx */ 4);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 1].at(/* levelCartIdx */ 5), /* levelCompressedIdx */ 5);

    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 1].at(/* levelCartIdx */ 6), /* levelCompressedIdx */ 2);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 1].at(/* levelCartIdx */ 7), /* levelCompressedIdx */ 3);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 1].at(/* levelCartIdx */ 10), /* levelCompressedIdx */ 6);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 1].at(/* levelCartIdx */ 11), /* levelCompressedIdx */ 7);

    // Check non-existing key level Cartesian values due to inactive parent cell
    BOOST_CHECK( levelCartToLevelCompressed[/* level */ 1].find(2) ==  levelCartToLevelCompressed[/* level */ 1].end());
    BOOST_CHECK( levelCartToLevelCompressed[/* level */ 1].find(3) ==  levelCartToLevelCompressed[/* level */ 1].end());
    BOOST_CHECK( levelCartToLevelCompressed[/* level */ 1].find(8) ==  levelCartToLevelCompressed[/* level */ 1].end());
    BOOST_CHECK( levelCartToLevelCompressed[/* level */ 1].find(9) ==  levelCartToLevelCompressed[/* level */ 1].end());
}

BOOST_AUTO_TEST_CASE(aFewMoreInactiveParentCells_serial)
{
    const std::string deckString = R"( RUNSPEC
  DIMENS
 -- NX NY NZ cells per x-,y-, and z-direction
     4 3 3 /
  GRID
  COORD -- grid corner coordinates (bounding box), 6*(NX +1)*(NY +1) values.
  0 0 0   0 0 3 -- bottom/top of pillar 0
  1 0 0   1 0 3
  2 0 0   2 0 3
  3 0 0   3 0 3
  4 0 0   4 0 3
  0 1 0   0 1 3
  1 1 0   1 1 3
  2 1 0   2 1 3
  3 1 0   3 1 3
  4 1 0   4 1 3
  0 2 0   0 2 3
  1 2 0   1 2 3
  2 2 0   2 2 3
  3 2 0   3 2 3
  4 2 0   4 2 3
  0 3 0   0 3 3
  1 3 0   1 3 3
  2 3 0   2 3 3
  3 3 0   3 3 3
  4 3 0   4 3 3 -- bottom/top of pillar 19
  /
  ZCORN
-- flat layers, no pinch-outs; NZ=3 → top and bottom per layer
  48*0   -- top layer k=0
  48*1   -- bottom layer k=0
  48*1   -- top layer k=1
  48*2   -- bottom layer k=1
  48*2   -- top layer k=2
  48*3   -- bottom layer k=2
  /
  ACTNUM
-- i=0..3, j=0..2
  1 1 1 1  -- layer k=0, j=0
  1 1 1 1  --            j=1
  0 1 1 0  --            j=2
  0 1 1 1  -- layer k=1, j=0
  0 1 0 0  --            j=1
  1 1 1 1  --            j=2
  1 1 1 1  -- layer k=2, j=0
  1 1 0 0  --            j=1
  1 1 0 1  --            j=2
  /
  PORO
  36*0.15
  /)";

    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckString,
                              /* cells_per_dim = */ {{3,3,3}, {3,3,3}},
                              /* startIJK_vec = */ {{0,0,1}, {2,2,2}},
                              /* endIJK_vec = */ {{3,2,2}, {4,3,3}},
                              /* lgr_name_vec = */ {"LGR1", "LGR2"});

    const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
    const auto levelCartToLevelCompressed = Opm::Lgr::levelCartesianToLevelCompressedMaps(grid, levelCartMapp);

    // LGR1 parent cell global ids = {10, 11, 13}
    // LGR2 parent cell global ids = {26}

    // Level zero grid               Level zero Cartesian indices (only active cells)
    // ------------------------
    //  k = 2 | 24  25   -  26 |        | 32 33  -  35 |
    //        | 22  23   -   - |        | 28 29  -   - |
    //        | 18  19  20  21 |        | 24 25  26 27 |
    // -------------------------        -----------------
    //  k = 1 | 14  15  16  17 |        | 20 21  22  23 |
    //        |  -  13   -   - |        | -  17   -   - |
    //        |  -  10  11  12 |        | -  13  14  15 |
    // -------------------------        -----------------
    // k = 0  | -   8   9   - |         | -   9  10   - |
    //        | 4   5   6   7 |         | 4   5   6   7 |
    //        | 0   1   2   3 |         | 0   1   2   3 |

    // Check the size of the map (27 active cells in a level zero grid of dimensions 4x3x3)
    BOOST_CHECK_EQUAL(levelCartToLevelCompressed[/* level */ 0].size(), /* level active cells count */ 27);

    // Check a few level Cartesian indices match the level active compressed indices
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 0].at(/* levelCartIdx */ 9), /* levelCompressedIdx */ 8);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 0].at(/* levelCartIdx */ 17), /* levelCompressedIdx */ 13);

    // Check non-existing key level Cartesian values due to inactive cell
    BOOST_CHECK( levelCartToLevelCompressed[/* level */ 0].find(30) ==  levelCartToLevelCompressed[/* level */ 0].end());

    // LGR1 parent cell global ids   |  - 13  - |         LGR1 dims = 9x6x3
    //                               |  - 10 11 |
    // LGR1 local/compressed indices                |     LGR1 level Cartesian indices                  LGR1 level Cartesian indices
    //                                              |     Note: INACTIVE PART IS ALSO ILLUSTRATED       ONLY ACTIVE CELLS
    // k = 2              | 78  79  80 |            |     153 154 155 |156 157 158 |159 160 161 |       -   -   - | 156 157 158 |  -   -   - |
    //          INACTIVE  | 75  76  77 | INACTIVE   |     144 145 146 |147 148 149 |150 151 152 |       -   -   - | 147 148 149 |  -   -   - |
    //                    | 72  73  74 |            |     135 136 137 |138 139 140 |141 142 143 |       -   -   - | 138 139 140 |  -   -   - |
    //         --------------------------------------     ---------------------------------------      ---------------------------------------
    //                    | 24  25  26 | 51  52  53 |     126 127 128 |129 130 131 |132 133 134 |       -   -   - |129 130 131 |132 133 134 |
    //          INACTIVE  | 21  22  23 | 48  49  50 |     117 118 119 |120 121 122 |123 124 125 |       -   -   - |120 121 122 |123 124 125 |
    //                    | 18  19  20 | 45  46  47 |     108 109 110 |111 112 113 |114 115 116 |       -   -   - |111 112 113 |114 115 116 |
    // ---------------------------------------------      ---------------------------------------      ---------------------------------------
    // k = 1              | 69  70  71 |            |      99 100 101 |102 103 104 |105 106 107 |       -   -   - |102 103 104 |  -   -   - |
    //          INACTIVE  | 66  67  68 | INACTIVE   |      90  91  92 | 93  94  95 | 96  97  98 |       -   -   - | 93  94  95 |  -   -   - |
    //                    | 63  64  65 |            |      81  82  83 | 84  85  86 | 87  88  89 |       -   -   - | 84  85  86 |  -   -   - |
    //         -------------------------------------      ---------------------------------------      ---------------------------------------
    //                    | 15  16  17 | 42  43  44 |      72  73  74 | 75  76  77 | 78  79  80 |       -   -   - | 75  76  77 | 78  79  80 |
    //          INACTIVE  | 12  13  14 | 39  40  41 |      63  64  65 | 66  67  68 | 69  70  71 |       -   -   - | 66  67  68 | 69  70  71 |
    //                    |  9  10  11 | 36  37  38 |      54  55  56 | 57  58  59 | 60  61  62 |       -   -   - | 57  58  59 | 60  61  62 |
    // ----------------------------------------------     ---------------------------------------     ---------------------------------------
    // k = 0              | 60  61  62 |            |      45  46  47 | 48  49  50 | 51  52  53 |       -   -   - | 48  49  50 |  -   -   - |
    //          INACTIVE  | 57  58  59 | INACTIVE   |      36  37  38 | 39  40  41 | 42  43  44 |       -   -   - | 39  40  41 |  -   -   - |
    //                    | 54  55  56 |            |      27  28  29 | 30  31  32 | 33  34  35 |       -   -   - | 30  31  32 |  -   -   - |
    //         --------------------------------------     ---------------------------------------     ---------------------------------------
    //                    |  6   7   8 | 33  34  35 |      18  19  20 | 21  22  23 | 24  25  26 |       -   -   - | 21  22  23 | 24  25  26 |
    //          INACTIVE  |  3   4   5 | 30  31  32 |       9  10  11 | 12  13  14 | 15  16  17 |       -   -   - | 12  13  14 | 15  16  17 |
    //                    |  0   1   2 | 27  28  29 |       0   1   2 |  3   4   5 |  6   7   8 |       -   -   - |  3   4   5 |  6   7   8 |
    // --------------------------------------------------------------------------------------------------------------------------------------

    // Check the size of the map for LGR1 (3 active parent cells, refined each into 3x3x3 children)
    BOOST_CHECK_EQUAL(levelCartToLevelCompressed[/* level */ 1].size(), /* level active cells count */ 3*(3*3*3));

    // Check a few  the level Cartesian indices match the level active compressed indices
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 1].at(/* levelCartIdx */ 13), /* levelCompressedIdx */ 4);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 1].at(/* levelCartIdx */ 25), /* levelCompressedIdx */ 34);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 1].at(/* levelCartIdx */ 40), /* levelCompressedIdx */ 58);

    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 1].at(/* levelCartIdx */ 60), /* levelCompressedIdx */ 36);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 1].at(/* levelCartIdx */ 67), /* levelCompressedIdx */ 13);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 1].at(/* levelCartIdx */ 95), /* levelCompressedIdx */ 68);

    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 1].at(/* levelCartIdx */ 122), /* levelCompressedIdx */ 23);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 1].at(/* levelCartIdx */ 123), /* levelCompressedIdx */ 48);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 1].at(/* levelCartIdx */ 148), /* levelCompressedIdx */ 76);

    // Check non-existing key level Cartesian values due to inactive parent cell
    BOOST_CHECK( levelCartToLevelCompressed[/* level */ 1].find(37) ==  levelCartToLevelCompressed[/* level */ 1].end());
    BOOST_CHECK( levelCartToLevelCompressed[/* level */ 1].find(106) ==  levelCartToLevelCompressed[/* level */ 1].end());
    BOOST_CHECK( levelCartToLevelCompressed[/* level */ 1].find(153) ==  levelCartToLevelCompressed[/* level */ 1].end());

    // LGR2 parent cell global ids   | -  26 |
    // LGR2 local/compressed indices   |     LGR2 level Cartesian indices            | LGR2 level Cartesian indices
    //                                 |     Note: INACTIVE PART IS ALSO ILLUSTRATED | ONLY ACTIVE cells
    // k = 2   |          | 24  25  26 |      48  49  50 | 51  52  53 |                 -   -   - | 51  52  53 |
    //         | INACTIVE | 21  22  23 |      42  43  44 | 45  46  47 |                 -   -   - | 45  46  47 |
    //         |          | 18  19  20 |      36  37  38 | 39  40  41 |                 -   -   - | 39  40  41 |
    // ---------------------------------     -------------------------                 -------------------------
    // k = 1   |          | 15  16  17 |      30  31  32 | 33  34  35 |                 -   -   - | 33  34  35 |
    //         | INACTIVE | 12  13  14 |      24  25  26 | 27  28  29 |                 -   -   - | 27  28  29 |
    //         |          |  9  10  11 |      18  19  20 | 21  22  23 |                 -   -   - | 21  22  23 |
    // ---------------------------------     -------------------------                 -------------------------
    // k = 0   |          |  6   7   8 |      12  13  14 | 15  16  17 |                 -   -   - | 15  16  17 |
    //         | INACTIVE |  3   4   5 |       6   7   8 |  9  10  11 |                 -   -   - |  9  10  11 |
    //         |          |  0   1   2 |       0   1   2 |  3   4   5 |                 -   -   - |  3   4   5 |
    // ---------------------------------     -------------------------

    // Check the size of the map for LGR2 (1 active parent cell, refined each into 3x3x3 children)
    BOOST_CHECK_EQUAL(levelCartToLevelCompressed[/* level */ 2].size(), /* level active cells count */ 1*(3*3*3));

    // Check a few level Cartesian indices match the level active compressed indices
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 2].at(/* levelCartIdx */ 3), /* levelCompressedIdx */ 0);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 2].at(/* levelCartIdx */ 10), /* levelCompressedIdx */ 4);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 2].at(/* levelCartIdx */ 17), /* levelCompressedIdx */ 8);

    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 2].at(/* levelCartIdx */ 23), /* levelCompressedIdx */ 11);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 2].at(/* levelCartIdx */ 28), /* levelCompressedIdx */ 13);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 2].at(/* levelCartIdx */ 33), /* levelCompressedIdx */ 15);

    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 2].at(/* levelCartIdx */ 40), /* levelCompressedIdx */ 19);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 2].at(/* levelCartIdx */ 46), /* levelCompressedIdx */ 22);
    BOOST_CHECK_EQUAL( levelCartToLevelCompressed[/* level */ 2].at(/* levelCartIdx */ 52), /* levelCompressedIdx */ 25);

    // Check non-existing key level Cartesian values due to inactive parent cells
    BOOST_CHECK( levelCartToLevelCompressed[/* level */ 2].find(7) ==  levelCartToLevelCompressed[/* level */ 2].end());
    BOOST_CHECK( levelCartToLevelCompressed[/* level */ 2].find(19) ==  levelCartToLevelCompressed[/* level */ 2].end());
    BOOST_CHECK( levelCartToLevelCompressed[/* level */ 2].find(49) ==  levelCartToLevelCompressed[/* level */ 2].end());
}
