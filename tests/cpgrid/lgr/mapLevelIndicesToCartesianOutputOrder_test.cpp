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

#define BOOST_TEST_MODULE MapLevelIndicesToCartesianOutputOrderTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>
#include <opm/grid/cpgrid/LgrOutputHelpers.hpp>
#include <tests/cpgrid/lgr/LgrChecks.hpp>

#include <algorithm> // for std::minmax_element, std::is_sorted
#include <string>
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

void checkOutputOrderLevelGrids(const Dune::CpGrid& grid,
                                const Opm::LevelCartesianIndexMapper<Dune::CpGrid>& levelCartMapp)
{
    for (int level = 1; level <= grid.maxLevel(); ++level) {

        const int lgr_cells = grid.levelGridView(level).size(0);
        if (lgr_cells) { // level can be empty in one process

            const auto toOutput = Opm::Lgr::mapLevelIndicesToCartesianOutputOrder(grid, levelCartMapp, level);

            // Create a vector with the cartesian indices sorted by element index.
            std::vector<int> simulatorContainer; // reordered linearized level Cartesian indices
            simulatorContainer.resize(lgr_cells);

            for (const auto& element : Dune::elements(grid.levelGridView(level))) {
                simulatorContainer[element.index()] = levelCartMapp.cartesianIndex(element.index(), level);
            }

            // Use toOutput to reorder simulatorContainer and test that the cartesian indices are now ascending.
            const auto outputContainer = Opm::Lgr::reorderForOutput(simulatorContainer, toOutput);

            const auto [minToOutput, maxToOutput] = std::minmax_element(toOutput.begin(), toOutput.end());
            BOOST_CHECK_EQUAL( *minToOutput, 0);
            BOOST_CHECK_EQUAL( *maxToOutput, lgr_cells-1 ); // shifted by -1 since indexing starts at 0

            BOOST_CHECK( std::is_sorted(outputContainer.begin(), outputContainer.end()) );
        }
    }
}

BOOST_AUTO_TEST_CASE(simpleTestReOrderLgr_serial)
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

    bool isSerial = grid.comm().size() == 1;
    if (isSerial) {
        // LGR1 parent cell global ids = {0, 1}, inactive cell in the middle "is not seen"
        // Level zero grid 3x1x1: k = 0 |  0   -   1
        const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
        // Level zero grid
        // --------------------
        // k = 0  0   -   1 |
        //-----------------------------------------------------------------------
        // LGR1 serial compressed idxs  | LGR1 Cartesian indices, []: inactive  |     LGR1 compressed indxs assumed by output:
        // ----------------------------------------------------------------------
        // k = 0     2  3 | INAC | 6  7 |       k = 0     6  7 | [8] [9] |10 11 |     k = 0     4 5 | - - | 6 7 |
        //           0  1 |      | 4  5 |                 0  1 | [2] [3] | 4  5 |               0 1 | - - | 2 3 |
        checkOutputOrderLevelGrids(grid, levelCartMapp);
    }
}

BOOST_AUTO_TEST_CASE(aFewMoreInactiveParentCellsTestCase)
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

    bool isParallel = grid.comm().size() > 1;
    if (isParallel) {
        grid.loadBalance(/*overlapLayers*/ 1,
                         /*partitionMethod*/ Dune::PartitionMethod::zoltanGoG,
                         /*imbalanceTol*/ 1.1,
                         /*level*/ 0);

        grid.addLgrsUpdateLeafView( /* cells_per_dim = */ {{3,3,3}, {3,3,3}},
                                    /* startIJK_vec = */ {{0,0,1}, {2,2,2}},
                                    /* endIJK_vec = */ {{3,2,2}, {4,3,3}},
                                    /* lgr_name_vec = */ {"LGR1", "LGR2"});
        grid.syncDistributedGlobalCellIds();
    }

    // LGR1 parent cell global ids = {10, 11, 13}
    // LGR2 parent cell global ids = {26}

    // Level zero grid
    // --------------------
    // k = 0  -   8   9   - |  k = 1 14  15  16  17 | k = 2 24  25   -  26 |
    //        4   5   6   7 |         -  13   -   - |       22  23   -   - |
    //        0   1   2   3 |         -  10  11  12 |       18  19  20  21 |

    // LGR1 parent cell global ids   |  - 13  - |
    //                               |  - 10 11 |
    // LGR1 local/compressed indices                |     LGR1 level Cartesian indices                  LGR1 compressed indxs assumed by output:
    //                                              |     Note: INACTIVE PART IS ALSO ILLUSTRATED
    // k = 2              | 78  79  80 |            |     153 154 155 |156 157 158 |159 160 161 |       -   -   - | 78  79  80 |  -   -   - |
    //          INACTIVE  | 75  76  77 | INACTIVE   |     144 145 146 |147 148 149 |150 151 152 |       -   -   - | 75  76  77 |  -   -   - |
    //                    | 72  73  74 |            |     135 136 137 |138 139 140 |141 142 143 |       -   -   - | 72  73  74 |  -   -   - |
    //         --------------------------------------     ---------------------------------------      ---------------------------------------
    //                    | 24  25  26 | 51  52  53 |     126 127 128 |129 130 131 |132 133 134 |       -   -   - | 66  67  68 | 69  70  71 |
    //          INACTIVE  | 21  22  23 | 48  49  50 |     117 118 119 |120 121 122 |123 124 125 |       -   -   - | 60  61  62 | 63  64  65 |
    //                    | 18  19  20 | 45  46  47 |     108 109 110 |111 112 113 |114 115 116 |       -   -   - | 54  55  56 | 57  58  59 |
    // ---------------------------------------------      ---------------------------------------      ---------------------------------------
    // k = 1              | 69  70  71 |            |      99 100 101 |102 103 104 |105 106 107 |       -   -   - | 51  52  53 |  -   -   - |
    //          INACTIVE  | 66  67  68 | INACTIVE   |      90  91  92 | 93  94  95 | 96  97  98 |       -   -   - | 48  49  50 |  -   -   - |
    //                    | 63  64  65 |            |      81  82  83 | 84  85  86 | 87  88  89 |       -   -   - | 45  46  47 |  -   -   - |
    //         -------------------------------------      ---------------------------------------      ---------------------------------------
    //                    | 15  16  17 | 42  43  44 |      72  73  74 | 75  76  77 | 78  79  80 |       -   -   - | 39  40  41 | 42  43  44 |
    //          INACTIVE  | 12  13  14 | 39  40  41 |      63  64  65 | 66  67  68 | 69  70  71 |       -   -   - | 33  34  35 | 36  37  38 |
    //                    |  9  10  11 | 36  37  38 |      54  55  56 | 57  58  59 | 60  61  62 |       -   -   - | 27  28  29 | 30  31  32 |
    // ----------------------------------------------     ---------------------------------------     ---------------------------------------
    // k = 0              | 60  61  62 |            |      45  46  47 | 48  49  50 | 51  52  53 |       -   -   - | 24  25  26 |  -   -   - |
    //          INACTIVE  | 57  58  59 | INACTIVE   |      36  37  38 | 39  40  41 | 42  43  44 |       -   -   - | 21  22  23 |  -   -   - |
    //                    | 54  55  56 |            |      27  28  29 | 30  31  32 | 33  34  35 |       -   -   - | 18  19  20 |  -   -   - |
    //         --------------------------------------     ---------------------------------------     ---------------------------------------
    //                    |  6   7   8 | 33  34  35 |      18  19  20 | 21  22  23 | 24  25  26 |       -   -   - | 12  13  14 | 15  16  17 |
    //          INACTIVE  |  3   4   5 | 30  31  32 |       9  10  11 | 12  13  14 | 15  16  17 |       -   -   - |  6   7   8 |  9  10  11 |
    //                    |  0   1   2 | 27  28  29 |       0   1   2 |  3   4   5 |  6   7   8 |       -   -   - |  0   1   2 |  3   4   5 |
    // --------------------------------------------------------------------------------------------------------------------------------------

    // LGR2 parent cell global ids   | -  26 |
    // LGR2 local/compressed indices   |     LGR2 level Cartesian indices
    //                                 |     Note: INACTIVE PART IS ALSO ILLUSTRATED | LGR2 compressed indxs assumed by output:
    // k = 2   |          | 24  25  26 |      48  49  50 | 51  52  53 |                 -   -   - | 25  26  27 |
    //         | INACTIVE | 21  22  23 |      42  43  44 | 45  46  47 |                 -   -   - | 22  23  24 |
    //         |          | 18  19  20 |      36  37  38 | 39  40  41 |                 -   -   - | 19  20  21 |
    // ---------------------------------     -------------------------                 -------------------------
    // k = 1   |          | 15  16  17 |      30  31  32 | 33  34  35 |                 -   -   - | 16  17  18 |
    //         | INACTIVE | 12  13  14 |      24  25  26 | 27  28  29 |                 -   -   - | 13  14  15 |
    //         |          |  9  10  11 |      18  19  20 | 21  22  23 |                 -   -   - | 10  11  12 |
    // ---------------------------------     -------------------------                 -------------------------
    // k = 0   |          |  6   7   8 |      12  13  14 | 15  16  17 |                 -   -   - |  6   7   8 |
    //         | INACTIVE |  3   4   5 |       6   7   8 |  9  10  11 |                 -   -   - |  3   4   5 |
    //         |          |  0   1   2 |       0   1   2 |  3   4   5 |                 -   -   - |  0   1   2 |
    // ---------------------------------     -------------------------

    const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
    checkOutputOrderLevelGrids(grid, levelCartMapp);
}
