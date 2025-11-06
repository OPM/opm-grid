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

#define BOOST_TEST_MODULE LgrsSharingFacesTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>

#include <tests/cpgrid/lgr/LgrChecks.hpp>

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

void checkRefinedFaceCountPerGrid(const Dune::CpGrid& grid,
                                  const std::vector<int>& expected_leafFaceCount_per_level)
{
    for (int level = 0; level <= grid.maxLevel(); ++level) {
        BOOST_CHECK_EQUAL( grid.currentData()[level]->numFaces(), expected_leafFaceCount_per_level[level]);
    }
    BOOST_CHECK_EQUAL( grid.currentData().back()->numFaces(), expected_leafFaceCount_per_level.back());
}

void checkExpectedFaceCentersAndAreas(const Dune::CpGrid& grid,
                                      const std::vector<std::vector<std::array<double,3>>> expected_leafFaceCenters_per_level,
                                      const std::vector<std::vector<double>>& expected_leafFaceAreas_per_level)
{
    // Check per refined level grid. Skip level 0
    for (int level = 1; level <= grid.maxLevel(); ++level) {
        for (const auto& element : Dune::elements(grid.levelGridView(level))) {
            if (element.isLeaf()) { // not needed for these test cases (no nested refinement)
                for (const auto& intersection : Dune::intersections(grid.levelGridView(level), element)) {
                    BOOST_CHECK( std::find( expected_leafFaceCenters_per_level[element.level()].begin(),
                                            expected_leafFaceCenters_per_level[element.level()].end(),
                                            std::array{intersection.geometry().center()[0],
                                                       intersection.geometry().center()[1],
                                                       intersection.geometry().center()[2]}) !=
                                 expected_leafFaceCenters_per_level[element.level()].end() );
                }
            }
        }
    }
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& intersection : Dune::intersections(grid.leafGridView(), element)) {

            std::array<double,3> center{intersection.geometry().center()[0],
                                        intersection.geometry().center()[1],
                                        intersection.geometry().center()[2]};

            BOOST_CHECK( std::find_if( expected_leafFaceCenters_per_level[element.level()].begin(),
                                       expected_leafFaceCenters_per_level[element.level()].end(),
                                       [&](const auto& c){return Opm::areClose(c, center);}) !=
                         expected_leafFaceCenters_per_level[element.level()].end() );

            // intersection.indexInInside() is
            // 0, 1 for I_FACE false, true respec.
            // 2, 3 for J_FACE false, true respec.
            // 4, 5 for K_FACE false, true respec.
            BOOST_CHECK_EQUAL( intersection.geometry().volume(),
                               expected_leafFaceAreas_per_level[element.level()][intersection.indexInInside()] );
        }
    }
}

void checkExistenceOfIntersectionBetweenAncestors(const Dune::CpGrid& grid)
{
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        int level = element.level();
        for (const auto& intersection : Dune::intersections(grid.leafGridView(), element)) {

            // Check that for an intersection between an elements from different levels
            // there is also an intersection between the origin cells ('oldest' ancestor belonging to level zero).
            // If the level is the same, then check existence of intersection between the father cells. 
            if(intersection.neighbor()) {
               
                bool same_level = intersection.inside().level() == intersection.outside().level();
             
                const auto& ancestorInside = (same_level && level)? intersection.inside().father() : intersection.inside().getOrigin();
                const auto& ancestorOutside = (same_level && level)? intersection.outside().father() : intersection.outside().getOrigin();

                  for (const auto& ancestorInIntersection : Dune::intersections(grid.levelGridView(ancestorInside.level()), ancestorInside)) {
                        bool found = false;
                        // Only check the intersections with the same tag and orientation
                        if (ancestorInIntersection.indexInInside() != intersection.indexInInside()) {
                            continue;
                        }
                        for (const auto& ancestorOutIntersection : Dune::intersections(grid.levelGridView(ancestorOutside.level()), ancestorOutside)) {
                            if (ancestorInIntersection.id() == ancestorOutIntersection.id()) {
                                found = true;
                                break;
                            }
                        }
                        BOOST_CHECK(found);
                    }
                }
            }
        }
    }


// Same for the 3 test cases "lgrsSharing[]Faces_easyCenterAndAreaCheck"
std::vector<double> expected_leafFaceAreas_level2 = {
    16, // I_FACE false
    16, // I_FACE true
    16, // J_FACE false
    16, // J_FACE true
    16, // K_FACE false
    16  // K_FACE true
};

BOOST_AUTO_TEST_CASE(lgrsSharingIFaces_easyCenterAndAreaCheck)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,1,1}, /* cell_sizes = */ {8.,8.,8.});
    // level 0 grid:  LGR1 parent cells {0,1}, refined into 1x2x2 children each
    //                LGR2 parent cell {2}, refined into 2x2x2 children
    //                Cell with index 3 is not involved on any refinement. It appears
    //                on the leaf grid with a total of 9 faces:
    //                - 4 I_FACES false oriented, 1 I_FACE true oriented,
    //                  1 J/K_FACE false, 1 J/K_FACE true.

    // Ensure NY and NZ (y- and z-direction subdivisions) are identical for LGR1 and LGR2
    // since I-face coincides between cells 1 (refined in level 1) and 2 (refined in level 2).
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{1,2,2}, {2,2,2}},
                               /* startIJK_vec = */ {{0,0,0}, {2,0,0}},
                               /* endIJK_vec = */ {{2,1,1}, {3,1,1}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2"});

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{1,2,2}, {2,2,2}},
                           /* lgr_name_vec = */ {"LGR1", "LGR2"},
                           /* isGlobalRefined = */ false);

    // Level 0 grid dimension 4x1x1 -> 5x1x1 I_FACES + 4x2x1 J_FACES + 4x1x2 K_FACES = 21 faces
    // Level 1 grid dimension 2x2x2 -> 3x2x2 I_FACES + 2x3x2 J_FACES + 2x2x3 K_FACES = 36 faces
    // Level 2 grid dimension 2x2x2 -> 3x2x2 I_FACES + 2x3x2 J_FACES + 2x2x3 K_FACES = 36 faces
    // Leaf faces: 36 + 36 - 4(I_FACE shared between LGR1 and LGR2) + 5 ("coarse faces") = 73
    checkRefinedFaceCountPerGrid(grid, /* expected_leafFaceCount_per_level = */ std::vector<int>{ 21, 36, 36, 73});

    // Cell with global id 3 is not involved on any refinement. It appears
    // on the leaf grid with a total of 9 faces
    std::vector<std::array<double,3>> expected_leafFaceCenters_level0 = {
        {24,2,2}, {24,2,6}, {24,6,2}, {24,6,6}, // refined I_FACE false
        {32,4,4}, // coarse I_FACE true
        {28,0,4}, // coarse J_FACE false
        {28,8,4}, // coarse J_FACE true
        {28,4,0}, // coarse K_FACE false
        {28,4,8}  // coarse K_FACE true
    };
    std::vector<double> expected_leafFaceAreas_level0 = {
        16, // (all 4) refined I_FACE false
        64, // coarse I_FACE true
        64, // coarse J_FACE false
        64, // coarse J_FACE true
        64, // coarse K_FACE false
        64  // coarse K_FACE true
    };

    std::vector<std::array<double,3>> expected_leafFaceCenters_level1 = {
        {0,2,2},  {0,2,6},  {0,6,2},  {0,6,6}, // I_FACES
        {8,2,2},  {8,2,6},  {8,6,2},  {8,6,6},
        {16,2,2}, {16,2,6}, {16,6,2}, {16,6,6},
        {4,0,2}, {12,0,2},  {4,0,6}, {12,0,6}, // J_FACES
        {4,4,2}, {12,4,2},  {4,4,6}, {12,4,6},
        {4,8,2}, {12,8,2},  {4,8,6}, {12,8,6},
        {4,2,0}, {12,2,0},  {4,6,0}, {12,6,0}, // K_FACES
        {4,2,4}, {12,2,4},  {4,6,4}, {12,6,4},
        {4,2,8}, {12,2,8},  {4,6,8}, {12,6,8}
    };
    std::vector<double> expected_leafFaceAreas_level1 = {
        16, // I_FACE false
        16, // I_FACE true
        32, // J_FACE false
        32, // J_FACE true
        32, // K_FACE false
        32  // K_FACE true
    };

    std::vector<std::array<double,3>> expected_leafFaceCenters_level2 = {
        {16,2,2}, {16,2,6}, {16,6,2}, {16,6,6}, // I_FACES
        {20,2,2}, {20,2,6}, {20,6,2}, {20,6,6},
        {24,2,2}, {24,2,6}, {24,6,2}, {24,6,6},
        {18,0,2}, {22,0,2}, {18,0,6}, {22,0,6}, // J_FACES
        {18,4,2}, {22,4,2}, {18,4,6}, {22,4,6},
        {18,8,2}, {22,8,2}, {18,8,6}, {22,8,6},
        {18,2,0}, {22,2,0}, {18,6,0}, {22,6,0}, // K_FACES
        {18,2,4}, {22,2,4}, {18,6,4}, {22,6,4},
        {18,2,8}, {22,2,8}, {18,6,8}, {22,6,8}
    };

    std::vector<std::vector<double>> expected_leafFaceAreas_per_level = {
        expected_leafFaceAreas_level0, expected_leafFaceAreas_level1, expected_leafFaceAreas_level2};
    std::vector<std::vector<std::array<double,3>>> expected_leafFaceCenters_per_level = {
        expected_leafFaceCenters_level0, expected_leafFaceCenters_level1, expected_leafFaceCenters_level2};

    checkExpectedFaceCentersAndAreas(grid, expected_leafFaceCenters_per_level, expected_leafFaceAreas_per_level);
    checkExistenceOfIntersectionBetweenAncestors(grid);
}

BOOST_AUTO_TEST_CASE(lgrsSharingJFaces_easyCenterAndAreaCheck)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {1,4,1}, /* cell_sizes = */ {8.,8.,8.});
    // level 0 grid:  LGR1 parent cells {0,1}, refined into 2x1x2 children each
    //                LGR2 parent cell {2}, refined into 2x2x2 children
    //                Cell with index 3 is not involved on any refinement. It appears
    //                on the leaf grid with a total of 9 faces:
    //                - 1 I/K_FACE false oriented, 1 I/K_FACE true oriented,
    //                  4 J_FACE false, 1 J_FACE true.

    // Ensure NX and NZ (x- and z-direction subdivisions) are identical for LGR1 and LGR2.
    // since J-face coincides between cells 1 (refined in level 1) and 2 (refined in level 2).
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,1,2}, {2,2,2}},
                               /* startIJK_vec = */ {{0,0,0}, {0,2,0}},
                               /* endIJK_vec = */ {{1,2,1}, {1,3,1}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2"});

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{2,1,2}, {2,2,2}},
                           /* lgr_name_vec = */ {"LGR1", "LGR2"},
                           /* isGlobalRefined = */ false);

    // Level 0 grid dimension 1x4x1 -> 2x4x1 I_FACES + 1x5x1 J_FACES + 1x4x2 K_FACES = 21 faces
    // Level 1 grid dimension 2x2x2 -> 3x2x2 I_FACES + 2x3x2 J_FACES + 2x2x3 K_FACES = 36 faces
    // Level 2 grid dimension 2x2x2 -> 3x2x2 I_FACES + 2x3x2 J_FACES + 2x2x3 K_FACES = 36 faces
    // Leaf faces: 36 + 36 - 4(J_FACE shared between LGR1 and LGR2) + 5 ("coarse faces") = 73
    checkRefinedFaceCountPerGrid(grid, /* expected_leafFaceCount_per_level = */ std::vector<int>{ 21, 36, 36, 73});

    // Cell with global id 3 is not involved on any refinement. It appears
    // on the leaf grid with a total of 9 faces
    std::vector<std::array<double,3>> expected_leafFaceCenters_level0 = {
        {0,28,4}, // coarse I_FACE false
        {8,28,4}, // coarse I_FACE true
        {2,24,2}, {2,24,6}, {6,24,2},{6,24,6}, // refined J_FACE false
        {4,32,4}, // coarse J_FACE true
        {4,28,0}, // coarse K_FACE false
        {4,28,8}, // coarse K_FACE true
    };
    std::vector<double> expected_leafFaceAreas_level0 = {
        64, // coarse I_FACE false
        64, // coarse I_FACE true
        16, // (all 4) refined J_FACE false
        64, // coarse J_FACE true
        64, // coarse K_FACE false
        64  // coarse K_FACE true
    };

    std::vector<std::array<double,3>> expected_leafFaceCenters_level1 = {
        {0,4,2},  {0,4,6}, {0,12,2}, {0,12,6}, // I_FACES
        {4,4,2},  {4,4,6}, {4,12,2}, {4,12,6},
        {8,4,2},  {8,4,6}, {8,12,2}, {8,12,6},
        {2,0,2},  {2,0,6},  {6,0,2},  {6,0,6}, // J_FACES
        {2,8,2},  {2,8,6},  {6,8,2},  {6,8,6},
        {2,16,2}, {2,16,6}, {6,16,2}, {6,16,6},
        {2,4,0},  {6,4,0}, {2,12,0}, {6,12,0}, // K_FACES
        {2,4,4},  {6,4,4}, {2,12,4}, {6,12,4},
        {2,4,8},  {6,4,8}, {2,12,8}, {6,12,8}
    };
    std::vector<double> expected_leafFaceAreas_level1 = {
        32, // I_FACE false
        32, // I_FACE true
        16, // J_FACE false
        16, // J_FACE true
        32, // K_FACE false
        32  // K_FACE true
    };

    std::vector<std::array<double,3>> expected_leafFaceCenters_level2 = {
        {0,18,2}, {0,22,2}, {0,18,6}, {0,22,6},
        {4,18,2}, {4,22,2}, {4,18,6}, {4,22,6},
        {8,18,2}, {8,22,2}, {8,18,6}, {8,22,6}, // I_FACES
        {2,16,2}, {2,16,6}, {6,16,2}, {6,16,6},
        {2,20,2}, {2,20,6}, {6,20,2}, {6,20,6},
        {2,24,2}, {2,24,6}, {6,24,2}, {6,24,6},  // J_FACES
        {2,18,0}, {2,22,0}, {6,18,0}, {6,22,0},
        {2,18,4}, {2,22,4}, {6,18,4}, {6,22,4},
        {2,18,8}, {2,22,8}, {6,18,8}, {6,22,8} // K_FACES
    };

    std::vector<std::vector<double>> expected_leafFaceAreas_per_level = {
        expected_leafFaceAreas_level0, expected_leafFaceAreas_level1, expected_leafFaceAreas_level2};
    std::vector<std::vector<std::array<double,3>>> expected_leafFaceCenters_per_level = {
        expected_leafFaceCenters_level0, expected_leafFaceCenters_level1, expected_leafFaceCenters_level2};

    checkExpectedFaceCentersAndAreas(grid, expected_leafFaceCenters_per_level, expected_leafFaceAreas_per_level);
    checkExistenceOfIntersectionBetweenAncestors(grid);
}

BOOST_AUTO_TEST_CASE(lgrsSharingKFaces_easyCenterAndAreaCheck)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {1,1,4}, /* cell_sizes = */ {8.,8.,8.});
    // level 0 grid:  LGR1 parent cells {0,1}, refined into 2x2x1 children each
    //                LGR2 parent cell {2}, refined into 2x2x2 children
    //                Cell with index 3 is not involved on any refinement. It appears
    //                on the leaf grid with a total of 9 faces:
    //                - 1 I/J_FACE false oriented, 1 I/J_FACE true oriented,
    //                  4 K_FACE false, 1 K_FACE true.

    // Ensure NX and NY (x- and y-direction subdivisions) are identical for LGR1 and LGR2
    // since K-face coincides between cells 1 (refined in level 1) and 2 (refined in level 2).
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,1}, {2,2,2}},
                               /* startIJK_vec = */ {{0,0,0}, {0,0,2}},
                               /* endIJK_vec = */ {{1,1,2}, {1,1,3}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2"});

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{2,2,1}, {2,2,2}},
                           /* lgr_name_vec = */ {"LGR1", "LGR2"},
                           /* isGlobalRefined = */ false);

    // Level 0 grid dimension 1x1x4 -> 2x1x4 I_FACES + 1x2x4 J_FACES + 1x1x5 K_FACES = 21 faces
    // Level 1 grid dimension 2x2x2 -> 3x2x2 I_FACES + 2x3x2 J_FACES + 2x2x3 K_FACES = 36 faces
    // Level 2 grid dimension 2x2x2 -> 3x2x2 I_FACES + 2x3x2 J_FACES + 2x2x3 K_FACES = 36 faces
    // Leaf faces: 36 + 36 - 4(K_FACE shared between LGR1 and LGR2) + 5 ("coarse faces") = 73
    checkRefinedFaceCountPerGrid(grid, /* expected_leafFaceCount_per_level = */ std::vector<int>{ 21, 36, 36, 73});

    // Cell with global id 3 is not involved on any refinement. It appears
    // on the leaf grid with a total of 9 faces
    std::vector<std::array<double,3>> expected_leafFaceCenters_level0 = {
        {0,4,28}, // coarse I_FACE false
        {8,4,28}, // coarse I_FACE true
        {4,0,28}, // coarse J_FACE false
        {4,8,28}, // coarse J_FACE true
        {2,2,24}, {2,6,24}, {6,2,24}, {6,6,24}, // refined K_FACE false
        {4,4,32}, // coarse K_FACE true
    };
    std::vector<double> expected_leafFaceAreas_level0 = {
        64, // coarse I_FACE false
        64, // coarse I_FACE true
        64, // coarse J_FACE false
        64, // coarse J_FACE true
        16, // (all 4) refined K_FACE false
        64  // coarse K_FACE true
    };

    std::vector<std::array<double,3>> expected_leafFaceCenters_level1 = {
        {0,2,4}, {0,6,4}, {0,2,12}, {0,6,12}, // I_FACES
        {4,2,4}, {4,6,4}, {4,2,12}, {4,6,12},
        {8,2,4}, {8,6,4}, {8,2,12}, {8,6,12},
        {2,0,4}, {2,0,12}, {6,0,4}, {6,0,12}, // J_FACES
        {2,4,4}, {2,4,12}, {6,4,4}, {6,4,12},
        {2,8,4}, {2,8,12}, {6,8,4}, {6,8,12},
        {2,2,0}, {2,6,0}, {6,2,0}, {6,6,0}, // K_FACES
        {2,2,8}, {2,6,8}, {6,2,8}, {6,6,8},
        {2,2,16}, {2,6,16}, {6,2,16}, {6,6,16}
    };
    std::vector<double> expected_leafFaceAreas_level1 = {
        32, // I_FACE false
        32, // I_FACE true
        32, // J_FACE false
        32, // J_FACE true
        16, // K_FACE false
        16  // K_FACE true
    };

    std::vector<std::array<double,3>> expected_leafFaceCenters_level2 = {
        {2,2,16}, {2,6,16},{6,2,16},{6,6,16},
        {2,2,20},{2,6,20}, {6,2,20},{6,6,20},
        {2,2,24},{2,6,24},{6,2,24},{6,6,24},
        {0,2,18}, {0,2,22}, {0,6,18}, {0,6,22},
        {4,2,18}, {4,2,22}, {4,6,18}, {4,6,22},
        {8,2,18}, {8,2,22}, {8,6,18}, {8,6,22},
        {2,0,18}, {2,0,22}, {6,0,18}, {6,0,22},
        {2,4,18}, {2,4,22}, {6,4,18}, {6,4,22},
        {2,8,18}, {2,8,22}, {6,8,18}, {6,8,22}
    };

    std::vector<std::vector<std::array<double,3>>> expected_leafFaceCenters_per_level = {
        expected_leafFaceCenters_level0, expected_leafFaceCenters_level1, expected_leafFaceCenters_level2};
    std::vector<std::vector<double>> expected_leafFaceAreas_per_level = {
        expected_leafFaceAreas_level0, expected_leafFaceAreas_level1, expected_leafFaceAreas_level2};

    checkExpectedFaceCentersAndAreas(grid, expected_leafFaceCenters_per_level, expected_leafFaceAreas_per_level);
    checkExistenceOfIntersectionBetweenAncestors(grid);
}

BOOST_AUTO_TEST_CASE(lgrsSharingIJKFaces)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.,1.,1.});
    //                          LGR1 parent cells    LGR2 parent cells   LGR3 parent cells   LGR4 parent cells
    // k = 2   32 33 34 35 |                                  33 34
    //         28 29 30 31 |          29 30                                                 28
    //         24 25 26 27 |          25 26
    // --------------------|    ------------------  ------------------  ------------------  ------------------
    // k = 1   20 21 22 23 |                                  21 22
    //         16 17 18 19 |          17 18                                                 16
    //         12 13 14 15 |          13 14
    // --------------------|    ------------------  ------------------  ------------------  ------------------
    // k = 0    8  9 10 11 |                                   9 10
    //          4  5  6  7 |                                                                 4
    //          0  1  2  3 |                                                 1  2
    //---------------------|
    // I-faces coincide between cells (16–17) and (28–29).
    // Ensure NY and NZ (y- and z-direction subdivisions) are identical for LGR1 and LGR4.
    // J-faces coincide between cells (17–21), (18–22), (29–33), and (30–34).
    // Ensure NX and NZ (x- and z-direction subdivisions) are identical for LGR1 and LGR2.
    // K-faces coincide between cells (1–13) and (2–14).
    // Ensure NX and NY (x- and y-direction subdivisions) are identical for LGR1 and LGR3.
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{4,1,2}, {4,3,2}, {4,1,3}, {5,1,2}},
                               /* startIJK_vec = */ {{1,0,1}, {1,2,0}, {1,0,0}, {0,1,0}},
                               /* endIJK_vec = */ {{3,2,3}, {3,3,3}, {3,1,1}, {1,2,3}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2", "LGR3", "LGR4"});

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{4,1,2}, {4,3,2}, {4,1,3}, {5,1,2}},
                           /* lgr_name_vec = */ {"LGR1", "LGR2", "LGR3", "LGR4"},
                           /* isGlobalRefined = */ false);
    
    checkExistenceOfIntersectionBetweenAncestors(grid);
}


BOOST_AUTO_TEST_CASE(lgrsSharingIFaces_nested)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {2,1,1}, /* cell_sizes = */ {1.,1.,1.});

    // I-face coincide between cells 1 (with children in level 1) and cell 2
    // (with children in level 2).
    // Ensure NY and NZ (y- and z-direction subdivisions) are identical for LGR1 and LGR2.
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{3,2,4}, {5,2,4}, {2,2,4}, {1,2,4}},
                               /* startIJK_vec = */ {{0,0,0}, {1,0,0}, {0,0,0}, {2,0,0}},
                               /* endIJK_vec = */ {{1,1,1}, {2,1,1}, {2,2,2}, {4,2,2}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2", "LGR3", "LGR4"},
                               /* lgr_parent_grid_name_vec = */ {"GLOBAL", "GLOBAL", "LGR1",  "LGR1"});

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{3,2,4}, {5,2,4}, {2,2,4}, {1,2,4}},
                           /* lgr_name_vec = */ {"LGR1", "LGR2", "LGR3", "LGR4"},
                           /* isGlobalRefined = */ false,
                           /* preRefineMaxLevel = */ 0,
                           /* isNested = */ true);
    
    checkExistenceOfIntersectionBetweenAncestors(grid);
}

BOOST_AUTO_TEST_CASE(lgrSharingJFaces_nested)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {1,2,1}, /* cell_sizes = */ {1.,1.,1.});

    // J-face coincides between cells (0–1).
    // Ensure NX and NZ (x- and z-direction subdivisions) are identical for LGR1 and LGR2.
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */        {{4,3,2},  {4,5,2}, {4,2,2}, {4,1,2}},
                               /* startIJK_vec = */             {{0,0,0},  {0,1,0}, {0,0,0}, {0,2,0}},
                               /* endIJK_vec = */               {{1,1,1},  {1,2,1}, {2,2,2}, {2,4,2}},
                               /* lgr_name_vec = */             { "LGR1",   "LGR2",  "LGR3", "LGR4"},
                               /* lgr_parent_grid_name_vec = */ {"GLOBAL", "GLOBAL", "LGR1",  "LGR1"});

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{4,3,2}, {4,5,2}, {4,2,2}, {4,1,2}},
                           /* lgr_name_vec = */ {"LGR1", "LGR2", "LGR3", "LGR4"},
                           /* isGlobalRefined = */ false,
                           /* preRefineMaxLevel = */ 0,
                           /* isNested = */ true);
    
    checkExistenceOfIntersectionBetweenAncestors(grid);
}

BOOST_AUTO_TEST_CASE(lgrSharingKFaces_nested)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {1,1,2}, /* cell_sizes = */ {1.,1.,1.});

    // K-face coincides between cells (0–1).
    // Ensure NX and NY (x- and y-direction subdivisions) are identical for LGR1 and LGR2.
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{4,2,3}, {4,2,5}, {4,2,2}, {4,2,1}},
                               /* startIJK_vec = */ {{0,0,0}, {0,0,1}, {0,0,0}, {0,0,2}},
                               /* endIJK_vec = */ {{1,1,1}, {1,1,2}, {2,2,2}, {2,2,4}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2", "LGR3", "LGR4"},
                                /* lgr_parent_grid_name_vec = */ {"GLOBAL", "GLOBAL", "LGR1",  "LGR1"});

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{4,2,3}, {4,2,5}, {4,2,2}, {4,2,1}},
                           /* lgr_name_vec = */ {"LGR1", "LGR2", "LGR3", "LGR4"},
                           /* isGlobalRefined = */ false,
                           /* preRefineMaxLevel = */ 0,
                           /* isNested = */ true);

    checkExistenceOfIntersectionBetweenAncestors(grid);
}
