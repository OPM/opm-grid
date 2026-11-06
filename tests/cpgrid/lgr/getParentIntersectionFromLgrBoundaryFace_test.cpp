/*
  Copyright 2026 Equinor ASA.

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

#define BOOST_TEST_MODULE GetParentIntersectionFromLgrBoundaryFaceTests
#include <boost/test/unit_test.hpp>


#include <tests/cpgrid/lgr/LgrChecks.hpp>

#include <stdexcept>
#include <set>
#include <string>

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

std::set<int> checkBothOrientations(const Dune::CpGrid& grid,
                                    int smallerLevel,
                                    int indexInInsideFromSmallerLevel,
                                    int indexInInsideFromLargerLevel)
{
    std::set<int> parentIntersectionToPointIndices{};
    
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& intersection : Dune::intersections(grid.leafGridView(), element)) {
            if (intersection.neighbor() && (intersection.inside().level() != intersection.outside().level())) {

                const auto& parentIntersection = grid.getParentIntersectionFromLgrBoundaryFace(intersection);

                const auto& parentInterToPoint = grid.currentData().front()->faceToPoint(parentIntersection.id());
                for (const auto& point : parentInterToPoint) {
                    parentIntersectionToPointIndices.insert(point);
                }
                

                if (intersection.inside().level() == smallerLevel)
                    BOOST_CHECK_EQUAL(parentIntersection.indexInInside(), indexInInsideFromSmallerLevel);
                else
                    BOOST_CHECK_EQUAL(parentIntersection.indexInInside(), indexInInsideFromLargerLevel);
            }
        }
    }
    return parentIntersectionToPointIndices;
}

// Recall that intersection.indexInInside() returns 0,1,2,3,4, or 5,
// mapping the face tag (I/J/K) and orientation (true/false) as follows:
// I false <-> 0
// I true  <-> 1
// J false <-> 2
// J true  <-> 3
// K true  <-> 4
// K false <-> 5

BOOST_AUTO_TEST_CASE(parentIntersection_I_FACE)
{
    const std::string deck_string =
        R"(
RUNSPEC
DIMENS
  2 1 1 /
GRID
DX
  2*8 /
DY
  2*4 /
DZ
  2*2 /
TOPS
  2*0 /
PORO
  2*0.15 /
)";
    // Grid: NX=2, NY=1, NZ=1  
    // Pillars = 3x2 = 6 [(NX+1) = 3 pillars in x, (NY+1) = 2 pillars in y]

    // Pillar indices (j=0):                   Pillar indices (j=1):
    //
    // k=1    1 -------- 3 -------- 5          k=1   7 -------- 9 -------- 11
    //        |          |          |                |          |          |
    //        |          |          |                |          |          |
    // k=0    0 -------- 2 -------- 4          k=0   6 -------- 8 -------- 10
    //       i=0        i=1        i=2              i=0        i=1        i=2

    // Elements in level zero grid share an I_FACE with vertex indices = {2,3,8,9}
    std::set<int> expectedCornerIndices{2,8,9,3};
    // From element with index 0, I_FACE true,   (this is intersection.indexInInside() == 1)
    // from element with index 1, I_FACE false.  (this is intersection.indexInInside() == 0)

    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid, deck_string,
                              /* cells_per_dim_vec = */ {{2,2,2}},
                              /* startIJK_vec = */      {{0,0,0}},
                              /* endIJK_vec = */        {{1,1,1}},
                              /* lgr_name_vec = */      {"LGR1"});

   

    // Element with index 0 is refined into 2x2x2 child cells. Some of those share their I_FACE true
    // with the un-refined/coarse element whose index in level zero is 1.

    // if intersection.inside() is the coarse cell,  parentIntersection indexInInside should be 0 (I-)
    // if intersection.inside() is the refined cell, parentIntersection indexInInside should be 1 (I+)
    const std::set<int> parentIntersectionToPoint = checkBothOrientations(grid,
                                                                          /* smallerLevel = */ 0,
                                                                          /* indexInInsideFromSmallerLevel = */ 0,
                                                                          /* indexInInsideFromLargerLevel = */ 1);
    BOOST_CHECK(parentIntersectionToPoint == expectedCornerIndices);
}

BOOST_AUTO_TEST_CASE(parentIntersection_J_FACE)
{
    const std::string deck_string =
        R"(
RUNSPEC
DIMENS
  1 2 1 /
GRID
DX
  2*8 /
DY
  2*4 /
DZ
  2*2 /
TOPS
  2*0 /
PORO
  2*0.15 /
)";
    // Grid: NX=1, NY=2, NZ=1  
    // Pillars = 2x3= 6 [(NX+1) = 2 pillars in x, (NY+1) = 3 pillars in y]

    // Pillar indices (j=0):        Pillar indices (j=1):        Pillar indices (j=2):   
    //
    // k=1    1 -------- 3          k=1   5 -------- 7          k=1   9---------11 
    //        |          |                |          |                |          |        
    //        |          |                |          |                |          |         
    // k=0    0 -------- 2          k=0   4 -------- 6          k=0   8 --------10 
    //       i=0        i=1              i=0        i=1              i=0        i=1       

    // Elements in level zero grid share an J_FACE with vertex indices = {4,6,5,7}
    std::set<int> expectedCornerIndices{4,6,5,7};
    // From element with index 0, J_FACE true,   (this is intersection.indexInInside() == 3)
    // from element with index 1, J_FACE false.  (this is intersection.indexInInside() == 2)

    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid, deck_string,
                              /* cells_per_dim_vec = */ {{2,2,2}},
                              /* startIJK_vec = */      {{0,0,0}},
                              /* endIJK_vec = */        {{1,1,1}},
                              /* lgr_name_vec = */      {"LGR1"});

    // Element with index 0 is refined into 2x2x2 child cells. Some of those share their J_FACE true
    // with the un-refined/coarse element whose index in level zero is 1.

    // if intersection.inside() is the coarse cell,  parentIntersection indexInInside should be 2 (J-)
    // if intersection.inside() is the refined cell, parentIntersection indexInInside should be 3 (J+)
    const auto parentIntersectionToPoint = checkBothOrientations(grid,
                                                                 /* smallerLevel = */ 0,
                                                                 /* indexInInsideFromSmallerLevel = */ 2,
                                                                 /* indexInInsideFromLargerLevel = */ 3);
    BOOST_CHECK(parentIntersectionToPoint == expectedCornerIndices);
}

BOOST_AUTO_TEST_CASE(parentIntersection_K_FACE)
{
    const std::string deck_string =
        R"(
RUNSPEC
DIMENS
  1 1 2 /
GRID
DX
  2*8 /
DY
  2*4 /
DZ
  2*2 /
TOPS
  2*0 /
PORO
  2*0.15 /
)";
    // Grid: NX=1, NY=1, NZ=2  
    // Pillars = 2x2= 4 [(NX+1) = 2 pillars in x, (NY+1) = 2 pillars in y]

    // Pillar indices (j=0):        Pillar indices (j=1):          
    //
    // k=2    2 -------- 5          k=1   8 -------- 11          
    //        |          |                |          |                       
    //        |          |                |          |     
    // k=1    1 -------- 4          k=1   7 -------- 10          
    //        |          |                |          |                       
    //        |          |                |          |                   
    // k=0    0 -------- 3          k=0   6 -------- 9          
    //       i=0        i=1              i=0        i=1

    // Elements in level zero grid share an K_FACE with vertex indices = {1,4,7,10}
    std::set<int> expectedCornerIndices{1,4,7,10};
    // From element with index 0, K_FACE true,   (this is intersection.indexInInside() == 5)
    // from element with index 1, K_FACE false.  (this is intersection.indexInInside() == 4)

    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid, deck_string,
                              /* cells_per_dim_vec = */ {{2,2,2}},
                              /* startIJK_vec = */      {{0,0,0}},
                              /* endIJK_vec = */        {{1,1,1}},
                              /* lgr_name_vec = */      {"LGR1"});

    // Element with index 0 is refined into 2x2x2 child cells. Some of those share their K_FACE true
    // with the un-refined/coarse element whose index in level zero is 1.

    // if intersection.inside() is the coarse cell,  parentIntersection indexInInside should be 4 (K-)
    // if intersection.inside() is the refined cell, parentIntersection indexInInside should be 5 (K+)
    const auto parentIntersectionToPoint = checkBothOrientations(grid,
                                                                 /* smallerLevel = */ 0,
                                                                 /* indexInInsideFromSmallerLevel= */ 4,
                                                                 /* indexInInsideFromLargerLevel = */ 5);
    BOOST_CHECK(parentIntersectionToPoint == expectedCornerIndices);
}

BOOST_AUTO_TEST_CASE(parentIntersectionBoundaryLgrs_X)
{
    const std::string deck_string =
        R"(
RUNSPEC
DIMENS
  2 1 1 /
GRID
DX
  2*8 /
DY
  2*4 /
DZ
  2*2 /
TOPS
  2*0 /
PORO
  2*0.15 /
)";
    // Grid: NX=2, NY=1, NZ=1  
    // Pillars = 3x2 = 6 [(NX+1) = 3 pillars in x, (NY+1) = 2 pillars in y]

    // Pillar indices (j=0):                   Pillar indices (j=1):
    //
    // k=1    1 -------- 3 -------- 5          k=1   7 -------- 9 -------- 11
    //        |          |          |                |          |          |
    //        |          |          |                |          |          |
    // k=0    0 -------- 2 -------- 4          k=0   6 -------- 8 -------- 10
    //       i=0        i=1        i=2              i=0        i=1        i=2

    // Elements in level zero grid share an I_FACE with vertex indices = {2,3,8,9}
    std::set<int> expectedCornerIndices{2,8,9,3};
    // From element with index 0, I_FACE true,   (this is intersection.indexInInside() == 1)
    // from element with index 1, I_FACE false.  (this is intersection.indexInInside() == 0)

    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid, deck_string,
                              /* cells_per_dim_vec = */ {{2,2,2}, {3,2,2}},
                              /* startIJK_vec = */      {{0,0,0}, {1,0,0}},
                              /* endIJK_vec = */        {{1,1,1}, {2,1,1}},
                              /* lgr_name_vec = */      {"LGR1", "LGR2"});

   

    // Element with index 0 is refined into 2x2x2 child cells. Some of those share their I_FACE true
    // with some refined cells from LGR2.

    // Element with index 1 is refined into 3x2x2 child cells. Some of those share their I_FACE false
    // with some refined cells from LGR1.

    // if intersection.inside() is in LGR1,  parentIntersection indexInInside should be 1 (I+)
    // if intersection.inside() is in LGR2, parentIntersection indexInInside should be 0 (I-)
    const std::set<int> parentIntersectionToPoint = checkBothOrientations(grid,
                                                                          /* smallerLevel = */ 1,
                                                                          /* indexInInsideromCoarse = */ 1,
                                                                          /* indexInInsideFromRefined = */ 0);
    BOOST_CHECK(parentIntersectionToPoint == expectedCornerIndices);
}

BOOST_AUTO_TEST_CASE(parentIntersectionAtBoundaryLgrs_Y)
{
    const std::string deck_string =
        R"(
RUNSPEC
DIMENS
  1 2 1 /
GRID
DX
  2*8 /
DY
  2*4 /
DZ
  2*2 /
TOPS
  2*0 /
PORO
  2*0.15 /
)";
    // Grid: NX=1, NY=2, NZ=1  
    // Pillars = 2x3= 6 [(NX+1) = 2 pillars in x, (NY+1) = 3 pillars in y]

    // Pillar indices (j=0):        Pillar indices (j=1):        Pillar indices (j=2):   
    //
    // k=1    1 -------- 3          k=1   5 -------- 7          k=1   9---------11 
    //        |          |                |          |                |          |        
    //        |          |                |          |                |          |         
    // k=0    0 -------- 2          k=0   4 -------- 6          k=0   8 --------10 
    //       i=0        i=1              i=0        i=1              i=0        i=1       

    // Elements in level zero grid share an J_FACE with vertex indices = {4,6,5,7}
    std::set<int> expectedCornerIndices{4,6,5,7};
    // From element with index 0, J_FACE true,   (this is intersection.indexInInside() == 3)
    // from element with index 1, J_FACE false.  (this is intersection.indexInInside() == 2)

    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid, deck_string,
                              /* cells_per_dim_vec = */ {{2,2,2}, {2,3,2}},
                              /* startIJK_vec = */      {{0,0,0}, {0,1,0}},
                              /* endIJK_vec = */        {{1,1,1}, {1,2,1}},
                              /* lgr_name_vec = */      {"LGR1",  "LGR2"});

    // Element with index 0 is refined into 2x2x2 child cells. Some of those share their J_FACE true
    // with some refined cells from LGR2.

    // Element with index 1 is refined into 2x3x2 child cells. Some of those share their J_FACE false
    // with some refined cells from LGR1.

    // if intersection.inside() is in LGR1,  parentIntersection indexInInside should be 3 (J+)
    // if intersection.inside() is in LGR2, parentIntersection indexInInside should be 2 (J-)
    const auto parentIntersectionToPoint = checkBothOrientations(grid,
                                                                 /* smallerLevel = */ 1,
                                                                 /* indexInInsideFromSmallerLevel = */ 3,
                                                                 /* indexInInsideFromLargerLevel = */ 2);
    BOOST_CHECK(parentIntersectionToPoint == expectedCornerIndices);
}

BOOST_AUTO_TEST_CASE(parentIntersectionAtBoundaryLgrs_Z)
{
    const std::string deck_string =
        R"(
RUNSPEC
DIMENS
  1 1 2 /
GRID
DX
  2*8 /
DY
  2*4 /
DZ
  2*2 /
TOPS
  2*0 /
PORO
  2*0.15 /
)";
    // Grid: NX=1, NY=1, NZ=2  
    // Pillars = 2x2= 4 [(NX+1) = 2 pillars in x, (NY+1) = 2 pillars in y]

    // Pillar indices (j=0):        Pillar indices (j=1):          
    //
    // k=2    2 -------- 5          k=1   8 -------- 11          
    //        |          |                |          |                       
    //        |          |                |          |     
    // k=1    1 -------- 4          k=1   7 -------- 10          
    //        |          |                |          |                       
    //        |          |                |          |                   
    // k=0    0 -------- 3          k=0   6 -------- 9          
    //       i=0        i=1              i=0        i=1

    // Elements in level zero grid share an K_FACE with vertex indices = {1,4,7,10}
    std::set<int> expectedCornerIndices{1,4,7,10};
    // From element with index 0, K_FACE true,   (this is intersection.indexInInside() == 5)
    // from element with index 1, K_FACE false.  (this is intersection.indexInInside() == 4)

    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid, deck_string,
                              /* cells_per_dim_vec = */ {{2,2,2}, {2,2,3}},
                              /* startIJK_vec = */      {{0,0,0}, {0,0,1}},
                              /* endIJK_vec = */        {{1,1,1}, {1,1,2}},
                              /* lgr_name_vec = */      {"LGR1",  "LGR2"});

    // Element with index 0 is refined into 2x2x2 child cells. Some of those share their K_FACE true
    // with some refined cells from LGR2.

    // Element with index 1 is refined into 2x2x3 child cells. Some of those share their K_FACE false
    // with some refined cells from LGR1.

    // if intersection.inside() is in LGR1, parentIntersection indexInInside should be 5 (K+)
    // if intersection.inside() is in LGR2, parentIntersection indexInInside should be 4 (K-)
    const auto parentIntersectionToPoint = checkBothOrientations(grid,
                                                                 /* smallerLevel = */ 1,
                                                                 /* indexInInsideFromSmallerLevel= */ 5,
                                                                 /* indexInInsideFromLargerLevel = */ 4);
    BOOST_CHECK(parentIntersectionToPoint == expectedCornerIndices);
}

// A vertically faulted LGR boundary. This is the configuration in which
// indexInInside() alone cannot identify the parent face, and in which the
// leaf face sits on the box's LOW side -- both of which a matching,
// unfaulted box boundary hides.
//
//        i=0 (coarse)      i=1 (refined, LGR1)
//   z=0  +-----------+
//        |  (0,0,0)  |- - - +-----------+  z=1
//   z=2  +-----------+      |  (1,0,0)  |
//        |  (0,0,1)  |      |           |
//   z=4  +-----------+ - - -+-----------+  z=3
//        |    ...    |      |    ...    |
//
// Each coarse cell overlaps two cells of the refined column, so the coarse
// cell has two level-0 I+ faces -- and each carries indexInInside() == 1.
BOOST_AUTO_TEST_CASE(parentIntersectionFaultedLowSideBoundary)
{
    // Columns offset by half a cell: column i=0 spans z 0..8, column i=1
    // spans z 1..9. TOPS is per cell, i fastest.
    const std::string deck_string =
        R"(
RUNSPEC
DIMENS
  2 1 4 /
GRID
CARFIN
'LGR1' 2 2 1 1 1 4 2 2 8 /
ENDFIN
DX
  8*8 /
DY
  8*4 /
DZ
  8*2 /
TOPS
  0 1  2 3  4 5  6 7 /
PORO
  8*0.15 /
)";

    // Only the Conforming backend refines across a fault.
    Dune::CpGrid grid;
    grid.setLgrBackend(Opm::Refinement::Backend::Conforming);
    Opm::createGridAndAddLgrs(grid, deck_string,
                              /* cells_per_dim_vec = */ {{2,2,2}},
                              /* startIJK_vec = */      {{1,0,0}},
                              /* endIJK_vec = */        {{2,1,4}},
                              /* lgr_name_vec = */      {"LGR1"});

    const auto& level0 = grid.levelGridView(0);
    int coarseToRefined = 0;
    int ambiguousSides = 0;
    int sideOnlySearchWouldMiss = 0;

    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& intersection : Dune::intersections(grid.leafGridView(), element)) {
            if (!intersection.neighbor() ||
                (intersection.inside().level() == intersection.outside().level()))
            {
                continue;
            }
            ++coarseToRefined;

            const bool insideIsCoarse = intersection.inside().level() == 0;

            // The refined column is at i=1, so the coarse cell sees the
            // interface on its I+ side and the refined cell on its I- side.
            // Getting this backwards is what a normal pointing down-axis does.
            BOOST_CHECK_EQUAL(intersection.indexInInside(), insideIsCoarse ? 1 : 0);
            BOOST_CHECK_EQUAL(intersection.indexInOutside(), insideIsCoarse ? 0 : 1);

            const auto& parentIntersection =
                grid.getParentIntersectionFromLgrBoundaryFace(intersection);

            // The parent face must be the level-0 face between the two cells'
            // level-0 ancestors, not merely some face on the same side.
            BOOST_CHECK(parentIntersection.neighbor());
            BOOST_CHECK_EQUAL(parentIntersection.inside().index(),
                              intersection.inside().getOrigin().index());
            BOOST_CHECK_EQUAL(parentIntersection.outside().index(),
                              intersection.outside().getOrigin().index());
            BOOST_CHECK_EQUAL(parentIntersection.indexInInside(),
                              intersection.indexInInside());

            // Matching on the side alone is genuinely ambiguous here, so the
            // cell-pair criterion above is doing real work: count the level-0
            // faces of this ancestor that share the leaf face's side, and note
            // whether the first of them - what a side-only search returns - is
            // the face we got.
            int sameSideCandidates = 0;
            int firstSameSideId = -1;
            for (const auto& candidate : Dune::intersections(level0, intersection.inside().getOrigin())) {
                if (candidate.indexInInside() != intersection.indexInInside()) {
                    continue;
                }
                if (firstSameSideId < 0) {
                    firstSameSideId = candidate.id();
                }
                ++sameSideCandidates;
            }
            BOOST_CHECK_GT(sameSideCandidates, 0);

            if (sameSideCandidates > 1) {
                ++ambiguousSides;
                if (firstSameSideId != parentIntersection.id()) {
                    ++sideOnlySearchWouldMiss;
                }
            }
        }
    }

    BOOST_CHECK_EQUAL(level0.size(0), 8);

    // Without any coarse<->refined face the loop above proves nothing.
    BOOST_CHECK_GT(coarseToRefined, 0);

    // The fault must actually produce cells with two faces on one side, and a
    // side-only search must actually pick the wrong one for some of them -
    // otherwise this deck is not testing what it claims to.
    BOOST_CHECK_GT(ambiguousSides, 0);
    BOOST_CHECK_GT(sideOnlySearchWouldMiss, 0);
}
