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

#define BOOST_TEST_MODULE ParentCellWithMoreThanOneFacePerTypeTests
#include <boost/test/unit_test.hpp>

#include <dune/common/fvector.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/LgrHelpers.hpp>

#include <tests/cpgrid/lgr/LgrChecks.hpp>


#include <string>

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

using Coordinate = Dune::FieldVector<double, 3>;

void checkNewVertices(const std::set<Coordinate,Opm::Lgr::FieldVectorLess>& collectedVertices,
                      const std::set<Coordinate,Opm::Lgr::FieldVectorLess>& expectedVertices)
{
    BOOST_CHECK_EQUAL( collectedVertices.size(), expectedVertices.size());

    for (const auto& expectedVertex : expectedVertices) {
        auto it = collectedVertices.find(expectedVertex);
        BOOST_CHECK(it != collectedVertices.end());
    }
}

void checkOverlapNewFace(const std::map<int,std::pair<int,std::vector<Coordinate>>>& overlapFaces,
                         int parentFaceIndex,
                         const std::vector<Coordinate>& expectedSortedNewFace)
{
    auto it = overlapFaces.find( parentFaceIndex );
    BOOST_CHECK(it != overlapFaces.end());

    const auto& [refinedFaceIDx, refinedFaceToCoord] = overlapFaces.at( parentFaceIndex );
    
    BOOST_CHECK_EQUAL(refinedFaceToCoord.size(), expectedSortedNewFace.size());

    for (std::size_t i = 0; i < expectedSortedNewFace.size(); ++i) {
        BOOST_CHECK_EQUAL(refinedFaceToCoord[i], expectedSortedNewFace[i]);
    }
}

BOOST_AUTO_TEST_CASE(parentCellWithMoreThanOne_I_FACE_true)
{
    // Level zero grid dims = 2x1x1
    //
    // cell 0
    // bottom face corners (0,0,0), (6,0,0), (0,6,0), (6,6,0)
    //    top face corners (0,0,8), (6,0,8), (0,6,8), (6,6,8)
    //
    // cell 1
    // bottom face corners (6,0,1), (12,0,1), (6,6,1),  (12,6,1)
    //    top face corners (6,0,9), (12,0,9), (12,6,9), (12,6,9)
    const std::string deckString =
        R"(RUNSPEC
DIMENS
 2 1 1 /

GRID

COORD
 0 0 0     0 0 9
 6 0 0     6 0 9
12 0 0    12 0 9

 0 6 0    0 6 9
 6 6 0    6 6 9
12 6 0   12 6 9
/

ZCORN
0 0 1 1  0 0 1 1
8 8 9 9  8 8 9 9
/

ACTNUM
2*1
/

PORO
2*0.15
/
)";

    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckString,
                              /* cells_per_dim_vec */ {{2,3,2}},
                              /* startIJK_vec */      {{0,0,0}},
                              /* endIJK_vec */        {{1,1,1}},
                              /* lgr_name_vec */      {"LGR1"});

    // LGR1 dimensions = {2,3,2}
    // LGR1 indices
    //
    // k = 1      |10    11|
    //            | 8     9|
    //            | 6     7|
    //            ----------
    // k = 0      | 4     5|
    //            | 2     3|
    //            | 0     1|
    //            ----------

    // Element 0 in level zero grid has two faces of type {I_FACE, true}
    //
    // Vertices of those faces lie on the plane x = 6    | After refinement, number of subdivisions in       LGR1 cell indices
    //                                                   | y- and z- directions:
    //              (6,0,8) ---------------- (6,6,8)     |  (6,0,8) --(6,2,8)-(6,4,8)--(6,6,8)               x-----x-----x-----x
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //                 |                      |          |     |         *       *        |                  |  7  *  9  *  11 |
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //                 |      face idx 2      |          |     |         *       *        |                  |     *     *     |
    //                 |                      |          |  (6,0,4) **(6,2,4)*(6,4,4)**(6,6,4)               x*****x*****x*****x
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //              (6,0,1) -----------------(6,6,1)     |  (6,0,1) --(?,?,?)-(?,?,?)--(6,6,1)               x- 1 -x- 3 -x- 5 -x
    //                 |      face idx 1      |          |     |         *       *        |                  |     *     *     |
    //              (6,0,0) -----------------(6,6,0)     |  (6,0,0) --(6,2,0)-(6,4,0)--(6,6,0)               x-----x-----x-----x
    //                                                   |
    //                                                   | The missing vertices are (6,2,1) and (6,4,1), appering in elements 1,3, or 5 in LGR1.
    //                                                   | In LGR1 element 1: (6,2,1)
    //                                                   | In LGR1 element 3: (6,2,1) and (6,4,1)
    //                                                   | In LGR1 element 5: (6,4,1)

    const auto& refinedGridData = *grid.currentData()[1];
    const auto& parentGridData = *grid.currentData()[0];
    const auto parentElem = Dune::cpgrid::Entity<0>(parentGridData, 0, true);

    std::map<Dune::FieldVector<double, 3>, int, Opm::Lgr::FieldVectorLess> vertexToIdx{};
    std::map<Dune::FieldVector<double, 3>, int, Opm::Lgr::FieldVectorLess> existingVtxInCoarseGridToItsIdx{};
    std::map<int, std::vector<std::vector<Dune::FieldVector<double,3>>>> allOverlapFaces{};
    std::vector<std::vector<std::pair<int,std::vector<Dune::FieldVector<double,3>>>>> vanishedRefinedFaceToNewRefinedFaces{};
    std::vector<std::vector<int>> refinedFaceIdx_to_parentFaceIdx{};

    vanishedRefinedFaceToNewRefinedFaces.resize(52);
    std::vector<bool> wrongRefinedCells(refinedGridData.size(0));


    std::vector<int> fullyContainedInParentFace{};
    fullyContainedInParentFace.resize(refinedGridData.numFaces());
                                   
    std::cout<< refinedGridData.size(3) << " vertices in LGR1 " << std::endl;
    BOOST_CHECK_EQUAL( refinedGridData.size(3), 42);
    // LGR1 dims 2x3x2 -> 3x4x3 vertices + 4 extra missing vertices  (6,0,1) --(6,2,1)-(6,4,1)--(6,6,1)
    
    std::cout<< refinedGridData.numFaces() << " faces in LGR1 " << std::endl;
    //BOOST_CHECK_EQUAL( refinedGridData.numFaces(), ); // 52-3+6 = 55 faces
      
    
    for (const auto& refinedElem : Dune::elements(grid.levelGridView(1))) {
     
        std::set<Coordinate,Opm::Lgr::FieldVectorLess> expectedVertices{};
        std::vector<Coordinate> expectedNewFaceInFace2{}; // {vertex '0', vertex '1', vertex '2', vertex '3'}
        std::vector<Coordinate> expectedNewFaceInFace1{};
        // Vertex order in I_FACE: 0->jk, 1-> (j+1)k, 2->(j+1)(k+1), 3->j(k+1)
        //
        //         j(k+1) <-'3' --------- '2'-> (j+1)(k+1)
        //                   |             |                        
        //                   |             |
        //             jk <-'0' --------- '1'-> (j+1)k

        if (refinedElem.index() ==  1){
            expectedVertices = {{6.,0.,1.}, {6., 2., 1.}};
            // this element has to have 7 faces: 1 I-,J-,J+,K-,K+, and 2 I+:
            //      (6,0,4) **(6,2,4)
            //         |         *        I_FACE, true with vertices (6,0,1),(6,2,1),(6,2,4),(6,0,4)
            //         |         *
            //      (6,0,1) --(6,2,1)
            //         |         *        I_FACE, true with vertices (6,0,0),(6,2,0),(6,2,1),(6,0,1)
            //      (6,0,0) --(6,2,0)
            expectedNewFaceInFace2 = {{6,0,1}, {6,2,1}, {6,2,4}, {6,0,4}}; 
            expectedNewFaceInFace1 = {{6,0,0}, {6,2,0}, {6,2,1}, {6,0,1}};

            // const auto& overlapFaces = allOverlapFaces.at( /* parent face index */ 2);

            // BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 2);

            // checkOverlapNewFace(overlapFaces, 2, expectedNewFaceInFace2);
            //checkOverlapNewFace(overlapFaces,  1, expectedNewFaceInFace1);
        }
        else if (refinedElem.index() ==  3){
            expectedVertices = {{6., 2., 1.}, {6., 4.,1.}};

            // this element has to have 7 faces: 1 I-,J-,J+,K-,K+, and 2 I+:
            //      (6,2,4) **(6,4,4)
            //         |         *        I_FACE, true with vertices (6,2,1),(6,4,1),(6,4,4),(6,2,4)
            //         |         *
            //      (6,2,1) --(6,4,1)
            //         |         *        I_FACE, true with vertices (6,2,0),(6,4,0),(6,4,1),(6,2,4)
            //      (6,2,0) --(6,4,0)
            expectedNewFaceInFace2 = {{6,2,1},{6,4,1},{6,4,4},{6,2,4}};
            expectedNewFaceInFace1 = {{6,2,0},{6,4,0},{6,4,1},{6,2,1}};

            //            BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 2);

            //checkOverlapNewFace(overlapFaces, 2, expectedNewFaceInFace2);
            // checkOverlapNewFace(overlapFaces, 1, expectedNewFaceInFace1);
        }
        else if (refinedElem.index() ==  5){
            expectedVertices = {{6., 4.,1.}, {6.,6.,1.}};

            // this element has to have 7 faces: 1 I-,J-,J+,K-,K+, and 2 I+:
            //      (6,4,4) **(6,6,4)
            //         |         *        I_FACE, true with vertices (6,4,1),(6,6,1),(6,6,4),(6,4,4)
            //         |         *
            //      (6,4,1) --(6,6,1)
            //         |         *        I_FACE, true with vertices (6,4,0),(6,6,0),(6,6,1),(6,4,1)
            //      (6,4,0) --(6,6,0)
            expectedNewFaceInFace2 = {{6,4,1},{6,6,1},{6,6,4},{6,4,4}};
            expectedNewFaceInFace1 = {{6,4,0},{6,6,0},{6,6,1},{6,4,1}};

            //            BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 2);

            //checkOverlapNewFace(overlapFaces, /* parent face index */ 2, expectedNewFaceInFace2);
            // checkOverlapNewFace(overlapFaces, /* parent face index */ 1, expectedNewFaceInFace1);
        }
        else if (refinedElem.index() ==  7) { 
            expectedNewFaceInFace2 = {{6,0,4}, {6,2,4}, {6,2,8}, {6,0,8}};
            //  BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 1);
            // checkOverlapNewFace(overlapFaces, /* parent face index */ 2, expectedNewFaceInFace2);
        }
        else if (refinedElem.index() ==  9) { 
            expectedNewFaceInFace2 = {{6,2,4}, {6,4,4}, {6,4,8}, {6,2,8}};
            //BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 1);
            //checkOverlapNewFace(overlapFaces, /* parent face index */ 2, expectedNewFaceInFace2);
        }
        else if (refinedElem.index() ==  11) { 
            expectedNewFaceInFace2 = {{6,4,4}, {6,6,4}, {6,6,8}, {6,4,8}};
            //BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 1);
            // checkOverlapNewFace(overlapFaces, /* parent face index */ 2, expectedNewFaceInFace2);
        }
        else {
            //BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 0);
        }
        // checkNewVertices(collectedVertices, expectedVertices); // expectedVertices is empty if refinedElem.index() != 1, 3, or 5
    }
}


BOOST_AUTO_TEST_CASE(parentCellWithMoreThanOne_I_FACE_false, *boost::unit_test::disabled())
{
    // Level zero grid dims = 2x1x1
    //
    // cell 0
    // bottom face corners (0,0,0), (6,0,0), (0,6,0), (6,6,0)
    //    top face corners (0,0,8), (6,0,8), (0,6,8), (6,6,8)
    //
    // cell 1
    // bottom face corners (6,0,1), (12,0,1), (6,6,1),  (12,6,1)
    //    top face corners (6,0,9), (12,0,9), (12,6,9), (12,6,9)
    const std::string deckString =
        R"(RUNSPEC
DIMENS
 2 1 1 /

GRID

COORD
 0 0 0     0 0 9
 6 0 0     6 0 9
12 0 0    12 0 9

 0 6 0    0 6 9
 6 6 0    6 6 9
12 6 0   12 6 9
/

ZCORN
0 0 1 1  0 0 1 1
8 8 9 9  8 8 9 9
/

ACTNUM
2*1
/

PORO
2*0.15
/
)";

    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckString,
                              /* cells_per_dim_vec */ {{2,3,2}},
                              /* startIJK_vec */      {{1,0,0}},
                              /* endIJK_vec */        {{2,1,1}},
                              /* lgr_name_vec */      {"LGR1"});

    // LGR1 dimensions = {2,3,2}
    // LGR1 indices
    //
    // k = 1      |10    11|
    //            | 8     9|
    //            | 6     7|
    //            ----------
    // k = 0      | 4     5|
    //            | 2     3|
    //            | 0     1|
    //            ----------

    // Element 1 in level zero grid has two faces of type {I_FACE, false}
    //
    // Vertices of those faces lie on the plane x = 6    | After refinement, number of subdivisions in       LGR1 cell indices
    //                                                   | y- and z- directions:
    //              (6,0,9) -----------------(6,6,9)     |  (6,0,9) --(6,2,9)-(6,4,9)--(6,6,9)               x*****x*****x*****x
    //                 |      face idx 3      |          |     |         *       *        |                  |     *     *     |
    //              (6,0,8) -----------------(6,6,8)     |  (6,0,8) --(?,?,?)-(?,?,?)--(6,6,8)               x- 6 -x- 8 -x-10--x
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //                 |      face idx 2      |          |  (6,0,5) **(6,2,5)*(6,4,5)**(6,6,5)               x*****x*****x*****x
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //                 |                      |          |     |         *       *        |                  |  0  *  2  *  4  |
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //              (6,0,1) ---------------- (6,6,1)     |  (6,0,1) --(6,2,1)-(6,4,1)--(6,6,1)               x-----x-----x-----x
    //                                                   |
    //                                                   | The missing vertices are (6,2,8) and (6,4,8), appering in elements 6,8, or 10 in LGR1.
    //                                                   | In LGR1 element  6: (6,2,8)
    //                                                   | In LGR1 element  8: (6,2,8) and (6,4,8)
    //                                                   | In LGR1 element 10: (6,4,8)

    const auto& refinedGridData = *grid.currentData()[1];
    const auto& parentGridData = *grid.currentData()[0];
    const auto parentElem = Dune::cpgrid::Entity<0>(parentGridData, 1, true);

    int numFaces = 52;

    std::map<Dune::FieldVector<double, 3>, int, Opm::Lgr::FieldVectorLess> vertexToIdx{};
    std::map<Dune::FieldVector<double, 3>, int, Opm::Lgr::FieldVectorLess> existingVtxInCoarseGridToItsIdx{};
    std::map<int, std::vector<std::vector<Dune::FieldVector<double,3>>>> allOverlapFaces{};
    std::vector<std::vector<std::pair<int,std::vector<Dune::FieldVector<double,3>>>>> vanishedRefinedFaceToNewRefinedFaces{};
    std::vector<std::vector<int>> refinedFaceIdx_to_parentFaceIdx{};

    std::vector<int> fullyContainedInParentFace{};
    fullyContainedInParentFace.resize(refinedGridData.numFaces());

    vanishedRefinedFaceToNewRefinedFaces.resize(52);
  

    for (const auto& refinedElem : Dune::elements(grid.levelGridView(1))) {

        const auto collectedVertices = Opm::Lgr::collectNewVertices<Coordinate>(refinedGridData,
                                                                                refinedElem,
                                                                                parentGridData,
                                                                                parentElem,
                                                                                vertexToIdx,
                                                                                existingVtxInCoarseGridToItsIdx,
                                                                                refinedFaceIdx_to_parentFaceIdx,
                                                                                allOverlapFaces,
                                                                                vanishedRefinedFaceToNewRefinedFaces,
                                                                                fullyContainedInParentFace);
    }

    for (const auto& refinedElem : Dune::elements(grid.levelGridView(1))) {

        std::set<Coordinate,Opm::Lgr::FieldVectorLess> expectedVertices{};
        std::vector<Coordinate> expectedNewFaceInFace3{};
        std::vector<Coordinate> expectedNewFaceInFace2{};
        // Vertex order in I_FACE: 0->jk, 1-> (j+1)k, 2->(j+1)(k+1), 3->j(k+1)
        //
        //         j(k+1) <-'3' --------- '2'-> (j+1)(k+1)
        //                   |             |                        
        //                   |             |
        //             jk <-'0' --------- '1'-> (j+1)k

        if (refinedElem.index() ==  6){
            expectedVertices = {{6.,0.,8.}, {6., 2., 8.}};

            // this element has to have 7 faces: 1 I+,J-,J+,K-,K+, and 2 I-:
            //      (6,0,9) --(6,2,9)
            //         |         *        I_FACE, false with vertices (6,0,8),(6,2,8),(6,2,9),(6,0,9)
            //      (6,0,8) **(6,2,8)
            //         |         *        I_FACE, false with vertices (6,0,5),(6,2,5),(6,2,8),(6,0,8)
            //         |         *
            //      (6,0,5) --(6,2,5)
            expectedNewFaceInFace3 = {{6,0,8},{6,2,8},{6,2,9},{6,0,9}};
            expectedNewFaceInFace2 = {{6,0,5},{6,2,5},{6,2,8},{6,0,8}};

            // BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 2);

            // checkOverlapNewFace(overlapFaces, /* parent face index */ 3, expectedNewFaceInFace3);
            //checkOverlapNewFace(overlapFaces, /* parent face index */ 2, expectedNewFaceInFace2);
        }
        else if (refinedElem.index() ==  8){
            expectedVertices = {{6., 2., 8.}, {6., 4.,8.}};

            // this element has to have 7 faces: 1 I+,J-,J+,K-,K+, and 2 I-:
            //      (6,2,9) --(6,4,9)
            //         |         *        I_FACE, false with vertices (6,2,8),(6,4,8),(6,4,9),(6,2,9)
            //      (6,2,8) **(6,4,8)
            //         |         *        I_FACE, false with vertices (6,2,5),(6,4,5),(6,4,8),(6,2,8)
            //         |         *
            //      (6,2,5) --(6,4,5)
            expectedNewFaceInFace3 = {{6,2,8},{6,4,8},{6,4,9},{6,2,9}};
            expectedNewFaceInFace2 = {{6,2,5},{6,4,5},{6,4,8},{6,2,8}};


            //            BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 2);

            // checkOverlapNewFace(overlapFaces, /* parent face index */ 3, expectedNewFaceInFace3);
            // checkOverlapNewFace(overlapFaces, /* parent face index */ 2, expectedNewFaceInFace2);
        }
        else if (refinedElem.index() ==  10){
            expectedVertices = {{6., 4.,8.}, {6.,6.,8.}};

            // this element has to have 7 faces: 1 I+,J-,J+,K-,K+, and 2 I-:
            //      (6,4,9) --(6,6,9)
            //         |         *        I_FACE, false with vertices (6,4,8),(6,6,8),(6,6,9),(6,4,9)
            //      (6,4,8) **(6,6,8)
            //         |         *        I_FACE, false with vertices (6,4,5),(6,6,5),(6,6,8),(6,4,8)
            //         |         *
            //      (6,4,5) --(6,6,5)
            expectedNewFaceInFace3 = {{6,4,8},{6,6,8},{6,6,9},{6,4,9}};
            expectedNewFaceInFace2 = {{6,4,5},{6,6,5},{6,6,8},{6,4,8}};

            // BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 2);

            // checkOverlapNewFace(overlapFaces, /* parent face index */ 3, expectedNewFaceInFace3);
            // checkOverlapNewFace(overlapFaces, /* parent face index */ 2, expectedNewFaceInFace2);
        }
        else if (refinedElem.index() ==  0) { 
            expectedNewFaceInFace2 = {{6,0,1}, {6,2,1}, {6,2,5}, {6,0,5}};
            // BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 1);
            // checkOverlapNewFace(overlapFaces, /* parent face index */ 2, expectedNewFaceInFace2);
        }
        else if (refinedElem.index() ==  2) { 
            expectedNewFaceInFace2 = {{6,2,1}, {6,4,1}, {6,4,5}, {6,2,5}};
            //  BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 1);
            // checkOverlapNewFace(overlapFaces, /* parent face index */ 2, expectedNewFaceInFace2);
        }
        else if (refinedElem.index() ==  4) { 
            expectedNewFaceInFace2 = {{6,4,1}, {6,6,1}, {6,6,5}, {6,4,5}};
            // BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 1);
            // checkOverlapNewFace(overlapFaces, /* parent face index */ 2, expectedNewFaceInFace2);
        }
        else {
            // BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 0);
        }
        // checkNewVertices(collectedVertices, expectedVertices); // expectedVertices empty if refinedElem.index() != 6,8, or 10
    }
}

BOOST_AUTO_TEST_CASE(parentCellWithMoreThanSixIntersections_J_FACE_true,*boost::unit_test::disabled())
{
    // Level zero grid dims = 1x2x1
    //
    // cell 0
    // bottom face corners (0,0,0), (6,0,0), (0,6,0), (6,6,0)
    //    top face corners (0,0,8), (6,0,8), (0,6,8), (6,6,8)
    //
    // cell 1
    // bottom face corners (0,6,1), (6,6,1), (0,12,1), (6,12,1)
    //    top face corners (0,6,9), (6,6,9), (0,12,9), (6,12,9)

    const std::string deckString =
        R"(RUNSPEC
DIMENS
 1 2 1 /

GRID

COORD
 0 0 0    0 0 9
 6 0 0    6 0 9

 0 6 0    0 6 9
 6 6 0    6 6 9

 0 12 0   0 12 9
 6 12 0   6 12 9
/

ZCORN
0 0 0 0  1 1 1 1
8 8 8 8  9 9 9 9
/

ACTNUM
2*1
/

PORO
2*0.15
/
)";

    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckString,
                              /* cells_per_dim_vec */ {{3,2,2}},
                              /* startIJK_vec */      {{0,0,0}},
                              /* endIJK_vec */        {{1,1,1}},
                              /* lgr_name_vec */      {"LGR1"});

    // LGR1 dimensions = {3,2,2}
    // LGR1 indices
    //
    // k = 1      | 9  10  11|
    //            | 6   7   8|
    //            ----------
    // k = 0      | 3   4   5|
    //            | 0   1   2|
    //            ------------


    // Element 0 in level zero grid has two faces of type {J_FACE, true}
    //
    // Vertices of those faces lie on the plane y = 6    | After refinement, number of subdivisions in       LGR1 cell indices
    //                                                   | y- and z- directions:
    //              (0,6,8) ---------------- (6,6,8)     |  (0,6,8) --(2,6,8)-(4,6,8)--(6,6,8)               x-----x-----x-----x
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //                 |                      |          |     |         *       *        |                  |  9  * 10  *  11 |
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //                 |       face idx 6     |          |     |         *       *        |                  |     *     *     |
    //                 |                      |          |  (0,6,4) **(2,6,4)*(4,6,4)**(6,6,4)               x*****x*****x*****x
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //              (0,6,1) -----------------(6,6,1)     |  (0,6,1) --(?,?,?)-(?,?,?)--(6,6,1)               x- 3 -x- 4 -x- 5 -x
    //                 |       face idx 5     |          |     |         *       *        |                  |     *     *     |
    //              (0,6,0) -----------------(6,6,0)     |  (0,6,0) --(2,6,0)-(4,6,0)--(6,6,0)               x-----x-----x-----x
    //                                                   |
    //                                                   | The missing vertices are (2,6,1) and (4,6,1), appering in elements 3,4 or 5 in LGR1.
    //                                                   | In LGR1 element 3: (2,6,1)
    //                                                   | In LGR1 element 4: (2,6,1) and (4,6,1)
    //                                                   | In LGR1 element 5: (4,6,1)

    const auto& refinedGridData = *grid.currentData()[1];
    const auto& parentGridData = *grid.currentData()[0];
    const auto parentElem = Dune::cpgrid::Entity<0>(parentGridData, 0, true);

    int numFaces = 52;
    
    std::map<Dune::FieldVector<double, 3>, int, Opm::Lgr::FieldVectorLess> vertexToIdx{};
    std::map<Dune::FieldVector<double, 3>, int, Opm::Lgr::FieldVectorLess> existingVtxInCoarseGridToItsIdx{};
    std::map<int, std::vector<std::vector<Dune::FieldVector<double,3>>>> allOverlapFaces{};
    std::vector<std::vector<std::pair<int,std::vector<Dune::FieldVector<double,3>>>>> vanishedRefinedFaceToNewRefinedFaces{};
    std::vector<std::vector<int>> refinedFaceIdx_to_parentFaceIdx{};

    std::vector<int> fullyContainedInParentFace{};
    fullyContainedInParentFace.resize(refinedGridData.numFaces());

    vanishedRefinedFaceToNewRefinedFaces.resize(52);

    for (const auto& refinedElem : Dune::elements(grid.levelGridView(1))) {

        const auto collectedVertices = Opm::Lgr::collectNewVertices<Coordinate>(refinedGridData,
                                                                                refinedElem,
                                                                                parentGridData,
                                                                                parentElem,
                                                                                vertexToIdx,
                                                                                existingVtxInCoarseGridToItsIdx,
                                                                                refinedFaceIdx_to_parentFaceIdx,
                                                                                allOverlapFaces,
                                                                                vanishedRefinedFaceToNewRefinedFaces,
                                                                                fullyContainedInParentFace);
    }

    for (const auto& refinedElem : Dune::elements(grid.levelGridView(1))) {

        std::set<Coordinate,Opm::Lgr::FieldVectorLess> expectedVertices{};
        std::vector<Coordinate> expectedNewFaceInFace6{};
        std::vector<Coordinate> expectedNewFaceInFace5{};
        // Vertex order in J_FACE: 0->(i+1)k, 1-> ik, 2->i(k+1), 3->(i+1)(k+1)
        //
        //         i(k+1) <-'2' --------- '3'-> (i+1)(k+1)
        //                   |             |                        
        //                   |             |
        //             ik <-'1' --------- '0'-> (i+1)k

        if (refinedElem.index() ==  3){
            expectedVertices = {{0.,6.,1.}, {2., 6., 1.}};

            // this element has to have 7 faces: 1 I-,I+,J-,K-,K+, and 2 J+:
            //      (0,6,4) **(2,6,4)
            //         |         *        J_FACE, true with vertices (0,6,1),(2,6,1),(2,6,4),(,6,4)
            //         |         *
            //      (0,6,1) --(2,6,1)
            //         |         *        J_FACE, true with vertices (0,6,0),(2,6,0),(2,6,1),(0,6,1)
            //      (0,6,0) --(2,6,0)
            expectedNewFaceInFace6 = {{2,6,1}, {0,6,1}, {0,6,4}, {2,6,4}}; 
            expectedNewFaceInFace5 = {{2,6,0}, {0,6,0}, {0,6,1}, {2,6,1}}; 

            //   BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 2);

            // checkOverlapNewFace(overlapFaces, /* parent face index */ 6, expectedNewFaceInFace6);
            /// checkOverlapNewFace(overlapFaces, /* parent face index */ 5, expectedNewFaceInFace5);
        }
        else if (refinedElem.index() == 4){
            expectedVertices = {{2., 6., 1.}, {4., 6., 1.}};

            // this element has to have 7 faces: 1 I-,I+,J-,K-,K+, and 2 J+:
            //      (2,6,4) **(4,6,4)
            //         |         *        J_FACE, true with vertices (2,6,1),(4,6,1),(4,6,4),(2,6,4)
            //         |         *
            //      (2,6,1) --(4,6,1)
            //         |         *        J_FACE, true with vertices (2,6,0),(4,6,0),(4,6,1),(2,6,1)
            //      (2,6,0) --(4,6,0)
            expectedNewFaceInFace6 = {{4,6,1}, {2,6,1}, {2,6,4}, {4,6,4}}; 
            expectedNewFaceInFace5 = {{4,6,0}, {2,6,0}, {2,6,1}, {4,6,1}};

            // BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 2);

            // checkOverlapNewFace(overlapFaces, /* parent face index */ 6, expectedNewFaceInFace6);
            // checkOverlapNewFace(overlapFaces, /* parent face index */ 5, expectedNewFaceInFace5);
        }
        else if (refinedElem.index() == 5){
            expectedVertices = {{4., 6.,1.}, {6.,6.,1.}};

            // this element has to have 7 faces: 1 I-,I+,J-,K-,K+, and 2 J+:
            //      (4,6,4) **(6,6,4)
            //         |         *        J_FACE, true with vertices (4,6,1),(6,6,1),(6,6,4),(4,6,4)
            //         |         *
            //      (4,6,1) --(6,6,1)
            //         |         *        J_FACE, true with vertices (4,6,0),(6,6,0),(6,6,1),(4,6,1)
            //      (4,6,0) --(6,6,0)
            expectedNewFaceInFace6 = {{6,6,1}, {4,6,1}, {4,6,4}, {6,6,4}}; 
            expectedNewFaceInFace5 = {{6,6,0}, {4,6,0}, {4,6,1}, {6,6,1}};

            // BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 2);

            // checkOverlapNewFace(overlapFaces, /* parent face index */ 6, expectedNewFaceInFace6);
            // checkOverlapNewFace(overlapFaces, /* parent face index */ 5, expectedNewFaceInFace5);
        }
        else if (refinedElem.index() ==  9) { 
            expectedNewFaceInFace6 = {{2,6,4}, {0,6,4}, {0,6,8}, {2,6,8}};
            // BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 1);
            // checkOverlapNewFace(overlapFaces, /* parent face index */ 6, expectedNewFaceInFace6);
        }
        else if (refinedElem.index() ==  10) { 
            expectedNewFaceInFace6 = {{4,6,4}, {2,6,4}, {2,6,8}, {4,6,8}};
            // BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 1);
            // checkOverlapNewFace(overlapFaces, /* parent face index */ 6, expectedNewFaceInFace6);
        }
        else if (refinedElem.index() ==  11) { 
            expectedNewFaceInFace6 = {{6,6,4}, {4,6,4}, {4,6,8}, {6,6,8}};
            // BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 1);
            // checkOverlapNewFace(overlapFaces, /* parent face index */ 6, expectedNewFaceInFace6);
        }
        else {
            //BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 0);
        }
        // checkNewVertices(collectedVertices, expectedVertices);  // expectedVertices is empty if refinedElem.index() != 3,4, or 5
    }
}


BOOST_AUTO_TEST_CASE(parentCellWithMoreThanSixIntersections_J_FACE_false, *boost::unit_test::disabled())
{
    // Level zero grid dims = 1x2x1
    //
    // cell 0
    // bottom face corners (0,0,0), (6,0,0), (0,6,0), (6,6,0)
    //    top face corners (0,0,8), (6,0,8), (0,6,8), (6,6,8)
    //
    // cell 1
    // bottom face corners (0,6,1), (6,6,1), (0,12,1), (6,12,1)
    //    top face corners (0,6,9), (6,6,9), (0,12,9), (6,12,9)

    const std::string deckString =
        R"(RUNSPEC
DIMENS
 1 2 1 /

GRID

COORD
 0 0 0    0 0 9
 6 0 0    6 0 9

 0 6 0    0 6 9
 6 6 0    6 6 9

 0 12 0   0 12 9
 6 12 0   6 12 9
/

ZCORN
0 0 0 0  1 1 1 1
8 8 8 8  9 9 9 9
/

ACTNUM
2*1
/

PORO
2*0.15
/
)";

    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckString,
                              /* cells_per_dim_vec */ {{3,2,2}},
                              /* startIJK_vec */      {{0,1,0}},
                              /* endIJK_vec */        {{1,2,1}},
                              /* lgr_name_vec */      {"LGR1"});

    // LGR1 dimensions = {3,2,2}
    // LGR1 indices
    //
    // k = 1      | 9  10  11|
    //            | 6   7   8|
    //            ----------
    // k = 0      | 3   4   5|
    //            | 0   1   2|
    //            ------------

    // Element 1 in level zero grid has two faces of type {J_FACE, false}
    //
    // Vertices of those faces lie on the plane x = 6    | After refinement, number of subdivisions in       LGR1 cell indices
    //                                                   | y- and z- directions:
    //              (0,6,9) -----------------(6,6,9)     |  (0,6,9) --(2,6,9)-(4,6,9)--(6,6,9)               x*****x*****x*****x
    //                 |      face idx 7      |          |     |         *       *        |                  |     *     *     |
    //              (0,6,8) -----------------(6,6,8)     |  (0,6,8) --(?,?,?)-(?,?,?)--(6,6,8)               x- 6 -x- 7 -x- 8--x
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //                 |      face idx 6      |          |  (0,6,5) **(2,6,5)*(4,6,5)**(6,6,5)               x*****x*****x*****x
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //                 |                      |          |     |         *       *        |                  |  0  *  1  *  2  |
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //              (0,6,1) ---------------- (6,6,1)     |  (0,6,1) --(2,6,1)-(4,6,1)--(6,6,1)               x-----x-----x-----x
    //                                                   |
    //                                                   | The missing vertices are (2,6,8) and (4,6,8), appering in elements 6,7, or 8 in LGR1.
    //                                                   | In LGR1 element  6: (2,6,8)
    //                                                   | In LGR1 element  7: (2,6,8) and (4,6,8)
    //                                                   | In LGR1 element  8: (4,6,8)

    const auto& refinedGridData = *grid.currentData()[1];
    const auto& parentGridData = *grid.currentData()[0];
    const auto parentElem = Dune::cpgrid::Entity<0>(parentGridData, 1, true);

    int numFaces = 52;
    
    std::map<Dune::FieldVector<double, 3>, int, Opm::Lgr::FieldVectorLess> vertexToIdx{};
    std::map<Dune::FieldVector<double, 3>, int, Opm::Lgr::FieldVectorLess> existingVtxInCoarseGridToItsIdx{};
    std::map<int, std::vector<std::vector<Dune::FieldVector<double,3>>>> allOverlapFaces{};
    std::vector<std::vector<std::pair<int,std::vector<Dune::FieldVector<double,3>>>>> vanishedRefinedFaceToNewRefinedFaces{};
    std::vector<std::vector<int>> refinedFaceIdx_to_parentFaceIdx{};

    std::vector<int> fullyContainedInParentFace{};
    fullyContainedInParentFace.resize(refinedGridData.numFaces());


    vanishedRefinedFaceToNewRefinedFaces.resize(52);
    std::vector<bool> wrongRefinedCells(refinedGridData.size(0));

    for (const auto& refinedElem : Dune::elements(grid.levelGridView(1))) {

        const auto collectedVertices = Opm::Lgr::collectNewVertices<Coordinate>(refinedGridData,
                                                                                refinedElem,
                                                                                parentGridData,
                                                                                parentElem,
                                                                                vertexToIdx,
                                                                                existingVtxInCoarseGridToItsIdx,
                                                                                refinedFaceIdx_to_parentFaceIdx,
                                                                                allOverlapFaces,
                                                                                vanishedRefinedFaceToNewRefinedFaces,
                                                                                fullyContainedInParentFace);
    }

    for (const auto& refinedElem : Dune::elements(grid.levelGridView(1))) {

        std::set<Coordinate,Opm::Lgr::FieldVectorLess> expectedVertices{};
        std::vector<Coordinate> expectedNewFaceInFace7{};
        std::vector<Coordinate> expectedNewFaceInFace6{};
        // Vertex order in J_FACE: 0->(i+1)k, 1-> ik, 2->i(k+1), 3->(i+1)(k+1)
        //
        //         i(k+1) <-'2' --------- '3'-> (i+1)(k+1)
        //                   |             |                        
        //                   |             |
        //             ik <-'1' --------- '0'-> (i+1)k

        if (refinedElem.index() ==  6){
            expectedVertices = {{0.,6.,8.}, {2., 6., 8.}};

            // this element has to have 7 faces: 1 I-,I+,J+,K-,K+, and 2 J-:
            //      (0,6,9) --(2,6,9)
            //         |         *        J_FACE, false with vertices (0,6,8),(2,6,8),(2,6,9),(0,6,9)
            //      (0,6,8) **(2,6,8)
            //         |         *        J_FACE, false with vertices (0,6,5),(2,6,5),(2,6,8),(0,6,8)
            //         |         *
            //      (0,6,5) --(2,6,5)
            expectedNewFaceInFace7 = {{2,6,8}, {0,6,8}, {0,6,9}, {2,6,9}}; 
            expectedNewFaceInFace6 = {{2,6,5}, {0,6,5}, {0,6,8}, {2,6,8}};

            //BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 2);

            // checkOverlapNewFace(overlapFaces, /* parent face index */ 7, expectedNewFaceInFace7);
            // checkOverlapNewFace(overlapFaces, /* parent face index */ 6, expectedNewFaceInFace6);
        }
        else if (refinedElem.index() == 7){
            expectedVertices = {{2., 6., 8.}, {4., 6., 8.}};

            // this element has to have 7 faces: 1 I-,I+,J+,K-,K+, and 2 J-:
            //      (2,6,9) --(4,6,9)
            //         |         *        J_FACE, false with vertices (2,6,8),(4,6,8),(4,6,9),(2,6,9)
            //      (2,6,8) **(4,6,8)
            //         |         *        J_FACE, false with vertices (2,6,5),(4,6,5),(4,6,8),(2,6,8)
            //         |         *
            //      (2,6,5) --(4,6,5)
            expectedNewFaceInFace7 = {{4,6,8}, {2,6,8}, {2,6,9}, {4,6,9}}; 
            expectedNewFaceInFace6 = {{4,6,5}, {2,6,5}, {2,6,8}, {4,6,8}};

            // BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 2);

            // checkOverlapNewFace(overlapFaces, /* parent face index */ 7, expectedNewFaceInFace7);
           // checkOverlapNewFace(overlapFaces, /* parent face index */ 6, expectedNewFaceInFace6);
        }
        else if (refinedElem.index() == 8){
            expectedVertices = {{4., 6., 8.}, {6.,6.,8.}};

            // this element has to have 7 faces: 1 I-,I+,J+,K-,K+, and 2 J-:
            //      (4,6,9) --(6,6,9)
            //         |         *        J_FACE, false with vertices (4,6,8),(6,6,8),(6,6,9),(4,6,9)
            //      (4,6,8) **(6,6,8)
            //         |         *        J_FACE, false with vertices (4,6,5),(6,6,5),(6,6,8),(4,6,8)
            //         |         *
            //      (4,6,5) --(6,6,5)
            expectedNewFaceInFace7 = {{6,6,8}, {4,6,8}, {4,6,9}, {6,6,9}};
            expectedNewFaceInFace6 = {{6,6,5}, {4,6,5}, {4,6,8}, {6,6,8}};

            //BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 2);

            // checkOverlapNewFace(overlapFaces, /* parent face index */ 7, expectedNewFaceInFace7);
            // checkOverlapNewFace(overlapFaces, /* parent face index */ 6, expectedNewFaceInFace6);
        }
        else if (refinedElem.index() ==  0) { 
            expectedNewFaceInFace6 = {{2,6,1}, {0,6,1}, {0,6,5}, {2,6,5}};
            //   BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 1);
            //  checkOverlapNewFace(overlapFaces, /* parent face index */ 6, expectedNewFaceInFace6);
        }
         else if (refinedElem.index() ==  1) { 
             expectedNewFaceInFace6 = {{4,6,1}, {2,6,1}, {2,6,5}, {4,6,5}};
             // BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 1);
             // checkOverlapNewFace(overlapFaces, /* parent face index */ 6, expectedNewFaceInFace6);
        }
         else if (refinedElem.index() ==  2) {
             expectedNewFaceInFace6 = {{6,6,1},{4,6,1}, {4,6,5}, {6,6,5}};
             // BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 1);
             // checkOverlapNewFace(overlapFaces, /* parent face index */ 6, expectedNewFaceInFace6);
        }
        else {
            // BOOST_CHECK_EQUAL( overlapFaces.size(), /*expectedNewFacesSize */ 0);
        }

        // checkNewVertices(collectedVertices, expectedVertices);  // expectedVertices is empty if refinedElem.index() != 3,4, or 5
    }
}
