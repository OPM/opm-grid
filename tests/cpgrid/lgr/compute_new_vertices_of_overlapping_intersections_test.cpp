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
#include <opm/grid/cpgrid/LgrFaultHelpers.hpp>
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


void checkFaceCountPerType(int repeatedFaceType,
                           int expectedRepeatedFaceTypeCount,
                           const std::array<std::vector<int>,6>& classifiedFaces)
{
    for (int faceType = 0; faceType < 6; ++faceType) {
        int expectedFaceTypeCount = (faceType == repeatedFaceType) ? expectedRepeatedFaceTypeCount : 1;
        BOOST_CHECK_EQUAL( classifiedFaces[faceType].size(), expectedFaceTypeCount);
    }
}

// The test cases consist in a grid with only 2 cells in level zero.
// Only one of them gets refined, the other appears on the leaf with
// same volume and center, but more faces.
void checkFaceCountInLeafCoarseElem(const Dune::CpGrid& grid,
                                    int expectedTotalFaceCount,
                                    int repeatedFaceType,
                                    int expectedRepeatedFaceTypeCount)
{
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        if (element.level() == 0) {

            const auto& cellToFace = grid.currentLeafData().cellToFace(element.index());
            BOOST_CHECK_EQUAL( cellToFace.size(), expectedTotalFaceCount);

            const auto classifiedFaces =  Opm::Lgr::groupFaceIndicesByType(grid.currentLeafData(),
                                                                           element);

            checkFaceCountPerType(repeatedFaceType,
                                  expectedRepeatedFaceTypeCount,
                                  classifiedFaces);
        }
    }
}

void checkFaceToCoord(const Dune::cpgrid::CpGridData& refinedGridData,
                      const std::vector<int>& selectedFaces,
                      const std::set<Coordinate,Opm::Lgr::FieldVectorLess>& expectedFaceToCoord1,
                      const std::set<Coordinate,Opm::Lgr::FieldVectorLess>& expectedFaceToCoord2)
{

    for (const auto& faceIdx : selectedFaces) {

        const auto& faceToPoint = refinedGridData.faceToPoint(faceIdx);
        std::set<Coordinate,Opm::Lgr::FieldVectorLess> faceToCoord{};
        for (const auto& point : faceToPoint) {
            const auto pointEntity = Dune::cpgrid::Entity<3>(refinedGridData, point, true);
            faceToCoord.insert(pointEntity.geometry().center());
        }
        BOOST_CHECK( (faceToCoord == expectedFaceToCoord1) || (faceToCoord == expectedFaceToCoord2) );
    }
}

void checkNewRefinedFaces(const Dune::CpGrid& grid,
                          const Dune::cpgrid::CpGridData& refinedGridData,
                          const std::vector<std::vector<std::set<Coordinate,Opm::Lgr::FieldVectorLess>>>& selectedFaceToCoord,
                          int repeatedFaceType)
{

    for (const auto& refinedElem : Dune::elements(grid.levelGridView(1))) {

        const auto& cellToFace = refinedGridData.cellToFace(refinedElem.index());
        const auto classifiedFaces = Opm::Lgr::groupFaceIndicesByType(refinedGridData, refinedElem);

        if (selectedFaceToCoord[refinedElem.index()].size() == 2) {

            BOOST_CHECK_EQUAL( cellToFace.size(), 7);

            checkFaceCountPerType(repeatedFaceType,
                                  /* expectedRepeatedFaceTypeCount = */ 2,
                                  classifiedFaces);

            checkFaceToCoord(refinedGridData,
                             classifiedFaces[repeatedFaceType],
                             selectedFaceToCoord[refinedElem.index()][0],
                             selectedFaceToCoord[refinedElem.index()][1]);

        }
        else if (selectedFaceToCoord[refinedElem.index()].size() == 1){

            BOOST_CHECK_EQUAL( cellToFace.size(), 6);

            checkFaceCountPerType(/* repeatedFaceType = */ -1, // invalid face type
                                  /* expectedRepeatedFaceTypeCount = */ 0,
                                  classifiedFaces);

            checkFaceToCoord(refinedGridData,
                             classifiedFaces[repeatedFaceType],
                             selectedFaceToCoord[refinedElem.index()][0],
                             selectedFaceToCoord[refinedElem.index()][0]); // same, here there is no repeaetd face
        }
        else {
            BOOST_CHECK(selectedFaceToCoord[refinedElem.index()].empty());
            BOOST_CHECK_EQUAL( cellToFace.size(), 6);

            checkFaceCountPerType(/* repeatedFaceType = */ -1, // invalid face type
                                  /* expectedRepeatedFaceTypeCount = */ 0,
                                  classifiedFaces);
        }
    }
}

void expectedCellCountRespectToCellToFaceSize(const Dune::CpGrid& grid,
                                              const std::vector<int>& expectedCount7Faces)
{
    for (int level = 0; level <= grid.maxLevel(); ++level) {
        int count_7faces_level = 0;
        int count_6faces_level = 0;
        for (const auto& element : Dune::elements(grid.levelGridView(level))) {
            int count = grid.currentData()[level]->cellToFace(element.index()).size();
            if (count == 7) {
                ++count_7faces_level;
            }
            else {
                BOOST_CHECK_EQUAL(count, 6);
                ++count_6faces_level;
            }
        }
        BOOST_CHECK_EQUAL(count_7faces_level, expectedCount7Faces[level]);
        BOOST_CHECK_EQUAL(count_6faces_level, grid.levelGridView(level).size(0) - count_7faces_level);
    }
}

void checkLevelVertices(const Dune::CpGrid& grid,
                        const std::vector<std::set<Coordinate,Opm::Lgr::FieldVectorLess>>& expectedLgrVertices)
{
    for (int level = 1; level <= grid.maxLevel(); ++level) {
        std::set<Coordinate,Opm::Lgr::FieldVectorLess> levelVertices{};
        for (const auto& vertex : Dune::vertices(grid.levelGridView(level))) {
            levelVertices.insert(vertex.geometry().center());
        }

        // Note: FieldVectorLess determines where each element lives in the ordered tree (from the set).
        //       The insertion order is discarded once the element is placed.
        //       Both sets use the same FieldVectorLess comparator.
        //       Check they have the same size
        BOOST_REQUIRE_EQUAL( expectedLgrVertices[level-1].size(), levelVertices.size());

        auto expected =  expectedLgrVertices[level-1].begin();
        auto actual = levelVertices.begin();
        //       Check every element is the same
        for (; expected != expectedLgrVertices[level-1].end(); ++expected, ++actual) {
            BOOST_CHECK_EQUAL_COLLECTIONS(expected->begin(), expected->end(),
                                          actual->begin(), actual->end());
        }
    }
}


void checkLeafVertices(const Dune::CpGrid& grid,
                       const std::set<Coordinate,Opm::Lgr::FieldVectorLess>& expectedVertices)
{
    std::set<Coordinate,Opm::Lgr::FieldVectorLess> leafVertices{};
    for (const auto& vertex : Dune::vertices(grid.leafGridView())) {
        leafVertices.insert(vertex.geometry().center());
    }
    BOOST_REQUIRE_EQUAL( expectedVertices.size(), leafVertices.size());

    auto expected =  expectedVertices.begin();
    auto actual = leafVertices.begin();

    for (; expected != expectedVertices.end(); ++expected, ++actual) {
        BOOST_CHECK_EQUAL_COLLECTIONS(expected->begin(), expected->end(),
                                      actual->begin(), actual->end());
    }
}

bool equalFaces(const std::vector<std::vector<Dune::FieldVector<double,3>>>& actual,
                const std::vector<std::vector<Dune::FieldVector<double,3>>>& expected)
{
    if (actual.size() != expected.size())
        return false;

    for (std::size_t i = 0; i < actual.size(); ++i) {
        if (actual[i].size() != expected[i].size())
            return false;

        auto aIt = actual[i].begin();
        auto eIt = expected[i].begin();

        while (aIt != actual[i].end()) {
            for (int d = 0; d < 3; ++d) // coordinate by coordinate
                if ((*aIt)[d] != (*eIt)[d])
                    return false;

            ++aIt;
            ++eIt;
        }
    }

    return true;
}

void checkFaces(const Dune::cpgrid::CpGridData& gridData,
                const std::vector<std::vector<Coordinate>>& expectedFaces)
{
    std::vector<std::vector<Dune::FieldVector<double,3>>> actualFaces{};
    actualFaces.resize(gridData.numFaces());

    for (int i = 0; i < gridData.numFaces(); ++i) {
        for (const auto& vertexIdx : gridData.faceToPoint(i)){
            actualFaces[i].push_back(Dune::cpgrid::Entity<3>(gridData, vertexIdx, true).geometry().center());
        }
    }

    BOOST_CHECK(equalFaces(actualFaces, expectedFaces));
}


BOOST_AUTO_TEST_CASE(parentCellWithMoreThanOne_I_FACE_trueOriented_nonTrivialOverlap)
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

    BOOST_CHECK_EQUAL(parentGridData.cellToFace(parentElem.index()).size(), 7);

    BOOST_CHECK_EQUAL( refinedGridData.size(3), 40);
    // LGR1 dims 2x3x2 -> 3x4x3 vertices + 4 extra missing vertices  (6,0,1), (6,2,1), (6,4,1), and (6,6,1).
    BOOST_CHECK_EQUAL( refinedGridData.numFaces(), 55);
    // LGR1 dims 2x3x2 -> 52 faces (before correction due to missing points)
    // 3 of those 52 faces vanished and give origin to 6 new faces: 52 - 3 + 6 = 55 faces

    const auto& leafGridData = grid.currentLeafData();
    BOOST_CHECK_EQUAL( leafGridData.size(3), 46); // 46 = 40 + 6 (40 in level 1 + 6 vertices from cell_to_point_ from coarse element)
    BOOST_CHECK_EQUAL( leafGridData.numFaces(), 61); // 61 = 55 + 6 (55 in level 1 + 6 other faces from coarse element)

    // Originally, the element not involved in refinement
    // had  7 faces. It's neihgbor in level zero
    // got refined and the I_FACE that they shared has been
    // replaced by 6 refined faces. Then, the leaf element has
    // one of each I+,J-,J+, K-, K+, and 1 coarse + 6 refined I-.
    checkFaceCountInLeafCoarseElem(grid,
                                   /* expectedTotalFaceCount = */ 12,
                                   /* repeatedFaceType = */ 0, // 0->I-
                                   /* expectedRepeatedFaceTypeCount = */ 7);

    expectedCellCountRespectToCellToFaceSize(grid,
                                             /* expectedCount7Faces = */ std::vector<int>{2, 3, 3}); // level 0, level 1, leaf grid

    // Collect the expected data to later on check
    std::vector<std::vector<std::set<Coordinate,Opm::Lgr::FieldVectorLess>>> selectedFaceToCoord{};
    selectedFaceToCoord.resize(grid.levelGridView(1).size(0));

    for (const auto& refinedElem : Dune::elements(grid.levelGridView(1))) {

        std::set<Coordinate,Opm::Lgr::FieldVectorLess> expectedNewFaceInFace2{}; // {vertex '0', vertex '1', vertex '2', vertex '3'}
        std::set<Coordinate,Opm::Lgr::FieldVectorLess> expectedNewFaceInFace1{};
        // Vertex order in I_FACE: 0->jk, 1-> (j+1)k, 2->(j+1)(k+1), 3->j(k+1)
        //
        //         j(k+1) <-'3' --------- '2'-> (j+1)(k+1)
        //                   |             |
        //                   |             |
        //             jk <-'0' --------- '1'-> (j+1)k

        if (refinedElem.index() ==  1){
            // this element has to have 7 faces: 1 I-,J-,J+,K-,K+, and 2 I+:
            //      (6,0,4) **(6,2,4)
            //         |         *        I_FACE, true with vertices (6,0,1),(6,2,1),(6,2,4),(6,0,4)
            //         |         *
            //      (6,0,1) --(6,2,1)
            //         |         *        I_FACE, true with vertices (6,0,0),(6,2,0),(6,2,1),(6,0,1)
            //      (6,0,0) --(6,2,0)
            expectedNewFaceInFace2 = {{6,0,1}, {6,2,1}, {6,2,4}, {6,0,4}};
            expectedNewFaceInFace1 = {{6,0,0}, {6,2,0}, {6,2,1}, {6,0,1}};

            selectedFaceToCoord[1].push_back(expectedNewFaceInFace2);
            selectedFaceToCoord[1].push_back(expectedNewFaceInFace1);
        }
        else if (refinedElem.index() ==  3){
            // this element has to have 7 faces: 1 I-,J-,J+,K-,K+, and 2 I+:
            //      (6,2,4) **(6,4,4)
            //         |         *        I_FACE, true with vertices (6,2,1),(6,4,1),(6,4,4),(6,2,4)
            //         |         *
            //      (6,2,1) --(6,4,1)
            //         |         *        I_FACE, true with vertices (6,2,0),(6,4,0),(6,4,1),(6,2,4)
            //      (6,2,0) --(6,4,0)
            expectedNewFaceInFace2 = {{6,2,1},{6,4,1},{6,4,4},{6,2,4}};
            expectedNewFaceInFace1 = {{6,2,0},{6,4,0},{6,4,1},{6,2,1}};

            selectedFaceToCoord[3].push_back(expectedNewFaceInFace2);
            selectedFaceToCoord[3].push_back(expectedNewFaceInFace1);
        }
        else if (refinedElem.index() ==  5){
            // this element has to have 7 faces: 1 I-,J-,J+,K-,K+, and 2 I+:
            //      (6,4,4) **(6,6,4)
            //         |         *        I_FACE, true with vertices (6,4,1),(6,6,1),(6,6,4),(6,4,4)
            //         |         *
            //      (6,4,1) --(6,6,1)
            //         |         *        I_FACE, true with vertices (6,4,0),(6,6,0),(6,6,1),(6,4,1)
            //      (6,4,0) --(6,6,0)
            expectedNewFaceInFace2 = {{6,4,1},{6,6,1},{6,6,4},{6,4,4}};
            expectedNewFaceInFace1 = {{6,4,0},{6,6,0},{6,6,1},{6,4,1}};

            selectedFaceToCoord[5].push_back(expectedNewFaceInFace2);
            selectedFaceToCoord[5].push_back(expectedNewFaceInFace1);
        }
        else if (refinedElem.index() ==  7) {
            expectedNewFaceInFace2 = {{6,0,4}, {6,2,4}, {6,2,8}, {6,0,8}};
            selectedFaceToCoord[7].push_back(expectedNewFaceInFace2);
        }
        else if (refinedElem.index() ==  9) {
            expectedNewFaceInFace2 = {{6,2,4}, {6,4,4}, {6,4,8}, {6,2,8}};
            selectedFaceToCoord[9].push_back(expectedNewFaceInFace2);
        }
        else if (refinedElem.index() ==  11) {
            expectedNewFaceInFace2 = {{6,4,4}, {6,6,4}, {6,6,8}, {6,4,8}};
            selectedFaceToCoord[11].push_back(expectedNewFaceInFace2);
        }
    }
    checkNewRefinedFaces(grid, refinedGridData,
                         selectedFaceToCoord, /* repeatedFaceType = */ 1); // 1-> I+

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{2,3,2}},
                           /* lgr_name_vec = */ {"LGR1"});
}

BOOST_AUTO_TEST_CASE(parentCellWithMoreThanOne_I_FACE_trueOriented_trivialOverlap)
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
                              /* cells_per_dim_vec */ {{2,3,8}},
                              /* startIJK_vec */      {{0,0,0}},
                              /* endIJK_vec */        {{1,1,1}},
                              /* lgr_name_vec */      {"LGR1"});



    // Element 0 in level zero grid has two faces of type {I_FACE, true}
    //
    // Vertices of those faces lie on the plane x = 6    | After refinement, number of subdivisions in
    //                                                   | y- and z- directions:
    //              (6,0,8) ---------------- (6,6,8)     |  (6,0,8) --(6,2,8)-(6,4,8)--(6,6,8)
    //                 |                      |          |     |         *       *        |
    //                 |                      |          |  (6,0,7) **(6,2,7)*(6,4,7)**(6,6,7)
    //                 |                      |          |     |         *       *        |
    //                 |                      |          |  (6,0,6) **(6,2,6)*(6,4,6)**(6,6,6)
    //                 |                      |          |     |         *       *        |
    //                 |      face idx 2      |          |  (6,0,5) **(6,2,5)*(6,4,5)**(6,6,5)
    //                 |                      |          |     |         *       *        |
    //                 |                      |          |  (6,0,4) **(6,2,4)*(6,4,4)**(6,6,4)
    //                 |                      |          |     |         *       *        |
    //                 |                      |          |  (6,0,3) **(6,2,3)*(6,4,3)**(6,6,3)
    //                 |                      |          |     |         *       *        |
    //                 |                      |          |  (6,0,2) **(6,2,2)*(6,4,2)**(6,6,2)
    //                 |                      |          |     |         *       *        |
    //              (6,0,1) -----------------(6,6,1)     |  (6,0,1) --(6,2,1)-(6,4,1)--(6,6,1)
    //                 |      face idx 1      |          |     |         *       *        |
    //              (6,0,0) -----------------(6,6,0)     |  (6,0,0) --(6,2,0)-(6,4,0)--(6,6,0)
    //                                                   |

    const auto& refinedGridData = *grid.currentData()[1];
    const auto& parentGridData = *grid.currentData()[0];
    const auto parentElem = Dune::cpgrid::Entity<0>(parentGridData, 0, true);

    BOOST_CHECK_EQUAL(parentGridData.cellToFace(parentElem.index()).size(), 7);

    BOOST_CHECK_EQUAL( refinedGridData.size(3), 108);
    // LGR1 dims 2x3x8 -> 3x4x9 (= 108) vertices (4 "missing" vertices  (6,0,1), (6,2,1), (6,4,1), and (6,6,1) exist).
    BOOST_CHECK_EQUAL( refinedGridData.numFaces(), 190);
    // LGR1 dims 2x3x8 -> (72+64+54=) 190 faces (before and after correction).


    // Originally, the element not involved in refinement
    // had  7 faces. It's neihgbor in level zero
    // got refined and the I_FACE that they shared has been
    // replaced by 6 refined faces. Then, the leaf element has
    // one of each I+,J-,J+, K-, K+, and 1 coarse + 21 refined I-.
    checkFaceCountInLeafCoarseElem(grid,
                                   /* expectedTotalFaceCount = */ 27,
                                   /* repeatedFaceType = */ 0, // 0->I-
                                   /* expoectedRepeatedFaceTypeCount = */ 22);

    expectedCellCountRespectToCellToFaceSize(grid,
                                             /* expectedCount7Faces = */ std::vector<int>{2, 0}); // level 0, level 1

    const auto& leafGridData = grid.currentLeafData();

    BOOST_CHECK_EQUAL( leafGridData.size(3), 114); // 114 = 108 + 6 (108 in level 1 + 6 vertices from cell_to_point_ from coarse element)
    BOOST_CHECK_EQUAL( leafGridData.numFaces(), 196); // 196 = 190 + 6 (190 in level 1 + 6 other faces from coarse element)

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{2,3,8}},
                           /* lgr_name_vec = */ {"LGR1"});
}


BOOST_AUTO_TEST_CASE(parentCellWithMoreThanOne_I_FACE_false)
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

    BOOST_CHECK_EQUAL(parentGridData.cellToFace(parentElem.index()).size(), 7);

    BOOST_CHECK_EQUAL( refinedGridData.size(3), 40);
    // LGR1 dims 2x3x2 -> 3x4x3 vertices + 4 extra missing vertices  (6,0,8), (6,2,8), (6,4,8), and (6,6,8).
    BOOST_CHECK_EQUAL( refinedGridData.numFaces(), 55);
    // LGR1 dims 2x3x2 -> 52 faces (before correction due to missing points)
    // 3 of those 52 faces vanished and give origin to 6 new faces: 52 - 3 + 6 = 55 faces

    // Originally, the element not involved in refinement
    // had  7 faces. It's neihgbor in level zero
    // got refined and the I_FACE that they shared has been
    // replaced by 6 refined faces. Then, the leaf element has
    // 1 of each I-,J-,J+, K-, K+, and 1 coarse + 6 refined I+.
    checkFaceCountInLeafCoarseElem(grid,
                                   /* expectedTotalFaceCount = */ 12,
                                   /* repeatedFaceType = */ 1, // 1->I+
                                   /* expoectedRepeatedFaceTypeCount = */ 7);

    expectedCellCountRespectToCellToFaceSize(grid,
                                             /* expectedCount7Faces = */ std::vector<int>{2, 3}); // level 0, level 1

    // Collect the expected data to later on check
    std::vector<std::vector<std::set<Coordinate,Opm::Lgr::FieldVectorLess>>> selectedFaceToCoord{};
    selectedFaceToCoord.resize(grid.levelGridView(1).size(0));

    for (const auto& refinedElem : Dune::elements(grid.levelGridView(1))) {

        std::set<Coordinate, Opm::Lgr::FieldVectorLess> expectedNewFaceInFace3{};
        std::set<Coordinate, Opm::Lgr::FieldVectorLess> expectedNewFaceInFace2{};
        // Vertex order in I_FACE: 0->jk, 1-> (j+1)k, 2->(j+1)(k+1), 3->j(k+1)
        //
        //         j(k+1) <-'3' --------- '2'-> (j+1)(k+1)
        //                   |             |
        //                   |             |
        //             jk <-'0' --------- '1'-> (j+1)k

        if (refinedElem.index() ==  6){
            // this element has to have 7 faces: 1 I+,J-,J+,K-,K+, and 2 I-:
            //      (6,0,9) --(6,2,9)
            //         |         *        I_FACE, false with vertices (6,0,8),(6,2,8),(6,2,9),(6,0,9)
            //      (6,0,8) **(6,2,8)
            //         |         *        I_FACE, false with vertices (6,0,5),(6,2,5),(6,2,8),(6,0,8)
            //         |         *
            //      (6,0,5) --(6,2,5)
            expectedNewFaceInFace3 = {{6,0,8},{6,2,8},{6,2,9},{6,0,9}};
            expectedNewFaceInFace2 = {{6,0,5},{6,2,5},{6,2,8},{6,0,8}};

            selectedFaceToCoord[6].push_back(expectedNewFaceInFace3);
            selectedFaceToCoord[6].push_back(expectedNewFaceInFace2);
        }
        else if (refinedElem.index() ==  8){
            // this element has to have 7 faces: 1 I+,J-,J+,K-,K+, and 2 I-:
            //      (6,2,9) --(6,4,9)
            //         |         *        I_FACE, false with vertices (6,2,8),(6,4,8),(6,4,9),(6,2,9)
            //      (6,2,8) **(6,4,8)
            //         |         *        I_FACE, false with vertices (6,2,5),(6,4,5),(6,4,8),(6,2,8)
            //         |         *
            //      (6,2,5) --(6,4,5)
            expectedNewFaceInFace3 = {{6,2,8},{6,4,8},{6,4,9},{6,2,9}};
            expectedNewFaceInFace2 = {{6,2,5},{6,4,5},{6,4,8},{6,2,8}};

            selectedFaceToCoord[8].push_back(expectedNewFaceInFace3);
            selectedFaceToCoord[8].push_back(expectedNewFaceInFace2);
        }
        else if (refinedElem.index() ==  10){
            // this element has to have 7 faces: 1 I+,J-,J+,K-,K+, and 2 I-:
            //      (6,4,9) --(6,6,9)
            //         |         *        I_FACE, false with vertices (6,4,8),(6,6,8),(6,6,9),(6,4,9)
            //      (6,4,8) **(6,6,8)
            //         |         *        I_FACE, false with vertices (6,4,5),(6,6,5),(6,6,8),(6,4,8)
            //         |         *
            //      (6,4,5) --(6,6,5)
            expectedNewFaceInFace3 = {{6,4,8},{6,6,8},{6,6,9},{6,4,9}};
            expectedNewFaceInFace2 = {{6,4,5},{6,6,5},{6,6,8},{6,4,8}};

            selectedFaceToCoord[10].push_back(expectedNewFaceInFace3);
            selectedFaceToCoord[10].push_back(expectedNewFaceInFace2);
        }
        else if (refinedElem.index() ==  0) {
            expectedNewFaceInFace2 = {{6,0,1}, {6,2,1}, {6,2,5}, {6,0,5}};
            selectedFaceToCoord[0].push_back(expectedNewFaceInFace2);
        }
        else if (refinedElem.index() ==  2) {
            expectedNewFaceInFace2 = {{6,2,1}, {6,4,1}, {6,4,5}, {6,2,5}};
            selectedFaceToCoord[2].push_back(expectedNewFaceInFace2);
        }
        else if (refinedElem.index() ==  4) {
            expectedNewFaceInFace2 = {{6,4,1}, {6,6,1}, {6,6,5}, {6,4,5}};
            selectedFaceToCoord[4].push_back(expectedNewFaceInFace2);
        }
    }
    checkNewRefinedFaces(grid, refinedGridData,
                         selectedFaceToCoord, /* repeatedFaceType = */ 0); // 0-> I-

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{2,3,2}},
                           /* lgr_name_vec = */ {"LGR1"});
}

BOOST_AUTO_TEST_CASE(parentCellWithMoreThanSixIntersections_J_FACE_true)
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

    BOOST_CHECK_EQUAL(parentGridData.cellToFace(parentElem.index()).size(), 7);

    BOOST_CHECK_EQUAL( refinedGridData.size(3), 40);
    // LGR1 dims 3x2x2 -> 4x3x3 vertices + 4 extra missing vertices  (0,6,1), (2,6,1), (4,6,1), and (6,6,1).
    BOOST_CHECK_EQUAL( refinedGridData.numFaces(), 55);
    // LGR1 dims 3x2x2 -> 52 faces (before correction due to missing points)
    // 3 of those 52 faces vanished and give origin to 6 new faces: 52 - 3 + 6 = 55 faces

    // Originally, the element not involved in refinement
    // had  7 faces. It's neihgbor in level zero
    // got refined and the I_FACE that they shared has been
    // replaced by 6 refined faces. Then, the leaf element has
    // one of each I-,I+,J+, K-, K+, and 1 coarse + 6 refined J-.
    checkFaceCountInLeafCoarseElem(grid,
                                   /* expectedTotalFaceCount = */ 12,
                                   /* repeatedFaceType = */ 2, // 2->J-
                                   /* expoectedRepeatedFaceTypeCount = */ 7);

    expectedCellCountRespectToCellToFaceSize(grid,
                                             /* expectedCount7Faces = */ std::vector<int>{2, 3}); // level 0, level 1

    // Collect the expected data to later on check
    std::vector<std::vector<std::set<Coordinate,Opm::Lgr::FieldVectorLess>>> selectedFaceToCoord{};
    selectedFaceToCoord.resize(grid.levelGridView(1).size(0));

    for (const auto& refinedElem : Dune::elements(grid.levelGridView(1))) {

        std::set<Coordinate, Opm::Lgr::FieldVectorLess> expectedNewFaceInFace6{};
        std::set<Coordinate, Opm::Lgr::FieldVectorLess> expectedNewFaceInFace5{};
        // Vertex order in J_FACE: 0->(i+1)k, 1-> ik, 2->i(k+1), 3->(i+1)(k+1)
        //
        //         i(k+1) <-'2' --------- '3'-> (i+1)(k+1)
        //                   |             |
        //                   |             |
        //             ik <-'1' --------- '0'-> (i+1)k

        if (refinedElem.index() ==  3){
            // this element has to have 7 faces: 1 I-,I+,J-,K-,K+, and 2 J+:
            //      (0,6,4) **(2,6,4)
            //         |         *        J_FACE, true with vertices (0,6,1),(2,6,1),(2,6,4),(,6,4)
            //         |         *
            //      (0,6,1) --(2,6,1)
            //         |         *        J_FACE, true with vertices (0,6,0),(2,6,0),(2,6,1),(0,6,1)
            //      (0,6,0) --(2,6,0)
            expectedNewFaceInFace6 = {{2,6,1}, {0,6,1}, {0,6,4}, {2,6,4}};
            expectedNewFaceInFace5 = {{2,6,0}, {0,6,0}, {0,6,1}, {2,6,1}};

            selectedFaceToCoord[3].push_back(expectedNewFaceInFace6);
            selectedFaceToCoord[3].push_back(expectedNewFaceInFace5);
        }
        else if (refinedElem.index() == 4){
            // this element has to have 7 faces: 1 I-,I+,J-,K-,K+, and 2 J+:
            //      (2,6,4) **(4,6,4)
            //         |         *        J_FACE, true with vertices (2,6,1),(4,6,1),(4,6,4),(2,6,4)
            //         |         *
            //      (2,6,1) --(4,6,1)
            //         |         *        J_FACE, true with vertices (2,6,0),(4,6,0),(4,6,1),(2,6,1)
            //      (2,6,0) --(4,6,0)
            expectedNewFaceInFace6 = {{4,6,1}, {2,6,1}, {2,6,4}, {4,6,4}};
            expectedNewFaceInFace5 = {{4,6,0}, {2,6,0}, {2,6,1}, {4,6,1}};

            selectedFaceToCoord[4].push_back(expectedNewFaceInFace6);
            selectedFaceToCoord[4].push_back(expectedNewFaceInFace5);
        }
        else if (refinedElem.index() == 5){
            // this element has to have 7 faces: 1 I-,I+,J-,K-,K+, and 2 J+:
            //      (4,6,4) **(6,6,4)
            //         |         *        J_FACE, true with vertices (4,6,1),(6,6,1),(6,6,4),(4,6,4)
            //         |         *
            //      (4,6,1) --(6,6,1)
            //         |         *        J_FACE, true with vertices (4,6,0),(6,6,0),(6,6,1),(4,6,1)
            //      (4,6,0) --(6,6,0)
            expectedNewFaceInFace6 = {{6,6,1}, {4,6,1}, {4,6,4}, {6,6,4}};
            expectedNewFaceInFace5 = {{6,6,0}, {4,6,0}, {4,6,1}, {6,6,1}};

            selectedFaceToCoord[5].push_back(expectedNewFaceInFace6);
            selectedFaceToCoord[5].push_back(expectedNewFaceInFace5);
        }
        else if (refinedElem.index() ==  9) {
            expectedNewFaceInFace6 = {{2,6,4}, {0,6,4}, {0,6,8}, {2,6,8}};
            selectedFaceToCoord[9].push_back(expectedNewFaceInFace6);
        }
        else if (refinedElem.index() ==  10) {
            expectedNewFaceInFace6 = {{4,6,4}, {2,6,4}, {2,6,8}, {4,6,8}};
            selectedFaceToCoord[10].push_back(expectedNewFaceInFace6);
        }
        else if (refinedElem.index() ==  11) {
            expectedNewFaceInFace6 = {{6,6,4}, {4,6,4}, {4,6,8}, {6,6,8}};
            selectedFaceToCoord[11].push_back(expectedNewFaceInFace6);
        }
    }
    checkNewRefinedFaces(grid, refinedGridData,
                         selectedFaceToCoord, /* repeatedFaceType = */ 3); // 3-> J+

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{3,2,2}},
                           /* lgr_name_vec = */ {"LGR1"});
}


BOOST_AUTO_TEST_CASE(parentCellWithMoreThanSixIntersections_J_FACE_false)
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

    BOOST_CHECK_EQUAL(parentGridData.cellToFace(parentElem.index()).size(), 7);

    BOOST_CHECK_EQUAL( refinedGridData.size(3), 40);
    // LGR1 dims 3x2x2 -> 4x3x3 vertices + 4 extra missing vertices  (0,6,8), (2,6,8), (4,6,8), and (6,6,8).
    BOOST_CHECK_EQUAL( refinedGridData.numFaces(), 55);
    // LGR1 dims 3x2x2 -> 52 faces (before correction due to missing points)
    // 3 of those 52 faces vanished and give origin to 6 new faces: 52 - 3 + 6 = 55 faces

    // Originally, the element not involved in refinement
    // had  7 faces. It's neihgbor in level zero
    // got refined and the I_FACE that they shared has been
    // replaced by 6 refined faces. Then, the leaf element has
    // one of each I-,I+,J-, K-, K+, and 1 coarse + 6 refined J+.
    checkFaceCountInLeafCoarseElem(grid,
                                   /* expectedTotalFaceCount = */ 12,
                                   /* repeatedFaceType = */ 3, // 3->J+
                                   /* expoectedRepeatedFaceTypeCount = */ 7);

    expectedCellCountRespectToCellToFaceSize(grid,
                                             /* expectedCount7Faces = */ std::vector<int>{2, 3}); // level 0, level 1

    // Collect the expected data to later on check
    std::vector<std::vector<std::set<Coordinate,Opm::Lgr::FieldVectorLess>>> selectedFaceToCoord{};
    selectedFaceToCoord.resize(grid.levelGridView(1).size(0));

    for (const auto& refinedElem : Dune::elements(grid.levelGridView(1))) {

        std::set<Coordinate, Opm::Lgr::FieldVectorLess> expectedNewFaceInFace7{};
        std::set<Coordinate, Opm::Lgr::FieldVectorLess> expectedNewFaceInFace6{};
        // Vertex order in J_FACE: 0->(i+1)k, 1-> ik, 2->i(k+1), 3->(i+1)(k+1)
        //
        //         i(k+1) <-'2' --------- '3'-> (i+1)(k+1)
        //                   |             |
        //                   |             |
        //             ik <-'1' --------- '0'-> (i+1)k

        if (refinedElem.index() ==  6) {
            // this element has to have 7 faces: 1 I-,I+,J+,K-,K+, and 2 J-:
            //      (0,6,9) --(2,6,9)
            //         |         *        J_FACE, false with vertices (0,6,8),(2,6,8),(2,6,9),(0,6,9)
            //      (0,6,8) **(2,6,8)
            //         |         *        J_FACE, false with vertices (0,6,5),(2,6,5),(2,6,8),(0,6,8)
            //         |         *
            //      (0,6,5) --(2,6,5)
            expectedNewFaceInFace7 = {{2,6,8}, {0,6,8}, {0,6,9}, {2,6,9}};
            expectedNewFaceInFace6 = {{2,6,5}, {0,6,5}, {0,6,8}, {2,6,8}};

            selectedFaceToCoord[6].push_back(expectedNewFaceInFace7);
            selectedFaceToCoord[6].push_back(expectedNewFaceInFace6);
        }
        else if (refinedElem.index() == 7) {
            // this element has to have 7 faces: 1 I-,I+,J+,K-,K+, and 2 J-:
            //      (2,6,9) --(4,6,9)
            //         |         *        J_FACE, false with vertices (2,6,8),(4,6,8),(4,6,9),(2,6,9)
            //      (2,6,8) **(4,6,8)
            //         |         *        J_FACE, false with vertices (2,6,5),(4,6,5),(4,6,8),(2,6,8)
            //         |         *
            //      (2,6,5) --(4,6,5)
            expectedNewFaceInFace7 = {{4,6,8}, {2,6,8}, {2,6,9}, {4,6,9}};
            expectedNewFaceInFace6 = {{4,6,5}, {2,6,5}, {2,6,8}, {4,6,8}};

            selectedFaceToCoord[7].push_back(expectedNewFaceInFace7);
            selectedFaceToCoord[7].push_back(expectedNewFaceInFace6);
        }
        else if (refinedElem.index() == 8) {
            // this element has to have 7 faces: 1 I-,I+,J+,K-,K+, and 2 J-:
            //      (4,6,9) --(6,6,9)
            //         |         *        J_FACE, false with vertices (4,6,8),(6,6,8),(6,6,9),(4,6,9)
            //      (4,6,8) **(6,6,8)
            //         |         *        J_FACE, false with vertices (4,6,5),(6,6,5),(6,6,8),(4,6,8)
            //         |         *
            //      (4,6,5) --(6,6,5)
            expectedNewFaceInFace7 = {{6,6,8}, {4,6,8}, {4,6,9}, {6,6,9}};
            expectedNewFaceInFace6 = {{6,6,5}, {4,6,5}, {4,6,8}, {6,6,8}};

            selectedFaceToCoord[8].push_back(expectedNewFaceInFace7);
            selectedFaceToCoord[8].push_back(expectedNewFaceInFace6);
        }
        else if (refinedElem.index() ==  0) {
            expectedNewFaceInFace6 = {{2,6,1}, {0,6,1}, {0,6,5}, {2,6,5}};
            selectedFaceToCoord[0].push_back(expectedNewFaceInFace6);
        }
        else if (refinedElem.index() ==  1) {
            expectedNewFaceInFace6 = {{4,6,1}, {2,6,1}, {2,6,5}, {4,6,5}};
            selectedFaceToCoord[1].push_back(expectedNewFaceInFace6);
        }
        else if (refinedElem.index() ==  2) {
            expectedNewFaceInFace6 = {{6,6,1},{4,6,1}, {4,6,5}, {6,6,5}};
            selectedFaceToCoord[2].push_back(expectedNewFaceInFace6);
        }
    }
    checkNewRefinedFaces(grid, refinedGridData,
                         selectedFaceToCoord, /* repeatedFaceType = */ 2); // 2-> J-

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{3,2,2}},
                           /* lgr_name_vec = */ {"LGR1"});
}


BOOST_AUTO_TEST_CASE(neighboringSingleCellRefinementsDifferentLgrs)
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
                              {{2,3,2}, {2,3,2}}, // cells_per_dim_vec
                              {{0,0,0}, {1,0,0}}, // startIJK_vec
                              {{1,1,1}, {2,1,1}}, // endIJK_vec
                              {"LGR1", "LGR2"});  // lgr_name_vec*/

    // Element 0 and element 1 in level zero grid share an I_FACE (with face index 2)
    //
    // Vertices of those faces lie on the plane x = 6    | After refinement, number of subdivisions in
    //                                                   | y- and z- directions:
    //
    //              (6,0,9) -----------------(6,6,9)     |  (6,0,9) --(6,2,9)-(6,4,9)--(6,6,9)
    //                 |      face idx 3      |          |     |         *       *        |
    //              (6,0,8) ---------------- (6,6,8)     |  (6,0,8) --(6,2,8)-(6,4,8)--(6,6,8)
    //                 |                      |          |     |         *       *        |
    //                 |                      |          |     |         *       *        |
    //                 |                      |          |     |         *       *        |
    //                 |      face idx 2      |          |  (6,0,5) **(6,2,5)*(6,4,5)**(6,6,5)
    //                 |                      |          |     |         *       *        |
    //                 |                      |          |  (6,0,4) **(6,2,4)*(6,4,4)**(6,6,4)
    //                 |                      |          |     |         *       *        |
    //                 |                      |          |     |         *       *        |
    //              (6,0,1) -----------------(6,6,1)     |  (6,0,1) --(6,2,1)-(6,4,1)--(6,6,1)
    //                 |      face idx 1      |          |     |         *       *        |
    //              (6,0,0) -----------------(6,6,0)     |  (6,0,0) --(6,2,0)-(6,4,0)--(6,6,0)
    //                                                   |

    const auto& parentGridData = *grid.currentData()[0];
    const auto parent0 = Dune::cpgrid::Entity<0>(parentGridData, 0, true);
    const auto parent1 = Dune::cpgrid::Entity<0>(parentGridData, 1, true);

    BOOST_CHECK_EQUAL(parentGridData.cellToFace(parent0.index()).size(), 7);
    BOOST_CHECK_EQUAL(parentGridData.cellToFace(parent1.index()).size(), 7);

    const auto& refinedGrid1 = *grid.currentData()[1];
    BOOST_CHECK_EQUAL( refinedGrid1.size(3), 44);
    // LGR1 dims 3x2x2 -> 4x3x3 vertices + 4 extra missing vertices  (0,6,8), (2,6,8), (4,6,8), and (6,6,8)......
    BOOST_CHECK_EQUAL( refinedGrid1.numFaces(), 58);

    const auto& refinedGrid2 = *grid.currentData()[2];
    BOOST_CHECK_EQUAL( refinedGrid2.size(3), 44);
    // LGR1 dims 3x2x2 -> 4x3x3 vertices + 4 extra missing vertices  (0,6,8), (2,6,8), (4,6,8), and (6,6,8)........
    BOOST_CHECK_EQUAL( refinedGrid2.numFaces(), 58);

    const auto& leafGrid = grid.currentLeafData();

    BOOST_CHECK_EQUAL( leafGrid.size(3), 72); // 44 + 44 -16 -> 72
    BOOST_CHECK_EQUAL( leafGrid.numFaces(), 107); // 58 + 58 - 9 ->107

    expectedCellCountRespectToCellToFaceSize(grid,
                                             /* expectedCount7Faces = */ std::vector<int>{2, 6, 6}); // level 0, level 1, level2

    // Leaf grid should have 12 elements with 7 faces (due to I_FACE shared between
    // parent cells) and 24-12 = 12 with 6 faces
    int count_7faces = 0;
    int count_6faces = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        int count = leafGrid.cellToFace(element.index()).size();
        if (count == 7) {
            ++count_7faces;
        }
        else {
            BOOST_CHECK_EQUAL(count, 6);
            ++count_6faces;
        }
    }
    BOOST_CHECK_EQUAL(count_7faces, 12);
    BOOST_CHECK_EQUAL(count_6faces, grid.leafGridView().size(0) - count_7faces);


    Opm::checkGridWithLgrs(grid,
                           {{2,3,2}, {2,3,2}}, // cells_per_dim_vec
                           {"LGR1", "LGR2"}); // lgr_name_vec

}

BOOST_AUTO_TEST_CASE(neighboringSingleCellRefinementsSameLgr)
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
                              /* endIJK_vec */        {{2,1,1}},
                              /* lgr_name_vec */      {"LGR1"});

    // Element 0 and element 1 in level zero grid share an I_FACE (with face index 2)
    //
    // Vertices of those faces lie on the plane x = 6    | After refinement, number of subdivisions in
    //                                                   | y- and z- directions:
    //
    //              (6,0,9) -----------------(6,6,9)     |  (6,0,9) --(6,2,9)-(6,4,9)--(6,6,9)
    //                 |      face idx 3      |          |     |         *       *        |
    //              (6,0,8) ---------------- (6,6,8)     |  (6,0,8) --(6,2,8)-(6,4,8)--(6,6,8)
    //                 |                      |          |     |         *       *        |
    //                 |                      |          |     |         *       *        |
    //                 |                      |          |     |         *       *        |
    //                 |      face idx 2      |          |  (6,0,5) **(6,2,5)*(6,4,5)**(6,6,5)
    //                 |                      |          |     |         *       *        |
    //                 |                      |          |  (6,0,4) **(6,2,4)*(6,4,4)**(6,6,4)
    //                 |                      |          |     |         *       *        |
    //                 |                      |          |     |         *       *        |
    //              (6,0,1) -----------------(6,6,1)     |  (6,0,1) --(6,2,1)-(6,4,1)--(6,6,1)
    //                 |      face idx 1      |          |     |         *       *        |
    //              (6,0,0) -----------------(6,6,0)     |  (6,0,0) --(6,2,0)-(6,4,0)--(6,6,0)
    //                                                   |
    const auto& parentGridData = *grid.currentData()[0];
    const auto parent0 = Dune::cpgrid::Entity<0>(parentGridData, 0, true);
    const auto parent1 = Dune::cpgrid::Entity<0>(parentGridData, 1, true);

    BOOST_CHECK_EQUAL(parentGridData.cellToFace(parent0.index()).size(), 7);
    BOOST_CHECK_EQUAL(parentGridData.cellToFace(parent1.index()).size(), 7);

    const auto& refinedGrid1 = *grid.currentData()[1];
    BOOST_CHECK_EQUAL( refinedGrid1.size(3), 72);
    BOOST_CHECK_EQUAL( refinedGrid1.numFaces(), 107);

    const auto& leafGrid = grid.currentLeafData();
    BOOST_CHECK_EQUAL( leafGrid.size(3), 72); // 44 + 44 -16 -> 72
    BOOST_CHECK_EQUAL( leafGrid.numFaces(), 107); // 58 + 58 - 9 ->107

    expectedCellCountRespectToCellToFaceSize(grid,
                                             /* expectedCount7Faces = */ std::vector<int>{2, 12}); // level 0, level 1

    // Leaf grid should have 12 elements with 7 faces (due to I_FACE shared between
    // parent cells) and 24-12 = 12 with 6 faces
    int count_7faces = 0;
    int count_6faces = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        int count = leafGrid.cellToFace(element.index()).size();
        if (count == 7) {
            ++count_7faces;
        }
        else {
            BOOST_CHECK_EQUAL(count, 6);
            ++count_6faces;
        }
    }
    BOOST_CHECK_EQUAL(count_7faces, 12);
    BOOST_CHECK_EQUAL(count_6faces, grid.leafGridView().size(0) - count_7faces);

    Opm::checkGridWithLgrs(grid,
                           {{2,3,2}}, // cells_per_dim_vec
                           {"LGR1"}); // lgr_name_vec
}



BOOST_AUTO_TEST_CASE(simpleDiffLgr)
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
                              /* cells_per_dim_vec */ {{1,2,1}, {1,2,1}},
                              /* startIJK_vec */      {{0,0,0}, {1,0,0}},
                              /* endIJK_vec */        {{1,1,1}, {2,1,1}},
                              /* lgr_name_vec */      {"LGR1", "LGR2"});

    // Element 0 and element 1 in level zero grid share an I_FACE (with face index 2)
    //
    // Vertices of those faces lie on the plane x = 6    | After refinement, number of subdivisions in
    //                                                   | y- and z- directions:
    //
    //              (6,0,9) -----------------(6,6,9)     |  (6,0,9) --(6,3,9)---(6,6,9)
    //                 |      face idx 3      |          |     |         *         |
    //              (6,0,8) ---------------- (6,6,8)     |  (6,0,8) --(6,3,9)---(6,6,8)
    //                 |                      |          |     |         *         |
    //                 |                      |          |     |         *         |
    //                 |                      |          |     |         *         |
    //                 |      face idx 2      |          |     |         *         |
    //                 |                      |          |     |         *         |
    //                 |                      |          |     |         *         |
    //                 |                      |          |     |         *         |
    //                 |                      |          |     |         *         |
    //              (6,0,1) -----------------(6,6,1)     |  (6,0,1) --(6,3,1)---(6,6,1)
    //                 |      face idx 1      |          |     |         *         |
    //              (6,0,0) -----------------(6,6,0)     |  (6,0,0) --(6,3,0)---(6,6,0)
    //                                                   |
    const auto& parentGridData = *grid.currentData()[0];
    const auto parent0 = Dune::cpgrid::Entity<0>(parentGridData, 0, true);
    const auto parent1 = Dune::cpgrid::Entity<0>(parentGridData, 1, true);

    BOOST_CHECK_EQUAL(parentGridData.cellToFace(parent0.index()).size(), 7);
    BOOST_CHECK_EQUAL(parentGridData.cellToFace(parent1.index()).size(), 7);

    const auto& refinedGrid1 = *grid.currentData()[1];
    BOOST_CHECK_EQUAL( refinedGrid1.size(3), 15);
    BOOST_CHECK_EQUAL( refinedGrid1.numFaces(), 13);

    const auto& refinedGrid2 = *grid.currentData()[2];
    BOOST_CHECK_EQUAL( refinedGrid2.size(3), 15);
    BOOST_CHECK_EQUAL( refinedGrid2.numFaces(), 13);

    const auto& leafGrid = grid.currentLeafData();
    BOOST_CHECK_EQUAL( leafGrid.size(3), 24);
    BOOST_CHECK_EQUAL( leafGrid.numFaces(), 24);


    const std::set<Coordinate,Opm::Lgr::FieldVectorLess> expectedLGR1Vertices = {
        {0.,0.,0.}, {0.,0.,8.},              // pillar x = 0,  y = 0
        {6.,0.,0.}, {6.,0.,1.}, {6.,0.,8.},  // pillar x = 6,  y = 0
        {0.,3.,0.}, {0.,3.,8.},              // pillar x = 0,  y = 3
        {6.,3.,0.}, {6.,3.,1.}, {6.,3.,8.},  // pillar x = 6,  y = 3
        {0.,6.,0.}, {0.,6.,8.},              // pillar x = 0,  y = 6
        {6.,6.,0.}, {6.,6.,1.}, {6.,6.,8.}   // pillar x = 6,  y = 6
    };

    const std::set<Coordinate,Opm::Lgr::FieldVectorLess> expectedLGR2Vertices = {
        {6.,0.,1.}, {6.,0.,8.}, {6.,0.,9.}, // pillar x = 6,  y = 0
        {12.,0.,1.}, {12.,0.,9.},           // pillar x = 12, y = 0
        {6.,3.,1.}, {6.,3.,8.}, {6.,3.,9.}, // pillar x = 6,  y = 3
        {12.,3.,1.}, {12.,3.,9.},           // pillar x = 12, y = 3
        {6.,6.,1.}, {6.,6.,8.}, {6.,6.,9.}, // pillar x = 6,  y = 6
        {12.,6.,1.}, {12.,6.,9.}            // pillar x = 12, y = 6
    };

    const std::set<Coordinate,Opm::Lgr::FieldVectorLess> expectedLeafVertices = {
        {0.,0.,0.}, {0.,0.,8.},                          // pillar x = 0,  y = 0
        {6.,0.,0.}, {6.,0.,1.}, {6.,0.,8.}, {6.,0.,9.},  // pillar x = 6,  y = 0
        {12.,0.,1.}, {12.,0.,9.},                        // pillar x = 12, y = 0
        {0.,3.,0.}, {0.,3.,8.},                          // pillar x = 0,  y = 3
        {6.,3.,0.}, {6.,3.,1.}, {6.,3.,8.}, {6.,3.,9.},  // pillar x = 6,  y = 3
        {12.,3.,1.}, {12.,3.,9.},                        // pillar x = 12, y = 3
        {0.,6.,0.}, {0.,6.,8.},                          // pillar x = 0,  y = 6
        {6.,6.,0.}, {6.,6.,1.}, {6.,6.,8.}, {6.,6.,9.},  // pillar x = 6,  y = 6
        {12.,6.,1.}, {12.,6.,9.}                         // pillar x = 12, y = 6
    };

    checkLevelVertices(grid, std::vector{expectedLGR1Vertices, expectedLGR2Vertices});
    checkLeafVertices(grid, expectedLeafVertices);

    const std::vector<std::vector<Coordinate>> expectedLGR1Faces = {
        {{0.,0.,0.}, {6.,0.,0.}, {6.,3.,0.}, {0.,3.,0.}},      // K_FACE   z = 0, face 0
        {{0.,3.,0.}, {6.,3.,0.}, {6.,6.,0.}, {0.,6.,0.}},      // K_FACE   z = 0, face 1
        {{0.,0.,8.}, {6.,0.,8.}, {6.,3.,8.}, {0.,3.,8.}},      // K_FACE   z = 8, face 2
        {{0.,3.,8.}, {6.,3.,8.}, {6.,6.,8.}, {0.,6.,8.}},      // K_FACE   z = 8, face 3
        {{0.,0.,0.}, {0.,3.,0.}, {0.,3.,8.}, {0.,0.,8.}},      // I_FACE   x = 0, face 4
        {{0.,3.,0.}, {0.,6.,0.}, {0.,6.,8.}, {0.,3.,8.}},      // I_FACE   x = 0, face 5
        {{6.,0.,0.}, {6.,3.,0.}, {6.,3.,1.}, {6.,0.,1.}},      // I_FACE   x = 6, face 6
        {{6.,0.,1.}, {6.,3.,1.}, {6.,3.,8.}, {6.,0.,8.}},      // I_FACE   x = 6, face 7
        {{6.,3.,0.}, {6.,6.,0.}, {6.,6.,1.}, {6.,3.,1.}},      // I_FACE   x = 6, face 8
        {{6.,3.,1.}, {6.,6.,1.}, {6.,6.,8.}, {6.,3.,8.}},      // I_FACE   x = 6, face 9
        {{0.,0.,0.}, {6.,0.,0.}, {6.,0.,8.}, {0.,0.,8.}},      // J_FACE   y = 0, face 10
        {{0.,3.,0.}, {6.,3.,0.}, {6.,3.,8.}, {0.,3.,8.}},      // J_FACE   y = 3, face 11
        {{0.,6.,0.}, {6.,6.,0.}, {6.,6.,8.}, {0.,6.,8.}},      // J_FACE   y = 6, face 12
    };

    const std::vector<std::vector<Coordinate>> expectedLGR2Faces = {
        {{6.,0.,1.}, {12.,0.,1.}, {12.,3.,1.}, {6.,3.,1.}},    // K_FACE   z = 1,  face 0
        {{6.,3.,1.}, {12.,3.,1.}, {12.,6.,1.}, {6.,6.,1.}},    // K_FACE   z = 1,  face 1
        {{6.,0.,9.}, {12.,0.,9.}, {12.,3.,9.}, {6.,3.,9.}},    // K_FACE   z = 9,  face 2
        {{6.,3.,9.}, {12.,3.,9.}, {12.,6.,9.}, {6.,6.,9.}},    // K_FACE   z = 9,  face 3
        {{6.,0.,1.}, {6.,3.,1.}, {6.,3.,8.}, {6.,0.,8.}},      // I_FACE   x = 6,  face 4
        {{6.,0.,8.}, {6.,3.,8.}, {6.,3.,9.}, {6.,0.,9.}},      // I_FACE   x = 6,  face 5
        {{6.,3.,1.}, {6.,6.,1.}, {6.,6.,8.}, {6.,3.,8.}},      // I_FACE   x = 6,  face 6
        {{6.,3.,8.}, {6.,6.,8.}, {6.,6.,9.}, {6.,3.,9.}},      // I_FACE   x = 6,  face 7
        {{12.,0.,1.}, {12.,3.,1.}, {12.,3.,9.}, {12.,0.,9.}},  // I_FACE   x = 12, face 8
        {{12.,3.,1.}, {12.,6.,1.}, {12.,6.,9.}, {12.,3.,9.}},  // I_FACE   x = 12, face 9
        {{6.,0.,1.}, {12.,0.,1.}, {12.,0.,9.}, {6.,0.,9.}},    // J_FACE   y = 0,  face 10
        {{6.,3.,1.}, {12.,3.,1.}, {12.,3.,9.}, {6.,3.,9.}},    // J_FACE   y = 3,  face 11
        {{6.,6.,1.}, {12.,6.,1.}, {12.,6.,9.}, {6.,6.,9.}},    // J_FACE   y = 6,  face 12
    };

    const std::vector<std::vector<Coordinate>> expectedLeafFaces = {
        {{0.,0.,0.}, {6.,0.,0.}, {6.,3.,0.}, {0.,3.,0.}},      // K_FACE   z = 0,  face 0
        {{0.,3.,0.}, {6.,3.,0.}, {6.,6.,0.}, {0.,6.,0.}},      // K_FACE   z = 0,  face 1
        {{0.,0.,8.}, {6.,0.,8.}, {6.,3.,8.}, {0.,3.,8.}},      // K_FACE   z = 8,  face 2
        {{0.,3.,8.}, {6.,3.,8.}, {6.,6.,8.}, {0.,6.,8.}},      // K_FACE   z = 8,  face 3
        {{0.,0.,0.}, {0.,3.,0.}, {0.,3.,8.}, {0.,0.,8.}},      // I_FACE   x = 0,  face 6
        {{0.,3.,0.}, {0.,6.,0.}, {0.,6.,8.}, {0.,3.,8.}},      // I_FACE   x = 0,  face 5
        {{6.,0.,0.}, {6.,3.,0.}, {6.,3.,1.}, {6.,0.,1.}},      // I_FACE   x = 6,  face 6
        {{6.,3.,0.}, {6.,6.,0.}, {6.,6.,1.}, {6.,3.,1.}},      // I_FACE   x = 6,  face 7
        {{0.,0.,0.}, {6.,0.,0.}, {6.,0.,8.}, {0.,0.,8.}},      // J_FACE   y = 0,  face 8
        {{0.,3.,0.}, {6.,3.,0.}, {6.,3.,8.}, {0.,3.,8.}},      // J_FACE   y = 3,  face 9
        {{0.,6.,0.}, {6.,6.,0.}, {6.,6.,8.}, {0.,6.,8.}},      // J_FACE   y = 6,  face 10
        {{6.,0.,1.}, {12.,0.,1.}, {12.,3.,1.}, {6.,3.,1.}},    // K_FACE   z = 1,  face 11
        {{6.,3.,1.}, {12.,3.,1.}, {12.,6.,1.}, {6.,6.,1.}},    // K_FACE   z = 1,  face 12
        {{6.,0.,9.}, {12.,0.,9.}, {12.,3.,9.}, {6.,3.,9.}},    // K_FACE   z = 9,  face 13
        {{6.,3.,9.}, {12.,3.,9.}, {12.,6.,9.}, {6.,6.,9.}},    // K_FACE   z = 9,  face 14
        {{6.,0.,1.}, {6.,3.,1.}, {6.,3.,8.}, {6.,0.,8.}},      // I_FACE   x = 6,  face 15
        {{6.,0.,8.}, {6.,3.,8.}, {6.,3.,9.}, {6.,0.,9.}},      // I_FACE   x = 6,  face 16
        {{6.,3.,1.}, {6.,6.,1.}, {6.,6.,8.}, {6.,3.,8.}},      // I_FACE   x = 6,  face 17
        {{6.,3.,8.}, {6.,6.,8.}, {6.,6.,9.}, {6.,3.,9.}},      // I_FACE   x = 6,  face 18
        {{12.,0.,1.}, {12.,3.,1.}, {12.,3.,9.}, {12.,0.,9.}},  // I_FACE   x = 12, face 19
        {{12.,3.,1.}, {12.,6.,1.}, {12.,6.,9.}, {12.,3.,9.}},  // I_FACE   x = 12, face 20
        {{6.,0.,1.}, {12.,0.,1.}, {12.,0.,9.}, {6.,0.,9.}},    // J_FACE   y = 0,  face 21
        {{6.,3.,1.}, {12.,3.,1.}, {12.,3.,9.}, {6.,3.,9.}},    // J_FACE   y = 3,  face 22
        {{6.,6.,1.}, {12.,6.,1.}, {12.,6.,9.}, {6.,6.,9.}},    // J_FACE   y = 6,  face 23
    };

    checkFaces(refinedGrid1, expectedLGR1Faces);
    checkFaces(refinedGrid2, expectedLGR2Faces);
    checkFaces(leafGrid, expectedLeafFaces);

    Opm::checkGridWithLgrs(grid,
                           {{1,2,1}, {1,2,1}}, // cells_per_dim_vec
                           {"LGR1", "LGR2"}); // lgr_name_vec
}


BOOST_AUTO_TEST_CASE(simpleSameLgr)
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
                              /* cells_per_dim_vec */ {{1,2,1}},
                              /* startIJK_vec */      {{0,0,0}},
                              /* endIJK_vec */        {{2,1,1}},
                              /* lgr_name_vec */      {"LGR1"});

    // Element 0 and element 1 in level zero grid share an I_FACE (with face index 2)
    //
    // Vertices of those faces lie on the plane x = 6    | After refinement, number of subdivisions in
    //                                                   | y- and z- directions:
    //
    //              (6,0,9) -----------------(6,6,9)     |  (6,0,9) --(6,3,9)---(6,6,9)
    //                 |      face idx 3      |          |     |         *         |
    //              (6,0,8) ---------------- (6,6,8)     |  (6,0,8) --(6,3,9)---(6,6,8)
    //                 |                      |          |     |         *         |
    //                 |                      |          |     |         *         |
    //                 |                      |          |     |         *         |
    //                 |      face idx 2      |          |     |         *         |
    //                 |                      |          |     |         *         |
    //                 |                      |          |     |         *         |
    //                 |                      |          |     |         *         |
    //                 |                      |          |     |         *         |
    //              (6,0,1) -----------------(6,6,1)     |  (6,0,1) --(6,3,1)---(6,6,1)
    //                 |      face idx 1      |          |     |         *         |
    //              (6,0,0) -----------------(6,6,0)     |  (6,0,0) --(6,3,0)---(6,6,0)
    //                                                   |
    const auto& parentGridData = *grid.currentData()[0];
    const auto parent0 = Dune::cpgrid::Entity<0>(parentGridData, 0, true);
    const auto parent1 = Dune::cpgrid::Entity<0>(parentGridData, 1, true);

    BOOST_CHECK_EQUAL(parentGridData.cellToFace(parent0.index()).size(), 7);
    BOOST_CHECK_EQUAL(parentGridData.cellToFace(parent1.index()).size(), 7);

    const auto& refinedGrid1 = *grid.currentData()[1];
    BOOST_CHECK_EQUAL( refinedGrid1.size(3), 24);
    BOOST_CHECK_EQUAL( refinedGrid1.numFaces(), 24);

    const auto& leafGrid = grid.currentLeafData();
    BOOST_CHECK_EQUAL( leafGrid.size(3), 24);
    BOOST_CHECK_EQUAL( leafGrid.numFaces(), 24);

    const std::set<Coordinate,Opm::Lgr::FieldVectorLess> expectedVertices = {
        {0.,0.,0.}, {0.,0.,8.},                          // pillar x = 0,  y = 0
        {6.,0.,0.}, {6.,0.,1.}, {6.,0.,8.}, {6.,0.,9.},  // pillar x = 6,  y = 0
        {12.,0.,1.}, {12.,0.,9.},                        // pillar x = 12, y = 0
        {0.,3.,0.}, {0.,3.,8.},                          // pillar x = 0,  y = 3
        {6.,3.,0.}, {6.,3.,1.}, {6.,3.,8.}, {6.,3.,9.},  // pillar x = 6,  y = 3
        {12.,3.,1.}, {12.,3.,9.},                        // pillar x = 12, y = 3
        {0.,6.,0.}, {0.,6.,8.},                          // pillar x = 0,  y = 6
        {6.,6.,0.}, {6.,6.,1.}, {6.,6.,8.}, {6.,6.,9.},  // pillar x = 6,  y = 6
        {12.,6.,1.}, {12.,6.,9.}                         // pillar x = 12, y = 6
    };

    checkLevelVertices(grid, std::vector{expectedVertices});
    checkLeafVertices(grid, expectedVertices);

    const std::vector<std::vector<Coordinate>> expectedFaces = {
        {{0.,0.,0.}, {6.,0.,0.}, {6.,3.,0.}, {0.,3.,0.}},      // K_FACE   z = 0,  face 0
        {{0.,3.,0.}, {6.,3.,0.}, {6.,6.,0.}, {0.,6.,0.}},      // K_FACE   z = 0,  face 1
        {{0.,0.,8.}, {6.,0.,8.}, {6.,3.,8.}, {0.,3.,8.}},      // K_FACE   z = 8,  face 2
        {{0.,3.,8.}, {6.,3.,8.}, {6.,6.,8.}, {0.,6.,8.}},      // K_FACE   z = 8,  face 3
        {{0.,0.,0.}, {0.,3.,0.}, {0.,3.,8.}, {0.,0.,8.}},      // I_FACE   x = 0,  face 4
        {{0.,3.,0.}, {0.,6.,0.}, {0.,6.,8.}, {0.,3.,8.}},      // I_FACE   x = 0,  face 5
        {{6.,0.,0.}, {6.,3.,0.}, {6.,3.,1.}, {6.,0.,1.}},      // I_FACE   x = 6,  face 6
        {{6.,3.,0.}, {6.,6.,0.}, {6.,6.,1.}, {6.,3.,1.}},      // I_FACE   x = 6,  face 7
        {{0.,0.,0.}, {6.,0.,0.}, {6.,0.,8.}, {0.,0.,8.}},      // J_FACE   y = 0,  face 8
        {{0.,3.,0.}, {6.,3.,0.}, {6.,3.,8.}, {0.,3.,8.}},      // J_FACE   y = 3,  face 9
        {{0.,6.,0.}, {6.,6.,0.}, {6.,6.,8.}, {0.,6.,8.}},      // J_FACE   y = 6,  face 10
        {{6.,0.,1.}, {12.,0.,1.}, {12.,3.,1.}, {6.,3.,1.}},    // K_FACE   z = 1,  face 11
        {{6.,3.,1.}, {12.,3.,1.}, {12.,6.,1.}, {6.,6.,1.}},    // K_FACE   z = 1,  face 12
        {{6.,0.,9.}, {12.,0.,9.}, {12.,3.,9.}, {6.,3.,9.}},    // K_FACE   z = 9,  face 13
        {{6.,3.,9.}, {12.,3.,9.}, {12.,6.,9.}, {6.,6.,9.}},    // K_FACE   z = 9,  face 14
        {{6.,0.,1.}, {6.,3.,1.}, {6.,3.,8.}, {6.,0.,8.}},      // I_FACE   x = 6,  face 15
        {{6.,0.,8.}, {6.,3.,8.}, {6.,3.,9.}, {6.,0.,9.}},      // I_FACE   x = 6,  face 16
        {{6.,3.,1.}, {6.,6.,1.}, {6.,6.,8.}, {6.,3.,8.}},      // I_FACE   x = 6,  face 17
        {{6.,3.,8.}, {6.,6.,8.}, {6.,6.,9.}, {6.,3.,9.}},      // I_FACE   x = 6,  face 18
        {{12.,0.,1.}, {12.,3.,1.}, {12.,3.,9.}, {12.,0.,9.}},  // I_FACE   x = 12, face 19
        {{12.,3.,1.}, {12.,6.,1.}, {12.,6.,9.}, {12.,3.,9.}},  // I_FACE   x = 12, face 20
        {{6.,0.,1.}, {12.,0.,1.}, {12.,0.,9.}, {6.,0.,9.}},    // J_FACE   y = 0,  face 21
        {{6.,3.,1.}, {12.,3.,1.}, {12.,3.,9.}, {6.,3.,9.}},    // J_FACE   y = 3,  face 22
        {{6.,6.,1.}, {12.,6.,1.}, {12.,6.,9.}, {6.,6.,9.}},    // J_FACE   y = 6,  face 23
    };

    checkFaces(refinedGrid1, expectedFaces);
    checkFaces(leafGrid, expectedFaces);

    Opm::checkGridWithLgrs(grid,
                           {{1,2,1}}, // cells_per_dim_vec
                           {"LGR1"}); //  lgr_name_vec
}
