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

#define BOOST_TEST_MODULE LgrWithFaultsTests
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

// Test that corner and face relationships are correctly updated when
// LGRs are refined in the presence of faults.
//
// The tests use three basic deck configurations: 2x1x1, 1x2x1, and 2x1x2.
// The grid cells are shifted such that every cell in the level-zero grid
// has seven faces. Normally, a cell has one face of each type
// (I-, I+, J-, J+, K-, K+), but faults can cause a cell to have multiple
// faces of the same type.
//
// These grids are refined in different ways to cover several cases,
// including refining both cells, refining only one cell, an LGR containing
// a fault, and multiple LGRs whose boundaries coincide with faults.

using Coordinate = Dune::FieldVector<double, 3>;

// In all the test cases, which are quite simple grids of dimensions
// 2x1x1, 1x2x1, or 2x1x2, the grid cells are shifted such that
// every cell in the level-zero grid has seven faces.
void checkFaultInLevelZeroGrid(const Dune::CpGrid& grid)
{
    const auto& parentGridData = *grid.currentData()[0];
    for (const auto& element : Dune::elements(grid.levelGridView(0))) {
        BOOST_CHECK_EQUAL(parentGridData.cellToFace(element.index()).size(), 7);
    }
}

// In a without-faults-grid of dimensions {nx, ny, nz},
// - the total amount of vertices is given by (nx+1)*(ny+1)*(nz+1)
// - the total amount of faces is given by ((nx+1)*ny*nz) + (nx*(ny+1)*nz) + (nx*ny*(nz+1))
// Since now the refinement is aware of
// (a) parent cell having more than one face per type
//     (e.g. level zero grid having cell 0 with two I+ faces and cell 1 with two I- faces)
// (b) neighboring refinements with faults (creating new vertices and faces)
// the total amount of vertices and faces in presence of faults may be larger than
// (nx+1)*(ny+1)*(nz+1) and ((nx+1)*ny*nz) + (nx*(ny+1)*nz) + (nx*ny*(nz+1)) respectively.
void checkVertexAndFaceCount(const Dune::cpgrid::CpGridData& gridData,
                             const std::array<int,3>& nxnynz,
                             int newVertexCount,
                             int vanishedFaceCount,
                             int newFaceCount)
{
    const auto& [nx,ny,nz] = nxnynz;
    BOOST_CHECK_EQUAL( gridData.size(3), (nx+1)*(ny+1)*(nz+1) + newVertexCount);
    BOOST_CHECK_EQUAL( gridData.numFaces(), ((nx+1)*ny*nz) + (nx*(ny+1)*nz) + (nx*ny*(nz+1)) - vanishedFaceCount + newFaceCount); 
}

void checkVertexAndFaceLeafCount(const Dune::CpGrid& grid,
                                 int sharedBetweenLevelsVertexCount,
                                 int sharedBetweenLevelsFaceCount)
{   int vertexCount = 0;
    int faceCount = 0;
    for (int level = 0; level <= grid.maxLevel(); ++level) {
        vertexCount += grid.currentData()[level]->size(3);
        faceCount += grid.currentData()[level]->numFaces();
    }
    
    BOOST_CHECK_EQUAL( grid.size(3), vertexCount - sharedBetweenLevelsVertexCount);    
    BOOST_CHECK_EQUAL( grid.numFaces(), faceCount - sharedBetweenLevelsFaceCount); 
}

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
                const std::vector<std::vector<Coordinate>>& expectedFaces,
                const std::vector<std::pair<int,int>>& expectedFaceToCell)
{
    std::vector<std::vector<Dune::FieldVector<double,3>>> actualFaces{};
    actualFaces.resize(gridData.numFaces());

    for (int i = 0; i < gridData.numFaces(); ++i) {
        
        const auto& faceToCell = gridData.faceToCell(i);
        BOOST_CHECK(faceToCell.size() <= 2);
        
        const auto& [cell1, cell2] = expectedFaceToCell[i];
        int expectedFaceToCellSize = (cell2 == -1)? 1 : 2;
        
        if (expectedFaceToCellSize>1) {
            BOOST_CHECK_EQUAL( faceToCell[0].index(), cell1);
            BOOST_CHECK_EQUAL( faceToCell[1].index(), cell2);
        }
        else {
            BOOST_CHECK_EQUAL( faceToCell[0].index(), cell1);
        }
        
        
        for (const auto& vertexIdx : gridData.faceToPoint(i)){
            actualFaces[i].push_back(Dune::cpgrid::Entity<3>(gridData, vertexIdx, true).geometry().center());
        }
    }

    BOOST_CHECK(equalFaces(actualFaces, expectedFaces));
}

void checkCellToFace(const Dune::CpGrid& grid, int level,
                     const std::vector<std::vector<std::vector<Coordinate>>>& expectedCellToFace)
{
    bool isLeafGrid = (level<0);
    int cell_count = isLeafGrid? grid.leafGridView().size(0) : grid.levelGridView(level).size(0);
    
    BOOST_CHECK_EQUAL(cell_count, expectedCellToFace.size());
    
    const auto& elements = isLeafGrid? Dune::elements(grid.leafGridView()) : Dune::elements(grid.levelGridView(level));
    const auto& gridData = isLeafGrid? grid.currentLeafData() : *grid.currentData()[level];
    
    for (const auto& element : elements) {
        const auto& cellToFace = gridData.cellToFace(element.index());
        
        std::vector<std::vector<Dune::FieldVector<double,3>>> actualCellToFace{};
        actualCellToFace.resize(cellToFace.size());

        BOOST_CHECK_EQUAL( cellToFace.size(), expectedCellToFace[element.index()].size());

        int count = 0;
        for (const auto& face : cellToFace) {
            for (const auto& vertexIdx : gridData.faceToPoint(face.index())){
                actualCellToFace[count].push_back(Dune::cpgrid::Entity<3>(gridData, vertexIdx, true).geometry().center());
            }
            ++count;
        }

        BOOST_CHECK(equalFaces(actualCellToFace, expectedCellToFace[element.index()]));
    }
}

// Level zero grid dims = 2x1x1
//
// cell 0
// bottom face corners (0,0,0), (6,0,0), (0,6,0), (6,6,0)
//    top face corners (0,0,8), (6,0,8), (0,6,8), (6,6,8)
//
// cell 1
// bottom face corners (6,0,1), (12,0,1), (6,6,1),  (12,6,1)
//    top face corners (6,0,9), (12,0,9), (12,6,9), (12,6,9)
const std::string deckTwoCellsInXDirGrid =
    R"(RUNSPEC
DIMENS
 2 1 1 /

GRID

COORD
 0 0 0     0 0 9
 6 0 0     6 0 9
12 0 0    12 0 9

 0 6 0     0 6 9
 6 6 0     6 6 9
12 6 0    12 6 9
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

// In the test cases "simpleRefinementDiffLgrs" and "simpleRefinementSameLgr"
// each parent cell gets refined into 1x2x1 children. In this way, it's pretty
// easy and not so long to provide the expected cell-to-face and face-to-cell
// relationships, as well as the vertices coordinates and the face vertices
// coordinates for level and leaf grids.
// The leaf grid view in both test cases coincide, therefore we introduce below
// a few "expected-leaf-grid" data that will be used for both test cases.

// Element 0 and element 1 in level zero grid share an I_FACE (with face index 2)
//
// Vertices of those faces lie on the plane x = 6    | After refinement, number of subdivisions in
//                                                   | y- and z- directions:
//
//              (6,0,9) -----------------(6,6,9)     |  (6,0,9) --(6,3,9)---(6,6,9)
//                 |      face idx 3      |          |     |         *         |
//              (6,0,8) ---------------- (6,6,8)     |  (6,0,8) --(6,3,8)---(6,6,8)
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

// "simpleRefinementDiffLgrs"
//
// Recall that for a level zero input grid of dimensions 2x1x1 that gets refined
// LGR1 and LGR2, with 1x2x1 children per cell,
//
// k = 0 |  cell 0 | cell 1 |  level zero grid
//
// Parent cell from level zero refined into 1x2x2 children:
// k = 0, j = 1  | cell 1 (LGR*) |    level 1 and 2 grid cell indices are equal. 
//               |---------------|
//        j = 0  | cell 0 (LGR*) |
//
// k = 0, j = 1  | cell 1 (from LGR1) | cell 3 (from LGR2)|    leaf grid cell indices
//               |--------------------|-------------------|
//        j = 0  | cell 0 (from LGR1) | cell 2 (from LGR2)|
//
//
//                LGR1                       | LGR2
// New vertices: (6,0,1),(6,3,1),(6,6,1)     | (6,0,8),(6,3,8),(6,6,8)
//                                           |
// New faces:                                |
//              (6,0,8) --(6,3,8)---(6,6,8)  |  (6,0,9) --(6,3,9)---(6,6,9)
//                 |         *         |     |     |         *         |  
//                 |         *         |     |  (6,0,8) --(6,3,8)---(6,6,8) 
//                 |         *         |     |     |         *         | 
//              (6,0,1) --(6,3,1)---(6,6,1)  |     |         *         | 
//                 |         *         |     |     |         *         | 
//              (6,0,0) --(6,3,0)---(6,6,0)  |  (6,0,1) --(6,3,1)---(6,6,1)

// "simpleRefinementSameLgr"
//
// Recall that for a level zero input grid of dimensions 2x1x1 that gets refined
// into one LGR1 with 1x2x1 children per cell,
//
// k = 0 |  cell 0 | cell 1 |  level zero grid
//
// Entire level zero grid gets refined, each parent cell has 1x2x1 children:
//
// k = 0, j = 1  | cell 1  | cell 3 |    level 1 and leaf grid cell indices
//               |---------|--------|
//        j = 0  | cell 0  | cell 2 |


// To avoid creating containers with the same content, here we provide a few
// cell-to-face and face-to-cell relationships that appear in both test cases.

const std::vector<std::vector<Coordinate>> expectedLeafCellToFace0 = {
    {{0.,0.,0.}, {6.,0.,0.}, {6.,3.,0.}, {0.,3.,0.}},      // K_FACE   z = 0, face 0
    {{0.,0.,8.}, {6.,0.,8.}, {6.,3.,8.}, {0.,3.,8.}},      // K_FACE   z = 8, face 1
    {{0.,0.,0.}, {0.,3.,0.}, {0.,3.,8.}, {0.,0.,8.}},      // I_FACE   x = 0, face 2
    {{6.,0.,0.}, {6.,3.,0.}, {6.,3.,1.}, {6.,0.,1.}},      // I_FACE   x = 6, face 3
    {{6.,0.,1.}, {6.,3.,1.}, {6.,3.,8.}, {6.,0.,8.}},      // I_FACE   x = 6, face 4
    {{0.,0.,0.}, {6.,0.,0.}, {6.,0.,8.}, {0.,0.,8.}},      // J_FACE   y = 0, face 5
    {{0.,3.,0.}, {6.,3.,0.}, {6.,3.,8.}, {0.,3.,8.}}       // J_FACE   y = 3, face 6
};
    
const std::vector<std::vector<Coordinate>> expectedLeafCellToFace1 = {
    {{0.,3.,0.}, {6.,3.,0.}, {6.,6.,0.}, {0.,6.,0.}},      // K_FACE   z = 0, face 0
    {{0.,3.,8.}, {6.,3.,8.}, {6.,6.,8.}, {0.,6.,8.}},      // K_FACE   z = 8, face 1
    {{0.,3.,0.}, {0.,6.,0.}, {0.,6.,8.}, {0.,3.,8.}},      // I_FACE   x = 0, face 2
    {{6.,3.,0.}, {6.,6.,0.}, {6.,6.,1.}, {6.,3.,1.}},      // I_FACE   x = 6, face 3
    {{6.,3.,1.}, {6.,6.,1.}, {6.,6.,8.}, {6.,3.,8.}},      // I_FACE   x = 6, face 4
    {{0.,3.,0.}, {6.,3.,0.}, {6.,3.,8.}, {0.,3.,8.}},      // J_FACE   y = 3, face 5
    {{0.,6.,0.}, {6.,6.,0.}, {6.,6.,8.}, {0.,6.,8.}}       // J_FACE   y = 6, face 6
};

const std::vector<std::vector<Coordinate>> expectedLeafCellToFace2 = {
    {{6.,0.,1.}, {12.,0.,1.}, {12.,3.,1.}, {6.,3.,1.}},    // K_FACE   z = 1,  face 0
    {{6.,0.,9.}, {12.,0.,9.}, {12.,3.,9.}, {6.,3.,9.}},    // K_FACE   z = 9,  face 1
    {{6.,0.,1.}, {6.,3.,1.}, {6.,3.,8.}, {6.,0.,8.}},      // I_FACE   x = 6,  face 2
    {{6.,0.,8.}, {6.,3.,8.}, {6.,3.,9.}, {6.,0.,9.}},      // I_FACE   x = 6,  face 3
    {{12.,0.,1.}, {12.,3.,1.}, {12.,3.,9.}, {12.,0.,9.}},  // I_FACE   x = 12, face 4
    {{6.,0.,1.}, {12.,0.,1.}, {12.,0.,9.}, {6.,0.,9.}},    // J_FACE   y = 0,  face 5
    {{6.,3.,1.}, {12.,3.,1.}, {12.,3.,9.}, {6.,3.,9.}},    // J_FACE   y = 3,  face 6
};
    
const std::vector<std::vector<Coordinate>> expectedLeafCellToFace3 = {
    {{6.,3.,1.}, {12.,3.,1.}, {12.,6.,1.}, {6.,6.,1.}},    // K_FACE   z = 1,  face 0
    {{6.,3.,9.}, {12.,3.,9.}, {12.,6.,9.}, {6.,6.,9.}},    // K_FACE   z = 9,  face 1
    {{6.,3.,1.}, {6.,6.,1.}, {6.,6.,8.}, {6.,3.,8.}},      // I_FACE   x = 6,  face 2
    {{6.,3.,8.}, {6.,6.,8.}, {6.,6.,9.}, {6.,3.,9.}},      // I_FACE   x = 6,  face 3
    {{12.,3.,1.}, {12.,6.,1.}, {12.,6.,9.}, {12.,3.,9.}},  // I_FACE   x = 12, face 4
    {{6.,3.,1.}, {12.,3.,1.}, {12.,3.,9.}, {6.,3.,9.}},    // J_FACE   y = 3,  face 5
    {{6.,6.,1.}, {12.,6.,1.}, {12.,6.,9.}, {6.,6.,9.}},    // J_FACE   y = 6,  face 6
};

// Leaf faces whose faceToCell has 2 elements
//
//     (6,0,8) --(6,3,9)---(6,6,8)
//        |         *         |
//        |         *         |
//        |         *         |         leaf face index 15 to cells = {0, 2}
//        |  face   *  face   |         leaf face index 17 to cells = {1, 3}
//        |   idx   *   idx   |
//        |   15    *   17    |
//        |         *         |
//        |         *         |
//     (6,0,1) --(6,3,1)---(6,6,1)

const std::vector<std::vector<Coordinate>> expectedLeafFaces = {
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

const std::vector<std::pair<int,int>> expectedLeafFaceToCell = {
    {0, -1},      // K_FACE   z = 0,  face 0  {{0.,0.,0.}, {6.,0.,0.}, {6.,3.,0.}, {0.,3.,0.}}
    {1, -1},      // K_FACE   z = 0,  face 1  {{0.,3.,0.}, {6.,3.,0.}, {6.,6.,0.}, {0.,6.,0.}}
    {0, -1},      // K_FACE   z = 8,  face 2  {{0.,0.,8.}, {6.,0.,8.}, {6.,3.,8.}, {0.,3.,8.}} 
    {1, -1},      // K_FACE   z = 8,  face 3  {{0.,3.,8.}, {6.,3.,8.}, {6.,6.,8.}, {0.,6.,8.}}
    {0, -1},      // I_FACE   x = 0,  face 4  {{0.,0.,0.}, {0.,3.,0.}, {0.,3.,8.}, {0.,0.,8.}}
    {1, -1},      // I_FACE   x = 0,  face 5  {{0.,3.,0.}, {0.,6.,0.}, {0.,6.,8.}, {0.,3.,8.}}
    {0, -1},      // I_FACE   x = 6,  face 6  {{6.,0.,0.}, {6.,3.,0.}, {6.,3.,1.}, {6.,0.,1.}}
    {1, -1},      // I_FACE   x = 6,  face 7  {{6.,3.,0.}, {6.,6.,0.}, {6.,6.,1.}, {6.,3.,1.}}
    {0, -1},      // J_FACE   y = 0,  face 8  {{0.,0.,0.}, {6.,0.,0.}, {6.,0.,8.}, {0.,0.,8.}}
    {0,  1},      // J_FACE   y = 3,  face 9  {{0.,3.,0.}, {6.,3.,0.}, {6.,3.,8.}, {0.,3.,8.}}
    {1, -1},      // J_FACE   y = 6,  face 10 {{0.,6.,0.}, {6.,6.,0.}, {6.,6.,8.}, {0.,6.,8.}}
    {2, -1},      // K_FACE   z = 1,  face 11 {{6.,0.,1.}, {12.,0.,1.}, {12.,3.,1.}, {6.,3.,1.}}
    {3, -1},      // K_FACE   z = 1,  face 12 {{6.,3.,1.}, {12.,3.,1.}, {12.,6.,1.}, {6.,6.,1.}}
    {2, -1},      // K_FACE   z = 9,  face 13 {{6.,0.,9.}, {12.,0.,9.}, {12.,3.,9.}, {6.,3.,9.}}
    {3, -1},      // K_FACE   z = 9,  face 14 {{6.,3.,9.}, {12.,3.,9.}, {12.,6.,9.}, {6.,6.,9.}}
    {0,  2},      // I_FACE   x = 6,  face 15 {{6.,0.,1.}, {6.,3.,1.}, {6.,3.,8.}, {6.,0.,8.}}
    {2, -1},      // I_FACE   x = 6,  face 16 {{6.,0.,8.}, {6.,3.,8.}, {6.,3.,9.}, {6.,0.,9.}}
    {1,  3},      // I_FACE   x = 6,  face 17 {{6.,3.,1.}, {6.,6.,1.}, {6.,6.,8.}, {6.,3.,8.}}
    {3, -1},      // I_FACE   x = 6,  face 18 {{6.,3.,8.}, {6.,6.,8.}, {6.,6.,9.}, {6.,3.,9.}}
    {2, -1},      // I_FACE   x = 12, face 19 {{12.,0.,1.}, {12.,3.,1.}, {12.,3.,9.}, {12.,0.,9.}}
    {3, -1},      // I_FACE   x = 12, face 20 {{12.,3.,1.}, {12.,6.,1.}, {12.,6.,9.}, {12.,3.,9.}}
    {2, -1},      // J_FACE   y = 0,  face 21 {{6.,0.,1.}, {12.,0.,1.}, {12.,0.,9.}, {6.,0.,9.}}
    {2,  3},      // J_FACE   y = 3,  face 22 {{6.,3.,1.}, {12.,3.,1.}, {12.,3.,9.}, {6.,3.,9.}}
    {3, -1},      // J_FACE   y = 6,  face 23 {{6.,6.,1.}, {12.,6.,1.}, {12.,6.,9.}, {6.,6.,9.}}
};

BOOST_AUTO_TEST_CASE(simpleRefinementDiffLgrs)
{
    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckTwoCellsInXDirGrid,
                              /* cells_per_dim_vec */ {{1,2,1}, {1,2,1}},
                              /* startIJK_vec */      {{0,0,0}, {1,0,0}},
                              /* endIJK_vec */        {{1,1,1}, {2,1,1}},
                              /* lgr_name_vec */      {"LGR1", "LGR2"});

    checkFaultInLevelZeroGrid(grid);

    // In a grid of dimensions {nx, ny, nz},
    // - the total amount of vertices is given by (nx+1)*(ny+1)*(nz+1)
    // - the total amount of faces is given by ((nx+1)*ny*nz) + (nx*(ny+1)*nz) + (nx*ny*(nz+1))
    // Since now the refinement is aware of
    // (a) parent cell having more than one face per type
    //     (like in this level zero grid having cell 0 with two I+ faces and cell 1 with two I- faces)
    // (b) neighboring refinements with faults (creating new vertices and faces)
    // the total amount of vertices and faces in presence of faults may be larger than
    // (nx+1)*(ny+1)*(nz+1) and ((nx+1)*ny*nz) + (nx*(ny+1)*nz) + (nx*ny*(nz+1)) respectively.

    const auto& refinedGrid1 = *grid.currentData()[1];
    const auto& refinedGrid2 = *grid.currentData()[2];
    // for LGR1 and LGR2, {nx,ny,nz} = {1,2,1}
    // (nx+1)*(ny+1)*(nz+1) + "new vertices" = (2x3x2) + 3 = 15
    // ((nx+1)*ny*nz) + (nx*(ny+1)*nz) + (nx*ny*(nz+1)) - vanished faces + new faces = (2x2x1) + (1x3x1) + (1x2x2) - 2 + 4 = 13
    checkVertexAndFaceCount(refinedGrid1, /* nxnynz = */ {1,2,1}, /* newVertexCount = */ 3,
                            /* vanishedFaceCount = */ 2, /* newFaceCount = */ 4);
    checkVertexAndFaceCount(refinedGrid2, /* nxnynz = */ {1,2,1}, /* newVertexCount = */ 3,
                            /* vanishedFaceCount = */ 2, /* newFaceCount = */ 4);
    
    const auto& leafGrid = grid.currentLeafData();
    // level 0 + level 1 + level 2 vertices - shared vertices =
    // 16 + 15 + 15 - 10 (shared level 0 and 1) - 10 (shared level 0 and 2) - 2 (shared level 1 and 2, not level 0)  = 24
    // level 0 + level 1 + level 2 faces - shared/vanished faces = 13 + 13 + 13 - 13 (vanished level 0) - 2 (shared level 1 and level 2)  = 24
    checkVertexAndFaceLeafCount(grid, /* sharedBetweenLevelsVertexCount = */ 22, /* sharedBetweenLevelsFaceCount = */ 15);

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

    checkCellToFace(grid, /* level = */ 1, std::vector{expectedLeafCellToFace0,  expectedLeafCellToFace1});
    checkCellToFace(grid, /* level = */ 2, std::vector{expectedLeafCellToFace2,  expectedLeafCellToFace3});
    checkCellToFace(grid, /* level = */ -1 /* leafGrid!*/,
                    std::vector{expectedLeafCellToFace0,  expectedLeafCellToFace1,
                                expectedLeafCellToFace2,  expectedLeafCellToFace3});

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
    const std::vector<std::pair<int,int>> expectedLGR1FaceToCell = {
        {0, -1}, // K_FACE   z = 0, face 0 {{0.,0.,0.}, {6.,0.,0.}, {6.,3.,0.}, {0.,3.,0.}},      
        {1, -1}, // K_FACE   z = 0, face 1 {{0.,3.,0.}, {6.,3.,0.}, {6.,6.,0.}, {0.,6.,0.}},      
        {0, -1}, // K_FACE   z = 8, face 2 {{0.,0.,8.}, {6.,0.,8.}, {6.,3.,8.}, {0.,3.,8.}},      
        {1, -1}, // K_FACE   z = 8, face 3 {{0.,3.,8.}, {6.,3.,8.}, {6.,6.,8.}, {0.,6.,8.}},     
        {0, -1}, // I_FACE   x = 0, face 4 {{0.,0.,0.}, {0.,3.,0.}, {0.,3.,8.}, {0.,0.,8.}},      
        {1, -1}, // I_FACE   x = 0, face 5 {{0.,3.,0.}, {0.,6.,0.}, {0.,6.,8.}, {0.,3.,8.}},    
        {0, -1}, // I_FACE   x = 6, face 6  {{6.,0.,0.}, {6.,3.,0.}, {6.,3.,1.}, {6.,0.,1.}},     
        {0, -1}, // I_FACE   x = 6, face 7 {{6.,0.,1.}, {6.,3.,1.}, {6.,3.,8.}, {6.,0.,8.}},     
        {1, -1}, // I_FACE   x = 6, face 8 {{6.,3.,0.}, {6.,6.,0.}, {6.,6.,1.}, {6.,3.,1.}},      
        {1, -1}, // I_FACE   x = 6, face 9 {{6.,3.,1.}, {6.,6.,1.}, {6.,6.,8.}, {6.,3.,8.}},     
        {0, -1}, // J_FACE   y = 0, face 10 {{0.,0.,0.}, {6.,0.,0.}, {6.,0.,8.}, {0.,0.,8.}},     
        {0,  1}, // J_FACE   y = 3, face 11  {{0.,3.,0.}, {6.,3.,0.}, {6.,3.,8.}, {0.,3.,8.}},     
        {1, -1}  // J_FACE   y = 6, face 12 {{0.,6.,0.}, {6.,6.,0.}, {6.,6.,8.}, {0.,6.,8.}},      
    };
    checkFaces(refinedGrid1, expectedLGR1Faces, expectedLGR1FaceToCell);
 
 
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
    const std::vector<std::pair<int,int>> expectedLGR2FaceToCell = {
        {0, -1}, // K_FACE   z = 1,  face 0  {{6.,0.,1.}, {12.,0.,1.}, {12.,3.,1.}, {6.,3.,1.}},   
        {1, -1}, // K_FACE   z = 1,  face 1  {{6.,3.,1.}, {12.,3.,1.}, {12.,6.,1.}, {6.,6.,1.}},    
        {0, -1}, // K_FACE   z = 9,  face 2  {{6.,0.,9.}, {12.,0.,9.}, {12.,3.,9.}, {6.,3.,9.}},   
        {1, -1}, // K_FACE   z = 9,  face 3  {{6.,3.,9.}, {12.,3.,9.}, {12.,6.,9.}, {6.,6.,9.}},   
        {0, -1}, // I_FACE   x = 6,  face 4  {{6.,0.,1.}, {6.,3.,1.}, {6.,3.,8.}, {6.,0.,8.}},      
        {0, -1}, // I_FACE   x = 6,  face 5  {{6.,0.,8.}, {6.,3.,8.}, {6.,3.,9.}, {6.,0.,9.}},     
        {1, -1}, // I_FACE   x = 6,  face 6  {{6.,3.,1.}, {6.,6.,1.}, {6.,6.,8.}, {6.,3.,8.}},      
        {1, -1}, // I_FACE   x = 6,  face 7  {{6.,3.,8.}, {6.,6.,8.}, {6.,6.,9.}, {6.,3.,9.}},     
        {0, -1}, // I_FACE   x = 12, face 8  {{12.,0.,1.}, {12.,3.,1.}, {12.,3.,9.}, {12.,0.,9.}}, 
        {1, -1}, // I_FACE   x = 12, face 9  {{12.,3.,1.}, {12.,6.,1.}, {12.,6.,9.}, {12.,3.,9.}}, 
        {0, -1}, // J_FACE   y = 0,  face 10 {{6.,0.,1.}, {12.,0.,1.}, {12.,0.,9.}, {6.,0.,9.}},  
        {0,  1}, // J_FACE   y = 3,  face 11 {{6.,3.,1.}, {12.,3.,1.}, {12.,3.,9.}, {6.,3.,9.}},    
        {1, -1}  // J_FACE   y = 6,  face 12 {{6.,6.,1.}, {12.,6.,1.}, {12.,6.,9.}, {6.,6.,9.}},    
    };
    checkFaces(refinedGrid2, expectedLGR2Faces, expectedLGR2FaceToCell);
    
    checkFaces(leafGrid, expectedLeafFaces, expectedLeafFaceToCell);

    Opm::checkGridWithLgrs(grid,
                           {{1,2,1}, {1,2,1}}, // cells_per_dim_vec
                           {"LGR1", "LGR2"});  // lgr_name_vec
}

BOOST_AUTO_TEST_CASE(simpleRefinementSameLgr)
{
    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckTwoCellsInXDirGrid,
                              /* cells_per_dim_vec */ {{1,2,1}},
                              /* startIJK_vec */      {{0,0,0}},
                              /* endIJK_vec */        {{2,1,1}},
                              /* lgr_name_vec */      {"LGR1"});

    const auto& refinedGrid1 = *grid.currentData()[1];
    // for LGR1, {nx,ny,nz} = {2,2,1}
    // (nx+1)*(ny+1)*(nz+1) + "new vertices" = (3x3x2) + 6 = 24
    // ((nx+1)*ny*nz) + (nx*(ny+1)*nz) + (nx*ny*(nz+1)) - vanished faces + new faces = (3x2x1) + (2x3x1) + (2x2x1) - 2 + 6 = 24
    checkVertexAndFaceCount(refinedGrid1, /* nxnynz = */ {2,2,1}, /* newVertexCount = */ 6,
                            /* vanishedFaceCount = */ 2, /* newFaceCount = */ 6);
    
    const auto& leafGrid = grid.currentLeafData();
    // level 0 + level 1 vertices - shared vertices = 16 + 24 - 16 (shared level 0 and 1) = 24
    // level 0 + level 1 faces - shared/vanished faces = 13 + 24 - 13 (vanished level 0) = 24
    checkVertexAndFaceLeafCount(grid, /* sharedBetweenLevelsVertexCount = */ 16, /* sharedBetweenLevelsFaceCount = */ 13);

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

    checkCellToFace(grid, /* level = */ 1, std::vector{expectedLeafCellToFace0,  expectedLeafCellToFace1,
                                                       expectedLeafCellToFace2,  expectedLeafCellToFace3});
    checkCellToFace(grid, /* level = */ -1 /* leafGrid!*/,
                    std::vector{expectedLeafCellToFace0,  expectedLeafCellToFace1,
                                expectedLeafCellToFace2,  expectedLeafCellToFace3});
    
    checkFaces(refinedGrid1, expectedLeafFaces, expectedLeafFaceToCell);
    checkFaces(leafGrid, expectedLeafFaces, expectedLeafFaceToCell);

    Opm::checkGridWithLgrs(grid,
                           {{1,2,1}}, // cells_per_dim_vec
                           {"LGR1"}); //  lgr_name_vec
}

BOOST_AUTO_TEST_CASE(parentCellWithMoreThanOne_I_FACE_trueOriented_nonTrivialOverlap)
{
    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckTwoCellsInXDirGrid,
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

    const auto& refinedGrid = *grid.currentData()[1];
    // for LGR1, {nx,ny,nz} = {2,3,2}
    // (nx+1)*(ny+1)*(nz+1) + "new vertices" = (3x4x3) + 4 = 36 + 4 = 40
    // new vertices  (6,0,1), (6,2,1), (6,4,1), and (6,6,1).
    // ((nx+1)*ny*nz) + (nx*(ny+1)*nz) + (nx*ny*(nz+1)) - vanished faces + new faces = (3x3x2) + (2x4x2) + (2x3x3) - 3 + 6 = 52 - 3 + 6 = 55
    // LGR1 dims 2x3x2 -> 52 faces (before correction due to missing points)
    // 3 of those 52 faces vanished and give origin to 6 new faces: 52 - 3 + 6 = 55 faces
    checkVertexAndFaceCount(refinedGrid, /* nxnynz = */ {2,3,2}, /* newVertexCount = */ 4,
                            /* vanishedFaceCount = */ 3, /* newFaceCount = */ 6);
   
    // level 0 + level 1 vertices - shared vertices = 16 + 40 - 10 (shared level 0 and 1) = 46
    // level 0 + level 1 faces - shared/vanished faces = 13 + 55 - 7 (vanished level 0) = 61
    checkVertexAndFaceLeafCount(grid, /* sharedBetweenLevelsVertexCount = */ 10, /* sharedBetweenLevelsFaceCount = */ 7);


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
    checkNewRefinedFaces(grid, refinedGrid,
                         selectedFaceToCoord, /* repeatedFaceType = */ 1); // 1-> I+

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{2,3,2}},
                           /* lgr_name_vec = */ {"LGR1"});
}

BOOST_AUTO_TEST_CASE(parentCellWithMoreThanOne_I_FACE_trueOriented_trivialOverlap)
{
    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckTwoCellsInXDirGrid,
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
    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckTwoCellsInXDirGrid,
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

BOOST_AUTO_TEST_CASE(neighboringSingleCellRefinementsDifferentLgrs)
{
    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckTwoCellsInXDirGrid,
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
    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckTwoCellsInXDirGrid,
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


// Level zero grid dims = 1x2x1
//
// cell 0
// bottom face corners (0,0,0), (6,0,0), (0,6,0), (6,6,0)
//    top face corners (0,0,8), (6,0,8), (0,6,8), (6,6,8)
//
// cell 1
// bottom face corners (0,6,1), (6,6,1), (0,12,1), (6,12,1)
//    top face corners (0,6,9), (6,6,9), (0,12,9), (6,12,9)

const std::string deckTwoCellsInYDirGrid =
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

BOOST_AUTO_TEST_CASE(parentCellWithMoreThanSixIntersections_J_FACE_true)
{
    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckTwoCellsInYDirGrid,
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
    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckTwoCellsInYDirGrid,
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

// Level zero grid dims = 2x2x1
//
// Cell 0 (x=0..6,   y=0..6)
// bottom face corners (0,0,0), (6,0,0), (0,6,1), (6,6,1)
//    top face corners (0,0,8), (6,0,8), (0,6,9), (6,6,9)
//
// Cell 1 (x=6..12,  y=0..6)
// bottom face corners (6,0,1), (12,0,1), (6,6,2), (12,6,2)
//    top face corners (6,0,9), (12,0,9), (6,6,10), (12,6,10)
//
// Cell 2 (x=0..6,   y=6..12)
// bottom face corners (0,6,0), (6,6,0), (0,12,1), (6,12,1)
//    top face corners (0,6,8), (6,6,8), (0,12,9), (6,12,9)
//
// Cell 3 (x=6..12,  y=6..12)
// bottom face corners (6,6,1), (12,6,1), (6,12,2), (12,12,2)
//    top face corners (6,6,9), (12,6,9), (6,12,10), (12,12,10)

const std::string deckTwoColumnsGrid =
    R"(RUNSPEC
DIMENS
 2 1 2 /

GRID

COORD
 0  0 0     0  0 17
 6  0 0     6  0 17
12  0 0    12  0 17

 0  6 0     0  6 17
 6  6 0     6  6 17
12  6 0    12  6 17
/

ZCORN
0 0 1 1  0 0 1 1
8 8 9 9  8 8 9 9

 8  8  9  9     8  8  9  9
16 16 17 17    16 16 17 17
/

ACTNUM
4*1
/

PORO
4*0.15
/
)";

BOOST_AUTO_TEST_CASE(faultBetweenTwoColumnsSameLgr)
{
    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckTwoColumnsGrid,
                              /* cells_per_dim_vec */{{1,2,1}},
                              /* startIJK_vec */      {{0,0,0}},
                              /* endIJK_vec */        {{2,1,2}},
                              /* lgr_name_vec */      {"LGR1"});

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{1,2,1}},
                           /* lgr_name_vec = */ {"LGR1"});

}

BOOST_AUTO_TEST_CASE(faultBetweenTwoColumnsDifferentVerticalLgrs)
{
    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckTwoColumnsGrid,
                              /* cells_per_dim_vec */{{1,2,1}, {1,2,1}},
                              /* startIJK_vec */      {{0,0,0}, {1,0,0}},
                              /* endIJK_vec */        {{1,1,2}, {2,1,2}},
                              /* lgr_name_vec */      {"LGR1", "LGR2"});

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{1,2,1}, {1,2,1}},
                           /* lgr_name_vec = */ {"LGR1", "LGR2"});

}

BOOST_AUTO_TEST_CASE(faultBetweenTwoColumnsDifferentHorizontalLgrs)
{
    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckTwoColumnsGrid,
                              /* cells_per_dim_vec */{{1,2,1}, {1,2,1}},
                              /* startIJK_vec */      {{0,0,0}, {0,0,1}},
                              /* endIJK_vec */        {{2,1,1}, {2,1,2}},
                              /* lgr_name_vec */      {"LGR1", "LGR2"});

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{1,2,1}, {1,2,1}},
                           /* lgr_name_vec = */ {"LGR1", "LGR2"});

}
