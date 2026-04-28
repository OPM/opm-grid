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

#define BOOST_TEST_MODULE FaultFullyInLgrTest
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

BOOST_AUTO_TEST_CASE(parentCellWithMoreThanSixIntersections_I_FACE)
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
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //                 |                      |          |  (6,0,4) **(6,2,4)*(6,4,4)**(6,6,4)               x*****x*****x*****x
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //              (6,0,1) -----------------(6,6,1)     |  (6,0,1) --(?,?,?)-(?,?,?)--(6,6,1)               x- 1 -x- 3 -x- 5 -x
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //              (6,0,0) -----------------(6,6,0)     |  (6,0,0) --(6,2,0)-(6,4,0)--(6,6,0)               x-----x-----x-----x
    //                                                   |
    //                                                   | The missing vertices are (6,2,1) and (6,4,1), appering in elements 1,3, or 5 in LGR1.
    //                                                   | In LGR1 element 1: (6,2,1)
    //                                                   | In LGR1 element 3: (6,2,1) and (6,4,1)
    //                                                   | In LGR1 element 5: (6,4,1)

    const auto& refinedGridData = *grid.currentData()[1];
    const auto& parentGridData = *grid.currentData()[0];
    const auto parentElem = Dune::cpgrid::Entity<0>(parentGridData, 0, true);

    using Coordinate = Dune::FieldVector<double, 3>;

    for (const auto& refinedElem : Dune::elements(grid.levelGridView(1))) {

        const auto collectedVertices = Opm::Lgr::collectNewVertices<Coordinate>(refinedGridData, refinedElem, parentGridData, parentElem);
        std::set<Coordinate,Opm::Lgr::FieldVectorLess> expectedVertices{};

        if (refinedElem.index() ==  1){
            expectedVertices = {{6., 2., 1.}};
            BOOST_CHECK_EQUAL( collectedVertices.size(), 1);
        }
        else if (refinedElem.index() ==  3){
            expectedVertices = {{6., 2., 1.}, {6., 4.,1.}};
            BOOST_CHECK_EQUAL( collectedVertices.size(), 2);
        }
        else if (refinedElem.index() ==  5){
            expectedVertices = {{6., 4.,1.}};
            BOOST_CHECK_EQUAL( collectedVertices.size(), 1);
        }

        for (const auto& expectedVertex : expectedVertices) { // empty if refinedElem.index() != 1, 3, or 5
            auto it = collectedVertices.find(expectedVertex);
            BOOST_CHECK(it != collectedVertices.end());
        }
    }
}

BOOST_AUTO_TEST_CASE(parentCellWithMoreThanSixIntersections_J_FACE)
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
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //                 |                      |          |  (0,6,4) **(2,6,4)*(4,6,4)**(6,6,4)               x*****x*****x*****x
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //              (0,6,1) -----------------(6,6,1)     |  (0,6,1) --(?,?,?)-(?,?,?)--(6,6,1)               x- 3 -x- 4 -x- 5 -x
    //                 |                      |          |     |         *       *        |                  |     *     *     |
    //              (0,6,0) -----------------(6,6,0)     |  (0,6,0) --(2,6,0)-(4,6,0)--(6,6,0)               x-----x-----x-----x
    //                                                   |
    //                                                   | The missing vertices are (2,6,1) and (4,6,1), appering in elements 3,4 or 5 in LGR1.
    //                                                   | In LGR1 element 3: (2,6,1)
    //                                                   | In LGR1 element 4: (2,6,1) and (4,6,1)
    //                                                   | In LGR1 element 5: (4,6,1)

    const auto& refinedGridData = *grid.currentData()[1];
    const auto& parentGridData = *grid.currentData()[0];
    const auto parentElem = Dune::cpgrid::Entity<0>(parentGridData, 0, true);

    using Coordinate = Dune::FieldVector<double, 3>;

    for (const auto& refinedElem : Dune::elements(grid.levelGridView(1))) {

        const auto collectedVertices = Opm::Lgr::collectNewVertices<Coordinate>(refinedGridData, refinedElem, parentGridData, parentElem);
        std::set<Coordinate,Opm::Lgr::FieldVectorLess> expectedVertices{};

        if (refinedElem.index() ==  3){
            expectedVertices = {{2., 6., 1.}};
            BOOST_CHECK_EQUAL( collectedVertices.size(), 1);
        }
        else if (refinedElem.index() == 4){
            expectedVertices = {{2., 6., 1.}, {4., 6., 1.}};
            BOOST_CHECK_EQUAL( collectedVertices.size(), 2);
        }
        else if (refinedElem.index() == 5){
            expectedVertices = {{4., 6.,1.}};
            BOOST_CHECK_EQUAL( collectedVertices.size(), 1);
        }

        for (const auto& expectedVertex : expectedVertices) { // empty if refinedElem.index() != 3,4, or 5
            auto it = collectedVertices.find(expectedVertex);
            BOOST_CHECK(it != collectedVertices.end());
        }
    }
}
