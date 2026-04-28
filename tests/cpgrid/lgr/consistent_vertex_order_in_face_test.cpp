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

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>

#include <tests/cpgrid/lgr/LgrChecks.hpp>

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

void checkConsistentVertexOrderInFaceToPoint(const Dune::cpgrid::CpGridData& gridData,
                                             const Dune::cpgrid::Entity<0>& element)
{
    const auto& cellToFace = gridData.cellToFace(element.index());

    for (const auto& face : cellToFace) {

        const auto& faceToPoint = gridData.faceToPoint(face.index());

        const auto& vertex0 =  Dune::cpgrid::Entity<3>(gridData, faceToPoint[0], true).geometry().center();
        const auto& vertex1 =  Dune::cpgrid::Entity<3>(gridData, faceToPoint[1], true).geometry().center();
        const auto& vertex2 =  Dune::cpgrid::Entity<3>(gridData, faceToPoint[2], true).geometry().center();
        const auto& vertex3 =  Dune::cpgrid::Entity<3>(gridData, faceToPoint[3], true).geometry().center();

        const auto& faceTag = gridData.faceTag(face.index());

        if (faceTag == 0) {
            BOOST_CHECK( vertex0[1] == vertex3[1]);
            BOOST_CHECK( vertex1[1] == vertex2[1]);

            BOOST_CHECK( vertex0[2] == vertex1[2]);
            BOOST_CHECK( vertex2[2] == vertex3[2]);
        }
        else if (faceTag == 1) {
            BOOST_CHECK( vertex0[0] == vertex3[0]);
            BOOST_CHECK( vertex1[0] == vertex2[0]);

            BOOST_CHECK( vertex0[2] == vertex1[2]);
            BOOST_CHECK( vertex2[2] == vertex3[2]);
        }
        else if (faceTag == 2) {
            BOOST_CHECK( vertex0[0] == vertex3[0]);
            BOOST_CHECK( vertex1[0] == vertex2[0]);

            BOOST_CHECK( vertex0[1] == vertex1[1]);
            BOOST_CHECK( vertex2[1] == vertex3[1]);
        }
    }
}

BOOST_AUTO_TEST_CASE(consistentVertexOrderInFaceToPointBeforeAndAfterRefinement)
{
    // Level zero grid dims = 1x1x1
    //
    // cell 0
    // bottom face corners (0,0,0), (2,0,0), (0,2,0), (2,2,0)
    //    top face corners (0,0,2), (2,0,2), (0,2,2), (2,2,2)
    const std::string deckString =
        R"(
RUNSPEC
DIMENS
 1 1 1 /

GRID

COORD
 0 0 0     0 0 2
 2 0 0     2 0 2

 0 2 0     0 2 2
 2 2 0     2 2 2
/

ZCORN
 0 0 0 0  2 2 2 2
/

ACTNUM
 1
/

PORO
0.15
/
)";

    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckString,
                              /* cells_per_dim_vec */ {{2,2,2}},
                              /* startIJK_vec */ {{0,0,0}},
                              /* endIJK_vec */ {{1,1,1}},
                              /* lgr_name_vec */ {"LGR1"});

    for (const auto& element : Dune::elements(grid.levelGridView(0))) {
        checkConsistentVertexOrderInFaceToPoint(*grid.currentData().front(), element);
    }
    
    for (const auto& element : Dune::elements(grid.levelGridView(1))) {
        checkConsistentVertexOrderInFaceToPoint(*grid.currentData()[1], element);
    }
    
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        checkConsistentVertexOrderInFaceToPoint(grid.currentLeafData(), element);
    }    
}

