/*
  Copyright 2023, 2025 Equinor ASA.

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

#define BOOST_TEST_MODULE CuboidShapeTest
#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/grid/CpGrid.hpp>

struct Fixture
{
    Fixture()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(m_argc, m_argv);
        Opm::OpmLog::setupSimpleDefaultLogging();
    }

    static int rank()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        return Dune::MPIHelper::instance(m_argc, m_argv).rank();
    }
};

BOOST_GLOBAL_FIXTURE(Fixture);

bool checkCuboidShape(const Dune::cpgrid::Entity<0>& element, const Dune::CpGrid& grid)
{
    const auto& levelZeroGrid = grid.currentData().front();

    const auto cellToPoint = levelZeroGrid->cellToPoint(element.index());
    // bottom face corners {0,1,2,3}, top face corners {4,5,6,7}

    // Compute 'cuboid' volume with corners: |corn[1]-corn[0]|x|corn[3]-corn[1]|x|corn[5]-corn[1]|
    std::vector<Dune::cpgrid::Geometry<0,3>::GlobalCoordinate> aFewCorners;
    aFewCorners.resize(4); // {'0', '1', '3', '5'}
    aFewCorners[0] = (*(levelZeroGrid-> getGeometry().geomVector(std::integral_constant<int,3>()))).get(cellToPoint[0]).center();
    aFewCorners[1] = (*(levelZeroGrid-> getGeometry().geomVector(std::integral_constant<int,3>()))).get(cellToPoint[1]).center();
    aFewCorners[2] = (*(levelZeroGrid -> getGeometry().geomVector(std::integral_constant<int,3>()))).get(cellToPoint[3]).center();
    aFewCorners[3] = (*(levelZeroGrid -> getGeometry().geomVector(std::integral_constant<int,3>()))).get(cellToPoint[5]).center();
    //  l = length. b = breadth. h = height.
    double  length, breadth, height;
    length = std::sqrt( ((aFewCorners[1][0] -aFewCorners[0][0])*(aFewCorners[1][0] -aFewCorners[0][0])) +
                        ((aFewCorners[1][1] -aFewCorners[0][1])*(aFewCorners[1][1] -aFewCorners[0][1])) +
                        ((aFewCorners[1][2] -aFewCorners[0][2])*(aFewCorners[1][2] -aFewCorners[0][2])));
    breadth = std::sqrt( ((aFewCorners[1][0] -aFewCorners[2][0])*(aFewCorners[1][0] -aFewCorners[2][0])) +
                         ((aFewCorners[1][1] -aFewCorners[2][1])*(aFewCorners[1][1] -aFewCorners[2][1])) +
                         ((aFewCorners[1][2] -aFewCorners[2][2])*(aFewCorners[1][2] -aFewCorners[2][2])));
    height = std::sqrt( ((aFewCorners[1][0] -aFewCorners[3][0])*(aFewCorners[1][0] -aFewCorners[3][0])) +
                        ((aFewCorners[1][1] -aFewCorners[3][1])*(aFewCorners[1][1] -aFewCorners[3][1])) +
                        ((aFewCorners[1][2] -aFewCorners[3][2])*(aFewCorners[1][2] -aFewCorners[3][2])));
    const double cuboidVolume = length*breadth*height;
    const auto cellVolume =  (*(levelZeroGrid-> getGeometry().geomVector(std::integral_constant<int,0>())))[Dune::cpgrid::EntityRep<0>(element.index(), true)].volume();
    
    return (std::abs(cuboidVolume - cellVolume) <  1e-6);
}

void createEclGridCpGrid_and_checkCuboidShape(const std::string& deckString)
{
    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);

    Dune::CpGrid grid;
    Opm::EclipseGrid eclGrid(deck);

    grid.processEclipseFormat(&eclGrid, nullptr, false, false, false);
    
    for (const auto& element : Dune::elements(grid.leafGridView()))
    {
        BOOST_CHECK(!checkCuboidShape(element, grid));
    }
}

BOOST_AUTO_TEST_CASE(nonCuboidCell)
{/*
   Cell corners:                     COORD
   0 {0, 0, 0}                       line 1: corners 0 and 4
   1 {1, 0, 0}                       line 2: corners 1 and 5
   2 {0, 1, 0}                       line 3: corners 2 and 6
   3 {1, 1, 0}                       line 4: corners 3 and 7
   4 {0, 0, 1}
   5 {1, 0, 0}  coincides with 1
   6 {0, 1, 1}
   7 {1, 1, 0}  coincides with 3 */
    const std::string deckString =
        R"(RUNSPEC
        DIMENS
        1  1  1 /
        GRID
        COORD
        0 0 0  0 0 1
        1 0 0  1 0 0
        0 1 0  0 1 1
        1 1 0  1 1 0
        /
        ZCORN
        6*0
        2*1
        /
        ACTNUM
        1*1
        /
        )";

    createEclGridCpGrid_and_checkCuboidShape(deckString);
}
