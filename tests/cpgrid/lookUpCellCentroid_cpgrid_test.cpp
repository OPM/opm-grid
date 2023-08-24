//===========================================================================
//
// File: lookUpCellCentroid_cpgrid_test.cpp
//
// Created: Wed 26.07.2023 16:00:00
//
// Author(s): Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2023 Equinor ASA.

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

#define BOOST_TEST_MODULE LookUpCellCentroidTest
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif

#include <dune/common/version.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/LookUpCellCentroid.hh>

#include <sstream>
#include <iostream>

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

void createEclGridCpGrid_and_checkCentroid(const std::string& deckString)
{
    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);

    Dune::CpGrid grid;
    Opm::EclipseGrid eclGrid(deck);

    grid.processEclipseFormat(&eclGrid, nullptr, false, false, false);

    const auto& leafGridView = grid.leafGridView();
    const Dune::CartesianIndexMapper<Dune::CpGrid> gridCartMapper(grid);

    const Opm::LookUpCellCentroid<Dune::CpGrid, Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>
        lookUpCellCentroid(leafGridView, gridCartMapper, &eclGrid);

    for (const auto& element: Dune::elements(leafGridView)){
        const auto& elemEclCentroid = eclGrid.getCellCenter(gridCartMapper.cartesianIndex(element.index()));
        const auto& elemCpGridEclCentroid_Entity = grid.getEclCentroid(element);
        const auto& elemCpGridEclCentroid_Index = grid.getEclCentroid(element.index());
        const auto& centroid = lookUpCellCentroid(element.index());
        for (int coord = 0; coord < 3; ++coord)
        {
            BOOST_CHECK_EQUAL(elemEclCentroid[coord], elemCpGridEclCentroid_Entity[coord]);
            BOOST_CHECK_EQUAL(elemEclCentroid[coord], elemCpGridEclCentroid_Index[coord]);
            BOOST_CHECK_EQUAL(elemEclCentroid[coord], centroid[coord]);
        }
    }
}


BOOST_AUTO_TEST_CASE(stringA)
{
    const std::string deckString =
        R"( RUNSPEC
        DIMENS
        1  1  5 /
        GRID
        COORD
        0 0 0
        0 0 1
        1 0 0
        1 0 1
        0 1 0
        0 1 1
        1 1 0
        1 1 1
        /
        ZCORN
        4*0
        8*1
        8*2
        8*3
        8*4
        4*5
        /
        ACTNUM
        5*1
        /
        PORO
        5*0.15
        /)";

    createEclGridCpGrid_and_checkCentroid(deckString);
}


BOOST_AUTO_TEST_CASE(stringB)
{
    const std::string deckString =
        R"(RUNSPEC
        DIMENS
        1  1  3 /
        GRID
        COORD
        0 0 0  0 0 3
        1 0 0  1 0 3
        0 1 0  0 1 3
        1 1 0  1 1 3
        /
        ZCORN
        4*0
        8*1
        8*1.1
        4*2.1
        /
        ACTNUM
        3*1
        /
        PINCH
        0.02   NOGAP   1*   1*
        /
        )";

    createEclGridCpGrid_and_checkCentroid(deckString);
}

BOOST_AUTO_TEST_CASE(stringC)
{/*
   Cell corners:      COORD
   0 {0, 0, 2}        line 1: corners 0 and 4
   1 {4, 1, 1}        line 2: corners 1 and 5
   2 {2, 5, 1}        line 3: corners 2 and 6
   3 {5, 4, 1}        line 4: corners 3 and 7
   4 {1, 1, 5}
   5 {4, 0, 4}
   6 {0, 4, 5}
   7 {4, 3, 4} */
    const std::string deckString =
        R"(RUNSPEC
        DIMENS
        1  1  1 /
        GRID
        COORD
        0 0 2  1 1 5
        4 1 1  4 0 4
        2 5 1  0 4 5
        5 4 1  4 3 4
        /
        ZCORN
        4*0
        4*1
        /
        ACTNUM
        1*1
        /
        )";

    createEclGridCpGrid_and_checkCentroid(deckString);
}

BOOST_AUTO_TEST_CASE(stringD)
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
        4*0
        4*1
        /
        ACTNUM
        1*1
        /
        )";

    createEclGridCpGrid_and_checkCentroid(deckString);
}

BOOST_AUTO_TEST_CASE(stringE)
{/*
   Cell corners:       COORD
   0  {0, 0, 0}        line 1: corners 0 and 4
   1  {2, 0, 0}        line 2: corners 1 and 5
   2  {0, 1, 0}        line 3: corners 2 and 6
   3  {2, 1, 0}        line 4: corners 3 and 7
   4  {1, 0, 1}
   5  {3, 0, 1}
   6  {1, 1, 1}
   7  {3, 1, 1} */
    const std::string deckString =
        R"(RUNSPEC
        DIMENS
        1  1  1 /
        GRID
        COORD
        0 0 0  1 0 1
        2 0 0  3 0 1
        0 1 0  1 1 1
        2 1 0  3 1 1
        /
        ZCORN
        4*0.5
        4*1
        /
        ACTNUM
        1*1
        /
        )";

    createEclGridCpGrid_and_checkCentroid(deckString);
}

