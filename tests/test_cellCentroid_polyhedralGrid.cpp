//===========================================================================
//
// File: test_cellCentroid_polyhedralGrid.cpp
//
// Created: Tue 01.08.2023 10:34:00
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

#define BOOST_TEST_MODULE CellCentroidPolyhedralGridTest
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
#include <opm/grid/polyhedralgrid.hh>
#include <opm/grid/LookUpCellCentroid.hh>

#include <dune/common/unused.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <opm/grid/cpgrid/dgfparser.hh>
#include <opm/grid/polyhedralgrid/dgfparser.hh>

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

void createEclGridPolyhedralGrid_and_checkCentroid(const std::string& deckString)
{
    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);
    Opm::EclipseGrid eclGrid(deck);

    std::vector<double> porv;
    Dune::PolyhedralGrid<3,3> grid(eclGrid, porv);


    const auto& leafGridView = grid.leafGridView();

    using GridView = std::remove_cv_t< typename std::remove_reference<decltype(leafGridView)>::type>;

    const Dune::CartesianIndexMapper<Dune::PolyhedralGrid<3,3>> gridCartMapper(grid);
    const Opm::LookUpCellCentroid<Dune::PolyhedralGrid<3,3>,GridView> lookUpCellCentroid(leafGridView, gridCartMapper, &eclGrid);
    // Mapper
    const Dune::MultipleCodimMultipleGeomTypeMapper<GridView> mapper(leafGridView, Dune::mcmgElementLayout());

    for (const auto& element: Dune::elements(leafGridView)){
        const int idx = mapper.index(element);
        const auto& cartIdx = gridCartMapper.cartesianIndex(idx);
        const auto& elemEclCentroid = eclGrid.getCellCenter(cartIdx);
        const auto& centroid = lookUpCellCentroid.getCentroidFromEclGrid(idx);
        for (int coord = 0; coord < 3; ++coord)
        {
            BOOST_CHECK_EQUAL(elemEclCentroid[coord], centroid[coord]);
        }
    }
}


BOOST_AUTO_TEST_CASE(PolyGridFromEcl)
{
#if HAVE_ECL_INPUT
    const char *deckString =
        "RUNSPEC\n"
        "METRIC\n"
        "DIMENS\n"
        "4 4 4 /\n"
        "GRID\n"
        "DXV\n"
        "4*1 /\n"
        "DYV\n"
        "4*1 /\n"
        "DZ\n"
        "16*1 /\n"
        "TOPS\n"
        "16*100.0 /\n";

    createEclGridPolyhedralGrid_and_checkCentroid(deckString);
#endif
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

    createEclGridPolyhedralGrid_and_checkCentroid(deckString);
}
