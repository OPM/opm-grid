//===========================================================================
//
// File: lgr_cartesian_idx_test.cpp
//
// Created: Sep 26 2024 09:15
//
// Author(s): Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================
/*
  Copyright 2024 Equinor ASA.

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

#define BOOST_TEST_MODULE LGRTests
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/DefaultGeometryPolicy.hpp>
#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/cpgrid/EntityRep.hpp>
#include <opm/grid/cpgrid/Geometry.hpp>
#include <opm/grid/cpgrid/CartesianIndexMapperCollection.hpp>
#include <opm/grid/LookUpData.hh>

#include <dune/common/version.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <cassert>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <map>

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

void checkGlobalCellLgr(Dune::CpGrid& grid)
{
    const Dune::CartesianIndexMapper<Dune::CpGrid> mapper{grid};
    const Opm::CartesianIndexMapperCollection<Dune::CpGrid> cartMappCollection{grid};

    for (const auto& element : elements(grid.leafGridView()))
    {
        // How to get the Cartesian Index of a cell on the leaf grid view.
        // 1. The "default" Cartesian Index of a cell on the leaf (which is a mixed grid with coarse and refined cells) is equal to the Cartesian Index of:
        //    (a) the equivalent cell in level zero (same center and volume, cell not involved in any refinement).
        //    (b) the parent cell in level zero when the leaf cell is a refined one (it has a equivalent cell that belongs to a refined level grid).
        //    This default Cartesian Index for leaf cells coincide with the globalCell values.
        const auto& cartesian_idx_from_leaf_elem =  mapper.cartesianIndex( element.index() );
        const auto& global_cell_idx_leaf = grid.globalCell()[element.index()]; // parent cell index when the cell is a refined one.
        BOOST_CHECK_EQUAL(cartesian_idx_from_leaf_elem, global_cell_idx_leaf);
        // global_ijk represents the ijk values of the equivalent cell on the level zero, or parent cell if the leaf cell is a refined one.
        // Notice that all the refined cells of a same parent cell will get the same global_cell_ value and the same global_ijk (since they
        // inherit the parent cell value).
        std::array<int,3> global_ijk = {0,0,0};
        mapper.cartesianCoordinate( element.index(), global_ijk);

        // How to get the Level Cartesian Index of a cell on the leaf grid view.
        // Each LGR can be seen as a Cartesian Grid itself, with its own (local) logical_cartesian_size and its own (local) Cartesian indices.
        // Given a leaf cell, via the CartesianIndexMapper and its method cartesianIndexLevel(...), we get the local-Cartesian-index.
        const auto& levelCartMapp = cartMappCollection.getLevelMapper(element.level());
        const auto& cartesian_idx_from_level_elem =  levelCartMapp.cartesianIndex( element.getLevelElem().index());
        const auto& global_cell_idx_level = grid.currentData()[element.level()]->globalCell()[element.getLevelElem().index()];
        BOOST_CHECK_EQUAL(cartesian_idx_from_level_elem, global_cell_idx_level);
        // local_ijk represents the ijk values of the equivalent cell on the level its was born.
        std::array<int,3> local_ijk = {0,0,0};
        levelCartMapp.cartesianCoordinate( element.getLevelElem().index(), local_ijk);

        // For leaf cells that were not involved in any refinement, global_ijk and local_ijk must coincide.
        if(element.level()==0)
        {
            BOOST_CHECK_EQUAL( global_ijk[0], local_ijk[0]);
            BOOST_CHECK_EQUAL( global_ijk[1], local_ijk[1]);
            BOOST_CHECK_EQUAL( global_ijk[2], local_ijk[2]);
        }
    }

    const auto& localCartesianIdxSets_to_leafIdx = grid.mapLocalCartesianIndexSetsToLeafIndexSet();
    const auto& leafIdx_to_localCartesianIdxSets = grid.mapLeafIndexSetToLocalCartesianIndexSets();
    
    for (int level = 0; level < grid.maxLevel(); ++level)
    {
        for (const auto& element : elements(grid.levelGridView(level)))
        {
            const auto& global_cell_level = grid.currentData()[element.level()]->globalCell()[element.index()];
            if(element.isLeaf()) {
                const auto& leaf_idx = localCartesianIdxSets_to_leafIdx[element.level()].at(global_cell_level);
                BOOST_CHECK_EQUAL( leafIdx_to_localCartesianIdxSets[leaf_idx][0], element.level());
                BOOST_CHECK_EQUAL( leafIdx_to_localCartesianIdxSets[leaf_idx][1], global_cell_level);
            }
            else {
                BOOST_CHECK_THROW( localCartesianIdxSets_to_leafIdx[element.level()].at(global_cell_level), std::out_of_range);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(refine_one_cell)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);

    const std::array<int, 3> cells_per_dim = {2,2,2};
    const std::array<int, 3> startIJK = {1,0,1};
    const std::array<int, 3> endIJK = {2,1,2};// Single Cell! (with index 13 in level zero grid), LGR dimensions 2x2x2
    const std::string lgr_name = {"LGR1"};
    grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

    checkGlobalCellLgr(grid);
}

BOOST_AUTO_TEST_CASE(three_lgrs)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);

    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}, {4,4,4}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,2}, {3,2,2}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,1,1}, {1,1,3}, {4,3,3}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3"};
    grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

    checkGlobalCellLgr(grid);
}

BOOST_AUTO_TEST_CASE(inactiveCells_in_lgrs)
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
        0
        1
        1
        1
        0
        /
        PORO
        5*0.15
        /)";

    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,3}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,2}, {1,1,5}};
    // LGR1 cell indices = {0,1}, LGR2 cell indices = {3,4}.
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};

    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);
    Opm::EclipseState es(deck);

    Dune::CpGrid grid;
    grid.processEclipseFormat(&es.getInputGrid(), &es, false, false, false);

    grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
    
    checkGlobalCellLgr(grid);
}
