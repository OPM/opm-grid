//===========================================================================
//
// File: coord_zcorn_lgr_test.cpp
//
// Created: Nov 5 10:05 2024
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

#define BOOST_TEST_MODULE CoordZcornLgrTest
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/LookUpData.hh>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/DefaultGeometryPolicy.hpp>
#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/cpgrid/EntityRep.hpp>
#include <opm/grid/cpgrid/Geometry.hpp>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>

#include <dune/common/version.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/Grid/LgrCollection.hpp>
#include <opm/input/eclipse/EclipseState/Grid/Carfin.hpp>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <sstream>

struct Fixture {
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

void
checkGlobalCellLgr(Dune::CpGrid& grid)
{
    const Dune::CartesianIndexMapper<Dune::CpGrid> mapper {grid};

    const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);

    for (const auto& element : elements(grid.leafGridView())) {
        // How to get the Cartesian Index of a cell on the leaf grid view.
        // 1. The "default" Cartesian Index of a cell on the leaf (which is a mixed grid with coarse and refined cells)
        // is equal to the Cartesian Index of:
        //    (a) the equivalent cell in level zero (same center and volume, cell not involved in any refinement).
        //    (b) the parent cell in level zero when the leaf cell is a refined one (it has a equivalent cell that
        //    belongs to a refined level grid). This default Cartesian Index for leaf cells coincide with the globalCell
        //    values.
        const auto& cartesian_idx_from_leaf_elem = mapper.cartesianIndex(element.index());
        const auto& global_cell_idx_leaf
            = grid.globalCell()[element.index()]; // parent cell index when the cell is a refined one.
        BOOST_CHECK_EQUAL(cartesian_idx_from_leaf_elem, global_cell_idx_leaf);
        // global_ijk represents the ijk values of the equivalent cell on the level zero, or parent cell if the leaf
        // cell is a refined one. Notice that all the refined cells of a same parent cell will get the same global_cell_
        // value and the same global_ijk (since they inherit the parent cell value).
        std::array<int, 3> global_ijk = {0, 0, 0};
        mapper.cartesianCoordinate(element.index(), global_ijk);

        // How to get the Level Cartesian Index of a cell on the leaf grid view.
        // Each LGR can be seen as a Cartesian Grid itself, with its own (local) logical_cartesian_size and its own
        // (local) Cartesian indices. Given a leaf cell, via the CartesianIndexMapper and its method
        // cartesianIndexLevel(...), we get the local-Cartesian-index.
        const auto& cartesian_idx_from_level_elem
            = levelCartMapp.cartesianIndex(element.getEquivLevelElem().index(), element.level());
        const auto& global_cell_idx_level
            = grid.currentData()[element.level()]->globalCell()[element.getEquivLevelElem().index()];
        BOOST_CHECK_EQUAL(cartesian_idx_from_level_elem, global_cell_idx_level);
        // local_ijk represents the ijk values of the equivalent cell on the level its was born.
        std::array<int, 3> local_ijk = {0, 0, 0};
        levelCartMapp.cartesianCoordinate(element.getEquivLevelElem().index(), local_ijk, element.level());

        // For leaf cells that were not involved in any refinement, global_ijk and local_ijk must coincide.
        if (element.level() == 0) {
            BOOST_CHECK_EQUAL(global_ijk[0], local_ijk[0]);
            BOOST_CHECK_EQUAL(global_ijk[1], local_ijk[1]);
            BOOST_CHECK_EQUAL(global_ijk[2], local_ijk[2]);
        }
    }

    const auto& localCartesianIdxSets_to_leafIdx = grid.mapLocalCartesianIndexSetsToLeafIndexSet();
    const auto& leafIdx_to_localCartesianIdxSets = grid.mapLeafIndexSetToLocalCartesianIndexSets();

    for (int level = 0; level < grid.maxLevel()+1; ++level) {
        for (const auto& element : elements(grid.levelGridView(level))) {
            const auto& global_cell_level = grid.currentData()[element.level()]->globalCell()[element.index()];
            if (element.isLeaf()) {
                const auto& leaf_idx = localCartesianIdxSets_to_leafIdx[element.level()].at(global_cell_level);
                BOOST_CHECK_EQUAL(leafIdx_to_localCartesianIdxSets[leaf_idx][0], element.level());
                BOOST_CHECK_EQUAL(leafIdx_to_localCartesianIdxSets[leaf_idx][1], global_cell_level);
            } else {
                BOOST_CHECK_THROW(localCartesianIdxSets_to_leafIdx[element.level()].at(global_cell_level),
                                  std::out_of_range);
            }
        }
    }


    for (int level = 1; level < grid.maxLevel()+1; ++level) {
        // Get COORD values for LGRs
        // COORD defines a set of coordinate lines or pillars for a reservoir grid via an array.
        // A total of 6 x (NX+1) x (NY+1) lines must be specified for each coordinate data set (or reservoir).
        // For Cartesian geometry, each line is defined by the (x, y, z) coordinates of two distinct points on the line. The
        // lines are entered with I cycling fastest then J.
        //
        // For each LGR, we collect the coordinates of the vertices. Attention: pillar ordering not guaranteed.
        std::vector<std::array<double,3>> xyz_coord(grid.currentData()[level]->size(3)); // NOT ORDERED IN A 'PILLAR' WAY
        for(const auto& point : vertices(grid.levelGridView(level))) {
            for(int i = 0; i<3; ++i) {
                xyz_coord[point.index()][i] = point.geometry().center()[i];
            }
            std::cout<< "level: " << level << " x: " << xyz_coord[point.index()][0] << " y: "<< xyz_coord[point.index()][1] << " z: "<< xyz_coord[point.index()][2] << std::endl;
        }

        // Get ZCORN values for LGRs
        // ZCORN defines the depth of each corner point of a grid block on the pillars defining the reservoir grid. A
        // total of 8 x NX x NY x NZ values are needed to fully define all the depths in the model. The depths
        // specifying the top of the first layer are entered first with one point for each pillar for each grid block. The
        // points are entered with the X axis cycling fastest. Next come the depths of the bottom of the first layer. The
        // top of layer two follows etc.
        //
        // Expected values in each LGR: 8*refined_cells_in_that_lgr
        // Each refined cell has its 8 corners stored in <refined_level_grid>.cell_to_point_[ refined element index ]
        // in the following order:
        //     BOTTOM FACE          TOP FACE
        //     2 --- 3              6 --- 7
        //    /     /              /     /
        //   0 --- 1              4 --- 5
        //
        std::vector<std::vector<std::array<double,3>>> xyz_for_zcorn(grid.currentData()[level]->size(0)); // size total amount of cells in the level refined grid.
        for (const auto& element : elements(grid.levelGridView(level))) {
            xyz_for_zcorn[element.index()].reserve(8);
            for (int point_idx = 0; point_idx < 8; ++point_idx) {
                const auto& point = element.subEntity<3>(point_idx);
                std::array<double,3> xyz = {0., 0., 0.};
                for(int i = 0; i<3; ++i) {
                    xyz[i] = point.geometry().center()[i];
                    std::cout<< "level: " << level << " x: " << xyz[0] << " y: "<< xyz[1] << " z: "<< xyz[2] << std::endl;
                }
                xyz_for_zcorn[element.index()].push_back(xyz);
            }
        }
    }
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
        CARFIN
       -- NAME I1-I2 J1-J2 K1-K2 NX NY NZ
       'LGR1'   1  1  1  1  1  2  2  2  4/
        ENDFIN

        CARFIN
        -- NAME I1-I2 J1-J2 K1-K2 NX NY NZ
        'LGR2'  1  1  1  1  4  5   3  3   6/
        ENDFIN

        PORO
        5*0.15
        /)";

    // Each parent cell from LGR1-block-to-be-refined will be refined into 2x2x2 child cells in each direction x-,y-,
    // and z- respect. Block of parent cells for LGR1 starts at {0,0,0} 'ijk' from the global CpGrid (or {1,1,1} in
    // deck). Block of parent cells for LGR1 ends at {1,1,2} 'ijk' from the global CpGrid (or {1,1,2} in deck). Block of
    // parent cell global indices for LGR1 = {0 (inactive -> no refinement takes place), 1 (active -> cell gets
    // refined)} LGR1 dim = {2,2,4} refined cells in each direction. (Block of parent cells dim = {1,1,2}).

    // Each parent cell from LGR2-block-to-be-refined will be refined into 3x3x3 child cells in each direction x-,y-,
    // and z- respect. Block of parent cells for LGR2 starts at {0,0,3} 'ijk' from the global CpGrid (or {1,1,4} in
    // deck) Block of parent cells for LGR2 ends at {1,1,5} 'ijk' from the global CpGrid (or {1,1,5} in deck). Block of
    // parent cell glboal indices for LGR2 = {3 (active -> cell gets refined), 4 (inactive -> no refinement takes
    // place)} LGR2 dim = {3,3,6} refined cells in each direction. (Block of parent cells dim = {1,1,2}).

    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);
    Opm::EclipseState es(deck);

    Dune::CpGrid grid;
    grid.processEclipseFormat(&es.getInputGrid(), &es, false, false, false); // CARFIN ignore!

    Dune::CpGrid lgr1_grid;
    Opm::LgrCollection lgrs = es.getLgrs();
    const auto& lgr1 = es.getLgrs().getLgr("LGR1");
    BOOST_CHECK_EQUAL(lgr1.NAME(), "LGR1");


    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,3}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,2}, {1,1,5}};
    // LGR1 cell indices = {0,1}, LGR2 cell indices = {3,4}.
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};

    grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

    checkGlobalCellLgr(grid);
}
