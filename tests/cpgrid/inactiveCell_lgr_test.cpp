//===========================================================================
//
// File: inactiveCells_lgr_test.cpp
//
// Created:  July 04 2024 11:06:00
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

#include <config.h>

#define NVERBOSE

#define BOOST_TEST_MODULE InactiveCellsLgr

#include <boost/test/unit_test.hpp>

#include <dune/common/version.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/grid/CpGrid.hpp>

#include <iostream>
#include <vector>
#include <utility>

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

void testInactiveCellsLgrs(const std::string& deckString,
                           const std::vector<std::array<int,3>>& cells_per_dim_vec,
                           const std::vector<std::array<int,3>>& startIJK_vec,
                           const std::vector<std::array<int,3>>& endIJK_vec,
                           const std::vector<std::string>& lgr_name_vec)
{
    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);
    Opm::EclipseState es(deck);

    Dune::CpGrid grid;
    grid.processEclipseFormat(&es.getInputGrid(), &es, false, false, false);


    grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

    auto& data = grid.data_;

    BOOST_CHECK(data.size() == startIJK_vec.size() + 2);
    BOOST_CHECK( (*data[0]).child_to_parent_cells_.empty());
    BOOST_CHECK(grid.lgr_names_["GLOBAL"] == 0);
    const auto& all_parent_cell_indices = (*data[0]).getPatchesCells(startIJK_vec, endIJK_vec);

    for (long unsigned int level = 1; level < startIJK_vec.size() +1; ++level) // only 1 when there is only 1 patch
    {
        BOOST_CHECK( (*data[level]).parent_to_children_cells_.empty());
        BOOST_CHECK(grid.lgr_names_[lgr_name_vec[level-1]] == static_cast<int>(level));

        const auto& patch_cells = (*data[0]).getPatchCells(startIJK_vec[level-1], endIJK_vec[level-1]);

        // GLOBAL grid
        for (int cell = 0; cell <  data[0]-> size(0); ++cell)
        {
            Dune::cpgrid::Entity<0> entity = Dune::cpgrid::Entity<0>((*grid.data_[0]), cell, true);
            BOOST_CHECK( entity.hasFather() == false);
            BOOST_CHECK_THROW(entity.father(), std::logic_error);
            BOOST_CHECK_THROW(entity.geometryInFather(), std::logic_error);
            BOOST_CHECK( entity.getOrigin() ==  entity);
            BOOST_CHECK( entity.getOrigin().level() == 0);
            auto it = entity.hbegin(grid.maxLevel());
            auto endIt = entity.hend(grid.maxLevel());
            const auto& [lgr, childrenList] = (*grid.data_[0]).parent_to_children_cells_[cell];
            if (data[0]->getMark(entity) == 0){
                BOOST_CHECK_EQUAL(lgr, -1);
                BOOST_CHECK(childrenList.empty());
                BOOST_CHECK( entity.isLeaf() == true);
                // If it == endIt, then entity.isLeaf() true (when dristibuted_data_ is empty)
                BOOST_CHECK( it == endIt);
            }
            else{
                BOOST_CHECK(lgr != -1);
                BOOST_CHECK(childrenList.size() > 1);
                // Auxiliary int to check amount of children
                double referenceElemOneParent_volume = 0.;
                std::array<double,3> referenceElem_entity_center = {0.,0.,0.}; // Expected {.5,.5,.5}
                for (const auto& child : childrenList) {
                    BOOST_CHECK( child != -1);
                    BOOST_CHECK( data[lgr]-> child_to_parent_cells_[child][0] == 0);
                    BOOST_CHECK( data[lgr]-> child_to_parent_cells_[child][1] == cell);

                    const auto& childElem =  Dune::cpgrid::Entity<0>(*data[lgr], child, true);
                    BOOST_CHECK(childElem.hasFather() == true);
                    BOOST_CHECK(childElem.level() == lgr);
                    referenceElemOneParent_volume += childElem.geometryInFather().volume();
                    for (int c = 0; c < 3; ++c)  {
                        referenceElem_entity_center[c] += (childElem.geometryInFather().center())[c];
                    }
                }
                BOOST_CHECK_EQUAL( entity.isLeaf(), false); // parent cells do not appear in the LeafView
                // Auxiliary int to check hierarchic iterator functionality
                double referenceElemOneParent_volume_it = 0.;
                std::array<double,3> referenceElem_entity_center_it = {0.,0.,0.}; // Expected {.5,.5,.5}
                // If it != endIt, then entity.isLeaf() false (when dristibuted_data_ is empty)
                BOOST_CHECK( it != endIt );
                for (; it != endIt; ++it)
                {
                    // Do something with the son available through it->
                    BOOST_CHECK(it ->hasFather() == true);
                    BOOST_CHECK(it ->level() == lgr);
                    referenceElemOneParent_volume_it += it-> geometryInFather().volume();
                    for (int c = 0; c < 3; ++c)
                    {
                        referenceElem_entity_center_it[c] += (it-> geometryInFather().center())[c];
                    }
                }
                for (int c = 0; c < 3; ++c)
                {
                    referenceElem_entity_center[c]
                        /= cells_per_dim_vec[lgr-1][0]*cells_per_dim_vec[lgr-1][1]*cells_per_dim_vec[lgr-1][2];
                    referenceElem_entity_center_it[c]
                        /= cells_per_dim_vec[lgr-1][0]*cells_per_dim_vec[lgr-1][1]*cells_per_dim_vec[lgr-1][2];
                }
                BOOST_CHECK_CLOSE(referenceElemOneParent_volume, 1, 1e-13);
                BOOST_CHECK_CLOSE(referenceElem_entity_center[0], .5, 1e-13);
                BOOST_CHECK_CLOSE(referenceElem_entity_center[1], .5, 1e-13);
                BOOST_CHECK_CLOSE(referenceElem_entity_center[2], .5, 1e-13);
                BOOST_CHECK_CLOSE(referenceElemOneParent_volume_it, 1, 1e-13);
                BOOST_CHECK_CLOSE(referenceElem_entity_center_it[0], .5, 1e-13);
                BOOST_CHECK_CLOSE(referenceElem_entity_center_it[1], .5, 1e-13);
                BOOST_CHECK_CLOSE(referenceElem_entity_center_it[2], .5, 1e-13);
            }
            BOOST_CHECK( entity.level() == 0);
        }

        // LGRs
        for (int cell = 0; cell <  data[level]-> size(0); ++cell)
        {
            Dune::cpgrid::Entity<0> entity = Dune::cpgrid::Entity<0>((*grid.data_[level]), cell, true);
            BOOST_CHECK( entity.hasFather() == true);
            BOOST_CHECK( entity.getOrigin() ==  entity.father());
            BOOST_CHECK( entity.index() == (data[level] -> global_cell_[entity.index()])); // global_cell_ = {0,1,..., total cells -1}
            BOOST_CHECK( entity.getOrigin().level() == 0);
            BOOST_CHECK_CLOSE(entity.geometryInFather().volume(),
                              1./(cells_per_dim_vec[level-1][0]*cells_per_dim_vec[level-1][1]*cells_per_dim_vec[level-1][2]), 1e-6);
            BOOST_CHECK(entity.father().level() == 0);
            // Check entity.father() has been marked for refinement (corresponding LGR parents)
            BOOST_CHECK_EQUAL( data[0]->getMark(entity.father()), 1);
            const auto& child_to_parent = (*data[level]).child_to_parent_cells_[cell];
            BOOST_CHECK_EQUAL( child_to_parent[0] == -1, false);
            BOOST_CHECK_EQUAL( child_to_parent[0] == 0, true);
            BOOST_CHECK_EQUAL( child_to_parent[1], entity.father().index());
            BOOST_CHECK( std::get<0>((*data[0]).parent_to_children_cells_[child_to_parent[1]]) == entity.level());
            BOOST_CHECK_EQUAL((std::find(std::get<1>((*data[0]).parent_to_children_cells_[child_to_parent[1]]).begin(),
                                         std::get<1>((*data[0]).parent_to_children_cells_[child_to_parent[1]]).end(),
                                         entity.index()) ==
                               std::get<1>((*data[0]).parent_to_children_cells_[child_to_parent[1]]).end()) , false);
            // Check amount of children cells of the parent cell
            BOOST_CHECK_EQUAL(std::get<1>((*data[0]).parent_to_children_cells_[child_to_parent[1]]).size(),
                              cells_per_dim_vec[level-1][0]*cells_per_dim_vec[level-1][1]*cells_per_dim_vec[level-1][2]);
            BOOST_CHECK( entity.level() == static_cast<int>(level));
            BOOST_CHECK( entity.isLeaf() == true);
            auto it = entity.hbegin(grid.maxLevel());
            auto endIt = entity.hend(grid.maxLevel());
            // If entity.isLeaf(), then it == endIt (when dristibuted_data_ is empty)
            BOOST_CHECK( it == endIt);
        }

        // LeafView faces
        for (int face = 0; face <  data[startIJK_vec.size()+1]-> face_to_cell_.size(); ++face)
        {
            const auto& faceToPoint =  (*data[startIJK_vec.size() +1]).face_to_point_[face];
            BOOST_CHECK(faceToPoint.size() == 4);
            for (int i = 0; i < 4; ++i) {
                BOOST_CHECK((*data[startIJK_vec.size() +1]).face_to_point_[face][i] != -1);
            }

            Dune::cpgrid::EntityRep<1> faceEntity(face, true);
            BOOST_CHECK((*data[startIJK_vec.size() +1]).face_to_cell_[faceEntity].size() < 3);
        }

        // LeafView
        for (int cell = 0; cell <  data[startIJK_vec.size()+1]-> size(0); ++cell)
        {
            BOOST_CHECK( data[startIJK_vec.size()+1] -> cell_to_point_[cell].size() == 8);
            for (int i = 0; i < 8; ++i)
            {
                BOOST_CHECK( data[startIJK_vec.size()+1] -> cell_to_point_[cell][i] != -1);
            }
            Dune::cpgrid::Entity<0> entity = Dune::cpgrid::Entity<0>((*grid.data_[startIJK_vec.size()+1]), cell, true);
            for (int i = 0; i < data[startIJK_vec.size()+1] -> cell_to_face_[entity].size(); ++i)
            {
                BOOST_CHECK( data[startIJK_vec.size()+1] -> cell_to_face_[entity][i].index() != -1);
            }
            const auto& child_to_parent = (*data[startIJK_vec.size()+1]).child_to_parent_cells_[cell];
            const auto& level_cellIdx = (*data[startIJK_vec.size()+1]).leaf_to_level_cells_[entity.index()];
            auto it = entity.hbegin(grid.maxLevel());
            auto endIt = entity.hend(grid.maxLevel());
            BOOST_CHECK(entity.isLeaf());
            BOOST_CHECK(it == endIt);
            if (entity.hasFather()){
                BOOST_CHECK_CLOSE(entity.geometryInFather().volume(),
                                  1./(cells_per_dim_vec[entity.level()-1][0]
                                      *cells_per_dim_vec[entity.level()-1][1]
                                      *cells_per_dim_vec[entity.level()-1][2]), 1e-6);
                BOOST_CHECK(entity.father().level() == 0);
                BOOST_CHECK_EQUAL( data[0]->getMark(entity.father()), 1);
                BOOST_CHECK(child_to_parent[0] != -1);
                BOOST_CHECK_EQUAL( child_to_parent[0] == 0, true);
                BOOST_CHECK_EQUAL( child_to_parent[1], entity.father().index());
                BOOST_CHECK( entity.father() == entity.getOrigin());
                BOOST_CHECK(  (data[startIJK_vec.size() +1] -> global_cell_[entity.index()]) ==
                              (data[0] -> global_cell_[entity.getOrigin().index()]) );
                BOOST_CHECK( entity.getOrigin().level() == 0);
                BOOST_CHECK( std::get<0>((*data[0]).parent_to_children_cells_[child_to_parent[1]]) == entity.level());
                BOOST_CHECK_EQUAL((std::find(std::get<1>((*data[0]).parent_to_children_cells_[child_to_parent[1]]).begin(),
                                             std::get<1>((*data[0]).parent_to_children_cells_[child_to_parent[1]]).end(),
                                             level_cellIdx[1]) ==
                                   std::get<1>((*data[0]).parent_to_children_cells_[child_to_parent[1]]).end()) , false);
                // Check amount of children cells of the parent cell
                BOOST_CHECK_EQUAL(std::get<1>((*data[0]).parent_to_children_cells_[child_to_parent[1]]).size(),
                                  cells_per_dim_vec[entity.level()-1][0]*
                                  cells_per_dim_vec[entity.level()-1][1]*cells_per_dim_vec[entity.level()-1][2]);
                BOOST_CHECK( entity.father().isLeaf() == false);
                BOOST_CHECK( (entity.level() > 0) || (entity.level() < static_cast<int>(startIJK_vec.size()) +1));
                BOOST_CHECK( level_cellIdx[0] == entity.level());
            }
            else{
                BOOST_CHECK_THROW(entity.father(), std::logic_error);
                BOOST_CHECK_THROW(entity.geometryInFather(), std::logic_error);
                BOOST_CHECK_EQUAL(child_to_parent[0], -1);
                BOOST_CHECK_EQUAL(child_to_parent[1], -1);
                BOOST_CHECK( level_cellIdx[0] == 0);
                BOOST_CHECK( std::get<0>((*data[0]).parent_to_children_cells_[level_cellIdx[1]]) == -1);
                BOOST_CHECK( std::get<1>((*data[0]).parent_to_children_cells_[level_cellIdx[1]]).empty());
                BOOST_CHECK( entity.level() == 0);
                // Get index of the cell in level 0
                const auto& entityOldIdx =  (*data[startIJK_vec.size()+1]).leaf_to_level_cells_[entity.index()][1];
                BOOST_CHECK( entity.getOrigin().index() == entityOldIdx);
                BOOST_CHECK( entity.getOrigin().level() == 0);
                // Get IJK of the old index
                std::array<int,3> entityOldIJK;
                (*data[0]).getIJK(entityOldIdx, entityOldIJK); // ijk
            }
        }
    } // end-level-for-loop

    BOOST_CHECK( static_cast<int>(startIJK_vec.size()) == grid.maxLevel());
    BOOST_CHECK( (*data[data.size()-1]).parent_to_children_cells_.empty());

    for (long unsigned int l = 0; l < startIJK_vec.size() +1; ++l) // level 0,1,2,... , last patch
    {
        const auto& view = grid.levelGridView(l);
        for (const auto& element: elements(view)){
            BOOST_CHECK_EQUAL(element.level(), l);
        }
    }

    std::set<int> allIds_set;
    std::vector<int> allIds_vec;
    allIds_vec.reserve(data.back()->size(0) + data.back()->size(3));
    for (const auto& element: elements(grid.leafGridView())){
        const auto& localId = data.back()->localIdSet().id(element);
        const auto& globalId = data.back()->globalIdSet().id(element);
        // In serial run, local and global id coincide:
        BOOST_CHECK_EQUAL(localId, globalId);
        allIds_set.insert(localId);
        allIds_vec.push_back(localId);
        // Check that the global_id_set_ptr_ has the correct id (id from the level where the entity was born).
        BOOST_CHECK_EQUAL( grid.globalIdSet().id(element), (*data[element.level()]).localIdSet().id(element.getEquivLevelElem()));
    }
    // Check injectivity of the map local_id_set_ (and, indirectly, global_id_set_) after adding cell ids.
    BOOST_CHECK( allIds_set.size() == allIds_vec.size());

    for (const auto& point: vertices(grid.leafGridView())){
        const auto& localId = data.back()->localIdSet().id(point);
        const auto& globalId = data.back()->globalIdSet().id(point);
        BOOST_CHECK_EQUAL(localId, globalId);
        allIds_set.insert(localId);
        allIds_vec.push_back(localId);
    }
    // Check injectivity of the map local_id_set_ (and, indirectly, global_id_set_) after adding point ids.
    BOOST_CHECK( allIds_set.size() == allIds_vec.size());

    // Local/Global id sets for level grids (level 0, 1, ..., maxLevel)
    for (int level = 0; level < grid.maxLevel() +1; ++level)
    {
        std::set<int> levelIds_set;
        std::vector<int> levelIds_vec;
        levelIds_vec.reserve(data[level]->size(0) + data[level]->size(3));

        for (const auto& element: elements(grid.levelGridView(level))){
            const auto& localId = data[level]->localIdSet().id(element);
            const auto& globalId = data[level]->globalIdSet().id(element);
            // In serial run, local and global id coincide:
            BOOST_CHECK_EQUAL(localId, globalId);
            levelIds_set.insert(localId);
            levelIds_vec.push_back(localId);
            // Search in the leaf grid view elements for the element with the same id, if it exists.
            if (auto itIsLeaf = std::find_if( elements(grid.leafGridView()).begin(),
                                              elements(grid.leafGridView()).end(),
                                              [localId, data](const Dune::cpgrid::Entity<0>& leafElem)
                                              { return (localId == data.back()->localIdSet().id(leafElem)); });
                itIsLeaf != elements(grid.leafGridView()).end()) {
                BOOST_CHECK( itIsLeaf->getEquivLevelElem() == element);
            }
            if (element.isLeaf()) { // Check that the id of a cell not involved in any further refinement appears on the IdSet of the leaf grid view.
                BOOST_CHECK( std::find(allIds_set.begin(), allIds_set.end(), localId) != allIds_set.end());
            }
            else { // Check that the id of a cell that vanished during refinement does not appear on the IdSet of the leaf grid view.
                BOOST_CHECK( std::find(allIds_set.begin(), allIds_set.end(), localId) == allIds_set.end());
            }
            const auto& idx = data[level]->indexSet().index(element);
            // In serial run, local and global id coincide:
            BOOST_CHECK_EQUAL(idx, element.index());
        }

        for (const auto& point : vertices(grid.levelGridView(level))) {
            const auto& localId = data[level]->localIdSet().id(point);
            const auto& globalId = data[level]->globalIdSet().id(point);
            BOOST_CHECK_EQUAL(localId, globalId);
            levelIds_set.insert(localId);
            levelIds_vec.push_back(localId);
            // Search in the leaf grid view elements for the element with the same id, if it exists.
            if (auto itIsLeaf = std::find_if( vertices(grid.leafGridView()).begin(),
                                              vertices(grid.leafGridView()).end(),
                                              [localId, data](const Dune::cpgrid::Entity<3>& leafPoint)
                                              { return (localId == data.back()->localIdSet().id(leafPoint)); });
                itIsLeaf != vertices(grid.leafGridView()).end()) {
                BOOST_CHECK( (*itIsLeaf).geometry().center() == point.geometry().center());
            }
        }
        // Check injectivity of the map local_id_set_ (and, indirectly, global_id_set_)
        BOOST_CHECK( levelIds_set.size() == levelIds_vec.size());
    }
}

BOOST_GLOBAL_FIXTURE(Fixture);


BOOST_AUTO_TEST_CASE(inactiveCells_in_at_least_one_lgr)
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
        1
        /
        PORO
        5*0.15
        /)";

    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}};
    // LGR1 cell indices: {0,1}, 0 INACTIVE CELL, 1 ACTIVE CELL
    // LGR2 cell indices: {4} ACTIVE CELL
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,4}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,2}, {1,1,5}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    testInactiveCellsLgrs(deckString, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
}

BOOST_AUTO_TEST_CASE(inactiveCells_outside_lgr)
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
        1
        1
        1
        0
        1
        /
        PORO
        5*0.15
        /)";

    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}, {4,4,4}};
    // LGR1 cell indices: {0}
    // LGR2 cell indices: {2}
    // LGR3 cell indices: {4}
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,2}, {0,0,4}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,1}, {1,1,3}, {1,1,5}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3"};
    testInactiveCellsLgrs(deckString, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
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
    testInactiveCellsLgrs(deckString, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
}



BOOST_AUTO_TEST_CASE(inactiveCells_in_and_outside_lgrs)
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
        1
        0
        0
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
    testInactiveCellsLgrs(deckString, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
}

BOOST_AUTO_TEST_CASE(inactiveCells_only_outside_lgrs)
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
        1
        1
        0
        0
        1
        /
        PORO
        5*0.15
        /)";

    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,4}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,2}, {1,1,5}};
    // LGR1 cell indices = {0,1}, LGR2 cell indices = {4}.
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    testInactiveCellsLgrs(deckString, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
}


BOOST_AUTO_TEST_CASE(warning_allInactiveCells_in_lgr)
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
        0
        0
        0
        1
        /
        PORO
        5*0.15
        /)";

    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,4}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,2}, {1,1,5}};
    // LGR1 cell indices = {0,1}, LGR2 cell indices = {4}.
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    testInactiveCellsLgrs(deckString, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
}


