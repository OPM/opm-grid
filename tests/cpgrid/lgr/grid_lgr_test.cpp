//===========================================================================
//
// File: grid_lgr_test.cpp
//
// Created: ?? 2023
//
// Author(s): Markus Blatt        <markus@dr-blatt.de>
//            Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================
/*
  Copyright 2022-2023 Equinor ASA.

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
#include <opm/grid/LookUpData.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <sstream>
#include <iostream>
#include <cstdlib>
#include <cmath>

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

#define CHECK_COORDINATES(c1, c2)                                       \
    for (int c = 0; c < 3; c++) {                                       \
        BOOST_TEST(c1[c] == c2[c], boost::test_tools::tolerance(1e-12)); \
    }


void check_refinedPatch_grid(const std::array<int,3> cells_per_dim,
                             const std::array<int,3> start_ijk,
                             const std::array<int,3> end_ijk,
                             const Dune::cpgrid::EntityVariable<Dune::cpgrid::Geometry<3, 3>,0> refined_cells,
                             const Dune::cpgrid::EntityVariable<Dune::cpgrid::Geometry<2,3>,1> refined_faces,
                             const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>> refined_corners)
{
    const std::array<int,3> patch_dim = {end_ijk[0]-start_ijk[0], end_ijk[1]-start_ijk[1], end_ijk[2]-start_ijk[2]};
    if ((patch_dim[0] == 0) || (patch_dim[1] == 0) || (patch_dim[2] == 0)) {
        OPM_THROW(std::logic_error, "Empty patch. Cannot convert patch into cell.");
    }
    // Check amount of refined faces.
    int count_faces = (cells_per_dim[0]*patch_dim[0]*cells_per_dim[1]*patch_dim[1]*((cells_per_dim[2]*patch_dim[2])+1))// 'bottom/top faces'
        +  (((cells_per_dim[0]*patch_dim[0])+1)*cells_per_dim[1]*patch_dim[1]*cells_per_dim[2]*patch_dim[2]) // 'front/back faces'
        + (cells_per_dim[0]*patch_dim[0]*((cells_per_dim[1]*patch_dim[1]) +1)*cells_per_dim[2]*patch_dim[2]);  // 'left/right faces'
    BOOST_CHECK_EQUAL(refined_faces.size(), count_faces);
    // Check amount of refined corners.
    int count_corners = ((cells_per_dim[0]*patch_dim[0])+1)*((cells_per_dim[1]*patch_dim[1])+1)*((cells_per_dim[2]*patch_dim[2])+1);
    BOOST_CHECK_EQUAL(refined_corners.size(), count_corners);

    int count_cells = cells_per_dim[0]*patch_dim[0]*cells_per_dim[1]*patch_dim[1]*cells_per_dim[2]*patch_dim[2];
    BOOST_CHECK_EQUAL(refined_cells.size(), count_cells);


}

void refinePatch_and_check(Dune::CpGrid& coarse_grid,
                           const std::vector<std::array<int,3>>& cells_per_dim_vec,
                           const std::vector<std::array<int,3>>& startIJK_vec,
                           const std::vector<std::array<int,3>>& endIJK_vec,
                           const std::vector<std::string>& lgr_name_vec)
{
    auto& data = coarse_grid.currentData();
    // Add LGRs and update grid.
    try
    {
        coarse_grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

        BOOST_CHECK(data.size() == startIJK_vec.size() + 2);
        BOOST_CHECK( (*data[0]).child_to_parent_cells_.empty());
        BOOST_CHECK(coarse_grid.getLgrNameToLevel().at("GLOBAL") == 0);
        const auto& all_parent_cell_indices = (*data[0]).getPatchesCells(startIJK_vec, endIJK_vec);

        // GLOBAL grid
        for (int cell = 0; cell <  data[0]-> size(0); ++cell)
        {
            Dune::cpgrid::Entity<0> entity = Dune::cpgrid::Entity<0>(*data[0], cell, true);
            BOOST_CHECK( entity.hasFather() == false);
            BOOST_CHECK_THROW(entity.father(), std::logic_error);
            BOOST_CHECK_THROW(entity.geometryInFather(), std::logic_error);
            BOOST_CHECK( entity.getOrigin() ==  entity);
            BOOST_CHECK( entity.getOrigin().level() == 0);

            auto it = entity.hbegin(coarse_grid.maxLevel()); // With element.level(), fails
            auto endIt = entity.hend(coarse_grid.maxLevel());
            const auto& [lgr, childrenList] = (*data[0]).parent_to_children_cells_[entity.index()];
            // If it == endIt, then entity.isLeaf() true (when dristibuted_data_ is empty)
            if (it == endIt){
                BOOST_CHECK_EQUAL(lgr, -1);
                BOOST_CHECK(childrenList.empty());
                BOOST_CHECK( entity.isLeaf() == true);
                BOOST_CHECK( entity.mightVanish() == false);
            }
            else{
                BOOST_CHECK(lgr != -1);
                // If it != endIt, then entity.isLeaf() false (when dristibuted_data_ is empty)
                BOOST_CHECK( entity.mightVanish() == true);
                BOOST_CHECK( entity.isNew() == false);
                BOOST_CHECK_EQUAL( entity.isLeaf(), false); // parent cells do not appear in the LeafView

                double referenceElemOneParent_volume = 0.;
                std::array<double,3> referenceElem_entity_center = {0.,0.,0.}; // Expected {.5,.5,.5}
                for (const auto& child : childrenList) {
                    BOOST_CHECK( child != -1);
                    BOOST_CHECK( data[lgr]-> child_to_parent_cells_[child][0] == 0);
                    BOOST_CHECK( data[lgr]-> child_to_parent_cells_[child][1] == entity.index());

                    const auto& childElem =  Dune::cpgrid::Entity<0>(*data[lgr], child, true);
                    const auto& childGeoInFather = childElem.geometryInFather();

                    BOOST_CHECK(childElem.hasFather() == true);
                    BOOST_CHECK(childElem.level() == lgr);
                    referenceElemOneParent_volume += childGeoInFather.volume();
                    for (int c = 0; c < 3; ++c)  {
                        referenceElem_entity_center[c] += childGeoInFather.center()[c];
                    }
                }

                double referenceElemOneParent_volume_it = 0.;
                std::array<double,3> referenceElem_entity_center_it = {0.,0.,0.}; // Expected {.5,.5,.5}
                for (; it != endIt; ++it)
                {
                    BOOST_CHECK(it ->hasFather() == true);
                    BOOST_CHECK(it ->level() == lgr);
                    const auto& geoInFather =  it-> geometryInFather();

                    referenceElemOneParent_volume_it += geoInFather.volume();
                    for (int c = 0; c < 3; ++c)
                    {
                        referenceElem_entity_center_it[c] += geoInFather.center()[c];
                    }
                }
                for (int c = 0; c < 3; ++c) {
                    referenceElem_entity_center[c] /= cells_per_dim_vec[lgr-1][0]*cells_per_dim_vec[lgr-1][1]*cells_per_dim_vec[lgr-1][2];
                    referenceElem_entity_center_it[c] /= cells_per_dim_vec[lgr-1][0]*cells_per_dim_vec[lgr-1][1]*cells_per_dim_vec[lgr-1][2];
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
        } // end-GLOBAL/LevelZero-grid


        for (long unsigned int level = 1; level < startIJK_vec.size() +1; ++level) // only 1 when there is only 1 patch
        {
            check_refinedPatch_grid(cells_per_dim_vec[level-1], startIJK_vec[level-1], endIJK_vec[level-1],
                                    (*data[level]).geometry_.template geomVector<0>(),
                                    (*data[level]).geometry_.template geomVector<1>(),
                                    (*data[level]).geometry_.template geomVector<3>());
            BOOST_CHECK( (*data[level]).parent_to_children_cells_.empty());
            BOOST_CHECK(coarse_grid.getLgrNameToLevel().at(lgr_name_vec[level-1]) == static_cast<int>(level));

            auto itMin = std::min_element((data[level] -> global_cell_).begin(),  (data[level] -> global_cell_).end());
            auto itMax = std::max_element((data[level] -> global_cell_).begin(),  (data[level] -> global_cell_).end());
            BOOST_CHECK_EQUAL( *itMin, 0);
            const auto& maxCartesianIdxLevel = data[level]->logical_cartesian_size_[0]*data[level]->logical_cartesian_size_[1]* data[level]->logical_cartesian_size_[2] -1;
            BOOST_CHECK_EQUAL( *itMax, maxCartesianIdxLevel);

            // LGRs
            for (int cell = 0; cell <  data[level]-> size(0); ++cell)
            {
                Dune::cpgrid::Entity<0> entity = Dune::cpgrid::Entity<0>(*data[level], cell, true);
                BOOST_CHECK( entity.hasFather() == true);
                BOOST_CHECK( entity.getOrigin() ==  entity.father());
                BOOST_CHECK( entity.getOrigin().level() == 0);
                BOOST_CHECK_CLOSE(entity.geometryInFather().volume(),
                                  1./(cells_per_dim_vec[level-1][0]*cells_per_dim_vec[level-1][1]*cells_per_dim_vec[level-1][2]), 1e-6);
                BOOST_CHECK(entity.father().level() == 0);
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
                auto it = entity.hbegin(coarse_grid.maxLevel());
                auto endIt = entity.hend(coarse_grid.maxLevel());
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
                Dune::cpgrid::Entity<0> entity = Dune::cpgrid::Entity<0>(*data[startIJK_vec.size()+1], cell, true);
                for (int i = 0; i < data[startIJK_vec.size()+1] -> cell_to_face_[entity].size(); ++i)
                {
                    BOOST_CHECK( data[startIJK_vec.size()+1] -> cell_to_face_[entity][i].index() != -1);
                }
                const auto& child_to_parent = (*data[startIJK_vec.size()+1]).child_to_parent_cells_[cell];
                const auto& level_cellIdx = (*data[startIJK_vec.size()+1]).leaf_to_level_cells_[entity.index()];
                auto it = entity.hbegin(coarse_grid.maxLevel());
                auto endIt = entity.hend(coarse_grid.maxLevel());
                BOOST_CHECK(entity.isLeaf());
                BOOST_CHECK(it == endIt);
                if (entity.hasFather()){
                    BOOST_CHECK_CLOSE(entity.geometryInFather().volume(),
                                      1./(cells_per_dim_vec[entity.level()-1][0]
                                          *cells_per_dim_vec[entity.level()-1][1]
                                          *cells_per_dim_vec[entity.level()-1][2]), 1e-6);
                    BOOST_CHECK(entity.father().level() == 0);
                    BOOST_CHECK_EQUAL( (std::find(all_parent_cell_indices.begin(), all_parent_cell_indices.end(),
                                                  entity.father().index()) == all_parent_cell_indices.end()), false);
                    BOOST_CHECK(child_to_parent[0] != -1);
                    BOOST_CHECK_EQUAL( child_to_parent[0] == 0, true);
                    BOOST_CHECK_EQUAL( child_to_parent[1], entity.father().index());
                    BOOST_CHECK( entity.father() == entity.getOrigin());
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
                    // Get the entity cell_to_face_ on the LeafView
                    const auto& leaf_cell_to_face =
                        (*data[startIJK_vec.size()+1]).cell_to_face_[Dune::cpgrid::EntityRep<0>(entity.index(), true)];
                    int face_count = 0;
                    bool touch_patch_onLeftFace = false; // false -> 0, true -> 1
                    bool touch_patch_onRightFace = false;
                    bool touch_patch_onFrontFace = false;
                    bool touch_patch_onBackFace = false;
                    bool touch_patch_onBottomFace = false;
                    bool touch_patch_onTopFace = false;
                    // CHECK LEFT FACE(S)
                    for (long unsigned int patch = 0; patch < startIJK_vec.size(); ++patch)
                    {
                        // LEFT Check if i == endIJK_vec[patch][0] (and j, k in the range) for some patch.
                        // If true, increase amount of faces by cells_per_dim_vec[patch][1]*cells_per_dim_vec[patch][2]
                        // Auxiliary bool
                        touch_patch_onLeftFace = touch_patch_onLeftFace ||
                            ((entityOldIJK[0] == endIJK_vec[patch][0]) &&
                             (entityOldIJK[1] >= startIJK_vec[patch][1]) && (entityOldIJK[1] < endIJK_vec[patch][1]) &&
                             (entityOldIJK[2] >= startIJK_vec[patch][2]) && (entityOldIJK[2] < endIJK_vec[patch][2]));
                        if (touch_patch_onLeftFace)
                        {
                            face_count += touch_patch_onLeftFace*(cells_per_dim_vec[patch][1]*cells_per_dim_vec[patch][2]);
                            break;
                        }
                    } // end patch-for-loop
                    // CHECK RIGHT FACE(S)
                    for (long unsigned int patch = 0; patch < startIJK_vec.size(); ++patch)
                    {
                        // RIGHT Check if i+1 == startIJK_vec[patch][0] (and j, k in the range) for some patch.
                        // If true, increase amount of faces by cells_per_dim_vec[patch][1]*cells_per_dim_vec[patch][2]
                        // Auxiliary bool
                        touch_patch_onRightFace = touch_patch_onRightFace ||
                            ((entityOldIJK[0]+1 == startIJK_vec[patch][0]) &&
                             (entityOldIJK[1] >= startIJK_vec[patch][1]) && (entityOldIJK[1] < endIJK_vec[patch][1]) &&
                             (entityOldIJK[2] >= startIJK_vec[patch][2]) && (entityOldIJK[2] < endIJK_vec[patch][2]));
                        if (touch_patch_onRightFace)
                        {
                            face_count += touch_patch_onRightFace*(cells_per_dim_vec[patch][1]*cells_per_dim_vec[patch][2]);
                            break;
                        }
                    } // end patch-for-loop
                    // CHECK FRONT FACE(S)
                    for (long unsigned int patch = 0; patch < startIJK_vec.size(); ++patch)
                    {
                        // FRONT Check if j == endIJK_vec[patch][1] (and i, k in the range) for some patch.
                        // If true, increase amount of faces by cells_per_dim_vec[patch][0]*cells_per_dim_vec[patch][2]
                        // Auxiliary bool
                        touch_patch_onFrontFace = touch_patch_onFrontFace ||
                            ((entityOldIJK[1] == endIJK_vec[patch][1]) &&
                             (entityOldIJK[0] >= startIJK_vec[patch][0]) && (entityOldIJK[0] < endIJK_vec[patch][0]) &&
                             (entityOldIJK[2] >= startIJK_vec[patch][2]) && (entityOldIJK[2] < endIJK_vec[patch][2]));
                        if (touch_patch_onFrontFace)
                        {
                            face_count +=  touch_patch_onFrontFace*(cells_per_dim_vec[patch][0]*cells_per_dim_vec[patch][2]);
                            break;
                        }
                    } // end patch-for-loop
                    // CHECK BACK FACE(S)
                    for (long unsigned int patch = 0; patch < startIJK_vec.size(); ++patch)
                    {
                        // BACK Check if j+1 == statIJK_vec[patch][1] (and i, k in the range) for some patch.
                        // If true, increase amount of faces by cells_per_dim_vec[patch][0]*cells_per_dim_vec[patch][2]
                        // Auxiliary bool
                        touch_patch_onBackFace = touch_patch_onBackFace ||
                            ((entityOldIJK[1]+1 == startIJK_vec[patch][1]) &&
                             (entityOldIJK[0] >= startIJK_vec[patch][0]) && (entityOldIJK[0] < endIJK_vec[patch][0]) &&
                             (entityOldIJK[2] >= startIJK_vec[patch][2]) && (entityOldIJK[2] < endIJK_vec[patch][2]));
                        if (touch_patch_onBackFace)
                        {
                            face_count += touch_patch_onBackFace*(cells_per_dim_vec[patch][0]*cells_per_dim_vec[patch][2]);
                            break;
                        }
                    } // end patch-for-loop
                    // CHECK BOTTOM FACE(S)
                    for (long unsigned int patch = 0; patch < startIJK_vec.size(); ++patch)
                    {
                        // BOTTOM Check if k == endIJK_vec[patch][2] (and i, j in the range) for some patch.
                        // If true, increase amount of faces by cells_per_dim_vec[patch][0]*cells_per_dim_vec[patch][1]
                        // Auxiliary bool
                        touch_patch_onBottomFace = touch_patch_onBottomFace ||
                            ((entityOldIJK[2] == endIJK_vec[patch][2]) &&
                             (entityOldIJK[0] >= startIJK_vec[patch][0]) && (entityOldIJK[0] < endIJK_vec[patch][0]) &&
                             (entityOldIJK[1] >= startIJK_vec[patch][1]) && (entityOldIJK[1] < endIJK_vec[patch][1]));
                        if (touch_patch_onBottomFace)
                        {
                            face_count += touch_patch_onBottomFace*(cells_per_dim_vec[patch][0]*cells_per_dim_vec[patch][1]);
                            break;
                        }
                    } // end patch-for-loop
                    // CHECK TOP FACE(S)
                    for (long unsigned int patch = 0; patch < startIJK_vec.size(); ++patch)
                    {
                        // TOP Check if k+1 == startIJK_vec[patch][2] (and i, j in the range) for some patch.
                        // If true, increase amount of faces by cells_per_dim_vec[patch][0]*cells_per_dim_vec[patch][2]
                        // Auxiliary bool
                        touch_patch_onTopFace = touch_patch_onTopFace ||
                            ((entityOldIJK[2]+1 == startIJK_vec[patch][2]) &&
                             (entityOldIJK[0] >= startIJK_vec[patch][0]) && (entityOldIJK[0] < endIJK_vec[patch][0]) &&
                             (entityOldIJK[1] >= startIJK_vec[patch][1]) && (entityOldIJK[1] < endIJK_vec[patch][1]));
                        if (touch_patch_onTopFace)
                        {
                            face_count += touch_patch_onTopFace*(cells_per_dim_vec[patch][0]*cells_per_dim_vec[patch][1]);
                            break;
                        }
                    } // end patch-for-loop
                    if (!touch_patch_onLeftFace)
                    {
                        face_count +=1;
                    }
                    if (!touch_patch_onRightFace)
                    {
                        face_count +=1;
                    }
                    if (!touch_patch_onFrontFace)
                    {
                        face_count +=1;
                    }
                    if (!touch_patch_onBackFace)
                    {
                        face_count +=1;
                    }
                    if (!touch_patch_onBottomFace)
                    {
                        face_count +=1;
                    }
                    if (!touch_patch_onTopFace)
                    {
                        face_count +=1;
                    }
                    BOOST_CHECK( leaf_cell_to_face.size() == face_count);
                } // end else
            }
        } // end-level-for-loop

        BOOST_CHECK( static_cast<int>(startIJK_vec.size()) == coarse_grid.maxLevel());
        BOOST_CHECK( (*data[data.size()-1]).parent_to_children_cells_.empty());

        auto itMin = std::min_element((data.back() -> global_cell_).begin(),  (data.back()-> global_cell_).end());
        auto itMax = std::max_element((data.back() -> global_cell_).begin(),  (data.back() -> global_cell_).end());
        BOOST_CHECK_EQUAL( *itMin, 0);
        const auto& maxCartesianIdx = coarse_grid.logicalCartesianSize()[0]*coarse_grid.logicalCartesianSize()[1]*coarse_grid.logicalCartesianSize()[2] -1;
        BOOST_CHECK_EQUAL( *itMax, maxCartesianIdx);

        for (long unsigned int l = 0; l < startIJK_vec.size() +1; ++l) // level 0,1,2,... , last patch
        {
            const auto& view = coarse_grid.levelGridView(l);
            for (const auto& element: elements(view)){
                BOOST_CHECK_EQUAL(element.level(), l);
            }
        }

        std::vector<int> leaf_to_parent_cell; // To store parent cell index, when leaf cell has a parent. Empty entry otherwise.
        leaf_to_parent_cell.resize(data[startIJK_vec.size()+1]-> size(0)); // Correct size.

        Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LeafGridView> leafMapper(coarse_grid.leafGridView(), Dune::mcmgElementLayout());
        Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LevelGridView> level0Mapper(coarse_grid.levelGridView(0), Dune::mcmgElementLayout());

        for (const auto& element: elements(coarse_grid.leafGridView())){
            BOOST_CHECK( ((element.level() >= 0) || (element.level() < static_cast<int>(startIJK_vec.size()) +1)));
            if (element.hasFather()) { // leaf_cell has a father!
                leaf_to_parent_cell[leafMapper.index(element)] = level0Mapper.index(element.father());
                const auto& parent_id = data[0]->localIdSet().id(element.father());
                BOOST_CHECK(element.index() == leafMapper.index(element));
                BOOST_CHECK(element.father().index() == leaf_to_parent_cell[element.index()]);
                BOOST_CHECK(element.father().index() == parent_id);
                BOOST_CHECK(element.father().index() == level0Mapper.index(element.father()));
            }
        }

        std::set<int> allIds_set;
        std::vector<int> allIds_vec;
        allIds_vec.reserve(data.back()->size(0) + data.back()->size(3));
        for (const auto& element: elements(coarse_grid.leafGridView())){
            const auto& localId = data.back()->localIdSet().id(element);
            const auto& globalId = data.back()->globalIdSet().id(element);
            // In serial run, local and global id coincide:
            BOOST_CHECK_EQUAL(localId, globalId);
            allIds_set.insert(localId);
            allIds_vec.push_back(localId);
            // Check that the global_id_set_ptr_ has the correct id (id from the level where the entity was born).
            BOOST_CHECK_EQUAL( coarse_grid.globalIdSet().id(element), data[element.level()]->localIdSet().id(element.getLevelElem()));
        }
        // Check injectivity of the map local_id_set_ (and, indirectly, global_id_set_) after adding cell ids.
        BOOST_CHECK( allIds_set.size() == allIds_vec.size());

        for (const auto& point : vertices(coarse_grid.leafGridView())){
            const auto& localId = data.back()->localIdSet().id(point);
            const auto& globalId = data.back()->globalIdSet().id(point);
            BOOST_CHECK_EQUAL(localId, globalId);
            allIds_set.insert(localId);
            allIds_vec.push_back(localId);
            // Check that the global_id_set_ptr_ has the correct id (id from the level where the entity was born).
        }
        // Check injectivity of the map local_id_set_ (and, indirectly, global_id_set_) after adding point ids.
        BOOST_CHECK( allIds_set.size() == allIds_vec.size());

        // Local/Global id sets for level grids (level 0, 1, ..., maxLevel)
        for (int level = 0; level < coarse_grid.maxLevel() +1; ++level)
        {
            std::set<int> levelIds_set;
            std::vector<int> levelIds_vec;
            levelIds_vec.reserve(data[level]->size(0) + data[level]->size(3));
            for (const auto& element: elements(coarse_grid.levelGridView(level))){
                const auto& localId = data[level]->localIdSet().id(element);
                const auto& globalId = data[level]->globalIdSet().id(element);
                // In serial run, local and global id coincide:
                BOOST_CHECK_EQUAL(localId, globalId);
                levelIds_set.insert(localId);
                levelIds_vec.push_back(localId);
                // Search in the leaf grid view elements for the element with the same id, if it exists.
                if (auto itIsLeaf = std::find_if( elements(coarse_grid.leafGridView()).begin(),
                                                  elements(coarse_grid.leafGridView()).end(),
                                                  [localId, data](const Dune::cpgrid::Entity<0>& leafElem)
                                                  { return (localId == data.back()->localIdSet().id(leafElem)); });
                    itIsLeaf != elements(coarse_grid.leafGridView()).end()) {
                    BOOST_CHECK( itIsLeaf->getLevelElem() == element);
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

            for (const auto& point : vertices(coarse_grid.levelGridView(level))) {
                const auto& localId = data[level]->localIdSet().id(point);
                const auto& globalId = data[level]->globalIdSet().id(point);
                BOOST_CHECK_EQUAL(localId, globalId);
                levelIds_set.insert(localId);
                levelIds_vec.push_back(localId);
                // Search in the leaf grid view elements for the element with the same id, if it exists.
                if (auto itIsLeaf = std::find_if( vertices(coarse_grid.leafGridView()).begin(),
                                                  vertices(coarse_grid.leafGridView()).end(),
                                                  [localId, data](const Dune::cpgrid::Entity<3>& leafPoint)
                                                  { return (localId == data.back()->localIdSet().id(leafPoint)); });
                    itIsLeaf != vertices(coarse_grid.leafGridView()).end()) {
                    BOOST_CHECK( (*itIsLeaf).geometry().center() == point.geometry().center() );
                }
            }
            // Check injectivity of the map local_id_set_ (and, indirectly, global_id_set_)
            BOOST_CHECK( levelIds_set.size() == levelIds_vec.size());
        }
    }
    catch (const std::exception& e) {
        std::cout << e.what() << "'\n";
    }
}


BOOST_AUTO_TEST_CASE(refine_patch_different_cell_sizes)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {2.0, 1.0, 4.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {2,2,2};
    const std::array<int, 3> startIJK = {1,0,1};
    const std::array<int, 3> endIJK = {3,2,3};  // patch_dim = {3-1, 2-0, 3-1} ={2,2,2}
    const std::string lgr_name = {"LGR1"};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    refinePatch_and_check(coarse_grid, {cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});
    BOOST_CHECK_EQUAL(coarse_grid.currentData()[0]->patchesShareFace({startIJK}, {endIJK}), false);
}

BOOST_AUTO_TEST_CASE(refine_patch)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {2,2,2};
    const std::array<int, 3> startIJK = {1,0,1};
    const std::array<int, 3> endIJK = {3,2,3};  // patch_dim = {3-1, 2-0, 3-1} ={2,2,2}
    const std::string lgr_name = {"LGR1"};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    refinePatch_and_check(coarse_grid, {cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});
    BOOST_CHECK_EQUAL(coarse_grid.currentData()[0]->patchesShareFace({startIJK}, {endIJK}), false);
}

BOOST_AUTO_TEST_CASE(refine_patch_one_cell)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {2,2,2};
    const std::array<int, 3> startIJK = {1,0,1};
    const std::array<int, 3> endIJK = {2,1,2};  // patch_dim = {2-1, 1-0, 2-1} ={1,1,1} -> Single Cell!
    const std::string lgr_name = {"LGR1"};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    refinePatch_and_check(coarse_grid, {cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});
    BOOST_CHECK_EQUAL(coarse_grid.currentData()[0]->patchesShareFace({startIJK}, {endIJK}), false);
}

BOOST_AUTO_TEST_CASE(lgrs_disjointPatches)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}, {4,4,4}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,2}, {3,2,2}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,1,1}, {1,1,3}, {4,3,3}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3"};
    refinePatch_and_check(coarse_grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
    BOOST_CHECK_EQUAL(coarse_grid.currentData()[0]->patchesShareFace(startIJK_vec, endIJK_vec), false);
}

BOOST_AUTO_TEST_CASE(lgrs_disjointPatchesB)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {3,2,0}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,2,1}, {4,3,3}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    refinePatch_and_check(coarse_grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
    BOOST_CHECK_EQUAL(coarse_grid.currentData()[0]->patchesShareFace(startIJK_vec, endIJK_vec), false);
}

BOOST_AUTO_TEST_CASE(patches_share_corner)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}, {4,4,4}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {1,1,1}, {3,2,2}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,1}, {2,2,2}, {4,3,3}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3"};
    refinePatch_and_check(coarse_grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
    BOOST_CHECK_EQUAL(coarse_grid.currentData()[0]->patchesShareFace(startIJK_vec, endIJK_vec), false);
}

BOOST_AUTO_TEST_CASE(patches_share_edge)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {2,2,2}, {2,2,2}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {2,0,1}, {3,2,2}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,1,1}, {3,1,2}, {4,3,3}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3"};
    refinePatch_and_check(coarse_grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
    BOOST_CHECK_EQUAL(coarse_grid.currentData()[0]->patchesShareFace(startIJK_vec, endIJK_vec), false);
}

BOOST_AUTO_TEST_CASE(patches_share_face)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {2,2,2}, {2,2,2}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {2,0,0}, {3,2,2}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,1,1}, {3,1,1}, {4,3,3}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3"};
    refinePatch_and_check(coarse_grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
}

BOOST_AUTO_TEST_CASE(patches_share_faceB)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{4,2,2}, {3,2,2}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {1,0,0}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,1}, {2,1,1}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    refinePatch_and_check(coarse_grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
}

BOOST_AUTO_TEST_CASE(patches_share_Iface_with_diff_cells_per_dim)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{4,2,5}, {3,2,2}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {1,0,0}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,1}, {2,1,1}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    refinePatch_and_check(coarse_grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
}

BOOST_AUTO_TEST_CASE(patches_share_Jface_with_diff_cells_per_dim)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{4,2,5}, {3,2,2}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,1,0}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,1}, {1,2,1}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    refinePatch_and_check(coarse_grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
}

BOOST_AUTO_TEST_CASE(patches_share_Kface_with_diff_cells_per_dim)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{4,2,5}, {3,2,2}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,1}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,1}, {1,1,2}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    refinePatch_and_check(coarse_grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
}


void check_global_refine(const Dune::CpGrid& refined_grid, const Dune::CpGrid& equiv_fine_grid)
{

    const auto& refined_data = refined_grid.currentData();
    const auto& equiv_data = equiv_fine_grid.currentData();

    BOOST_CHECK(refined_data.size()==3);

    const auto& refined_leaf = *refined_data[2];
    const auto& equiv_leaf = *equiv_data[0];

    // Check the container sizes
    BOOST_CHECK_EQUAL(refined_leaf.face_to_cell_.size(), equiv_leaf.face_to_cell_.size());
    BOOST_CHECK_EQUAL(refined_leaf.face_to_point_.size(), equiv_leaf.face_to_point_.size());
    BOOST_CHECK_EQUAL(refined_leaf.face_normals_.size(), equiv_leaf.face_normals_.size());
    BOOST_CHECK_EQUAL(refined_leaf.face_tag_.size(), equiv_leaf.face_tag_.size());
    BOOST_CHECK_EQUAL(refined_leaf.cell_to_point_.size(), equiv_leaf.cell_to_point_.size());
    BOOST_CHECK_EQUAL(refined_leaf.cell_to_face_.size(), equiv_leaf.cell_to_face_.size());
    BOOST_CHECK_EQUAL(refined_leaf.geomVector<3>().size(), equiv_leaf.geomVector<3>().size());

    BOOST_CHECK_EQUAL(refined_grid.size(3), equiv_fine_grid.size(3));
    BOOST_CHECK_EQUAL(refined_grid.size(0), equiv_fine_grid.size(0));

    for(const auto& point: refined_leaf.geomVector<3>())
    {
        auto equiv_point_iter = equiv_leaf.geomVector<3>().begin();
        while ((equiv_point_iter != equiv_leaf.geomVector<3>().end()) && (point.center() != equiv_point_iter->center())) {
            ++equiv_point_iter;
        }
        CHECK_COORDINATES(point.center(), equiv_point_iter->center());
        for(const auto& coord: point.center())
            BOOST_TEST(std::isfinite(coord));

    }
    for(const auto& cell: refined_leaf.geomVector<3>())
    {
        auto equiv_cell_iter = equiv_leaf.geomVector<3>().begin();
        while ((equiv_cell_iter != equiv_leaf.geomVector<3>().end()) && (cell.center() != equiv_cell_iter->center())) {
            ++equiv_cell_iter;
        }
        CHECK_COORDINATES(cell.center(), equiv_cell_iter->center());
        for(const auto& coord: cell.center())
            BOOST_TEST(std::isfinite(coord));
        BOOST_CHECK_CLOSE(cell.volume(), equiv_cell_iter->volume(), 1e-6);
    }

    /////
    const auto& grid_view = refined_grid.leafGridView();
    const auto& equiv_grid_view = equiv_fine_grid.leafGridView();


    for(const auto& element: elements(grid_view))
    {
        BOOST_CHECK( element.getOrigin().level() == 0);
        auto equiv_element_iter = equiv_grid_view.begin<0>();
        bool closedCenter =  (std::abs(element.geometry().center()[0] - equiv_element_iter->geometry().center()[0]) < 1e-12) &&
            (std::abs(element.geometry().center()[1] - equiv_element_iter->geometry().center()[1]) < 1e-12) &&
            (std::abs(element.geometry().center()[2] - equiv_element_iter->geometry().center()[2])< 1e-12);

        while ((equiv_element_iter != equiv_grid_view.end<0>()) && (!closedCenter)) {
            ++equiv_element_iter;
            closedCenter = (std::abs(element.geometry().center()[0] - equiv_element_iter->geometry().center()[0]) < 1e-12) &&
                (std::abs(element.geometry().center()[1] - equiv_element_iter->geometry().center()[1]) < 1e-12) &&
                (std::abs(element.geometry().center()[2] - equiv_element_iter->geometry().center()[2])< 1e-12);
        }
        for(const auto& intersection: intersections(grid_view, element)) {
            // find matching intersection (needed as ordering is allowed to be different
            bool matching_intersection_found = false;
            for(auto& intersection_match: intersections(equiv_grid_view, *equiv_element_iter)) {
                if(intersection_match.indexInInside() == intersection.indexInInside()) {
                    BOOST_CHECK(intersection_match.neighbor() == intersection.neighbor());

                    if(intersection.neighbor()) {
                        BOOST_CHECK(intersection_match.indexInOutside() == intersection.indexInOutside());
                    }

                    CHECK_COORDINATES(intersection_match.centerUnitOuterNormal(), intersection.centerUnitOuterNormal());
                    const auto& geom_match = intersection_match.geometry();
                    BOOST_TEST(0.0 == 1e-11, boost::test_tools::tolerance(1e-8));
                    const auto& geom =  intersection.geometry();
                    BOOST_CHECK_CLOSE(geom_match.volume(), geom.volume(), 1e-13);
                    CHECK_COORDINATES(geom_match.center(), geom.center());
                    BOOST_CHECK(geom_match.corners() == geom.corners());

                    decltype(geom.corner(0)) sum_match{}, sum{};

                    for(int cor = 0; cor < geom.corners(); ++cor) {
                        sum += geom.corner(cor);
                        sum_match += geom_match.corner(1);
                    }
                    CHECK_COORDINATES(sum, sum_match);
                    matching_intersection_found = true;
                    break;
                }
            } // end-for-loop-intersection_match
            BOOST_CHECK(matching_intersection_found);
        }
    }
}


BOOST_AUTO_TEST_CASE(global_refine)
{
    // Create a 4x3x3 grid with length 4x3x3
    // and refine each cells into 4 children cells
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {2,2,2};
    const std::array<int, 3> startIJK = {0,0,0};
    const std::array<int, 3> endIJK = {4,3,3};
    const std::string lgr_name = "LGR";
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    coarse_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

    // Create a 8x6x6 grid with length 4x3x3
    Dune::CpGrid fine_grid;
    const std::array<double, 3> fine_cell_sizes = {0.5, 0.5, 0.5};
    const std::array<int, 3> fine_grid_dim = {8,6,6};
    fine_grid.createCartesian(fine_grid_dim, fine_cell_sizes);

    check_global_refine(coarse_grid, fine_grid);
}


BOOST_AUTO_TEST_CASE(global_norefine)
{
    // Create a 4x3x3 grid with length 4x3x3
    // and refine each cells into 4 children cells
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {1,1,1};
    const std::array<int, 3> startIJK = {0,0,0};
    const std::array<int, 3> endIJK = {4,3,3};
    const std::string lgr_name = "LGR";
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    coarse_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

    // Create a 4x3x3 grid with length 4x3x3
    Dune::CpGrid fine_grid;
    const std::array<double, 3> fine_cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> fine_grid_dim = {4,3,3};
    fine_grid.createCartesian(fine_grid_dim, fine_cell_sizes);

    check_global_refine(coarse_grid, fine_grid);
}
