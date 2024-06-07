//===========================================================================
//
// File: adapt_cpgrid_test.cpp
//
// Created: Apr 22 2024 16:45
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

#define BOOST_TEST_MODULE AdaptTests
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

void markAndAdapt_check(Dune::CpGrid& coarse_grid,
                        const std::array<int,3>& cells_per_dim,
                        const std::vector<int>& markedCells,
                        Dune::CpGrid& other_grid,
                        bool isBlockShape,
                        bool hasBeenRefinedAtLeastOnce,
                        bool isGlobalRefinement)
{
    const int startingGridIdx = coarse_grid.chooseData().size() -1; // size before calling adapt

    std::vector<int> assignRefinedLevel(coarse_grid.chooseData()[startingGridIdx]->size(0));

    for (const auto& elemIdx : markedCells)
    {
        const auto& elem =  Dune::cpgrid::Entity<0>(*(coarse_grid.chooseData()[startingGridIdx]), elemIdx, true);
        coarse_grid.mark(1, elem);
        coarse_grid.getMark(elem);
        assignRefinedLevel[elemIdx] = coarse_grid.maxLevel() + 1;
        BOOST_CHECK( coarse_grid.getMark(elem) == 1);
        BOOST_CHECK( elem.mightVanish() == true);
    }
    bool preAdapt = coarse_grid.preAdapt();
    const auto& data = coarse_grid.chooseData();
    if(preAdapt) {
        coarse_grid.adapt({cells_per_dim}, assignRefinedLevel, {"LGR"});
        coarse_grid.postAdapt();
        BOOST_CHECK(static_cast<int>(data.size()) == startingGridIdx+3);
        const auto& adapted_leaf = *data[startingGridIdx+2];

        if(isBlockShape) {
            const auto& blockRefinement_data = other_grid.chooseData();
            const auto& blockRefinement_leaf = *blockRefinement_data.back();

            // Check the container sizes
            BOOST_CHECK_EQUAL(adapted_leaf.geomVector<3>().size(), blockRefinement_leaf.geomVector<3>().size());
            BOOST_CHECK_EQUAL(adapted_leaf.face_to_cell_.size(), blockRefinement_leaf.face_to_cell_.size());
            BOOST_CHECK_EQUAL(adapted_leaf.face_to_point_.size(), blockRefinement_leaf.face_to_point_.size());
            BOOST_CHECK_EQUAL(adapted_leaf.face_normals_.size(), blockRefinement_leaf.face_normals_.size());
            BOOST_CHECK_EQUAL(adapted_leaf.face_tag_.size(), blockRefinement_leaf.face_tag_.size());
            BOOST_CHECK_EQUAL(adapted_leaf.cell_to_point_.size(), blockRefinement_leaf.cell_to_point_.size());
            BOOST_CHECK_EQUAL(adapted_leaf.cell_to_face_.size(), blockRefinement_leaf.cell_to_face_.size());
            BOOST_CHECK_EQUAL(coarse_grid.size(3), other_grid.size(3));
            BOOST_CHECK_EQUAL(coarse_grid.size(0), other_grid.size(0));
            /** Checking amount of corners and cells in the refined grid fails fails. To be fixed.*/
            /*if(!hasBeenRefinedAtLeastOnce) {
              BOOST_CHECK_EQUAL(coarse_grid.size(1,0), other_grid.size(1,0)); // equal amount of cells in level 1
              BOOST_CHECK_EQUAL(coarse_grid.size(1,3), other_grid.size(1,3)); // equal amount of corners in level 1
              }*/

            for(const auto& point: adapted_leaf.geomVector<3>()){
                auto equiv_point_iter = blockRefinement_leaf.geomVector<3>().begin();
                while (point.center() != equiv_point_iter->center()) {
                    ++equiv_point_iter;
                }
                CHECK_COORDINATES(point.center(), equiv_point_iter->center());
                for(const auto& coord: point.center())
                    BOOST_TEST(std::isfinite(coord));

            }
            for(const auto& cell: adapted_leaf.geomVector<3>()) {
                auto equiv_cell_iter = blockRefinement_leaf.geomVector<3>().begin();
                while (cell.center() != equiv_cell_iter->center()) {
                    ++equiv_cell_iter;
                }
                CHECK_COORDINATES(cell.center(), equiv_cell_iter->center());
                for(const auto& coord: cell.center())
                    BOOST_TEST(std::isfinite(coord));
                BOOST_CHECK_CLOSE(cell.volume(), equiv_cell_iter->volume(), 1e-6);
            }

            /////  THE FOLLOWING CODE WOULD FIT FOR TESTING GLOBAL REFINEMENT
            if (isGlobalRefinement) {
                const auto& grid_view = coarse_grid.leafGridView();
                const auto& equiv_grid_view = other_grid.leafGridView();


                for(const auto& element: elements(grid_view)) {
                    BOOST_CHECK( element.getOrigin().level() == 0);
                    auto equiv_element_iter = equiv_grid_view.begin<0>();
                    bool closedCenter =  (std::abs(element.geometry().center()[0] - equiv_element_iter->geometry().center()[0]) < 1e-12) &&
                        (std::abs(element.geometry().center()[1] - equiv_element_iter->geometry().center()[1]) < 1e-12) &&
                        (std::abs(element.geometry().center()[2] - equiv_element_iter->geometry().center()[2])< 1e-12);

                    while (!closedCenter) {
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
                                bool closedGeomCenter =  (std::abs(geom_match.center()[0] - geom.center()[0]) < 1e-12) &&
                                    (std::abs(geom_match.center()[1] - geom.center()[1]) < 1e-12) &&
                                    (std::abs(geom_match.center()[2] - geom.center()[2])< 1e-12);
                                if (!closedGeomCenter) {
                                    break; // Check next intersection_match
                                }
                                BOOST_CHECK_CLOSE(geom_match.volume(), geom.volume(), 1e-6);
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
            } // end-if-isGlobalRefinement
        } // end-if-isBlockShape

        const auto& grid_view = coarse_grid.leafGridView();

        const auto& preAdapt_view = coarse_grid.levelGridView(startingGridIdx);
        // Note: preAdapt grid in level "startingGridIdx"
        //       refined grid in level  "startingGridIdx+1"
        //       adapted grid in level  "startingGridIdx+2"
        Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LeafGridView> adaptMapper(grid_view, Dune::mcmgElementLayout());
        Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LevelGridView> preAdaptMapper(preAdapt_view, Dune::mcmgElementLayout());
        const auto& adapt_idSet = adapted_leaf.local_id_set_;
        const auto& preAdapt_idSet = (*data[startingGridIdx+1]).local_id_set_;

        for(const auto& element: elements(grid_view)) {
            // postAdapt() has been called, therefore every element gets marked with 0
            BOOST_CHECK( coarse_grid.getMark(element) == 0);
            BOOST_CHECK(  adapted_leaf.cell_to_point_[element.index()].size() == 8);
            for (int i = 0; i < 8; ++i) {
                BOOST_CHECK(  adapted_leaf.cell_to_point_[element.index()][i] != -1);
            }
            for (int i = 0; i <  adapted_leaf.cell_to_face_[element].size(); ++i) {
                BOOST_CHECK(  adapted_leaf.cell_to_face_[element][i].index() != -1);
            }
            const auto& child_to_parent =  adapted_leaf.child_to_parent_cells_[element.index()];
            const auto& level_cellIdx =  adapted_leaf.leaf_to_level_cells_[element.index()];
            auto it = element.hbegin(level_cellIdx[0] /*level*/);//coarse_grid.maxLevel());
            auto endIt = element.hend(level_cellIdx[0] /*level*/); //coarse_grid.maxLevel());
            BOOST_CHECK(element.isLeaf());
            BOOST_CHECK(it == endIt);
            if (element.hasFather()){
                BOOST_CHECK( element.isNew() == true);
                BOOST_CHECK_CLOSE(element.geometryInFather().volume(), 1./(cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]), 1e-6);
                if (hasBeenRefinedAtLeastOnce){
                    BOOST_CHECK(element.father().level() <= startingGridIdx); // Remove? Not that useful...
                    BOOST_CHECK( element.getOrigin().level() <= startingGridIdx);
                }
                else {
                    BOOST_CHECK(element.father().level() == 0);
                    BOOST_CHECK( element.getOrigin().level() == 0);  // To do: check Entity::getOrigin()
                }
                BOOST_CHECK_EQUAL( (std::find(markedCells.begin(), markedCells.end(), element.father().index()) == markedCells.end()), false);
                BOOST_CHECK(child_to_parent[0] != -1);
                BOOST_CHECK_EQUAL( child_to_parent[0] == startingGridIdx, true);
                BOOST_CHECK_EQUAL( child_to_parent[1], element.father().index());
                BOOST_CHECK( element.father() == element.getOrigin());
                BOOST_CHECK(  ( adapted_leaf.global_cell_[element.index()]) == (data[startingGridIdx]->global_cell_[element.getOrigin().index()]) );
                BOOST_CHECK( std::get<0>(data[startingGridIdx]->parent_to_children_cells_[child_to_parent[1]]) == element.level());
                // Check amount of children cells of the parent cell
                BOOST_CHECK_EQUAL(std::get<1>(data[startingGridIdx]->parent_to_children_cells_[child_to_parent[1]]).size(),
                                  cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]);
                BOOST_CHECK( element.father().isLeaf() == false);
                BOOST_CHECK( (element.level() > 0) || (element.level() < coarse_grid.maxLevel() +1));
                BOOST_CHECK( level_cellIdx[0] == element.level());
                //
                const auto& id = (*adapt_idSet).id(element);
                const auto& parent_id = (*preAdapt_idSet).id(element.father());
                BOOST_CHECK(element.index() == id);
                BOOST_CHECK(element.index() == adaptMapper.index(element));
                BOOST_CHECK(element.father().index() == parent_id);
                BOOST_CHECK(element.father().index() == preAdaptMapper.index(element.father()));
            }
            else{
                BOOST_CHECK_THROW(element.father(), std::logic_error);
                BOOST_CHECK_THROW(element.geometryInFather(), std::logic_error);
                BOOST_CHECK_EQUAL(child_to_parent[0], -1);
                BOOST_CHECK_EQUAL(child_to_parent[1], -1);
                if (hasBeenRefinedAtLeastOnce){
                    BOOST_CHECK( level_cellIdx[0] == startingGridIdx); // Equal when the grid has been refined only once. Remove this check?
                    BOOST_CHECK( element.level() == startingGridIdx);
                    BOOST_CHECK( element.getOrigin().level() <= startingGridIdx);
                }
                else  {
                    BOOST_CHECK( level_cellIdx[0] == 0);
                    BOOST_CHECK( element.level() == 0);
                    BOOST_CHECK( element.getOrigin().level() == 0);
                }
                BOOST_CHECK( std::get<0>(data[startingGridIdx]-> parent_to_children_cells_[level_cellIdx[1]]) == -1);
                BOOST_CHECK( std::get<1>(data[startingGridIdx]->parent_to_children_cells_[level_cellIdx[1]]).empty());
                // Get index of the cell in level 0
                const auto& entityOldIdx =   adapted_leaf.leaf_to_level_cells_[element.index()][1];
                BOOST_CHECK( element.getOrigin().index() == entityOldIdx);
                BOOST_CHECK( element.getOrigin().level() == 0);
                BOOST_CHECK( element.isNew() == false);
            }
            BOOST_CHECK( element.mightVanish() == false); // marks get rewrtitten and set to 0 via postAdapt call
        } // end-element-for-loop

        // Some checks on the preAdapt grid
        for(const auto& element: elements(preAdapt_view)) {
            if (!hasBeenRefinedAtLeastOnce){
                BOOST_CHECK( element.hasFather() == false);
                BOOST_CHECK_THROW(element.father(), std::logic_error);
                BOOST_CHECK_THROW(element.geometryInFather(), std::logic_error);
                BOOST_CHECK( element.getOrigin() ==  element);
                BOOST_CHECK( element.getOrigin().level() == startingGridIdx);
                BOOST_CHECK( element.isNew() == false);
            }
            auto it = element.hbegin(coarse_grid.maxLevel()); // With element.level(), fails
            auto endIt = element.hend(coarse_grid.maxLevel());
            const auto& [lgr, childrenList] = (*data[startingGridIdx]).parent_to_children_cells_[element.index()];
            if (std::find(markedCells.begin(), markedCells.end(), element.index()) == markedCells.end()){
                BOOST_CHECK_EQUAL(lgr, -1);
                BOOST_CHECK(childrenList.empty());
                BOOST_CHECK( element.isLeaf() == true);
                // If it == endIt, then entity.isLeaf() true (when dristibuted_data_ is empty)
                BOOST_CHECK( it == endIt);
                BOOST_CHECK( element.mightVanish() == false);
            }
            else{
                BOOST_CHECK(lgr != -1);
                BOOST_CHECK(static_cast<int>(childrenList.size()) == cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]);
                // If it != endIt, then entity.isLeaf() false (when dristibuted_data_ is empty)
                BOOST_CHECK_EQUAL( it == endIt, false);
                BOOST_CHECK( element.mightVanish() == true);
                BOOST_CHECK( element.isNew() == false);
                BOOST_CHECK_EQUAL( element.isLeaf(), false); // parent cells do not appear in the LeafView
                // Auxiliary int to check amount of children
                double referenceElemOneParent_volume = 0.;
                std::array<double,3> referenceElem_entity_center = {0.,0.,0.}; // Expected {.5,.5,.5}
                for (const auto& child : childrenList) {
                    BOOST_CHECK( child != -1);
                    BOOST_CHECK( data[startingGridIdx+1]-> child_to_parent_cells_[child][0] == startingGridIdx);  //
                    BOOST_CHECK( data[startingGridIdx+1]-> child_to_parent_cells_[child][1] == element.index());

                    const auto& childElem =  Dune::cpgrid::Entity<0>(*data[startingGridIdx+1], child, true);
                    BOOST_CHECK(childElem.hasFather() == true);
                    BOOST_CHECK(childElem.level() == lgr);
                    referenceElemOneParent_volume += childElem.geometryInFather().volume();
                    for (int c = 0; c < 3; ++c)  {
                        referenceElem_entity_center[c] += (childElem.geometryInFather().center())[c];
                    }
                }
                /// Auxiliary int to check amount of children
                double referenceElemOneParent_volume_it = 0.;
                std::array<double,3> referenceElem_entity_center_it = {0.,0.,0.}; // Expected {.5,.5,.5}
                for (; it != endIt; ++it)
                {
                    BOOST_CHECK(it ->hasFather() == true);
                    BOOST_CHECK(it ->level() == lgr);
                    referenceElemOneParent_volume_it += it-> geometryInFather().volume();
                    for (int c = 0; c < 3; ++c)
                    {
                        referenceElem_entity_center_it[c] += (it-> geometryInFather().center())[c];
                    }
                }
                for (int c = 0; c < 3; ++c) {
                    referenceElem_entity_center[c] /= cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2];
                    referenceElem_entity_center_it[c] /= cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2];
                }
                BOOST_CHECK_CLOSE(referenceElemOneParent_volume, 1, 1e-6);
                BOOST_CHECK_CLOSE(referenceElem_entity_center[0], .5, 1e-6);
                BOOST_CHECK_CLOSE(referenceElem_entity_center[1], .5, 1e-6);
                BOOST_CHECK_CLOSE(referenceElem_entity_center[2], .5, 1e-6);

                BOOST_CHECK_CLOSE(referenceElemOneParent_volume_it, 1, 1e-6);
                BOOST_CHECK_CLOSE(referenceElem_entity_center_it[0], .5, 1e-6);
                BOOST_CHECK_CLOSE(referenceElem_entity_center_it[1], .5, 1e-6);
                BOOST_CHECK_CLOSE(referenceElem_entity_center_it[2], .5, 1e-6);
            }
            BOOST_CHECK( element.level() == 0);
        } // end-preAdaptElements-for-loop
    } // end-if-preAdapt
}

BOOST_AUTO_TEST_CASE(doNothing)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {2,2,2};
    std::vector<int> markedCells;
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, false, false, false);
}

BOOST_AUTO_TEST_CASE(globalRefinement)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {2,2,2};
    std::vector<int> markedCells(36);
    std::iota(markedCells.begin(), markedCells.end(), 0);
    coarse_grid.createCartesian(grid_dim, cell_sizes);

    // Create other grid for comparison
    Dune::CpGrid other_grid;
    other_grid.createCartesian(grid_dim, cell_sizes);

    const std::array<int, 3> startIJK = {0,0,0};
    const std::array<int, 3> endIJK = {4,3,3};
    const std::string lgr_name = {"LGR1"};
    other_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, other_grid, true, false, true);
}

BOOST_AUTO_TEST_CASE(mark2consequtiveCells)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {2,2,2};
    coarse_grid.createCartesian(grid_dim, cell_sizes);


    // Create other grid for comparison
    Dune::CpGrid other_grid;
    other_grid.createCartesian(grid_dim, cell_sizes);

    const std::array<int, 3> startIJK = {2,0,0};
    const std::array<int, 3> endIJK = {4,1,1};  // -> marked elements 2 and 3
    const std::string lgr_name = {"LGR1"};
    other_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

    std::vector<int> markedCells = {2,3};
    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, other_grid, true, false, false);
}

BOOST_AUTO_TEST_CASE(mark2InteriorConsequtiveCells)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {2,2,2};
    coarse_grid.createCartesian(grid_dim, cell_sizes);


    // Create other grid for comparison
    Dune::CpGrid other_grid;
    other_grid.createCartesian(grid_dim, cell_sizes);

    const std::array<int, 3> startIJK = {1,1,1};
    const std::array<int, 3> endIJK = {3,2,2};  // -> marked elements 17 and 18
    const std::string lgr_name = {"LGR1"};
    other_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

    std::vector<int> markedCells = {17,18};
    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, other_grid, true, false, false);
}

BOOST_AUTO_TEST_CASE(markNonBlockShapeCells)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {2,2,2};
    std::vector<int> markedCells = {0}; //,1,2,5,13};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, false, false, false);
}


BOOST_AUTO_TEST_CASE(markNonBlockShapeCells_II)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {2,3,4};
    std::vector<int> markedCells = {1,4,6,9,17,22,28,32,33};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, false, false, false);
}

/*BOOST_AUTO_TEST_CASE(adaptFromAMixedGrid)
  {
  // Create a grid
  Dune::CpGrid coarse_grid;
  const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
  const std::array<int, 3> grid_dim = {4,3,3};
  const std::array<int, 3> cells_per_dim = {2,2,2};
  coarse_grid.createCartesian(grid_dim, cell_sizes);


  const std::array<int, 3> startIJK = {1,1,1};
  const std::array<int, 3> endIJK = {3,2,2};  // -> marked elements 17 and 18
  const std::string lgr_name = {"LGR1"};
  coarse_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

  std::vector<int> markedCells = {0,1}; // coarse cells
  markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, true, true, false);
  }

  BOOST_AUTO_TEST_CASE(adaptFromAMixedGridRefinedCell)
  {
  // Create a grid
  Dune::CpGrid coarse_grid;
  const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
  const std::array<int, 3> grid_dim = {4,3,3};
  const std::array<int, 3> cells_per_dim = {2,2,2};
  coarse_grid.createCartesian(grid_dim, cell_sizes);


  const std::array<int, 3> startIJK = {1,1,1};
  const std::array<int, 3> endIJK = {3,2,2};  // -> marked elements 17 and 18
  const std::string lgr_name = {"LGR1"};
  coarse_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

  std::vector<int> markedCells = {34}; // {34, 35} refined cells (with parent cell 17) -> check conversion of corners between neighboring lgrs definition!
  // Error:  Cannot convert corner index from one LGR to its neighboring LGR
  // 42 refined cell with parent cell 18.
  markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, true, true, false);
  }


  BOOST_AUTO_TEST_CASE(adaptFromAMixedGridMixedCells)
  {
  // Create a grid
  Dune::CpGrid coarse_grid;
  const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
  const std::array<int, 3> grid_dim = {4,3,3};
  const std::array<int, 3> cells_per_dim = {2,2,2};
  coarse_grid.createCartesian(grid_dim, cell_sizes);


  const std::array<int, 3> startIJK = {1,1,1};
  const std::array<int, 3> endIJK = {3,2,2};  // -> marked elements 17 and 18
  const std::string lgr_name = {"LGR1"};
  coarse_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

  std::vector<int> markedCells = {2,3,34}; // {34, 41} refined cells (with parent cell 17),
  // 42, 49 refined cell with parent cell 18-> check conversion of corners between neighboring lgrs definition!
  markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, false, true, false);
  }*/
