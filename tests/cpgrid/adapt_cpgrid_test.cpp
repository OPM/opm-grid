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

void markAndAdapt_check(Dune::CpGrid& grid,
                        const std::array<int,3>& cells_per_dim,
                        const std::vector<int>& markedCells,
                        Dune::CpGrid& other_grid,
                        bool isBlockShape,
                        bool hasBeenRefinedAtLeastOnce,
                        bool isGlobalRefinement)
{
    const int startingGridIdx = grid.currentData().size() -1; // size before calling adapt
    std::vector<int> assignRefinedLevel(grid.currentData()[startingGridIdx]->size(0));

    for (const auto& elemIdx : markedCells)
    {
        const auto& elem =  Dune::cpgrid::Entity<0>(*(grid.currentData()[startingGridIdx]), elemIdx, true);
        grid.mark(1, elem);
        assignRefinedLevel[elemIdx] = grid.maxLevel() + 1;
        BOOST_CHECK( grid.getMark(elem) == 1);
        BOOST_CHECK( elem.mightVanish() == true);
    }

    const auto& data = grid.currentData();

    if(grid.preAdapt()) {

        grid.adapt({cells_per_dim}, assignRefinedLevel, {"LGR"+std::to_string(grid.maxLevel() +1)});
        grid.postAdapt();

        BOOST_CHECK(static_cast<int>(data.size()) == grid.maxLevel() +2);

        const auto& adapted_leaf = *data[grid.maxLevel()+1];
        const auto& grid_view = grid.leafGridView();

        if(isBlockShape) { // For a mixed grid that gets refined a second time, isBlockShape == false, even though the marked elements form a block.
            const auto& blockRefinement_data = other_grid.currentData();
            const auto& blockRefinement_leaf = *blockRefinement_data.back();

            // Check the container sizes
            BOOST_CHECK_EQUAL(adapted_leaf.numFaces(), blockRefinement_leaf.numFaces());
            BOOST_CHECK_EQUAL(grid.size(3), other_grid.size(3));
            BOOST_CHECK_EQUAL(grid.size(0), other_grid.size(0));
            BOOST_CHECK_EQUAL(grid.size(1,0), other_grid.size(1,0)); // equal amount of cells in level 1
            BOOST_CHECK_EQUAL(grid.size(1,3), other_grid.size(1,3)); // equal amount of corners in level 1


            for(const auto& point: adapted_leaf.geomVector<3>()){
                auto equiv_point_iter = blockRefinement_leaf.geomVector<3>().begin();
                while ((equiv_point_iter != blockRefinement_leaf.geomVector<3>().end()) && (point.center() != equiv_point_iter->center())) {
                    ++equiv_point_iter;
                }
                CHECK_COORDINATES(point.center(), equiv_point_iter->center());
                for(const auto& coord: point.center())
                    BOOST_TEST(std::isfinite(coord));

            }
            for(const auto& cell: adapted_leaf.geomVector<3>()) {
                auto equiv_cell_iter = blockRefinement_leaf.geomVector<3>().begin();
                while ((equiv_cell_iter != blockRefinement_leaf.geomVector<3>().end()) && (cell.center() != equiv_cell_iter->center())) {
                    ++equiv_cell_iter;
                }
                CHECK_COORDINATES(cell.center(), equiv_cell_iter->center());
                for(const auto& coord: cell.center())
                    BOOST_TEST(std::isfinite(coord));
                BOOST_CHECK_CLOSE(cell.volume(), equiv_cell_iter->volume(), 1e-24);
            }

            /////  THE FOLLOWING CODE WOULD FIT FOR TESTING GLOBAL REFINEMENT
            if (isGlobalRefinement) {
                const auto& equiv_grid_view = other_grid.leafGridView();

                for(const auto& element: elements(grid_view)) {
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

        Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LeafGridView> adaptMapper(grid_view, Dune::mcmgElementLayout());

        for(const auto& element: elements(grid_view)) {
            BOOST_CHECK( grid.getMark(element) == 0); // postAdapt() has been called, therefore every element gets marked with 0

            // Check element has 8 different corners
            BOOST_CHECK( element.geometry().corners() == 8); // Geometry::corners() always return 8.
            for (int i = 0; i < 8; ++i){
                for (int j = i+1; j < 8; ++j){
                    BOOST_CHECK( element.geometry().corner(i) != element.geometry().corner(j));
                }
            }
            // Check intersections of element have valid ids.
            for (const auto& intersection : intersections(grid_view, element))
            {
                BOOST_CHECK( intersection.id() > -1);
            }

            auto it = element.hbegin(grid.maxLevel());
            auto endIt = element.hend(grid.maxLevel());
            BOOST_CHECK(element.isLeaf());
            BOOST_CHECK(it == endIt);
            if (element.hasFather()){
                BOOST_CHECK( element.isNew() == true);
                BOOST_CHECK_CLOSE(element.geometryInFather().volume(), 1./(cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]), 1e-6);
                if (hasBeenRefinedAtLeastOnce){
                    BOOST_CHECK( element.father().level() <= startingGridIdx);
                    BOOST_CHECK( element.getOrigin().level() <= startingGridIdx);
                    BOOST_CHECK_EQUAL( data[element.father().level()] ->getMark(element.father()), 1);
                }
                else {
                    BOOST_CHECK( element.father().level() == 0);
                    BOOST_CHECK( element.getOrigin().level() == 0);
                    BOOST_CHECK_EQUAL( (std::find(markedCells.begin(), markedCells.end(), element.father().index()) == markedCells.end()), false);
                }
                BOOST_CHECK( element.father() == element.getOrigin());
                if(!hasBeenRefinedAtLeastOnce)
                {
                    BOOST_CHECK( element.father().isLeaf() == false);
                }
                BOOST_CHECK( (element.level() > 0) || (element.level() < grid.maxLevel() +1));
                BOOST_CHECK(element.index() == adaptMapper.index(element));
                BOOST_CHECK(element.index() == adaptMapper.index(element));
                /** Not ideal to define this for each element. Remove?*/
                const auto& preAdapt_view = grid.levelGridView(element.father().level());
                Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LevelGridView> preAdaptMapper(preAdapt_view, Dune::mcmgElementLayout());
                BOOST_CHECK(element.father().index() == preAdaptMapper.index(element.father()));
            }
            else{
                BOOST_CHECK_THROW(element.father(), std::logic_error);
                BOOST_CHECK_THROW(element.geometryInFather(), std::logic_error);
                if (hasBeenRefinedAtLeastOnce){
                    BOOST_CHECK( element.level() <= startingGridIdx);
                    BOOST_CHECK( element.getOrigin().level() <= startingGridIdx);
                }
                else  {
                    BOOST_CHECK( element.level() == 0);
                }
                BOOST_CHECK( element.getOrigin().level() == 0);
                BOOST_CHECK( element.isNew() == false);
            }
            BOOST_CHECK( element.mightVanish() == false); // marks get rewrtitten and set to 0 via postAdapt call
        } // end-element-for-loop



        if (startingGridIdx == 0) {
            const auto& preAdapt_view = grid.levelGridView(startingGridIdx);
            Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LevelGridView> preAdaptMapper(preAdapt_view, Dune::mcmgElementLayout());
            // Some checks on the preAdapt grid
            for(const auto& element: elements(preAdapt_view)) {

                BOOST_CHECK( element.hasFather() == false);
                BOOST_CHECK_THROW(element.father(), std::logic_error);
                BOOST_CHECK_THROW(element.geometryInFather(), std::logic_error);
                BOOST_CHECK( element.getOrigin() ==  element);
                BOOST_CHECK( element.getOrigin().level() == startingGridIdx);
                BOOST_CHECK( element.isNew() == false);
                auto it = element.hbegin(grid.maxLevel()); // With element.level(), fails
                auto endIt = element.hend(grid.maxLevel());
                if (std::find(markedCells.begin(), markedCells.end(), element.index()) == markedCells.end()){
                    BOOST_CHECK( element.isLeaf() == true);
                    // If it == endIt, then entity.isLeaf() true (when dristibuted_data_ is empty)
                    BOOST_CHECK( it == endIt);
                    BOOST_CHECK( element.mightVanish() == false);
                }
                else{
                    // If it != endIt, then entity.isLeaf() false (when dristibuted_data_ is empty)
                    BOOST_CHECK_EQUAL( it == endIt, false);
                    BOOST_CHECK( element.mightVanish() == true);
                    BOOST_CHECK( element.isNew() == false);
                    BOOST_CHECK_EQUAL( element.isLeaf(), false); // parent cells do not appear in the LeafView

                    /// Auxiliary int to check amount of children
                    double referenceElemOneParent_volume_it = 0.;
                    std::array<double,3> referenceElem_entity_center_it = {0.,0.,0.}; // Expected {.5,.5,.5}
                    for (; it != endIt; ++it)
                    {
                        BOOST_CHECK(it ->hasFather() == true);
                        referenceElemOneParent_volume_it += it-> geometryInFather().volume();
                        for (int c = 0; c < 3; ++c)
                        {
                            referenceElem_entity_center_it[c] += (it-> geometryInFather().center())[c];
                        }
                    }
                    for (int c = 0; c < 3; ++c) {
                        referenceElem_entity_center_it[c] /= cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2];
                    }
                    BOOST_CHECK_CLOSE(referenceElemOneParent_volume_it, 1, 1e-13);
                    BOOST_CHECK_CLOSE(referenceElem_entity_center_it[0], .5, 1e-13);
                    BOOST_CHECK_CLOSE(referenceElem_entity_center_it[1], .5, 1e-13);
                    BOOST_CHECK_CLOSE(referenceElem_entity_center_it[2], .5, 1e-13);
                }
                BOOST_CHECK( element.level() == 0);
            } // end-preAdaptElements-for-loop
        } // end-startingGridIdx==0

        std::set<int> allIds_set;
        std::vector<int> allIds_vec;
        allIds_vec.reserve(data.back()->size(0) + data.back()->size(3));
        for (const auto& element: elements(grid_view)){
            const auto& localId = data.back()->localIdSet().id(element);
            const auto& globalId = data.back()->globalIdSet().id(element);
            // In serial run, local and global id coincide:
            BOOST_CHECK_EQUAL(localId, globalId);
            allIds_set.insert(localId);
            allIds_vec.push_back(localId);
            // Check that the global_id_set_ptr_ has the correct id (id from the level where the entity was born).
            BOOST_CHECK_EQUAL( grid.globalIdSet().id(element), data[element.level()]->localIdSet().id(element.getEquivLevelElem()));
        }
        // Check injectivity of the map local_id_set_ (and, indirectly, global_id_set_) after adding cell ids.
        BOOST_CHECK( allIds_set.size() == allIds_vec.size());

        for (const auto& point: vertices(grid_view)){
            const auto& localId = data.back()->localIdSet().id(point);
            const auto& globalId = data.back()->globalIdSet().id(point);
            BOOST_CHECK_EQUAL(localId, globalId);
            allIds_set.insert(localId);
            allIds_vec.push_back(localId);
        }
        // Check injectivity of the map local_id_set_ (and, indirectly, global_id_set_) after adding point ids.
        BOOST_CHECK( allIds_set.size() == allIds_vec.size());
        // CpGrid supports only elements (cells) and vertices (corners). Total amount of ids for the leaf grid view should coincide
        // with the total amount of cells and corners on the leaf grid view.
        BOOST_CHECK( static_cast<int>(allIds_set.size()) == (data.back()->size(0) + data.back()->size(3)));


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
                // The following check is commented even though all the test cases pass it. However, runnning this file
                // with it (uncommented) takes ~2.5 minutes.
                // Search in the leaf grid view elements for the element with the same id, if it exists.
                /*if (auto itIsLeaf = std::find_if( elements(coarse_grid.leafGridView()).begin(),
                  elements(coarse_grid.leafGridView()).end(),
                  [localId, data](const Dune::cpgrid::Entity<0>& leafElem)
                  { return (localId == data.back()->localIdSet().id(leafElem)); });
                  itIsLeaf != elements(coarse_grid.leafGridView()).end()) {
                  BOOST_CHECK( itIsLeaf->getEquivLevelElem() == element);
                  }*/
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
                // The following check is commented even though all the test cases pass it. However, runnning this file
                // with it (uncommented) takes ~2.5 minutes.
                /* // Search in the leaf grid view elements for the element with the same id, if it exists.
                   if (auto itIsLeaf = std::find_if( vertices(coarse_grid.leafGridView()).begin(),
                   vertices(coarse_grid.leafGridView()).end(),
                   [localId, data](const Dune::cpgrid::Entity<3>& leafPoint)
                   { return (localId == data.back()->localIdSet().id(leafPoint)); });
                   itIsLeaf != vertices(coarse_grid.leafGridView()).end()) {
                   BOOST_CHECK( (*itIsLeaf).geometry().center() == point.geometry().center() );
                   }*/
            }
            // Check injectivity of the map local_id_set_ (and, indirectly, global_id_set_)
            BOOST_CHECK( levelIds_set.size() == levelIds_vec.size());
            // CpGrid supports only elements (cells) and vertices (corners). Total amount of ids for each level grid should coincide
            // with the total amount of cells and corners on that level grid.
            BOOST_CHECK( static_cast<int>(levelIds_set.size()) == (data[level]->size(0) + data[level]->size(3)));
        }
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
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, true, false, true);
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

    // We set isBlockShape as false, even though global-refinement implies refinement of a block of cells.
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, other_grid, false, false, true);
}


BOOST_AUTO_TEST_CASE(doNothing_calling_globalRefine)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {2,2,2};
    std::vector<int> markedCells;
    coarse_grid.createCartesian(grid_dim, cell_sizes);

    // Create another grid
    Dune::CpGrid other_grid;
    other_grid.createCartesian(grid_dim, cell_sizes);
    other_grid.globalRefine(0);

    // We set isBlockShape as false, even though global-refinement implies refinement of a block of cells.
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, other_grid, true, false, true);
}

BOOST_AUTO_TEST_CASE(globalRefinement_calling_globalRefine)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);

    std::vector<int> markedCells(36);
    std::iota(markedCells.begin(), markedCells.end(), 0);

    // Create other grid for comparison
    Dune::CpGrid other_grid;
    other_grid.createCartesian(grid_dim, cell_sizes);
    other_grid.globalRefine(1);

    const std::array<int, 3> cells_per_dim = {2,2,2};

    // We set isBlockShape as false, even though global-refinement implies refinement of a block of cells.
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, other_grid, false, false, true);
}

BOOST_AUTO_TEST_CASE(calling_globalRefine_with_2)
{
    // Create a grid
    Dune::CpGrid equiv_fine_grid;
    const std::array<double, 3> cell_sizes = {0.5, 0.5, 0.5};
    const std::array<int, 3> grid_dim = {8,6,6};
    equiv_fine_grid.createCartesian(grid_dim, cell_sizes);

    std::vector<int> markedCells(288);
    std::iota(markedCells.begin(), markedCells.end(), 0);

    // Create other grid for comparison
    Dune::CpGrid other_grid;
    const std::array<double, 3> other_cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> other_grid_dim = {4,3,3};
    other_grid.createCartesian(other_grid_dim, other_cell_sizes);
    other_grid.globalRefine(2);

    const std::array<int, 3> cells_per_dim = {2,2,2};

    // We set isBlockShape as false, even though global-refinement implies refinement of a block of cells.
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    markAndAdapt_check(equiv_fine_grid, cells_per_dim, markedCells, other_grid, false, false, true);
}

BOOST_AUTO_TEST_CASE(calling_globalRefine_with_3)
{
    // Create a grid
    Dune::CpGrid equiv_fine_grid;
    const std::array<double, 3> cell_sizes = {0.25, 0.25, 0.25};
    const std::array<int, 3> grid_dim = {16,12,12};
    equiv_fine_grid.createCartesian(grid_dim, cell_sizes);

    std::vector<int> markedCells(2304);
    std::iota(markedCells.begin(), markedCells.end(), 0);

    // Create other grid for comparison
    Dune::CpGrid other_grid;
    const std::array<double, 3> other_cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> other_grid_dim = {4,3,3};
    other_grid.createCartesian(other_grid_dim, other_cell_sizes);
    other_grid.globalRefine(3);

    const std::array<int, 3> cells_per_dim = {2,2,2};

    // We set isBlockShape as false, even though global-refinement implies refinement of a block of cells.
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    markAndAdapt_check(equiv_fine_grid, cells_per_dim, markedCells, other_grid, false, false, true);
}

BOOST_AUTO_TEST_CASE(throw_globalRefine_with_negative_int)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {0.25, 0.25, 0.25};
    const std::array<int, 3> grid_dim = {16,12,12};
    grid.createCartesian(grid_dim, cell_sizes);

    BOOST_CHECK_THROW(grid.globalRefine(-5), std::logic_error);
}

BOOST_AUTO_TEST_CASE(throw_globalRefine_of_mixed_grid)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);

    const std::array<int, 3> cells_per_dim = {2,2,2};
    const std::array<int, 3> startIJK = {2,0,0};
    const std::array<int, 3> endIJK = {4,1,1};  // -> marked elements 2 and 3
    const std::string lgr_name = {"LGR1"};
    grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

    BOOST_CHECK_THROW(grid.globalRefine(1), std::logic_error);
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
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
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
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
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
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
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
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, false, false, false);
}

BOOST_AUTO_TEST_CASE(markNonBlockCells_compareAdapt)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {3,3,3};
    std::vector<int> markedCells = {1,4,6,9,17,22,28,32,33};
    coarse_grid.createCartesian(grid_dim, cell_sizes);

    // Create a grid
    Dune::CpGrid other_grid;
    other_grid.createCartesian(grid_dim, cell_sizes);
    for (const auto& elemIdx : markedCells)
    {
        const auto& elem =  Dune::cpgrid::Entity<0>(*(other_grid.currentData()[0]), elemIdx, true);
        other_grid.mark(1, elem);
    }
    other_grid.preAdapt();
    other_grid.adapt();
    other_grid.postAdapt();

    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, other_grid, false, false, false);
}

BOOST_AUTO_TEST_CASE(callAdaptMultipleTimes)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);

    const std::array<int, 3> cells_per_dim = {2,2,2};
    std::vector<int> markedCells1 = {1,4,6};
    std::vector<int> markedCells2 = {38, 43}; // Equivalent cells to level 0 cells with indices {17,22};
    std::vector<int> markedCells3 = {63, 67}; // Equivalent cells to level 0 cells with indices {28,32};

    std::vector<int> markedCells = {1,4,6,17,22,28,32};

    // Create a grid
    Dune::CpGrid other_grid;
    other_grid.createCartesian(grid_dim, cell_sizes);
    for (const auto& elemIdx : markedCells1)
    {
        const auto& elem =  Dune::cpgrid::Entity<0>(*(other_grid.currentData()[0]), elemIdx, true);
        other_grid.mark(1, elem);
    }
    other_grid.preAdapt();
    other_grid.adapt();
    other_grid.postAdapt();

    for (const auto& elemIdx : markedCells2)
    {
        const auto& elem =  Dune::cpgrid::Entity<0>(*(other_grid.currentData().back()), elemIdx, true);
        other_grid.mark(1, elem);
    }
    other_grid.preAdapt();
    other_grid.adapt();
    other_grid.postAdapt();

    for (const auto& elemIdx : markedCells3)
    {
        const auto& elem =  Dune::cpgrid::Entity<0>(*(other_grid.currentData().back()), elemIdx, true);
        other_grid.mark(1, elem);
    }
    other_grid.preAdapt();
    other_grid.adapt();
    other_grid.postAdapt();

    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, other_grid, false, false, false);
}

BOOST_AUTO_TEST_CASE(refineCoarseCells_in_mixedGrid) {
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {2,2,2};
    coarse_grid.createCartesian(grid_dim, cell_sizes);

    // LGR1 marked element with elemIdx = 3, refined into 8 children cells with leaf indices 3,...,10.
    const std::array<int, 3> startIJK = {3,0,0};
    const std::array<int, 3> endIJK = {4,1,1};
    const std::string lgr_name = {"LGR1"};
    coarse_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

    std::vector<int> markedCells = {0,1,11,15}; // coarse cells (in level 0 grid, this cell has index 8)
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, false, true, false);
}

BOOST_AUTO_TEST_CASE(refineInteriorRefinedCells_in_mixedGrid) {
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {3,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);

    // LGR1 marked elements with elemIdx = 17 and 18, refined into 27 children cells each,
    // with leaf indices 17+0,...,17+26 = 43 (children {level 0, cell index 17}),44,...,70 (children {level 0, cell index 18}).
    const std::array<int, 3> startIJK = {1,1,1};
    const std::array<int, 3> endIJK = {3,2,2};
    const std::string lgr_name = {"LGR1"};
    coarse_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

    // Cells 30, 31, 56, and 57 are refined cells, located in the interior of the refined-level-grid-1 (lgr 1 / level 1).
    // Therefore, cell_to_face_ for all of them has size 6. (Their faces have all 2 refined neigboring cells - (not one coarse cell, and one refined)).
    std::vector<int> markedCells = {30,31, 56,57};
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, false, true, false);
}


BOOST_AUTO_TEST_CASE(refineMixedCells_in_mixedGrid) {
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {3,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);

    // LGR1 marked elements with elemIdx = 17 and 18, refined into 27 children cells each,
    // with leaf indices 17+0,...,17+26 = 43 (children {level 0, cell index 17}),44,...,70 (children {level 0, cell index 18}).
    const std::array<int, 3> startIJK = {1,1,1};
    const std::array<int, 3> endIJK = {3,2,2};
    const std::string lgr_name = {"LGR1"};
    coarse_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

    // - Cells 30, 31, 56, and 57 are refined cells, located in the interior of the refined-level-grid-1 (lgr 1 / level 1).
    // Therefore, cell_to_face_ for all of them has size 6. (Their faces have all 2 refined neigboring cells - (not one coarse cell, and one refined)).
    // - Cells 0,1,2,12, and 15 are coarse cells, not touching the boundary of the LGR1 (cells 12 and 15 do share corners with LGR1 but do not share
    // any face. Therefore, the faces of cells 0,1,2,12,and 15 have all 1 or 2 neighboring coarse cells).
    std::vector<int> markedCells = {0,1,2,12,15,30,31,56,57};
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, false, true, false);
}


BOOST_AUTO_TEST_CASE(refineMixedCells_in_multiLevelGrid) {
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {3,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);

    // LGR1: element with elemIdx = 17, refined into 27 children cells with leaf indices 17+0,...,17+26 = 43 (children {level 0, cell index 17}).
    // LGR2: element with elemIdx = 18, refined into 27 children cells with leaf indices 44,...,70 (children {level 0, cell index 18}).
    const std::vector<std::array<int, 3>> startIJK_vec = {{1,1,1}, {2,1,1}};
    const std::vector<std::array<int, 3>> endIJK_vec = {{2,2,2}, {3,2,2}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    coarse_grid.addLgrsUpdateLeafView({cells_per_dim, cells_per_dim}, startIJK_vec, endIJK_vec, lgr_name_vec);

    // - Cells 30, 31, 56, and 57 are refined cells, located in the interior of the refined-level-grid-1 (lgr 1 / level 1).
    // Therefore, cell_to_face_ for all of them has size 6. (Their faces have all 2 refined neigboring cells - (not one coarse cell, and one refined)).
    // - Cells 0,1,2,12, and 15 are coarse cells, not touching the boundary of the LGR1 (cells 12 and 15 do share corners with LGR1 but do not share
    // any face. Therefore, the faces of cells 0,1,2,12,and 15 have all 1 or 2 neighboring coarse cells).
    std::vector<int> markedCells = {0,1,2,12,15,30,31,56,57};
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, false, true, false);
}


BOOST_AUTO_TEST_CASE(refineMixedCells_in_mixedGrid_II)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {3,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);

    // LGR1 marked elements with elemIdx = 17 and 18, refined into 27 children cells each,
    // with leaf indices 17+0,...,17+26 = 43 (children {level 0, cell index 17}),44,...,70 (children {level 0, cell index 18}).
    const std::array<int, 3> startIJK = {1,1,1};
    const std::array<int, 3> endIJK = {3,2,2};
    const std::string lgr_name = {"LGR1"};
    coarse_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

    // - Cells 72 (in level 0 with cell index 20) and 84 (in level 0 with cell index 32) are coarse cells,
    // sharing one K_FACE, that do not share faces with LGR1 (they do share corners).
    // - Cells 25,34,43 are refined cells, children of {level 0, cell index 17}, forming a collum.
    // - Cells 50,59,68 are refined cells, children of {level 0, cell index 18}, forming a collum.
    // Cells 25 and 50, 34 and 59, 43 and 68, share a face (the collums are next to each other).
    std::vector<int> markedCells = {25,34,43,50,59,68, 72, 84};
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, true, true, false);
}

BOOST_AUTO_TEST_CASE(cellTouchesLgrBoundary_throw)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim = {3,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);

    // LGR1 marked elements with elemIdx = 17 and 18, refined into 27 children cells each,
    // with leaf indices 17+0,...,17+26 = 43 (children {level 0, cell index 17}),44,...,70 (children {level 0, cell index 18}).
    const std::array<int, 3> startIJK = {1,1,1};
    const std::array<int, 3> endIJK = {3,2,2};
    const std::string lgr_name = {"LGR1"};
    coarse_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

    // - Coarse cells touching the LGR1 on its boundary.
    // Cell 5, 14, 16, 71, 73, and 82, touching the bottom, front, left, right, back, and the top of LGR1, respectively.
    std::vector<int> markedCells = {5,14,16,71,73,82};
    // The last three bool arguments represent: isBlockShape, hasBeenRefinedAtLeastOnce, isGlobalRefinement.
    BOOST_CHECK_THROW(markAndAdapt_check(coarse_grid, cells_per_dim, markedCells, coarse_grid, true, true, false), std::logic_error);
}
