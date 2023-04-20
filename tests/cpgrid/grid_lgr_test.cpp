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
    int count_faces = (cells_per_dim[0]*patch_dim[0]*cells_per_dim[1]*patch_dim[1]*((cells_per_dim[2]*patch_dim[2])+1)) // 'bottom/top faces'
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
                           const std::array<int, 3>& cells_per_dim,
                           const std::array<int,3>& start_ijk,
                           const std::array<int,3>& end_ijk)
{
    
    // Call createGridWithLgr()
    auto& data = coarse_grid.data_;
    coarse_grid.addLgrUpdateLeafView(cells_per_dim, start_ijk, end_ijk);
    BOOST_CHECK(data.size()==3);
    check_refinedPatch_grid(cells_per_dim, start_ijk, end_ijk,
                            (*coarse_grid.data_[1]).geometry_.template geomVector<0>(),
                            (*coarse_grid.data_[1]).geometry_.template geomVector<1>(),
                            (*coarse_grid.data_[1]).geometry_.template geomVector<3>());

    BOOST_CHECK( (*coarse_grid.data_[0]).child_to_parent_cells_.empty());
    BOOST_CHECK( (*coarse_grid.data_[1]).parent_to_children_cells_.empty());
    BOOST_CHECK( (*coarse_grid.data_[2]).parent_to_children_cells_.empty());
    const auto& [patch_corners, patch_faces, patch_cells] = (*coarse_grid.data_[0]).getPatchGeomIndices(start_ijk, end_ijk);
    double referenceElem_cellifiedPatch_volume = 0.; // Expected volume = patch_dim[0]*...[1]*...[2] (= patch_cells.size())
    for (int cell = 0; cell<  data[0]-> size(0); ++cell)
    {
        Dune::cpgrid::Entity<0> entity = Dune::cpgrid::Entity<0>((*coarse_grid.data_[0]), cell, true);
        BOOST_CHECK( entity.hasFather() == false);
        BOOST_CHECK_THROW(entity.father(), std::logic_error);
        BOOST_CHECK_THROW(entity.geometryInFather(), std::logic_error);
        auto it = entity.hbegin(coarse_grid.maxLevel());
        auto endIt = entity.hend(coarse_grid.maxLevel());
        const auto& parent_to_children = (*coarse_grid.data_[0]).parent_to_children_cells_[entity.index()];
        if (std::find(patch_cells.begin(), patch_cells.end(), cell) == patch_cells.end()){
            BOOST_CHECK_EQUAL(std::get<0>(parent_to_children), -1);
            BOOST_CHECK_EQUAL(std::get<1>(parent_to_children)[0], -1);
            BOOST_CHECK( entity.isLeaf() == true);   
            // If it == endIt, then entity.isLeaf() true (when dristibuted_data_ is empty)
            BOOST_CHECK( it == endIt);
        }
        else{
            BOOST_CHECK(std::get<1>(parent_to_children).size() > 0);
            BOOST_CHECK(std::get<0>(parent_to_children) != -1);
            BOOST_CHECK( entity.isLeaf() == false);
            // If it != endIt, then entity.isLeaf() false (when dristibuted_data_ is empty)
            BOOST_CHECK_EQUAL( it == endIt, false);
            // Auxiliary int to check amount of children
            int total_children = 0;
            for (; it != endIt; ++it)
            {
                // Do something with the son available through it->
                BOOST_CHECK(it ->hasFather() == true);
                BOOST_CHECK(it ->level() == 1);
                std::cout << it->index() << '\n';
                total_children += 1;
            }
            std::cout << "Entity Index: " << entity.index() << " has " << total_children << " children" << '\n';
            double referenceElem_entity_volume = 0.;
            std::array<double,3> referenceElem_entity_center = {0.,0.,0.}; // Expected {.5,.5,.5}
            for (const auto& child : std::get<1>(parent_to_children))
            {
                Dune::cpgrid::Entity<0> child_entity = Dune::cpgrid::Entity<0>(*data[1], child, true);
                BOOST_CHECK( child_entity.hasFather() == true);
                const auto& child_entity_geomInFather = child_entity.geometryInFather();
                referenceElem_entity_volume += child_entity_geomInFather.volume();
                referenceElem_cellifiedPatch_volume += child_entity_geomInFather.volume();
                for (int c = 0; c < 3; ++c)
                {
                    referenceElem_entity_center[c] += child_entity_geomInFather.center()[c];
                }
            }
            for (int c = 0; c < 3; ++c)
            {
                referenceElem_entity_center[c] /= cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2];
            }
            BOOST_CHECK( referenceElem_entity_volume == 1);
            BOOST_CHECK( referenceElem_entity_center[0] == 0.5);
            BOOST_CHECK( referenceElem_entity_center[1] == 0.5);
            BOOST_CHECK( referenceElem_entity_center[2] == 0.5);
        }
        BOOST_CHECK( entity.level() == 0);
    }
    BOOST_CHECK( referenceElem_cellifiedPatch_volume == patch_cells.size());
   



    
    for (int cell = 0; cell<  data[1]-> size(0); ++cell)
    {
        Dune::cpgrid::Entity<0> entity = Dune::cpgrid::Entity<0>((*coarse_grid.data_[1]), cell, true);
        BOOST_CHECK( entity.hasFather() == true);
        const auto& entity_geomInFather = entity.geometryInFather();
        BOOST_CHECK(entity_geomInFather.volume() == 1./((*data[1]).parent_to_children_cells_dim_[0]*
                                                        (*data[1]).parent_to_children_cells_dim_[1]*
                                                        (*data[1]).parent_to_children_cells_dim_[2]));
        BOOST_CHECK(entity.father().level() == 0);
        BOOST_CHECK_EQUAL( (std::find(patch_cells.begin(), patch_cells.end(), entity.father().index()) == patch_cells.end()), false);
        const auto& child_to_parent = (*coarse_grid.data_[1]).child_to_parent_cells_[cell];
        BOOST_CHECK_EQUAL( child_to_parent[0] == -1, false);
        BOOST_CHECK( entity.level() == 1);
        BOOST_CHECK( entity.level() == coarse_grid.maxLevel());
        BOOST_CHECK( entity.isLeaf() == true);
    }

    for (int cell = 0; cell<  data[2]-> size(0); ++cell)
    {
        Dune::cpgrid::Entity<0> entity = Dune::cpgrid::Entity<0>((*coarse_grid.data_[2]), cell, true);
        const auto& child_to_parent = (*coarse_grid.data_[2]).child_to_parent_cells_[cell];
        if (entity.hasFather()){
            const auto& entity_geomInFather = entity.geometryInFather();
            BOOST_CHECK(entity_geomInFather.volume() == 1./((*data[1]).parent_to_children_cells_dim_[0]*
                               (*data[1]).parent_to_children_cells_dim_[1]*(*data[1]).parent_to_children_cells_dim_[2]));
             BOOST_CHECK(entity.father().level() == 0);
             BOOST_CHECK_EQUAL( (std::find(patch_cells.begin(), patch_cells.end(), entity.father().index()) == patch_cells.end()), false);
             BOOST_CHECK(!(child_to_parent[0] == -1));
             BOOST_CHECK_EQUAL( child_to_parent[1], entity.father().index());
             BOOST_CHECK( entity.father().isLeaf() == false);
             BOOST_CHECK( entity.level() == 1);
             const auto& parent_to_children = (*coarse_grid.data_[0]).parent_to_children_cells_[entity.father().index()];
             double referenceElem_entity_volume = 0.;
             std::array<double,3> referenceElem_entity_center = {0.,0.,0.}; // Expected {.5,.5,.5}
             for (const auto& child : std::get<1>(parent_to_children))
             {
                 Dune::cpgrid::Entity<0> child_entity = Dune::cpgrid::Entity<0>(*data[1], child, true);
                 BOOST_CHECK( child_entity.hasFather() == true);
                 const auto& child_entity_geomInFather = child_entity.geometryInFather();
                 referenceElem_entity_volume += child_entity_geomInFather.volume();
                 referenceElem_cellifiedPatch_volume += child_entity_geomInFather.volume();
                 for (int c = 0; c < 3; ++c)
                 {
                     referenceElem_entity_center[c] += child_entity_geomInFather.center()[c];
                 }
             }
             for (int c = 0; c < 3; ++c)
             {
                 referenceElem_entity_center[c] /= cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2];
             }
             BOOST_CHECK( referenceElem_entity_volume == 1);
             BOOST_CHECK( referenceElem_entity_center[0] == 0.5);
             BOOST_CHECK( referenceElem_entity_center[1] == 0.5);
             BOOST_CHECK( referenceElem_entity_center[2] == 0.5);
        }
        else{
            BOOST_CHECK_THROW(entity.father(), std::logic_error);
            BOOST_CHECK_THROW(entity.geometryInFather(), std::logic_error);
            BOOST_CHECK_EQUAL( child_to_parent[0], -1);
            BOOST_CHECK( entity.level() == 0);
        }
        BOOST_CHECK( entity.isLeaf() == true);
    }

    for (int l = 0; l < 2; ++l)
    {
        const auto& view = coarse_grid.levelGridView(l);
        for (const auto& element: elements(view)){
            BOOST_CHECK_EQUAL(element.level(), l);
        }
    }

    const auto& leaf_view = coarse_grid.leafGridView();
    for (const auto& element: elements(leaf_view)){
         BOOST_CHECK( ((element.level() == 0) || (element.level() == 1)));
    }

} 

BOOST_AUTO_TEST_CASE(refine_patch)
{
     // Create a grid
    Dune::CpGrid coarse_grid;
    std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim_patch = {2,2,2};   
    std::array<int, 3> start_ijk = {1,0,1};
    std::array<int, 3> end_ijk = {3,2,3};  // then patch_dim = {3-1, 2-0, 3-1} ={2,2,2}
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    refinePatch_and_check(coarse_grid,
                          cells_per_dim_patch,
                          start_ijk, end_ijk);  
}

#define CHECK_COORDINATES(c1, c2)                                       \
    for (int c = 0; c < 3; c++) {                                        \
        BOOST_TEST(c1[c] == c2[c], boost::test_tools::tolerance(1e-12)); \
    }


void check_global_refine(const Dune::CpGrid& refined_grid, const Dune::CpGrid& equiv_fine_grid)
{

    const auto& refined_data = refined_grid.data_;
    const auto& equiv_data = equiv_fine_grid.data_;

    BOOST_CHECK(refined_data.size()==3);

    const auto& refined_leaf = *refined_data[2];
    const auto& equiv_leaf = *equiv_data[0];

    // Check the container sizes
    BOOST_CHECK_EQUAL(refined_leaf.face_to_cell_.size(), equiv_leaf.face_to_cell_.size());
    BOOST_CHECK_EQUAL(refined_leaf.face_to_point_.size(), equiv_leaf.face_to_point_.size());
    BOOST_CHECK_EQUAL(refined_leaf.cell_to_point_.size(), equiv_leaf.cell_to_point_.size());
    BOOST_CHECK_EQUAL(refined_leaf.face_normals_.size(), equiv_leaf.face_normals_.size());

    // Check that the points (ordering/coordinates) matches
    auto i = 0u;
    auto equiv_point_iter = equiv_leaf.geomVector<3>().begin();
    for(const auto& point: refined_leaf.geomVector<3>())
    {
        CHECK_COORDINATES(point.center(), equiv_point_iter->center());
        std::cout <<"point "<<i++<<": ";
        for(const auto& coord: point.center())
            std::cout << coord<< " ";
            std::cout <<std::endl;
        for(const auto& coord: point.center())
        BOOST_TEST(std::isfinite(coord)); 
        ++equiv_point_iter;
    }
    auto j = 0u;
    auto equiv_cell_iter = equiv_leaf.geomVector<3>().begin();
    for(const auto& cell: refined_leaf.geomVector<3>())
    {
        CHECK_COORDINATES(cell.center(), equiv_cell_iter->center());
        std::cout <<"cell "<<j++<<": ";
        for(const auto& coord: cell.center())
            std::cout << coord<< " ";
            std::cout <<std::endl;
        for(const auto& coord: cell.center())
        BOOST_TEST(std::isfinite(coord));
        BOOST_CHECK_CLOSE(cell.volume(), equiv_cell_iter->volume(), 1e-6);
        ++equiv_cell_iter;
    }

    /////
    const auto& grid_view = refined_grid.leafGridView();
    const auto& equiv_grid_view = equiv_fine_grid.leafGridView();

    auto equiv_element_iter = equiv_grid_view.begin<0>();
    for(const auto& element: elements(grid_view))
    {
        for(const auto& intersection: intersections(grid_view, element))
        {
            // find matching intersection (needed as ordering is allowed to be different
            bool matching_intersection_found = false;
            for(auto& intersection_match: intersections(equiv_grid_view, *equiv_element_iter))
            {
                if(intersection_match.indexInInside() == intersection.indexInInside())
                {
                    BOOST_CHECK(intersection_match.neighbor() == intersection.neighbor());

                    if(intersection.neighbor())
                    {
                        BOOST_CHECK(intersection_match.indexInOutside() == intersection.indexInOutside());
                    }

                    CHECK_COORDINATES(intersection_match.centerUnitOuterNormal(), intersection.centerUnitOuterNormal());
                    const auto& geom_match = intersection_match.geometry();
                    BOOST_TEST(0 == 1e-11, boost::test_tools::tolerance(1e-8));
                    const auto& geom =  intersection.geometry();
                    BOOST_CHECK_CLOSE(geom_match.volume(), geom.volume(), 1e-6);
                    CHECK_COORDINATES(geom_match.center(), geom.center());
                    BOOST_CHECK(geom_match.corners() == geom.corners());

                    decltype(geom.corner(0)) sum_match{}, sum{};

                    for(int cor = 0; cor < geom.corners(); ++cor)
                    {
                        sum += geom.corner(cor);
                        sum_match += geom_match.corner(1);
                    }
                    CHECK_COORDINATES(sum, sum_match);
                    matching_intersection_found = true;
                    break;
                }
            }
            BOOST_CHECK(matching_intersection_found);
        }
        ++equiv_element_iter;
    }
    /////
}


BOOST_AUTO_TEST_CASE(global_refine)
{
    // Create a 4x3x3 grid with length 4x3x3
    // and refine each cells into 4 children cells
    Dune::CpGrid coarse_grid;
    std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    std::array<int, 3> grid_dim = {4,3,3};
    std::array<int, 3> cells_per_dim_patch = {2,2,2};
    std::array<int, 3> start_ijk = {0,0,0};
    std::array<int, 3> end_ijk = {4,3,3};  
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    coarse_grid.addLgrUpdateLeafView(cells_per_dim_patch, start_ijk, end_ijk);

    // Create a 8x6x6 grid with length 4x3x3
    Dune::CpGrid fine_grid;
    std::array<double, 3> fine_cell_sizes = {0.5, 0.5, 0.5};
    std::array<int, 3> fine_grid_dim = {8,6,6};
    fine_grid.createCartesian(fine_grid_dim, fine_cell_sizes);

    check_global_refine(coarse_grid, fine_grid);
}


BOOST_AUTO_TEST_CASE(global_norefine)
{
    // Create a 4x3x3 grid with length 4x3x3
    // and refine each cells into 4 children cells
    Dune::CpGrid coarse_grid;
    std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    std::array<int, 3> grid_dim = {4,3,3};
    std::array<int, 3> cells_per_dim_patch = {1,1,1};
    std::array<int, 3> start_ijk = {0,0,0};
    std::array<int, 3> end_ijk = {4,3,3};  
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    coarse_grid.addLgrUpdateLeafView(cells_per_dim_patch, start_ijk, end_ijk);

    // Create a 8x6x6 grid with length 4x3x3
    Dune::CpGrid fine_grid;
    std::array<double, 3> fine_cell_sizes = {1.0, 1.0, 1.0};
    std::array<int, 3> fine_grid_dim = {4,3,3};
    fine_grid.createCartesian(fine_grid_dim, fine_cell_sizes);

    check_global_refine(coarse_grid, fine_grid);
}
