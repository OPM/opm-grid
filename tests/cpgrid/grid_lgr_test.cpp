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
                           const std::vector<std::array<int,3>>& cells_per_dim_vec,
                           const std::vector<std::array<int,3>>& startIJK_vec,
                           const std::vector<std::array<int,3>>& endIJK_vec)
{
    auto& data = coarse_grid.data_;
    // Add LGRs and update grid.
    const bool are_disjoint = (*data[0]).disjointPatches(startIJK_vec, endIJK_vec);
    if (are_disjoint){
        coarse_grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec);

        BOOST_CHECK(data.size() == startIJK_vec.size() + 2);
        BOOST_CHECK( (*data[0]).child_to_parent_cells_.empty());

        const auto& all_parent_cell_indices = (*data[0]).getPatchesCells(startIJK_vec, endIJK_vec);

        for (long unsigned int level = 1; level < startIJK_vec.size() +1; ++level) // only 1 when there is only 1 patch
        {
            check_refinedPatch_grid(cells_per_dim_vec[level-1], startIJK_vec[level-1], endIJK_vec[level-1],
                                    (*data[level]).geometry_.template geomVector<0>(),
                                    (*data[level]).geometry_.template geomVector<1>(),
                                    (*data[level]).geometry_.template geomVector<3>());
            BOOST_CHECK( (*data[level]).parent_to_children_cells_.empty());

            const auto& patch_cells = (*data[0]).getPatchCells(startIJK_vec[level-1], endIJK_vec[level-1]);

            // GLOBAL grid
            for (int cell = 0; cell<  data[0]-> size(0); ++cell)
            {
                Dune::cpgrid::Entity<0> entity = Dune::cpgrid::Entity<0>((*coarse_grid.data_[0]), cell, true);
                BOOST_CHECK( entity.hasFather() == false);
                BOOST_CHECK_THROW(entity.father(), std::logic_error);
                const auto& parent_to_children = (*coarse_grid.data_[0]).parent_to_children_cells_[cell];
                if (std::find(all_parent_cell_indices.begin(), all_parent_cell_indices.end(), cell) == all_parent_cell_indices.end()){
                    BOOST_CHECK_EQUAL(std::get<0>(parent_to_children), -1);
                    BOOST_CHECK(std::get<1>(parent_to_children).empty());
                    BOOST_CHECK( entity.isLeaf() == true);
                }
                else{
                    BOOST_CHECK( std::get<0>(parent_to_children) != -1);
                    BOOST_CHECK( std::get<1>(parent_to_children).size() > 1);
                    for (const auto& child : std::get<1>(parent_to_children)){
                        BOOST_CHECK( child != -1);
                        BOOST_CHECK( data[std::get<0>(parent_to_children)] -> child_to_parent_cells_[child][0] == 0);
                        BOOST_CHECK( data[std::get<0>(parent_to_children)] -> child_to_parent_cells_[child][1] == cell);
                    }
                    BOOST_CHECK_EQUAL( entity.isLeaf(), false); // parent cells do not appear in the LeafView
                }
                BOOST_CHECK( entity.level() == 0);
            }
            // LGRs
            for (int cell = 0; cell<  data[level]-> size(0); ++cell)
            {
                Dune::cpgrid::Entity<0> entity = Dune::cpgrid::Entity<0>((*coarse_grid.data_[level]), cell, true);
                BOOST_CHECK( entity.hasFather() == true);
                BOOST_CHECK(entity.father().level() == 0);
                // Check entity.father().index() belongs to the patch_cells (corresponding LGR parents)
                BOOST_CHECK_EQUAL((std::find(patch_cells.begin(), patch_cells.end(),
                                             entity.father().index()) == patch_cells.end()) , false);
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
            }
            // LeafView
            for (int cell = 0; cell <  data[startIJK_vec.size()+1]-> size(0); ++cell)
            {
                Dune::cpgrid::Entity<0> entity = Dune::cpgrid::Entity<0>((*coarse_grid.data_[startIJK_vec.size()+1]), cell, true);
                const auto& child_to_parent = (*data[startIJK_vec.size()+1]).child_to_parent_cells_[cell];
                const auto& level_cellIdx = (*data[startIJK_vec.size()+1]).leaf_to_level_cells_[entity.index()];
                if (entity.hasFather()){
                    BOOST_CHECK(entity.father().level() == 0);
                    BOOST_CHECK_EQUAL( (std::find(all_parent_cell_indices.begin(), all_parent_cell_indices.end(),
                                                  entity.father().index()) == all_parent_cell_indices.end()), false);
                    BOOST_CHECK(child_to_parent[0] != -1);
                    BOOST_CHECK_EQUAL( child_to_parent[0] == 0, true);
                    BOOST_CHECK_EQUAL( child_to_parent[1], entity.father().index());
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
                    BOOST_CHECK_EQUAL(child_to_parent[0], -1);
                    BOOST_CHECK_EQUAL(child_to_parent[1], -1);
                    BOOST_CHECK( level_cellIdx[0] == 0);
                    BOOST_CHECK( std::get<0>((*data[0]).parent_to_children_cells_[level_cellIdx[1]]) == -1);
                    BOOST_CHECK( std::get<1>((*data[0]).parent_to_children_cells_[level_cellIdx[1]]).empty());
                    BOOST_CHECK( entity.level() == 0);
                    // Get index of the cell in level 0
                    const auto& entityOldIdx =  (*data[startIJK_vec.size()+1]).leaf_to_level_cells_[entity.index()][1];
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
                    for (long unsigned int patch = 0; patch < startIJK_vec.size(); ++patch)
                    {
                        // LEFT Check if i == endIJK_vec[patch][0] (and j, k in the range) for some patch.
                        // If true, increase amount of faces by cells_per_dim_vec[patch][1]*cells_per_dim_vec[patch][2]
                        // RIGHT Check if i+1 == startIJK_vec[patch][0] (and j, k in the range) for some patch.
                        // If true, increase amount of faces by cells_per_dim_vec[patch][1]*cells_per_dim_vec[patch][2]
                        // Auxiliary bools
                        touch_patch_onLeftFace = touch_patch_onLeftFace ||
                            ((entityOldIJK[0] == endIJK_vec[patch][0]) &&
                             (entityOldIJK[1] >= startIJK_vec[patch][1]) && (entityOldIJK[1] < endIJK_vec[patch][1]) &&
                             (entityOldIJK[2] >= startIJK_vec[patch][2]) && (entityOldIJK[2] < endIJK_vec[patch][2]));
                        touch_patch_onRightFace = touch_patch_onRightFace ||
                            ((entityOldIJK[0]+1 == startIJK_vec[patch][0]) &&
                             (entityOldIJK[1] >= startIJK_vec[patch][1]) && (entityOldIJK[1] < endIJK_vec[patch][1]) &&
                             (entityOldIJK[2] >= startIJK_vec[patch][2]) && (entityOldIJK[2] < endIJK_vec[patch][2]));
                        if (touch_patch_onLeftFace || touch_patch_onRightFace) // true when at lesat one is true
                        {
                            face_count += (touch_patch_onLeftFace + touch_patch_onRightFace)*
                                (cells_per_dim_vec[patch][1]*cells_per_dim_vec[patch][2]);
                            // when both are true, we summ twice, otherwise, we sum once.
                            break;
                        }
                    } // end patch-for-loop
                    for (long unsigned int patch = 0; patch < startIJK_vec.size(); ++patch)
                    {
                        // FRONT Check if j == endIJK_vec[patch][1] (and i, k in the range) for some patch.
                        // If true, increase amount of faces by cells_per_dim_vec[patch][0]*cells_per_dim_vec[patch][2]
                        // BACK Check if j+1 == statIJK_vec[patch][1] (and i, k in the range) for some patch.
                        // If true, increase amount of faces by cells_per_dim_vec[patch][0]*cells_per_dim_vec[patch][2]
                        // Auxiliary bools
                        touch_patch_onFrontFace = touch_patch_onFrontFace ||
                            ((entityOldIJK[1] == endIJK_vec[patch][1]) &&
                             (entityOldIJK[0] >= startIJK_vec[patch][0]) && (entityOldIJK[0] < endIJK_vec[patch][0]) &&
                             (entityOldIJK[2] >= startIJK_vec[patch][2]) && (entityOldIJK[2] < endIJK_vec[patch][2]));
                        touch_patch_onBackFace = touch_patch_onBackFace ||
                            ((entityOldIJK[1]+1 == startIJK_vec[patch][1]) &&
                             (entityOldIJK[0] >= startIJK_vec[patch][0]) && (entityOldIJK[0] < endIJK_vec[patch][0]) &&
                             (entityOldIJK[2] >= startIJK_vec[patch][2]) && (entityOldIJK[2] < endIJK_vec[patch][2]));
                        if (touch_patch_onFrontFace || touch_patch_onBackFace)
                        {
                            face_count +=  (touch_patch_onFrontFace + touch_patch_onBackFace)*
                                (cells_per_dim_vec[patch][0]*cells_per_dim_vec[patch][2]);
                            break;
                        }
                    } // end patch-for-loop
                    for (long unsigned int patch = 0; patch < startIJK_vec.size(); ++patch)
                    {
                        // BOTTOM Check if k == endIJK_vec[patch][2] (and i, j in the range) for some patch.
                        // If true, increase amount of faces by cells_per_dim_vec[patch][0]*cells_per_dim_vec[patch][1]
                        // TOP Check if k+1 == startIJK_vec[patch][2] (and i, j in the range) for some patch.
                        // If true, increase amount of faces by cells_per_dim_vec[patch][0]*cells_per_dim_vec[patch][2]
                        // Auxiliary bools
                        touch_patch_onBottomFace = touch_patch_onBottomFace ||
                            ((entityOldIJK[2] == endIJK_vec[patch][2]) &&
                             (entityOldIJK[0] >= startIJK_vec[patch][0]) && (entityOldIJK[0] < endIJK_vec[patch][0]) &&
                             (entityOldIJK[1] >= startIJK_vec[patch][1]) && (entityOldIJK[1] < endIJK_vec[patch][1]));
                        touch_patch_onTopFace = touch_patch_onTopFace ||
                            ((entityOldIJK[2]+1 == startIJK_vec[patch][2]) &&
                             (entityOldIJK[0] >= startIJK_vec[patch][0]) && (entityOldIJK[0] < endIJK_vec[patch][0]) &&
                             (entityOldIJK[1] >= startIJK_vec[patch][1]) && (entityOldIJK[1] < endIJK_vec[patch][1]));
                        if (touch_patch_onBottomFace || touch_patch_onTopFace)
                        {
                            face_count +=  (touch_patch_onBottomFace + touch_patch_onTopFace)*
                                (cells_per_dim_vec[patch][0]*cells_per_dim_vec[patch][1]);
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
                    std::cout << "Entity old index: " << entityOldIdx << "Entity leaf idx: " << entity.index() << '\n';
                    std::cout << "Leaf cell_to_face_.size(): " << leaf_cell_to_face.size() << '\n';
                    std::cout << "Face count: " << face_count << '\n';
                    std::cout << "touch_patch_onLeftFace: " << touch_patch_onLeftFace << '\n';
                    std::cout << "touch_patch_onRightFace: " << touch_patch_onRightFace << '\n';
                    std::cout << "touch_patch_onFrontFace: " << touch_patch_onFrontFace << '\n';
                    std::cout << "touch_patch_onBackFace: " << touch_patch_onBackFace << '\n';
                    std::cout << "touch_patch_onBottomFace: " << touch_patch_onBottomFace << '\n';
                    std::cout << "touch_patch_onTopFace: " << touch_patch_onTopFace << '\n';
                    BOOST_CHECK( leaf_cell_to_face.size() == face_count);
                } // end else
                BOOST_CHECK( entity.isLeaf() == true);
            }
        }

        BOOST_CHECK( static_cast<int>(startIJK_vec.size()) == coarse_grid.maxLevel());
        BOOST_CHECK( (*data[data.size()-1]).parent_to_children_cells_.empty());

        for (long unsigned int l = 0; l < startIJK_vec.size() +1; ++l) // level 0,1,2,... , last patch
        {
            const auto& view = coarse_grid.levelGridView(l);
            for (const auto& element: elements(view)){
                BOOST_CHECK_EQUAL(element.level(), l);
            }
        }

        const auto& leaf_view = coarse_grid.leafGridView();
        for (const auto& element: elements(leaf_view)){
            BOOST_CHECK( ((element.level() >= 0) || (element.level() < static_cast<int>(startIJK_vec.size()) +1)));
        }
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
    std::array<int, 3> end_ijk = {3,2,3};  // patch_dim = {3-1, 2-0, 3-1} ={2,2,2}
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    refinePatch_and_check(coarse_grid, {cells_per_dim_patch}, {start_ijk}, {end_ijk});
}

BOOST_AUTO_TEST_CASE(refine_patch_one_cell)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    std::array<int, 3> grid_dim = {4,3,3};
    const std::array<int, 3> cells_per_dim_patch = {2,2,2};
    std::array<int, 3> start_ijk = {1,0,1};
    std::array<int, 3> end_ijk = {2,1,2};  // patch_dim = {2-1, 1-0, 2-1} ={1,1,1} -> Single Cell!
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    refinePatch_and_check(coarse_grid, {cells_per_dim_patch}, {start_ijk}, {end_ijk});
    // Choose a cell touching the boundary of the patch and check amount of faces
}

BOOST_AUTO_TEST_CASE(lgrs_disjointPatches)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}, {4,4,4}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,2}, {3,2,2}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,1,1}, {1,1,3}, {4,3,3}};
    refinePatch_and_check(coarse_grid, cells_per_dim_vec, startIJK_vec, endIJK_vec);
}

BOOST_AUTO_TEST_CASE(lgrs_disjointPatchesB)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {3,2,0}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,2,1}, {4,3,3}};
    refinePatch_and_check(coarse_grid, cells_per_dim_vec, startIJK_vec, endIJK_vec);
}

BOOST_AUTO_TEST_CASE(patches_share_corner)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}, {4,4,4}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {1,1,1}, {3,2,2}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,1}, {2,2,2}, {4,3,3}};
    BOOST_CHECK_THROW(coarse_grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec), std::logic_error);
    std::cout << "Patches are NOT disjoint" << "\n";
}

BOOST_AUTO_TEST_CASE(patches_share_corners)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {2,2,2}, {2,2,2}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {2,0,1}, {3,2,2}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,1,1}, {3,1,1}, {4,3,3}};
    BOOST_CHECK_THROW(coarse_grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec), std::logic_error);
    std::cout << "Patches are NOT disjoint" << "\n";
}

BOOST_AUTO_TEST_CASE(pathces_share_face)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {2,2,2}, {2,2,2}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {2,0,0}, {3,2,2}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,1,1}, {3,1,1}, {4,3,3}};
    BOOST_CHECK_THROW(coarse_grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec), std::logic_error);
    std::cout << "Patches are NOT disjoint" << "\n";
}

BOOST_AUTO_TEST_CASE(pathces_share_faceB)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {2,2,2}, {2,2,2}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,1}, {1,1,2}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,2,1}, {3,2,1}, {4,3,3}};
    BOOST_CHECK_THROW(coarse_grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec), std::logic_error);
    std::cout << "Patches are NOT disjoint" << "\n";
}

#define CHECK_COORDINATES(c1, c2)                                       \
    for (int c = 0; c < 3; c++) {                                       \
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
    coarse_grid.addLgrsUpdateLeafView({cells_per_dim_patch}, {start_ijk}, {end_ijk});

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
    coarse_grid.addLgrsUpdateLeafView({cells_per_dim_patch}, {start_ijk}, {end_ijk});

    // Create a 8x6x6 grid with length 4x3x3
    Dune::CpGrid fine_grid;
    std::array<double, 3> fine_cell_sizes = {1.0, 1.0, 1.0};
    std::array<int, 3> fine_grid_dim = {4,3,3};
    fine_grid.createCartesian(fine_grid_dim, fine_cell_sizes);

    check_global_refine(coarse_grid, fine_grid);
}
