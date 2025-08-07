//===========================================================================
//
// File: addLgrsOnDistributedGrid_test.cpp
//
// Created: Monday 5 August 2024 14:44
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

#define BOOST_TEST_MODULE LGRsOnDistributedGridTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/common/CommunicationUtils.hpp>
#include <opm/grid/CpGrid.hpp>
#include <tests/cpgrid/lgr/LgrChecks.hpp>

#include <array>
#include <string>
#include <utility>
#include <vector>


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

std::vector<int> createTestCartesianGridAndParts(Dune::CpGrid& grid)
{
    grid.createCartesian({4,3,3} /*grid_dim*/, {1.0, 1.0, 1.0} /*cell_sizes*/);

    std::vector<int> parts(36);
    std::vector<std::vector<int>> cells_per_rank = { {0,1,4,5,8,9,16,20,21},
                                                     {12,13,17,24,25,28,29,32,33},
                                                     {2,3,6,7,10,11,18,22,23},
                                                     {14,15,19,26,27,30,31,34,35} };
    for (int rank = 0; rank < 4; ++rank) {
        for (const auto& elemIdx : cells_per_rank[rank]) {
            parts[elemIdx] = rank;
        }
    }
    return parts;
}

BOOST_AUTO_TEST_CASE(fullyInteriorLgrsHaveUniqueVertexGlobalIds)
{
    Dune::CpGrid grid;
    auto parts  = createTestCartesianGridAndParts(grid);

    if(grid.comm().size()>1)
    {
        grid.loadBalance(parts); // ownerFirst = false, addCornerCells = false, overlapLayerSize =1
        // It's not necessary to change the default values since the LGRs are fully interior.
        // Fully interior LGRs: each LGR is surrounded by other interior cells.
        const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}, {4,4,4}, {2,2,2}};
        const std::vector<std::array<int,3>> startIJK_vec = {{0,1,0}, {0,0,2}, {3,2,0}, {3,0,2}};
        const std::vector<std::array<int,3>> endIJK_vec = {{1,3,1}, {1,1,3}, {4,3,1}, {4,2,3}};
        const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3", "LGR4"};
        grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

        Opm::checkGridWithLgrs(grid, cells_per_dim_vec, lgr_name_vec, /* isGlobalRefined = */ false);

        // LGR1 dim 2x4x2 -> 3x5x3 = 45 points
        // LGR2 dim 3x3x3 -> 4x4x4 = 64 points
        // LGR3 dim 4x4x4 -> 5x5x5 = 125 points
        // LGR4 dim 2x4x2 -> 3x5x3 = 45 points
        const std::vector<int>& expected_vertex_ids_per_lgr = { 45, 64, 125, 45};
        // Total global ids in leaf grid view for points: 80 + 33 + 56 + 117 + 33 = 319
        Opm::checkExpectedVertexGlobalIdsCount(grid, expected_vertex_ids_per_lgr, 319 /*leaf_expected_vertex_ids*/);
    }
}


BOOST_AUTO_TEST_CASE(interiorLgrWithOverlapNeighborHasUniqueVertexGlobalIdsIfAddCornerCellsIsTrue)
{
    Dune::CpGrid grid;
    auto parts  = createTestCartesianGridAndParts(grid);

    if(grid.comm().size()>1)
    {
        grid.loadBalance(parts, false /*ownerFirst*/, true /*addCornerCells*/); // overlapLayerSize = 1
        // Set addCornersCells to true to achieve unique vertices global ids for the leaf and refined level grids.

        // NOT fully interior LGRs.
        // In rank 2, interior cell '7' has an overlap neighboring cell '19' (cell '19' is owned by rank 3).
        const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}, {4,4,4}, {2,2,2}};
        const std::vector<std::array<int,3>> startIJK_vec = {{0,1,0}, {0,0,2}, {3,1,0}, {3,0,2}};
        const std::vector<std::array<int,3>> endIJK_vec = {{1,3,1}, {1,1,3}, {4,2,1}, {4,2,3}};
        const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3", "LGR4"};
        grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

        Opm::checkGridWithLgrs(grid, cells_per_dim_vec, lgr_name_vec, /* isGlobalRefined = */ false);

        // LGR1 dim 2x4x2 -> 3x5x3 = 45 points
        // LGR2 dim 3x3x3 -> 4x4x4 = 64 points
        // LGR3 dim 4x4x4 -> 5x5x5 = 125 points
        // LGR4 dim 2x4x2 -> 3x5x3 = 45 points
        const std::vector<int>& expected_vertex_ids_per_lgr = { 45, 64, 125, 45};
        // Total global ids in leaf grid view for points: 80 + 33 + 56 + 117 + 33 = 319
        Opm::checkExpectedVertexGlobalIdsCount(grid, expected_vertex_ids_per_lgr, 319 /*leaf_expected_vertex_ids*/);
    }
}

BOOST_AUTO_TEST_CASE(distributedLgrHasUniqueVertexGlobalIdsIfAddCornerCellsTrue)
{
    Dune::CpGrid grid;
    auto parts  = createTestCartesianGridAndParts(grid);

    if(grid.comm().size()>1)
    {
        grid.loadBalance(parts, false, true); // ownerFirst = false, addCornerCells = true, overlapLayerSize = 1
        // Set addCornersCells to true to achieve unique vertices global ids for the leaf and refined level grids.

        const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}};
        const std::vector<std::array<int,3>> startIJK_vec = {{0,2,0}};
        const std::vector<std::array<int,3>> endIJK_vec = {{3,3,1}};
        const std::vector<std::string> lgr_name_vec = {"LGR1"};
        // LGR1 element indices = 8,9 (rank 0), 10 (rank 2). Total 24 refined cells, 63 points (63-16 = 47 with new global id).
        grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

        Opm::checkGridWithLgrs(grid, cells_per_dim_vec, lgr_name_vec, /* isGlobalRefined = */ false);

        // Check global id is not duplicated for points for each LGR
        // LGR1 dim 6x2x2 -> 7x3x3 = 63 points
        // Total global ids in leaf grid view for points: 80 + (63 - 16) = 127
        Opm::checkExpectedVertexGlobalIdsCount(grid, {63} /*expected_vertex_ids_per_lgr*/, 127 /*leaf_expected_vertex_ids*/);
    }
}

BOOST_AUTO_TEST_CASE(distributedLgrFailsVertexGlobalIdsUniquenessWithAddCornerCellsTrue)
{
    Dune::CpGrid grid;
    auto parts   = createTestCartesianGridAndParts(grid);

    if(grid.comm().size()>1)
    {
        grid.loadBalance(parts, false /*ownerFirst*/, true /*addCornerCells*/); // overlapLayerSize = 1
        // Setting addCornersCells to true is not enough to get unique vertices global ids for the leaf and refined level grids.

        const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}};
        const std::vector<std::array<int,3>> startIJK_vec = {{1,0,0}};
        const std::vector<std::array<int,3>> endIJK_vec = {{3,2,2}};
        const std::vector<std::string> lgr_name_vec = {"LGR1"};
        // LGR1 element indices = {1,2,5,6,13,14,17,18} where
        // 1,5 are interior in rank 0, 13,17 in rank 1, 2,6,18 in rank 2, 14 in rank 3.
        grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
        Opm::checkGridWithLgrs(grid, cells_per_dim_vec, lgr_name_vec, /* isGlobalRefined = */ false);

        // Check global id is not duplicated for points for each LGR
        // LGR1 dim 4x4x4 -> 5x5x5 = 125 points
        const int expected_vertex_ids = 125;
        // Remark: the LGR is distributed in ALL processes BUT
        // - P0, P1, P2 see the entire LGR (due to the argument addCornerCells = true).
        // However, P3 does NOT see cell 5, but sees all the other parent cells {1,2,6,13,14,17,18},
        // In P3, {1,6,11,17,21} are (overlap) added corner cells, and {2,3,7,13, 18,22,,25,29,33} are overlap
        // cells coming from overlapLayerSize = 1.
        // P3 creates +7 duplicated ids on the unseen-in-P3 cell 5
        const int additional_unnecesary_ids = 7;
        // Total global ids in leaf grid view for points: 80 + (125 - 27) = 178
        Opm::checkExpectedVertexGlobalIdsCount(grid,
                                               {expected_vertex_ids + additional_unnecesary_ids},
                                               178 /*leaf_expected_ertex_ids*/ + additional_unnecesary_ids);
    }
}

BOOST_AUTO_TEST_CASE(callAdaptWithArgsIsEquivalentToCallAddLgrsUpdateLeafGridViewOnDistributedCoarseGrid)
{
    Dune::CpGrid grid;
    auto parts = createTestCartesianGridAndParts(grid);

    if(grid.comm().size()>1)
    {
        grid.loadBalance(parts, false /*ownerFirst*/, true /*addCornerCells*/); // overlapLayerSize = 1
        // Setting addCornersCells to true is not enough to get unique vertices global ids for the leaf and refined level grids.

        const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}};
        const std::vector<std::array<int,3>> startIJK_vec = {{1,0,0}};
        const std::vector<std::array<int,3>> endIJK_vec = {{3,2,2}};

        // 1,5 are interior in rank 0, 13,17 in rank 1, 2,6,18 in rank 2, 14 in rank 3.
        const std::vector<int>& marked_elemIdx = {1,2,5,6,13,14,17,18};
        std::vector<int> assignRefinedLevel(grid.currentData().front()->size(0));
        for (const auto& idx : marked_elemIdx)
            assignRefinedLevel[idx] = 1;

        grid.adapt(cells_per_dim_vec,
                   assignRefinedLevel,
                   {"LGR1"} /*lgr_name_vec*/,
                   startIJK_vec,
                   endIJK_vec);

        Opm::checkGridWithLgrs(grid, cells_per_dim_vec, {"LGR1"}, /* isGlobalRefined = */ false);
    }
}

BOOST_AUTO_TEST_CASE(callAdaptWithoutArgsOnDistributedCoarseGrid)
{
    Dune::CpGrid grid;
    auto parts = createTestCartesianGridAndParts(grid);

    if(grid.comm().size()>1)
    {
        grid.loadBalance(parts);

        const auto& leafGridView = grid.leafGridView();
        // Mark all elements -> 'indirect' global refinement
        for (const auto& element : elements(leafGridView)){
            grid.mark(1, element);
        }
        grid.preAdapt();
        grid.adapt();
        grid.postAdapt();

        Opm::checkGridWithLgrs(grid, {{2,2,2}} /*cells_per_dim_vec*/, {"GR1"} /*lgr_name_vec (GR: GLOBAL REFINEMENT)*/,  /* isGlobalRefined = */ true);
    }
}

BOOST_AUTO_TEST_CASE(callGlobalRefineOnceOnDistributedCoarseGrid)
{
    Dune::CpGrid grid;
    grid.createCartesian({4,3,3} /*grid_dim*/, {1.0, 1.0, 1.0} /*cell_sizes*/);

    if(grid.comm().size()>1)
    {
        grid.loadBalance();
        grid.globalRefine(1);

        Opm::checkGridWithLgrs(grid,  {{2,2,2}} /*cells_per_dim_vec*/, {"GR1"} /*lgr_name_vec (GR: GLOBAL REFINEMENT)*/,  /* isGlobalRefined = */ true);
    }
}
