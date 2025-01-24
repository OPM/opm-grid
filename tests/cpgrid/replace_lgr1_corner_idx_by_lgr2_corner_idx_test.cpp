//===========================================================================
//
// File: replace_corner_idx_neigh_cells_test.cpp
//
// Created: Friday 24.01.2025 11:55:00
//
// Author(s): Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================
/*
  Copyright 2025 Equinor ASA.

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

#define BOOST_TEST_MODULE ReplaceCornerIdxNeighCellsTests
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif


#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>

#include <array>


struct Fixture
{
    Fixture()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(m_argc, m_argv);
        Opm::OpmLog::setupSimpleDefaultLogging();
    }
};

BOOST_GLOBAL_FIXTURE(Fixture);

/* To add LGRs in a CpGrid, each marked element for refinement gets refined into a
   single-cell-refinement. To store new born refined entities [points or faces]
   only once, we detect equivalent entities via their indices in each single-cell-refinement.

   To test replaceLgr1CornerIdxByLgr2CornerIdx(...), we create two grids, each of them
   containing only one cell. Add different LGRs in each grid, in different cases to
   create all possible escenarios (the two coarse cells sharing I_FACE, J_FACE or K_FACE).

   To avoid friend declarations, we 'copy' here replaceLgr1CornerIdxByLgr2CornerIdx
   (and getRefinedCornerIJK).

   Why single-cell-refinements instead of adding LGRs in one grid?:
   replaceLgr1CornerIdxByLgr2CornerIdx is a private method in CpGrid, meant to be
   used between independent/un-related single-cell-refinements (stored in CpGridData objects).
   Therefore, creating a grid, adding LGRs to it, and checking corner indices afterwards in
   refined cells on bouandary of the LGRs is not the escenario where this method should be tested.
*/

std::array<int,3> getRefinedCornerIJK(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr)
{
    // Order defined in Geometry::refine
    //  (j*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (i*(cells_per_dim[2]+1)) + k
    std::array<int,3> ijk;
    ijk[2] = cornerIdxInLgr % (cells_per_dim[2] +1);
    cornerIdxInLgr -= ijk[2];
    cornerIdxInLgr /= (cells_per_dim[2] +1);
    ijk[0] = cornerIdxInLgr % (cells_per_dim[0]+1);
    cornerIdxInLgr -=ijk[0];
    ijk[1] = cornerIdxInLgr / (cells_per_dim[0]+1);
    return ijk;
}

// Mimic CpGrid::replaceLgr1CornerIdxByLgr2CornerIdx
int replaceLgr1CornerIdxByLgr2CornerIdx(const std::array<int,3>& cells_per_dim_lgr1,
                                        int cornerIdxLgr1,
                                        const std::array<int,3>& cells_per_dim_lgr2)
{
    const auto& ijkLgr1 = getRefinedCornerIJK(cells_per_dim_lgr1, cornerIdxLgr1);
    // Order defined in Geometry::refine
    // (j*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (i*(cells_per_dim[2]+1)) + k

    // On a parallel run, no symmetry between neighboring elements should be assumed. Therefore, all the six cases
    // (i = 0, cells_per_dim[0], j = 0, cells_per_dim[1], and k = 0, cells_per_dim[2]) have to be taken into account.
    // On a serial run, it would be enough to consider i = cells_per_dim[0], j = cells_per_dim[1], and k = cells_per_dim[2].
    // To cover all possible scenarios, serial and parallel, we consider the six cases.

    if (ijkLgr1[0] == cells_per_dim_lgr1[0]) { // same j, k, but i = 0
        return   (ijkLgr1[1]*(cells_per_dim_lgr2[0]+1)*(cells_per_dim_lgr2[2]+1)) + ijkLgr1[2];
    }
    if (ijkLgr1[1] == cells_per_dim_lgr1[1]) { // same i,k, but j = 0
        return  (ijkLgr1[0]*(cells_per_dim_lgr2[2]+1)) + ijkLgr1[2];
    }
    if (ijkLgr1[2] == cells_per_dim_lgr1[2]) { // same i,j, but k = 0
        return  (ijkLgr1[1]*(cells_per_dim_lgr2[0]+1)*(cells_per_dim_lgr2[2]+1)) + (ijkLgr1[0]*(cells_per_dim_lgr2[2]+1));
    }
    if (ijkLgr1[0] == 0) { // same j,k, but i = cells_per_dim[0]
        return   (ijkLgr1[1]*(cells_per_dim_lgr2[0]+1)*(cells_per_dim_lgr2[2]+1)) + (cells_per_dim_lgr2[0]*(cells_per_dim_lgr2[2]+1))+ ijkLgr1[2];
    }
    if (ijkLgr1[1] == 0) { // same i,k, but j = cells_per_dim[1]
        return  (cells_per_dim_lgr2[1]*(cells_per_dim_lgr2[0]+1)*(cells_per_dim_lgr2[2]+1)) + (ijkLgr1[0]*(cells_per_dim_lgr2[2]+1)) + ijkLgr1[2];
    }
    if (ijkLgr1[2] == 0) { // same i,j, but k = cells_per_dim[2]
        return  (ijkLgr1[1]*(cells_per_dim_lgr2[0]+1)*(cells_per_dim_lgr2[2]+1)) + (ijkLgr1[0]*(cells_per_dim_lgr2[2]+1)) + cells_per_dim_lgr2[2];
    }
    else {
        const auto& message = "Cannot convert corner index from one LGR to its neighboring LGR.";
        OPM_THROW(std::logic_error, message);
    }
}

const std::shared_ptr<Dune::cpgrid::CpGridData> createSingleCellGridAndRefine(const std::array<int,3>& lgr_dim)
{
    // Create two grids, one single cell in each grid.
    Dune::CpGrid lgr;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int,3> coarse_grid_dim = {1,1,1};
    lgr.createCartesian(coarse_grid_dim, cell_sizes);

    // Single-cell-refinement for the only cell contained in lgr grid.
    const auto& [lgr_ptr,
                 lgr_parentCorners_to_equivalentRefinedCorners,
                 lgr_parentFace_to_itsRefinedFaces,
                 lgr_parentCell_to_itsRefinedCells,
                 lgr_refinedFace_to_itsParentFace,
                 lgr_refinedCell_to_itsParentCell]
        = lgr.currentData().back()->refineSingleCell(lgr_dim, 0);
    return lgr_ptr;
}

void check4CornersOfSharedRefinedFaceBetweenSingleCellRefinements(const Dune::cpgrid::Entity<0>& elemLgrA,
                                                                  const Dune::cpgrid::Entity<0>& elemLgrB,
                                                                  const std::unordered_map<int, int>& lgrA_to_lgrB,
                                                                  const std::array<int,3>& lgrA_dim,
                                                                  const std::array<int,3>& lgrB_dim)
{
    for (const auto& [idx_in_cell_lgrA, idx_in_cell_lgrB] : lgrA_to_lgrB)
    {
        const auto& corn_lgrA = elemLgrA.subEntity<3>( idx_in_cell_lgrA ).index();
        const auto& corn_lgrB = elemLgrB.subEntity<3>( idx_in_cell_lgrB ).index();

        BOOST_CHECK_EQUAL( replaceLgr1CornerIdxByLgr2CornerIdx(lgrA_dim, corn_lgrA, lgrB_dim), corn_lgrB);
        BOOST_CHECK_EQUAL( replaceLgr1CornerIdxByLgr2CornerIdx(lgrB_dim, corn_lgrB, lgrA_dim), corn_lgrA);
    }
}


BOOST_AUTO_TEST_CASE(neighboring_singleCellRefinements_x)
{
    // lgr1 and lgr2 grids mimick single-cell-refinements sharing I_FACE: | 0 | 1 |

    // Refine grids lgr1 and lgr2.
    const std::array<int, 3> lgr1_dim = {3,3,3};
    const std::array<int, 3> lgr2_dim = {4,3,3};
    // Number of subbivisions in y- and z- direction must coincide, when cells
    // from different LGRs share I_FACEs. In x-direction, they can differ.

    const auto& lgr1_ptr = createSingleCellGridAndRefine(lgr1_dim);
    const auto& lgr2_ptr = createSingleCellGridAndRefine(lgr2_dim);

    // Recall how the 8 corners of a cell are stored in cell_to_point_
    //     6 -- 7    TOP FACE
    //    /    /
    //   4 -- 5
    //     2 -- 3    BOTTOM FACE
    //    /    /
    //   0 -- 1
    std::unordered_map<int, int> lgrLeft_to_lgrRight = {{1,0}, {3,2}, {5,4}, {7,6}};

    // Illustration of cells on the boundary between LGR1 and LGR2, for k=1,
    // if LGR1 refines cell 0 and LGR2 refined cell 1.
    //      LGR1   LGR2
    //  15 16 17 | 20 21 22 23
    //  12 13 14 | 16 17 18 19
    //   9 10 11 | 12 13 14 15

    // elem 14 from lgr1 corners 1,2,5,7
    // should coincide with
    // elem 16 from lgr2 corners 0,3,4,6 respectively.
    const auto& elem14_lgr1 =  Dune::cpgrid::Entity<0>(*(lgr1_ptr), 14, true);
    const auto& elem16_lgr2 =  Dune::cpgrid::Entity<0>(*(lgr2_ptr), 16, true);

    check4CornersOfSharedRefinedFaceBetweenSingleCellRefinements(elem14_lgr1,
                                                                 elem16_lgr2,
                                                                 lgrLeft_to_lgrRight,
                                                                 lgr1_dim,
                                                                 lgr2_dim);

    // Illustration of cells on the boundary between LGR1 and LGR2, for k=1,
    // if LGR2 refines cell 0 and LGR1 refines cell 1.
    //         LGR2   LGR1
    //  20 21 22 23 | 15 16 17
    //  16 17 18 19 | 12 13 14
    //  12 13 14 15 |  9 10 11

    // elem 12 from lgr1 corners 1,2,5,7
    // should coincide with
    // elem 19 from lgr2 corners 0,3,4,6 respectively.
    const auto& elem19_lgr2 =  Dune::cpgrid::Entity<0>(*(lgr2_ptr), 19, true);
    const auto& elem12_lgr1 =  Dune::cpgrid::Entity<0>(*(lgr1_ptr), 12, true);

    check4CornersOfSharedRefinedFaceBetweenSingleCellRefinements(elem19_lgr2,
                                                                 elem12_lgr1,
                                                                 lgrLeft_to_lgrRight,
                                                                 lgr2_dim,
                                                                 lgr1_dim);
}

BOOST_AUTO_TEST_CASE(neighboring_singleCellRefinements_y)
{
    // lgr1 and lgr2 mimick a coarse grid having 2 cells sharing J_FACE:
    //   / 1 /
    //   ---   shared J_FACE
    //  / 0 /
    //  ---

    // Refine grids lgr1 and lgr2.
    const std::array<int, 3> lgr1_dim = {3,3,3};
    const std::array<int, 3> lgr2_dim = {3,4,3};
    // Number of subbivisions in x- and z- direction must coincide, when cells
    // from different LGRs share J_FACEs. In y-direction, they can differ.
    const auto& lgr1_ptr = createSingleCellGridAndRefine(lgr1_dim);
    const auto& lgr2_ptr = createSingleCellGridAndRefine(lgr2_dim);

    // Recall how the 8 corners of a cell are stored in cell_to_point_
    //     6 -- 7    TOP FACE
    //    /    /
    //   4 -- 5
    //     2 -- 3    BOTTOM FACE
    //    /    /
    //   0 -- 1
    std::unordered_map<int, int> lgrBack_to_lgrFront = {{2,0}, {3,1}, {6,4}, {7,5}};

    // Illustration of cells on the boundary between LGR1 and LGR2, for k=1,
    // if LGR1 refines cell 0 and LGR2 refined cell 1.
    //        21 22 23  LGR2
    //       18 19 20
    //      15 16 17
    //     12 13 14
    //    ---------
    //    15 16 17      LGR1
    //   12 13 14
    //  9 10 11


    // elem 16 from lgr1 corners 0,1,4,5
    // should coincide with
    // elem 13 from lgr2 corners 2,3,6,7 respectively.
    const auto& elem16_lgr1 =  Dune::cpgrid::Entity<0>(*(lgr1_ptr), 16, true);
    const auto& elem13_lgr2 =  Dune::cpgrid::Entity<0>(*(lgr2_ptr), 13, true);

    check4CornersOfSharedRefinedFaceBetweenSingleCellRefinements(elem16_lgr1,
                                                                 elem13_lgr2,
                                                                 lgrBack_to_lgrFront,
                                                                 lgr1_dim,
                                                                 lgr2_dim);

    // Illustration of cells on the boundary between LGR1 and LGR2, for k=1,
    // if LGR1 refines cell 1 and LGR2 refined cell 0.
    //        15 16 17  LGR1
    //       12 13 14
    //      9 10 11
    //     -------
    //    21 22 23      LGR2
    //   18 19 20
    //  15 16 17
    // 12 13 14

    // elem 10 from lgr1 corners 0,1,4,5
    // should coincide with
    // elem 22 from lgr2 corners 2,3,6,7 respectively.
    const auto& elem22_lgr2 =  Dune::cpgrid::Entity<0>(*(lgr2_ptr), 22, true);
    const auto& elem10_lgr1 =  Dune::cpgrid::Entity<0>(*(lgr1_ptr), 10, true);

    check4CornersOfSharedRefinedFaceBetweenSingleCellRefinements(elem22_lgr2,
                                                                 elem10_lgr1,
                                                                 lgrBack_to_lgrFront,
                                                                 lgr2_dim,
                                                                 lgr1_dim);
}

BOOST_AUTO_TEST_CASE(neighboring_singleCellRefinements_z)
{
    // lgr1 and lgr2  mimick a coarse grid has 2 cells sharing K_FACE:
    //   / 1 /
    //   ---   shared K_FACE
    //  / 0 /
    //  ---

    // Refine grids lgr1 and lgr2.
    const std::array<int, 3> lgr1_dim = {3,3,3};
    const std::array<int, 3> lgr2_dim = {3,3,4};
    // Number of subbivisions in x- and y- direction must coincide, when cells
    // from different LGRs share K_FACEs. In z-direction, they can differ.
    const auto& lgr1_ptr = createSingleCellGridAndRefine(lgr1_dim);
    const auto& lgr2_ptr = createSingleCellGridAndRefine(lgr2_dim);

    // Recall how the 8 corners of a cell are stored in cell_to_point_
    //     6 -- 7    TOP FACE
    //    /    /
    //   4 -- 5
    //     2 -- 3    BOTTOM FACE
    //    /    /
    //   0 -- 1
    std::unordered_map<int, int> lgrTop_to_lgrBottom = {{4,0}, {5,1}, {6,2}, {7,3}};

    // Illustration of cells on the boundary between LGR1 and LGR2,
    // if LGR1 refines cell 1 and LGR2 refined cell 0.
    //         6  7  8  LGR2 BOTTOM
    //        3  4  5
    //       0  1  2
    //      ------
    //     24 25 26  LGR1 TOP
    //    21 22 23
    //   18 19 20

    // elem 22 from lgr1 corners 4,5,6,7
    // should coincide with
    // elem 4 from lgr2 corners 0,1,2,3 respectively.
    const auto& elem22_lgr1 =  Dune::cpgrid::Entity<0>(*(lgr1_ptr), 22, true);
    const auto& elem4_lgr2 =  Dune::cpgrid::Entity<0>(*(lgr2_ptr), 4, true);

    check4CornersOfSharedRefinedFaceBetweenSingleCellRefinements(elem22_lgr1,
                                                                 elem4_lgr2,
                                                                 lgrTop_to_lgrBottom,
                                                                 lgr1_dim,
                                                                 lgr2_dim);

    // Illustration of cells on the boundary between LGR1 and LGR2,
    // if LGR1 refines cell 1 and LGR2 refined cell 0.
    //         6  7  8  LGR1 BOTTOM
    //        3  4  5
    //       0  1  2
    //     -------
    //    33 34 35      LGR2 TOP
    //   30 31 32
    //  27 28 29


    // elem 4 from lgr1 corners 4,5,6,7
    // should coincide with
    // elem 31 from lgr2 corners 0,1,2,3 respectively.
    const auto& elem31_lgr2 =  Dune::cpgrid::Entity<0>(*(lgr2_ptr), 31, true);
    const auto& elem4_lgr1 =  Dune::cpgrid::Entity<0>(*(lgr1_ptr), 4, true);

    check4CornersOfSharedRefinedFaceBetweenSingleCellRefinements(elem31_lgr2,
                                                                 elem4_lgr1,
                                                                 lgrTop_to_lgrBottom,
                                                                 lgr2_dim,
                                                                 lgr1_dim);
}
