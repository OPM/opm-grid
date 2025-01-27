//===========================================================================
//
// File: replace_lgr1_face_idx_by_lgr2_face_idx_test.cpp
//
// Created: Monday 27.01.2025 10:32:00
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

#define BOOST_TEST_MODULE ReplaceLgr1FaceIdxByLgr2FaceIdxTests
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif


#include <opm/grid/CpGrid.hpp>

#include <array>
#include <memory>
#include <stdexcept>


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

   To test replaceLgr1FaceIdxByLgr2FaceIdx(...), we create two grids, each of them
   containing only one cell. Add different LGRs in each grid, in different cases to
   create all possible escenarios (the two coarse cells sharing I_FACE, J_FACE or K_FACE).

   To avoid friend declarations, we 'copy' here replaceLgr1FaceIdxByLgr2FaceIdx
   (and getRefinedFaceIJK).

   Why single-cell-refinements instead of adding LGRs in one grid?:
   replaceLgr1FaceIdxByLgr2FaceIdx is a private method in CpGrid, meant to be
   used between independent/un-related single-cell-refinements (stored in CpGridData objects).
   Therefore, creating a grid, adding LGRs to it, and checking face indices afterwards in
   refined cells on bouandary of the LGRs is not the escenario where this method should be tested.
*/


std::array<int,3> getRefinedFaceIJK(const std::array<int,3>& cells_per_dim,
                                    int faceIdxInLgr,
                                    const std::shared_ptr<Dune::cpgrid::CpGridData>& elemLgr_ptr)
{
    // Order defined in Geometry::refine
    // K_FACES  (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i
    // I_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1))
    //           + (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j
    // J_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2] +1))
    //                    + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2])
    //                    + (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2]) + k
    const auto& i_faces =  (cells_per_dim[0] +1)*cells_per_dim[1]*cells_per_dim[2];
    const auto& j_faces =  cells_per_dim[0]*(cells_per_dim[1]+1)*cells_per_dim[2];
    const auto& k_faces =  cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1);

    if (faceIdxInLgr >= i_faces + j_faces + k_faces) {
        OPM_THROW(std::logic_error, "Invalid face index from single-cell-refinement.\n");
    }

    const auto& faceEntity =  Dune::cpgrid::EntityRep<1>(faceIdxInLgr, true);
    const auto& faceTag = elemLgr_ptr ->faceTag(faceEntity);
    std::array<int,3> ijk;
    switch (faceTag) {
    case I_FACE:
        faceIdxInLgr -= (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1));
        // faceIdxInLgr =  (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j
        ijk[1] = faceIdxInLgr % cells_per_dim[1];
        faceIdxInLgr -= ijk[1]; // (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1])
        faceIdxInLgr /= cells_per_dim[1]; // (i*cells_per_dim[2]) + k
        ijk[2] = faceIdxInLgr % cells_per_dim[2];
        faceIdxInLgr -=ijk[2]; // i*cells_per_dim[2]
        ijk[0] = faceIdxInLgr / cells_per_dim[2];
        break;
    case J_FACE:
        faceIdxInLgr -=  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2] +1))
            + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2]);
        // faceIdxInLgr =  (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2]) + k
        ijk[2] = faceIdxInLgr % cells_per_dim[2];
        faceIdxInLgr -= ijk[2]; // (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2])
        faceIdxInLgr /= cells_per_dim[2]; // (j*cells_per_dim[0]) + i
        ijk[0] = faceIdxInLgr % cells_per_dim[0];
        faceIdxInLgr -=ijk[0]; // j*cells_per_dim[0]
        ijk[1] = faceIdxInLgr / cells_per_dim[0];
        break;
    case K_FACE:
        //  (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i
        ijk[0] = faceIdxInLgr % cells_per_dim[0];
        faceIdxInLgr -= ijk[0]; // (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0])
        faceIdxInLgr /= cells_per_dim[0]; // (k*cells_per_dim[1]) + j
        ijk[1] = faceIdxInLgr % cells_per_dim[1];
        faceIdxInLgr -=ijk[1]; // k*cells_per_dim[1]
        ijk[2] = faceIdxInLgr / cells_per_dim[1];
        break;
    default:
        OPM_THROW(std::logic_error, "FaceTag is not I, J, or K!");
    }
    return ijk;
}


int replaceLgr1FaceIdxByLgr2FaceIdx(const std::array<int,3>& cells_per_dim_lgr1,
                                    int faceIdxInLgr1,
                                    const std::shared_ptr<Dune::cpgrid::CpGridData>& elemLgr1_ptr,
                                    const std::array<int,3>& cells_per_dim_lgr2)
{
    const auto& ijkLgr1 = getRefinedFaceIJK(cells_per_dim_lgr1, faceIdxInLgr1, elemLgr1_ptr);
    // lgr1 represents an element index < lgr2 (neighboring cells sharing a face with lgr1-element)
    // Order defined in Geometry::refine
    // K_FACES (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i
    // I_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1))
    //           + (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j
    // J_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2] +1))
    //                    + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2])
    //                    + (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2]) + k
    const int& kFacesLgr2 = cells_per_dim_lgr2[0]*cells_per_dim_lgr2[1]*(cells_per_dim_lgr2[2]+1);
    const int& iFacesLgr2 = ((cells_per_dim_lgr2[0]+1)*cells_per_dim_lgr2[1]*cells_per_dim_lgr2[2]);

    const auto& face_lgr1 =  Dune::cpgrid::EntityRep<1>(faceIdxInLgr1, true);
    const auto& face_tag = elemLgr1_ptr-> faceTag(face_lgr1);

    if (face_tag == I_FACE) {
        assert( cells_per_dim_lgr1[1] == cells_per_dim_lgr2[1]);
        assert( cells_per_dim_lgr1[2] == cells_per_dim_lgr2[2]);
        if (ijkLgr1[0] == cells_per_dim_lgr1[0]) { // same j,k, but i = 0
            return  kFacesLgr2 + (ijkLgr1[2]*cells_per_dim_lgr2[1]) + ijkLgr1[1];
        }
        else { // same j,k, but i = cells_per_dim[0]
            return  kFacesLgr2 + (cells_per_dim_lgr2[0]*cells_per_dim_lgr2[1]*cells_per_dim_lgr2[2])
                + (ijkLgr1[2]*cells_per_dim_lgr2[1]) + ijkLgr1[1];
        }
    }
    if (face_tag == J_FACE) {
        assert( cells_per_dim_lgr1[0] == cells_per_dim_lgr2[0]);
        assert( cells_per_dim_lgr1[2] == cells_per_dim_lgr2[2]);
        if (ijkLgr1[1] == cells_per_dim_lgr1[1]) { // same i,k, but j = 0
            return kFacesLgr2 + iFacesLgr2 + (ijkLgr1[0]*cells_per_dim_lgr2[2]) + ijkLgr1[2];
        }
        else { // same i,k, but j = cells_per_dim[1]
            return kFacesLgr2 + iFacesLgr2 + (cells_per_dim_lgr2[1]*cells_per_dim_lgr2[0]*cells_per_dim_lgr2[2])
                + (ijkLgr1[0]*cells_per_dim_lgr2[2]) + ijkLgr1[2];
        }
    }
    if (face_tag == K_FACE) {
        assert( cells_per_dim_lgr1[0] == cells_per_dim_lgr2[0]);
        assert( cells_per_dim_lgr1[1] == cells_per_dim_lgr2[1]);
        if (ijkLgr1[2] == cells_per_dim_lgr1[2]) { // same i,j, but k = 0
            return  (ijkLgr1[1]*cells_per_dim_lgr2[0]) + ijkLgr1[0];
        }
        else{ // same i, j, but k = cells_per_dim[2]
            return  (cells_per_dim_lgr2[2]*cells_per_dim_lgr2[0]*cells_per_dim_lgr2[1])
                + (ijkLgr1[1]*cells_per_dim_lgr2[0]) + ijkLgr1[0];
        }
    }
    OPM_THROW(std::logic_error,  "Cannot convert face index from one LGR to its neighboring LGR.");
}


const std::shared_ptr<Dune::cpgrid::CpGridData> createSingleCellGridAndRefine(const std::array<int,3>& lgr_dim)
{
    Dune::CpGrid lgr;
    // Create a single-cell-grid
    const std::array<double,3>& cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int,3>& coarse_grid_dim = {1,1,1};
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


BOOST_AUTO_TEST_CASE(neighboring_singleCellRefinements_x)
{
    // lgr1 and lgr2 grids mimick single-cell-refinements sharing I_FACE: | 0 | 1 |

    // Refine grids lgr1 and lgr2.
    const std::array<int, 3>& lgr1_dim = {3,3,3};
    const std::array<int, 3>& lgr2_dim = {4,3,3};
    // Number of subbivisions in y- and z- direction must coincide, when cells
    // from different LGRs share I_FACEs. In x-direction, they can differ.

    const auto& lgr1_ptr = createSingleCellGridAndRefine(lgr1_dim);
    const auto& lgr2_ptr = createSingleCellGridAndRefine(lgr2_dim);

    // Illustration of cells on the boundary between LGR1 and LGR2, for k=1,
    // if LGR1 refines cell 0 and LGR2 refined cell 1.
    //      LGR1   LGR2
    //  15 16 17 | 20 21 22 23
    //  12 13 14 | 16 17 18 19
    //   9 10 11 | 12 13 14 15

    // elem 14 from lgr1 face {I true}
    // should coincide with
    // elem 16 from lgr2 face {I false}

    // Recall how the 6 faces of a cell are stored in cell_to_face_
    // in Geometry::refine
    // cell_to_face_ [ element ] = {{K false}, {J false}, {I false}, {I true}, {J true}, {K true}}.
    const auto& faceTrue_lgr1 = lgr1_ptr-> cellFace( /*elemIdx*/ 14, /*local face index*/ 3); // I true for elem14_lgr1
    const auto& faceFalse_lgr2 = lgr2_ptr -> cellFace( /*elemIdx*/ 16, /*local face index*/ 2); // I false for elem16_lgr2

    BOOST_CHECK_EQUAL( replaceLgr1FaceIdxByLgr2FaceIdx(lgr1_dim, faceTrue_lgr1, lgr1_ptr, lgr2_dim), faceFalse_lgr2);
    BOOST_CHECK_EQUAL( replaceLgr1FaceIdxByLgr2FaceIdx(lgr2_dim, faceFalse_lgr2, lgr2_ptr, lgr1_dim), faceTrue_lgr1);

    // Illustration of cells on the boundary between LGR1 and LGR2, for k=1,
    // if LGR2 refines cell 0 and LGR1 refines cell 1.
    //         LGR2   LGR1
    //  20 21 22 23 | 15 16 17
    //  16 17 18 19 | 12 13 14
    //  12 13 14 15 |  9 10 11

    // elem 12 from lgr1 face {I true}
    // should coincide with
    // elem 19 from lgr2 face {I false}

    const auto& faceTrue_lgr2 = lgr2_ptr-> cellFace( /*elemIdx*/ 19, /*local face index*/ 3);; // I true for elem19_lgr2
    const auto& faceFalse_lgr1 = lgr1_ptr -> cellFace( /*elemIdx*/ 12, /*local face index*/ 2); // I false for elem12_lgr1

    BOOST_CHECK_EQUAL( replaceLgr1FaceIdxByLgr2FaceIdx(lgr2_dim, faceTrue_lgr2, lgr2_ptr, lgr1_dim), faceFalse_lgr1);
    BOOST_CHECK_EQUAL( replaceLgr1FaceIdxByLgr2FaceIdx(lgr1_dim, faceFalse_lgr1, lgr1_ptr, lgr2_dim), faceTrue_lgr2);

    // lgr1 has 108 faces (with indices 0, ..., 107).
    // lgr2 has 141 faces (with indices 0, ..., 140).
    const auto& non_existing_face = 141; // non exisitng face index for both lgrs.
    BOOST_CHECK_THROW( replaceLgr1FaceIdxByLgr2FaceIdx(lgr1_dim, non_existing_face, lgr1_ptr, lgr2_dim), std::logic_error);
    BOOST_CHECK_THROW( replaceLgr1FaceIdxByLgr2FaceIdx(lgr2_dim, non_existing_face, lgr2_ptr, lgr1_dim), std::logic_error);
}

BOOST_AUTO_TEST_CASE(neighboring_singleCellRefinements_y)
{
    // lgr1 and lgr2 mimick a coarse grid having 2 cells sharing J_FACE:
    //   / 1 /
    //   ---   shared J_FACE
    //  / 0 /
    //  ---

    // Refine grids lgr1 and lgr2.
    const std::array<int, 3>& lgr1_dim = {3,3,3};
    const std::array<int, 3>& lgr2_dim = {3,4,3};
    // Number of subbivisions in x- and z- direction must coincide, when cells
    // from different LGRs share J_FACEs. In y-direction, they can differ.
    const auto& lgr1_ptr = createSingleCellGridAndRefine(lgr1_dim);
    const auto& lgr2_ptr = createSingleCellGridAndRefine(lgr2_dim);

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

    // elem 16 from lgr1 face {J true}
    // should coincide with
    // elem 13 from lgr2 face {J false}

    // Recall how the 6 faces of a cell are stored in cell_to_face_
    // in Geometry::refine
    // cell_to_face_ [ element ] = {{K false}, {J false}, {I false}, {I true}, {J true}, {K true}}.
    const auto& faceTrue_lgr1 = lgr1_ptr-> cellFace( /*elemIdx*/ 16, /*local face index*/ 4); // J true for elem16_lgr1
    const auto& faceFalse_lgr2 = lgr2_ptr -> cellFace( /*elemIdx*/ 13, /*local face index*/ 1); // J false for elem13_lgr2

    BOOST_CHECK_EQUAL( replaceLgr1FaceIdxByLgr2FaceIdx(lgr1_dim, faceTrue_lgr1, lgr1_ptr, lgr2_dim), faceFalse_lgr2);
    BOOST_CHECK_EQUAL( replaceLgr1FaceIdxByLgr2FaceIdx(lgr2_dim, faceFalse_lgr2, lgr2_ptr, lgr1_dim), faceTrue_lgr1);


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

    // elem 10 from lgr1 face {J false}
    // should coincide with
    // elem 22 from lgr2 face {J true}

    const auto& faceTrue_lgr2 = lgr2_ptr-> cellFace( /*elemIdx*/ 22, /*local face index*/ 4); // J true for elem22_lgr2
    const auto& faceFalse_lgr1 = lgr1_ptr -> cellFace( /*elemIdx*/ 10, /*local face index*/ 1); // J false for elem10_lgr1

    BOOST_CHECK_EQUAL( replaceLgr1FaceIdxByLgr2FaceIdx(lgr2_dim, faceTrue_lgr2, lgr2_ptr, lgr1_dim), faceFalse_lgr1);
    BOOST_CHECK_EQUAL( replaceLgr1FaceIdxByLgr2FaceIdx(lgr1_dim, faceFalse_lgr1, lgr1_ptr, lgr2_dim), faceTrue_lgr2);

}

BOOST_AUTO_TEST_CASE(neighboring_singleCellRefinements_z)
{
    // lgr1 and lgr2  mimick a coarse grid has 2 cells sharing K_FACE:
    //   / 1 /
    //   ---   shared K_FACE
    //  / 0 /
    //  ---

    // Refine grids lgr1 and lgr2.
    const std::array<int, 3>& lgr1_dim = {3,3,3};
    const std::array<int, 3>& lgr2_dim = {3,3,4};
    // Number of subbivisions in x- and y- direction must coincide, when cells
    // from different LGRs share K_FACEs. In z-direction, they can differ.
    const auto& lgr1_ptr = createSingleCellGridAndRefine(lgr1_dim);
    const auto& lgr2_ptr = createSingleCellGridAndRefine(lgr2_dim);

    // Illustration of cells on the boundary between LGR1 and LGR2,
    // if LGR1 refines cell 1 and LGR2 refined cell 0.
    //         6  7  8  LGR2 BOTTOM
    //        3  4  5
    //       0  1  2
    //      ------
    //     24 25 26  LGR1 TOP
    //    21 22 23
    //   18 19 20

    // elem 22 from lgr1 {K true}
    // should coincide with
    // elem 4 from lgr2 {K false}

    // Recall how the 6 faces of a cell are stored in cell_to_face_
    // in Geometry::refine
    // cell_to_face_ [ element ] = {{K false}, {J false}, {I false}, {I true}, {J true}, {K true}}.
    const auto& faceTrue_lgr1 = lgr1_ptr-> cellFace( /*elemIdx*/ 22, /*local face index*/ 5); // J true for elem22_lgr1
    const auto& faceFalse_lgr2 = lgr2_ptr -> cellFace( /*elemIdx*/ 4, /*local face index*/ 0); // K false for elem4_lgr2

    BOOST_CHECK_EQUAL( replaceLgr1FaceIdxByLgr2FaceIdx(lgr1_dim, faceTrue_lgr1, lgr1_ptr, lgr2_dim), faceFalse_lgr2);
    BOOST_CHECK_EQUAL( replaceLgr1FaceIdxByLgr2FaceIdx(lgr2_dim, faceFalse_lgr2, lgr2_ptr, lgr1_dim), faceTrue_lgr1);


    // Illustration of cells on the boundary between LGR1 and LGR2,
    // if LGR1 refines cell 1 and LGR2 refined cell 0.
    //         6  7  8  LGR1 BOTTOM
    //        3  4  5
    //       0  1  2
    //     -------
    //    33 34 35      LGR2 TOP
    //   30 31 32
    //  27 28 29


    // elem 4 from lgr1 face {K false}
    // should coincide with
    // elem 31 from lgr2 {K true}

    const auto& faceTrue_lgr2 = lgr2_ptr-> cellFace( /*elemIdx*/ 31, /*local face index*/ 5);; // K true for elem31_lgr2
    const auto& faceFalse_lgr1 = lgr1_ptr -> cellFace( /*elemIdx*/ 4, /*local face index*/ 0); // K false for elem4_lgr1

    BOOST_CHECK_EQUAL( replaceLgr1FaceIdxByLgr2FaceIdx(lgr2_dim, faceTrue_lgr2, lgr2_ptr, lgr1_dim), faceFalse_lgr1);
    BOOST_CHECK_EQUAL( replaceLgr1FaceIdxByLgr2FaceIdx(lgr1_dim, faceFalse_lgr1, lgr1_ptr, lgr2_dim), faceTrue_lgr2);
}

