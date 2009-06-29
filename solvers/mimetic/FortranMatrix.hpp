//===========================================================================
//
// File: FortranMatrix.hpp
//
// Created: Fri Jun 19 09:51:00 2009
//
// Author(s): Bård Skaflestad     <bard.skaflestad@sintef.no>
//            Atgeirr F Rasmussen <atgeirr@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009 SINTEF ICT, Applied Mathematics.
  Copyright 2009 Statoil ASA.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPENRS_FORTRANMATRIX_HEADER
#define OPENRS_FORTRANMATRIX_HEADER

#include <algorithm>    // For std::fill_n().
#include <vector>

#include <dune/mimetic/solvers/fortran.hpp>
#include <dune/mimetic/solvers/blas_lapack.hpp>

namespace Dune {
    template<typename T,
             bool ownsItsData> class FortranMatrix;

    template<typename T>
    class FortranMatrix<T,false> {
    public:
        FortranMatrix(int rows, int cols, T* data)
            : rows_(rows), cols_(cols), data_(data)
        {}
        int      numRows()          const { return rows_;     }
        int      numCols()          const { return cols_;     }
        int      leadingDimension() const { return numRows(); }

        T&       operator()(int row, int col)
        {
            return data_[row + col*rows_];
        }
        const T& operator()(int row, int col) const
        {
            return data_[row + col*rows_];
        }

        // Caveat emptor [Don't access in C order].
        const T* data() const { return data_; }
        T*       data()       { return data_; }
    private:
        int rows_, cols_;
        T   *data_;
    };

    template<typename T>
    class FortranMatrix<T,true> {
    public:
        FortranMatrix(int rows, int cols)
            : rows_(rows), cols_(cols), data(rows * cols, T(0))
        {}
        int      numRows()          const { return rows_;     }
        int      numCols()          const { return cols_;     }
        int      leadingDimension() const { return numRows(); }

        T&       operator()(int row, int col)
        {
            return data_[row + col*rows_];
        }
        const T& operator()(int row, int col) const
        {
            return data_[row + col*rows_];
        }

        // Caveat emptor [Don't access in C order].
        const T* data() const { return &data_[0]; }
        T*       data()       { return &data_[0]; }
    private:
        int            rows_, cols_;
        std::vector<T> data_;
    };



    template<typename T>
    void zero(FortranMatrix<T>& A)
    {
        std::fill_n(A.data(), A.numRows() * A.numCols(), T(0.0));
    }


    // A <- orth(A)
    template<typename T>
    int orthogonalizeColumns(FortranMatrix<T>& A   ,
                             std::vector<T>&   tau ,
                             std::vector<T>&   work);
    {
        int info = 0;

        // Generate factorization.  Matrix Q is implicitly represented.
        Dune::BLAS_LAPACK::GEQRF(A.numRows(), A.numCols()               ,
                                 A.data()   , A.leadingDimension()      ,
                                 &tau[0]    , &work[0], int(work.size()),
                                 info);

        if (info == 0)
        {
            // QR factorization successfully generated--extract
            // explict representation of orthogonal matrix Q.
            Dune::BLAS_LAPACK::ORGQR(A.numRows(), A.numCols(), A.numCols()  ,
                                     A.data()   , A.leadingDimension()      ,
                                     &tau[0]    , &work[0], int(work.size()),
                                     info);
        }

        return info;
    }

    // C <- a1*A*A' + a2*C
    // Assumes T is an arithmetic (floating point) type, and that C==C'.
    template<typename T>
    void symmetricUpdate(const T&                a1,
                         const FortranMatrix<T>& A ,
                         const T&                a2,
                         FortranMatrix<T>&       C )
    {
        Dune::BLAS_LAPACK::SYRK("Upper"     , "No transpose"      ,
                                C.numRows() , A.numCols()         ,
                                a1, A.data(), A.leadingDimension(),
                                a2, C.data(), C.leadingDimension());
    }


    // B <- A*B*A'
    // Assumes T is an arithmetic (floating point) type, and that A==A'.
    template<typename T>
    void symmetricUpdate(const FortranMatrix<T>& A, FortranMatrix<T>& B)
    {
        // B <- A*B
        Dune::BLAS_LAPACK::TRMM("Left" , "Upper", "No transpose", "Non-unit",
                                B.numRows(), B.numCols(), T(1.0),
                                A.data(), A.leadingDimension(),
                                B.data(), B.leadingDimension());

        // B <- B*A (== A * B_orig * A)
        Dune::BLAS_LAPACK::TRMM("Right", "Upper", "No transpose", "Non-unit",
                                B.numRows(), B.numCols(), T(1.0),
                                A.data(), A.leadingDimension(),
                                B.data(), B.leadingDimension());
    }


    template<typename T, bool transA, bool transB>
    void matMulAdd(const T&                a1,
                   const FortranMatrix<T>& A ,
                   const FortranMatrix<T>& B ,
                   const T&                a2,
                   FortranMatrix<T>&       C);


    template<typename T>
    void matMulAdd<T,false,true>(const T&                a1,
                                 const FortranMatrix<T>& A ,
                                 const FortranMatrix<T>& B ,
                                 const T&                a2,
                                 FortranMatrix<T>&       C)
    {
        ASSERT(A.numRows() == C.numRows());
        ASSERT(B.numRows() == C.numCols());
        ASSERT(A.numCols() == B.numCols());

        GEMM("No Transpose", "Transpose", A.numRows(), B.numCols(), A.numCols(),
             a1, A.data(), A.leadingDimension(), B.data(), B.leadingDimension(),
             a2, C.data(), C.leadingDimension());
    }


} // namespace Dune
#endif // OPENRS_FORTRANMATRIX_HEADER
