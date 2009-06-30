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
#include <ostream>
#include <vector>

#include <dune/grid/common/ErrorMacros.hpp>

#include <dune/solvers/mimetic/fortran.hpp>
#include <dune/solvers/mimetic/blas_lapack.hpp>

namespace Dune {
    template<typename T,
             bool ownsItsData> class FortranMatrix;

    template<typename T>
    class FortranMatrix<T,false> {
    public:
        typedef T value_type;

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
        typedef T value_type;

        FortranMatrix(int rows, int cols)
            : rows_(rows), cols_(cols), data_(rows * cols, T(0))
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



    template<class Matrix>
    void zero(Matrix& A)
    {
        std::fill_n(A.data(), A.numRows() * A.numCols(),
                    typename Matrix::value_type(0.0));
    }


    // A <- orth(A)
    template<class Matrix>
    int orthogonalizeColumns(Matrix& A)
    {
        static std::vector<typename Matrix::value_type> tau;
        static std::vector<typename Matrix::value_type> work;

        if (tau .size() <      A.numCols()) tau .resize(     A.numCols());
        if (work.size() < 64 * A.numRows()) work.resize(64 * A.numRows());  // 64 from ILAENV

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
    template<class Matrix>
    void symmetricUpdate(const typename Matrix::value_type& a1,
                         const Matrix&                      A ,
                         const typename Matrix::value_type& a2,
                         Matrix&                            C )
    {
        Dune::BLAS_LAPACK::SYRK("Upper"     , "No transpose"      ,
                                C.numRows() , A.numCols()         ,
                                a1, A.data(), A.leadingDimension(),
                                a2, C.data(), C.leadingDimension());

        // Account for SYRK (in this case) only updating the upper
        // leading n-by-n submatrix.
        //
        for (int j = 0; j < C.numCols(); ++j) {
            for (int i = j+1; i < C.numRows(); ++i) {
                C(i,j) = C(j,i);
            }
        }
    }


    // B <- A*B*A'
    // Assumes T is an arithmetic (floating point) type, and that A==A'.
    template<class Matrix>
    void symmetricUpdate(const Matrix& A, Matrix& B)
    {
        // B <- A*B
        Dune::BLAS_LAPACK::TRMM("Left" , "Upper", "No transpose", "Non-unit",
                                B.numRows(), B.numCols(),
                                typename Matrix::value_type(1.0),
                                A.data(), A.leadingDimension(),
                                B.data(), B.leadingDimension());

        // B <- B*A (== A * B_orig * A)
        Dune::BLAS_LAPACK::TRMM("Right", "Upper", "No transpose", "Non-unit",
                                B.numRows(), B.numCols(),
                                typename Matrix::value_type(1.0),
                                A.data(), A.leadingDimension(),
                                B.data(), B.leadingDimension());

        // Account for TRMM (in this case) only updating the upper
        // leading n-by-n submatrix.
        //
        for (int j = 0; j < B.numCols(); ++j) {
            for (int i = j+1; i < B.numRows(); ++i) {
                B(i,j) = B(j,i);
            }
        }
    }


    template<class Matrix>
    void matMulAdd_NN(const typename Matrix::value_type& a1,
                      const Matrix&                      A ,
                      const Matrix&                      B ,
                      const typename Matrix::value_type& a2,
                      Matrix&                            C)
    {
        ASSERT(A.numRows() == C.numRows());
        ASSERT(A.numCols() == B.numRows());
        ASSERT(B.numCols() == C.numCols());

        Dune::BLAS_LAPACK::GEMM("No Transpose", "No Transpose",
                                A.numRows(), B.numCols(), A.numCols(),
                                a1, A.data(), A.leadingDimension(),
                                B.data(), B.leadingDimension(),
                                a2, C.data(), C.leadingDimension());
    }


    template<class Matrix>
    void matMulAdd_NT(const typename Matrix::value_type& a1,
                      const Matrix&                      A ,
                      const Matrix&                      B ,
                      const typename Matrix::value_type& a2,
                      Matrix&                            C)
    {
        ASSERT(A.numRows() == C.numRows());
        ASSERT(B.numRows() == C.numCols());
        ASSERT(A.numCols() == B.numCols());

        Dune::BLAS_LAPACK::GEMM("No Transpose", "Transpose",
                                A.numRows(), B.numRows(), A.numCols(),
                                a1, A.data(), A.leadingDimension(),
                                B.data(), B.leadingDimension(),
                                a2, C.data(), C.leadingDimension());
    }


    template<class charT, class traits, class T, bool ownsItsData>
    std::basic_ostream<charT,traits>&
    operator<<(std::basic_ostream<charT,traits>& os, const FortranMatrix<T, ownsItsData>& A)
    {
        for (int i = 0; i < A.numRows(); ++i) {
            for (int j = 0; j < A.numCols(); ++j)
                os << A(i,j) << ' ';
            os << '\n';
        }

        return os;
    }


} // namespace Dune
#endif // OPENRS_FORTRANMATRIX_HEADER
