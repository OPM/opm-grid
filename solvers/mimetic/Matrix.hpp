//===========================================================================
//
// File: Matrix.hpp
//
// Created: Tue Jun 30 11:25:46 2009
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

#ifndef OPENRS_MATRIX_HEADER
#define OPENRS_MATRIX_HEADER

#include <algorithm>
#include <ostream>
#include <vector>

#include <dune/grid/common/ErrorMacros.hpp>

#include <dune/solvers/mimetic/fortran.hpp>
#include <dune/solvers/mimetic/blas_lapack.hpp>

namespace Dune {

    template<typename T>
    class OwnData {
    public :
        OwnData(int sz, const T* data)
        {
            if (data) {
                data_.assign(data, data + sz);
            } else {
                data_.resize(sz);
            }
        }

        typedef T ValueType;

        ValueType&       operator[](int i)       { return data_[i]; }
        const ValueType& operator[](int i) const { return data_[i]; }
    private:
        std::vector<T> data_;
    };


    template<typename T>
    class SharedData {
    public:
        SharedData(int sz, T* data)
            : sz_(sz), data_(data)
        {
            ASSERT(data_ != 0);
        }

        typedef  T ValueType;

        ValueType&       operator[](int i)       { return data_[i]; }
        const ValueType& operator[](int i) const { return data_[i]; }
    private:
        int sz_;
        T*  data_;
    };


    template<typename T, template<typename> class StoragePolicy>
    class CMatrix : private StoragePolicy<T> {
    public:
        CMatrix(int rows, int cols, T* data)
            : StoragePolicy<T>(rows * cols, data)
            , rows_(rows), cols_(cols)
        {}

        typedef typename StoragePolicy<T>::ValueType ValueType;

        int      numRows()          const { return rows_;     }
        int      numCols()          const { return cols_;     }
        int      leadingDimension() const { return numCols(); }

        ValueType&       operator()(int row, int col)
        {
            return this->operator[](idx(row, col));
        }
        const ValueType& operator()(int row, int col) const
        {
            return this->operator[](idx(row, col));
        }

        ValueType*       data()       { return &this->operator[](0); }
        const ValueType* data() const { return &this->operator[](0); }
    private:
        int rows_, cols_;

        int idx(int row, int col) const
        {
            ASSERT ((0 <= row) && (row < numRows()));
            ASSERT ((0 <= col) && (col < numCols()));

            return row*numCols() + col;
        }
    };


    template<typename T, template<typename> class StoragePolicy>
    class FortranMatrix : private StoragePolicy<T> {
    public:
        FortranMatrix(int rows, int cols, T* data)
            : StoragePolicy<T>(rows * cols, data)
            , rows_(rows), cols_(cols)
        {}

        typedef typename StoragePolicy<T>::ValueType ValueType;

        int      numRows()          const { return rows_;     }
        int      numCols()          const { return cols_;     }
        int      leadingDimension() const { return numRows(); }

        T&       operator()(int row, int col)
        {
            return this->operator[](idx(row, col));
        }
        const T& operator()(int row, int col) const
        {
            return this->operator[](idx(row, col));
        }

        T*       data()       { return &this->operator[](0); }
        const T* data() const { return &this->operator[](0); }
    private:
        int rows_, cols_;

        int idx(int row, int col) const
        {
            ASSERT ((0 <= row) && (row < numRows()));
            ASSERT ((0 <= col) && (col < numCols()));

            return row + col*numRows();
        }
    };



    template<class Matrix>
    void zero(Matrix& A)
    {
        std::fill_n(A.data(), A.numRows() * A.numCols(),
                    typename Matrix::ValueType(0.0));
    }


    template<typename T, template<typename> class StoragePolicy>
    T trace(const FortranMatrix<T,StoragePolicy>& A)
    {
        T ret = 0.0;
        for (int i = 0; i < std::min(A.numRows(), A.numCols()); ++i) {
            ret += A(i,i);
        }
        return ret;
    }


    template<typename T, template<typename> class StoragePolicy>
    T trace(const CMatrix<T,StoragePolicy>& A)
    {
        const FortranMatrix<T,SharedData> At(A.numCols(), A.numRows(),
                                             const_cast<T*>(A.data()));
        return trace(At);
    }


    // A <- orth(A)
    template<typename T, template<typename> class StoragePolicy>
    int orthogonalizeColumns(FortranMatrix<T,StoragePolicy>& A)
    {
        typedef typename FortranMatrix<T,StoragePolicy>::ValueType ValueType;

        static std::vector<ValueType> tau;
        static std::vector<ValueType> work;

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
    template<typename T, template<typename> class StoragePolicy>
    void symmetricUpdate(const T&                              a1,
                         const FortranMatrix<T,StoragePolicy>& A ,
                         const T&                              a2,
                         FortranMatrix<T,StoragePolicy>&       C )
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
    template<typename T, template<typename> class StoragePolicy>
    void symmetricUpdate(const FortranMatrix<T,StoragePolicy>& A,
                         FortranMatrix<T,StoragePolicy>&       B)
    {
        typedef typename FortranMatrix<T,StoragePolicy>::ValueType ValueType;

        // B <- A*B
        Dune::BLAS_LAPACK::TRMM("Left" , "Upper", "No transpose", "Non-unit",
                                B.numRows(), B.numCols(), ValueType(1.0),
                                A.data(), A.leadingDimension(),
                                B.data(), B.leadingDimension());

        // B <- B*A (== A * B_orig * A)
        Dune::BLAS_LAPACK::TRMM("Right", "Upper", "No transpose", "Non-unit",
                                B.numRows(), B.numCols(), ValueType(1.0),
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


    template<typename                 T  ,
             template<typename> class SP1,
             template<typename> class SP2,
             template<typename> class SP3>
    void matMulAdd_NN(const T&                    a1,
                      const FortranMatrix<T,SP1>& A ,
                      const FortranMatrix<T,SP2>& B ,
                      const T&                    a2,
                      FortranMatrix<T,SP3>&       C)
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


    template<typename                 T  ,
             template<typename> class SP1,
             template<typename> class SP2,
             template<typename> class SP3>
    void matMulAdd_NT(const T&                    a1,
                      const FortranMatrix<T,SP1>& A ,
                      const FortranMatrix<T,SP2>& B ,
                      const T&                    a2,
                      FortranMatrix<T,SP3>&       C)
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


    template<typename                 T  ,
             template<typename> class SP1,
             template<typename> class SP2,
             template<typename> class SP3>
    void matMulAdd_NN(const T&                    a1,
                      const FortranMatrix<T,SP1>& A ,
                      const CMatrix<T,SP2>&       B ,
                      const T&                    a2,
                      FortranMatrix<T,SP3>&       C)
    {
        typedef typename CMatrix<T,SP2>::ValueType  ValueType;
        typedef FortranMatrix<ValueType,SharedData> FMat;

        const FMat Bt(B.numCols(), B.numRows(),
                      const_cast<ValueType*>(B.data()));

        matMulAdd_NT(a1, A, Bt, a2, C);
    }


    template<class charT, class traits,
             typename                 T,
             template<typename> class StoragePolicy>
    std::basic_ostream<charT,traits>&
    operator<<(std::basic_ostream<charT,traits>& os,
               const FortranMatrix<T,StoragePolicy>& A)
    {
        for (int i = 0; i < A.numRows(); ++i) {
            for (int j = 0; j < A.numCols(); ++j)
                os << A(i,j) << ' ';
            os << '\n';
        }

        return os;
    }


    template<class charT, class traits,
             typename                 T,
             template<typename> class StoragePolicy>
    std::basic_ostream<charT,traits>&
    operator<<(std::basic_ostream<charT,traits>& os,
               const CMatrix<T,StoragePolicy>& A)
    {
        for (int i = 0; i < A.numRows(); ++i) {
            for (int j = 0; j < A.numCols(); ++j)
                os << A(i,j) << ' ';
            os << '\n';
        }

        return os;
    }
} // namespace Dune
#endif // OPENRS_MATRIX_HEADER
