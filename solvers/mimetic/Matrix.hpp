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
#include <boost/bind.hpp>

#include <dune/common/fvector.hh>

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

        T&       operator[](int i)       { return data_[i]; }
        const T& operator[](int i) const { return data_[i]; }

        int size() const { return data_.size(); }
    private:
        std::vector<T> data_;
    };


    template<typename T>
    class SharedData {
    public:
        SharedData(int sz, T* data)
            : sz_(sz), data_(data)
        {
            ASSERT((sz == 0) == (data == 0));
        }

        T&       operator[](int i)       { return data_[i]; }
        const T& operator[](int i) const { return data_[i]; }

        int size() const { return sz_; }
    private:
        int sz_;
        T*  data_;
    };


    template<typename T>
    class ImmutableSharedData {
    public:
        ImmutableSharedData(int sz, const T* data)
            : sz_(sz), data_(data)
        {
            ASSERT(data_ != 0);
        }

        const T& operator[](int i) const { return data_[i]; }

        int size() const { return sz_; }
    private:
        int sz_;
        const T*  data_;
    };


    class OrderingBase {
    public:
        OrderingBase()
            : rows_(0), cols_(0)
        {}

        OrderingBase(int rows, int cols)
            : rows_(rows), cols_(cols)
        {}

        int numRows() const { return rows_; }
        int numCols() const { return cols_; }

    private:
        int rows_, cols_;
    };


    class COrdering : public OrderingBase {
    public:
        COrdering()
            : OrderingBase()
        {}

        COrdering(int rows, int cols)
            : OrderingBase(rows, cols)
        {}

        int leadingDimension() const { return numCols(); }

        int idx(int row, int col) const
        {
            ASSERT ((0 <= row) && (row < numRows()));
            ASSERT ((0 <= col) && (col < numCols()));

            return row*numCols() + col;
        }
    };


    class FortranOrdering : public OrderingBase {
    public:
        FortranOrdering()
            : OrderingBase()
        {}

        FortranOrdering(int rows, int cols)
            : OrderingBase(rows, cols)
        {}

        int leadingDimension() const { return numRows(); }

        int idx(int row, int col) const
        {
            ASSERT ((0 <= row) && (row < numRows()));
            ASSERT ((0 <= col) && (col < numCols()));

            return row + col*numRows();
        }
    };


    template<typename                 T,
             template<typename> class StoragePolicy,
             class                    OrderingPolicy>
    class FullMatrix : private StoragePolicy<T>,
                       private OrderingPolicy
    {
    public:
        FullMatrix()
            : StoragePolicy<T>(0, 0),
              OrderingPolicy()
        {}

        template <typename DataPointer>
        FullMatrix(int rows, int cols, DataPointer data)
            : StoragePolicy<T>(rows * cols, data),
              OrderingPolicy(rows, cols)
        {}

        template <template<typename> class OtherSP>
        explicit FullMatrix(const FullMatrix<T, OtherSP, OrderingPolicy>& m)
            : StoragePolicy<T>(m.numRows()*m.numCols(), m.data()),
              OrderingPolicy(m.numRows(), m.numCols())
        {
        }

        template <template<typename> class OtherSP>
        void operator+= (const FullMatrix<T, OtherSP, OrderingPolicy>& m)
        {
            ASSERT(numRows() == m.numRows() && numCols() == m.numCols());
            std::transform(data(), data() + this->size(),
                           m.data(), data(), std::plus<T>());
        }

        void operator*= (const T& scalar)
        {
            std::transform(data(), data() + this->size(),
                           data(), boost::bind(std::multiplies<T>(), _1, scalar));
        }

        typedef T value_type;

        using OrderingPolicy::numRows;
        using OrderingPolicy::numCols;
        using OrderingPolicy::leadingDimension;

        value_type&       operator()(int row, int col)
        {
            return this->operator[](this->idx(row, col));
        }
        const value_type& operator()(int row, int col) const
        {
            return this->operator[](this->idx(row, col));
        }

        value_type*       data()       { return &this->operator[](0); }
        const value_type* data() const { return &this->operator[](0); }
    };


    // Convenience typedefs
    typedef FullMatrix<double, OwnData,             COrdering>        OwnCMatrix;
    typedef FullMatrix<double, SharedData,          COrdering>        SharedCMatrix;
    typedef FullMatrix<double, ImmutableSharedData, COrdering>        ImmutableCMatrix;


    typedef FullMatrix<double, OwnData,             FortranOrdering>  OwnFortranMatrix;
    typedef FullMatrix<double, SharedData,          FortranOrdering>  SharedFortranMatrix;
    typedef FullMatrix<double, ImmutableSharedData, FortranOrdering>  ImmutableFortranMatrix;


#if 0
    template<typename T, template<typename> class StoragePolicy>
    class CMatrix : private StoragePolicy<T> {
    public:
        CMatrix()
            : StoragePolicy<T>(0, 0),
              rows_(0), cols_(0)
        {}

        template <typename DataPointer>
        CMatrix(int rows, int cols, DataPointer data)
            : StoragePolicy<T>(rows * cols, data),
              rows_(rows), cols_(cols)
        {}

        template <class OtherMatrixType>
        CMatrix(const OtherMatrixType& m)
            : StoragePolicy<T>(m.numRows()*m.numCols(), m.data()),
              rows_(m.numRows()), cols_(m.numCols())
        {
        }

        template <template<typename> class OtherStoragePolicy>
        void operator+= (const CMatrix<T,OtherStoragePolicy>& m)
        {
            ASSERT(numRows() == m.numRows() && numCols() == m.numCols());
            std::transform(data(), data() + this->size(),
                           m.data(), data(), std::plus<T>());
        }

        void operator*= (const T& scalar)
        {
            std::transform(data(), data() + this->size(),
                           data(), boost::bind(std::multiplies<T>(), _1, scalar));
        }

        typedef T value_type;

        int      numRows()          const { return rows_;     }
        int      numCols()          const { return cols_;     }
        int      leadingDimension() const { return numCols(); }

        value_type&       operator()(int row, int col)
        {
            return this->operator[](idx(row, col));
        }
        const value_type& operator()(int row, int col) const
        {
            return this->operator[](idx(row, col));
        }

        value_type*       data()       { return &this->operator[](0); }
        const value_type* data() const { return &this->operator[](0); }
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
        template<typename DataPointer>
        FortranMatrix(int rows, int cols, DataPointer data)
            : StoragePolicy<T>(rows * cols, data)
            , rows_(rows), cols_(cols)
        {}

        typedef T value_type;

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
#endif


    template<class Matrix>
    void zero(Matrix& A)
    {
        std::fill_n(A.data(), A.numRows() * A.numCols(),
                    typename Matrix::value_type(0.0));
    }


    template<class Matrix>
    typename Matrix::value_type
    trace(const Matrix& A)
    {
        typename Matrix::value_type ret(0);

        for (int i = 0; i < std::min(A.numRows(), A.numCols()); ++i) {
            ret += A(i,i);
        }
        return ret;
    }


    template<class Matrix, int rows>
    FieldVector<typename Matrix::value_type, rows>
    prod(const Matrix& A, const FieldVector<typename Matrix::value_type,rows>& x)
    {
        const int cols = rows;
        ASSERT (A.numRows() == rows);
        ASSERT (A.numCols() == cols);

        FieldVector<typename Matrix::value_type, rows> res(0.0);
        for (int c = 0; c < cols; ++c) {
            for (int r = 0; r < rows; ++r) {
                res[r] += A(r, c)*x[c];
            }
        }
        return res;
    }


    // A <- orth(A)
    template<typename T, template<typename> class StoragePolicy>
    int orthogonalizeColumns(FullMatrix<T,StoragePolicy,FortranOrdering>& A)
    {
        typedef typename FullMatrix<T,StoragePolicy,FortranOrdering>::value_type value_type;

        static std::vector<value_type> tau;
        static std::vector<value_type> work;

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
    void symmetricUpdate(const T&                                           a1,
                         const FullMatrix<T,StoragePolicy,FortranOrdering>& A ,
                         const T&                                           a2,
                         FullMatrix<T,StoragePolicy,FortranOrdering>&       C )
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
    void symmetricUpdate(const FullMatrix<T,StoragePolicy,FortranOrdering>& A,
                         FullMatrix<T,StoragePolicy,FortranOrdering>&       B)
    {
        typedef typename FullMatrix<T,StoragePolicy,FortranOrdering>::value_type value_type;

        // B <- A*B
        Dune::BLAS_LAPACK::TRMM("Left" , "Upper", "No transpose", "Non-unit",
                                B.numRows(), B.numCols(), value_type(1.0),
                                A.data(), A.leadingDimension(),
                                B.data(), B.leadingDimension());

        // B <- B*A (== A * B_orig * A)
        Dune::BLAS_LAPACK::TRMM("Right", "Upper", "No transpose", "Non-unit",
                                B.numRows(), B.numCols(), value_type(1.0),
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
    void matMulAdd_NN(const T&                             a1,
                      const FullMatrix<T,SP1,FortranOrdering>& A ,
                      const FullMatrix<T,SP2,FortranOrdering>& B ,
                      const T&                             a2,
                      FullMatrix<T,SP3,FortranOrdering>&       C)
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
    void matMulAdd_NT(const T&                                 a1,
                      const FullMatrix<T,SP1,FortranOrdering>& A ,
                      const FullMatrix<T,SP2,FortranOrdering>& B ,
                      const T&                                 a2,
                      FullMatrix<T,SP3,FortranOrdering>&       C)
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
    void matMulAdd_NN(const T&                                 a1,
                      const FullMatrix<T,SP1,FortranOrdering>& A ,
                      const FullMatrix<T,SP2,COrdering>&       B ,
                      const T&                                 a2,
                      FullMatrix<T,SP3,FortranOrdering>&       C)
    {
        typedef typename FullMatrix<T,SP2,COrdering>::value_type           value_type;
        typedef FullMatrix<value_type,ImmutableSharedData,FortranOrdering> FMat;

        const FMat Bt(B.numCols(), B.numRows(), B.data());

        matMulAdd_NT(a1, A, Bt, a2, C);
    }


    template<class charT, class traits,
             typename T, template<typename> class SP, class OP>
    std::basic_ostream<charT,traits>&
    operator<<(std::basic_ostream<charT,traits>& os,
               const FullMatrix<T,SP,OP>& A)
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
