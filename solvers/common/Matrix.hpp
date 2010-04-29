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
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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
#include <dune/common/ErrorMacros.hpp>

#include <dune/solvers/common/fortran.hpp>
#include <dune/solvers/common/blas_lapack.hpp>

namespace Dune {

    // ----------------------------------------------------------------------
    // FullMatrix storage policies.
    //

    /// @brief
    ///    FullMatrix StoragePolicy which provides object owning
    ///    semantics.
    ///
    /// @tparam T
    ///    Element type of the FullMatrix.  Often @code T @endcode is
    ///    an alias for @code double @endcode.
    template<typename T>
    class OwnData {
    public:
        /// @brief Storage element access.
        ///
        /// @param [in] i
        ///    Linear element index.
        ///
        /// @return
        ///    Storage element at index @code i @endcode.
        T&       operator[](int i)       { return data_[i]; }
        const T& operator[](int i) const { return data_[i]; }


        /// @brief Data size query.
        ///
        /// @return Number of elements in storage array.
        int size() const { return data_.size(); }

        /// @brief Direct access to all data.
        ///
        /// @return Pointer to first element of storage array.
        T*       data()       { return &data_[0]; }
        const T* data() const { return &data_[0]; }

    protected:
        /// @brief Constructor.
        ///
        /// @param [in] sz
        ///    Number of elements in FullMatrix storage array.
        ///
        /// @param [in] data
        ///    Initial data vector.  If non-NULL, must contain @code
        ///    sz @endcode elements which will be assigned to a
        ///    freshly allocated storage array.  If NULL, a @code sz
        ///    @endcode element all-zero storage array will be
        ///    constructed.
        OwnData(int sz, const T* data)
        {
            if (data) {
                data_.assign(data, data + sz);
            } else {
                data_.resize(sz);
            }
        }

    private:
        std::vector<T> data_;
    };

    /// @brief
    ///    FullMatrix StoragePolicy which provides object sharing
    ///    semantics.
    ///
    /// @tparam T
    ///    Element type of the FullMatrix.  Often @code T @endcode is
    ///    an alias for @code double @endcode.
    template<typename T>
    class SharedData {
    public:
        /// @brief Storage element access.
        ///
        /// @param [in] i
        ///    Linear element index.
        ///
        /// @return
        ///    Storage element at index @code i @endcode.
        T&       operator[](int i)       { return data_[i]; }
        const T& operator[](int i) const { return data_[i]; }

        /// @brief Data size query.
        ///
        /// @return Number of elements in storage array.
        int size() const { return sz_; }

        /// @brief Direct access to all data.
        ///
        /// @return Pointer to first element of storage array.
        T*       data()       { return data_; }
        const T* data() const { return data_; }

    protected:
        /// @brief Constructor.
        ///
        /// @param [in] sz
        ///    Number of elements in FullMatrix storage array.
        ///
        /// @param [in] data
        ///    Initial data vector.  If non-NULL, must point to a @code
        ///    sz @endcode-element data vector.  If NULL, @code sz
        ///    @endcode must be zero as well.
        SharedData(int sz, T* data)
            : sz_(sz), data_(data)
        {
            ASSERT ((sz == 0) == (data == 0));
        }

    private:
        int sz_;
        T*  data_;
    };

    /// @brief
    ///    FullMatrix StoragePolicy which provides immutable object
    ///    sharing semantics.
    ///
    /// @tparam T
    ///    Element type of the FullMatrix.  Often @code T @endcode is
    ///    an alias for @code double @endcode.
    template<typename T>
    class ImmutableSharedData {
    public:
        /// @brief Storage element access.
        ///
        /// @param [in] i
        ///    Linear element index.
        ///
        /// @return
        ///    Storage element at index @code i @endcode.
        const T& operator[](int i) const { return data_[i]; }

        /// @brief Data size query.
        ///
        /// @return Number of elements in storage array.
        int size() const { return sz_; }

        /// @brief Direct access to all data.
        ///
        /// @return Pointer to first element of storage array.
        const T* data() const { return data_; }

    protected:
        /// @brief Constructor.
        ///
        /// @param [in] sz
        ///    Number of elements in FullMatrix storage array.
        ///
        /// @param [in] data
        ///    Initial data vector.  Must be non-NULL and point to a
        ///    @code sz @endcode-element data vector.
        ImmutableSharedData(int sz, const T* data)
            : sz_(sz), data_(data)
        {
            ASSERT (data_ != 0);
        }

    private:
        int sz_;
        const T*  data_;
    };





    // ----------------------------------------------------------------------
    // FullMatrix ordering policies.
    //

    /// @brief
    ///    FullMatrix OrderingPolicy base class.
    class OrderingBase {
    public:
        /// @brief
        ///    Retrieve the number of matrix rows.
        ///
        /// @return
        ///    Number of matrix rows.
        int numRows() const { return rows_; }

        /// @brief
        ///    Retrieve the number of matrix columns.
        ///
        /// @return
        ///    Number of matrix columns.
        int numCols() const { return cols_; }

    protected:
        /// @brief
        ///    Default constructor (yields 0-by-0 matrix).
        OrderingBase()
            : rows_(0), cols_(0)
        {}

        /// @brief
        ///    Constructor for matrix of non-zero size.
        ///
        /// @param [in] rows
        ///    Number of matrix rows.
        ///
        /// @param [in] cols
        ///    Number of matrix columns.
        OrderingBase(int rows, int cols)
            : rows_(rows), cols_(cols)
        {}

    private:
        int rows_, cols_;
    };


    /// @brief
    ///    FullMatrix C-ordering policy (column index cycling the most
    ///    rapidly).
    class COrdering : public OrderingBase {
    public:
        /// @brief
        ///    Retrieve the (BLAS/LAPACK) leading dimension of the
        ///    matrix storage array.
        ///
        /// @return
        ///    Leading dimension (i.e., the number of columns for a
        ///    C-ordered matrix) of the matrix storage array.
        int leadingDimension() const { return numCols(); }

        /// @brief
        ///    Retrieve the linear index (into storage array) of the
        ///    element at specific row/column position of the current
        ///    matrix.
        ///
        /// @param [in] row
        ///    Row position of the required element
        ///
        /// @param [in] col
        ///    Column position of the required element
        ///
        /// @return
        ///    Linear index of (row,col) pair.
        int idx(int row, int col) const
        {
            ASSERT ((0 <= row) && (row < numRows()));
            ASSERT ((0 <= col) && (col < numCols()));

            return row*numCols() + col;
        }

    protected:
        /// @brief Default constructor.
        COrdering()
            : OrderingBase()
        {}

        /// @brief
        ///    Constructor for C-ordered matrix of non-zero size.
        ///
        /// @param [in] rows
        ///    Number of matrix rows.
        ///
        /// @param [in] cols
        ///    Number of matrix columns.
        COrdering(int rows, int cols)
            : OrderingBase(rows, cols)
        {}
    };


    /// @brief
    ///    FullMatrix Fortran ordering policy (row index cycling the
    ///    most rapidly).
    class FortranOrdering : public OrderingBase {
    public:
        /// @brief
        ///    Retrieve the (BLAS/LAPACK) leading dimension of the
        ///    matrix storage array.
        ///
        /// @return
        ///    Leading dimension (i.e., the number of rows for a
        ///    Fortran ordered matrix) of the matrix storage array.
        int leadingDimension() const { return numRows(); }

        /// @brief
        ///    Retrieve the linear index (into storage array) of the
        ///    element at specific row/column position of the current
        ///    matrix.
        ///
        /// @param [in] row
        ///    Row position of the required element
        ///
        /// @param [in] col
        ///    Column position of the required element
        ///
        /// @return
        ///    Linear index of (row,col) pair.
        int idx(int row, int col) const
        {
            ASSERT ((0 <= row) && (row < numRows()));
            ASSERT ((0 <= col) && (col < numCols()));

            return row + col*numRows();
        }

    protected:
        /// @brief Default constructor.
        FortranOrdering()
            : OrderingBase()
        {}

        /// @brief
        ///    Constructor for Fortran ordered matrix of non-zero size.
        ///
        /// @param [in] rows
        ///    Number of matrix rows.
        ///
        /// @param [in] cols
        ///    Number of matrix columns.
        FortranOrdering(int rows, int cols)
            : OrderingBase(rows, cols)
        {}
    };




    // ----------------------------------------------------------------------
    // Class FullMatrix.
    //


    /// @brief
    ///    Dynamically sized m-by-n matrix with general element
    ///    storage (in a linear array) and element ordering.
    ///
    /// @tparam T
    ///    Element type.
    ///
    /// @tparam StoragePolicy
    ///    How to organize/store the matrix elements.  Parametrized on
    ///    the element type.
    ///
    /// @tparam OrderingPolicy
    ///    How to order the m-by-n matrix elements.  Typically
    ///    'COrdering' or 'FortranOrdering'.
    template<typename                 T,
             template<typename> class StoragePolicy,
             class                    OrderingPolicy>
    class FullMatrix : private StoragePolicy<T>,
                       private OrderingPolicy
    {
    public:
        /// @brief
        ///    Default constructor.
        FullMatrix()
            : StoragePolicy<T>(0, 0),
              OrderingPolicy()
        {}

        /// @brief
        ///    Constructor.
        ///
        /// @tparam DataPointer
        ///    Type representing a pointer to data with which to
        ///    initalize the matrix elements.  Assumed to support
        ///    ordinary dereferencing and pointer arithmetic.
        ///
        /// @param [in] rows
        ///    Number of matrix rows.
        ///
        /// @param [in] cols.
        ///    Number of matrix columns.
        ///
        /// @param [in] data.
        ///    Initial matrix data.  Interpretation of this data is
        ///    dependent upon the @code StoragePolicy @endcode.
        template <typename DataPointer>
        FullMatrix(int rows, int cols, DataPointer data)
            : StoragePolicy<T>(rows * cols, data),
              OrderingPolicy(rows, cols)
        {}

        /// @brief
        ///    Copy constructor.
        ///
        /// @tparam OtherSP
        ///    Storage policy of other matrix.
        ///
        /// @param [in] m
        ///    Constructor right hand side.
        template <template<typename> class OtherSP>
        explicit FullMatrix(const FullMatrix<T, OtherSP, OrderingPolicy>& m)
            : StoragePolicy<T>(m.numRows()*m.numCols(), m.data()),
              OrderingPolicy(m.numRows(), m.numCols())
        {
        }

        /// @brief
        ///    Assignment operator.
        ///
        /// @tparam OtherSP
        ///    Storage policy of other matrix.
        ///
        /// @tparam OtherOP
        ///    Ordering policy of other matrix.
        ///
        /// @param [in] m
        ///    Assignment right hand side.
        template <template<typename> class OtherSP, class OtherOP>
        FullMatrix& operator=(const FullMatrix<T, OtherSP, OtherOP>& m)
        {
            ASSERT(numRows() == m.numRows());
            ASSERT(numCols() == m.numCols());
            for (int r = 0; r < numRows(); ++r) {
                for (int c = 0; c < numCols(); ++c) {
                    this->operator()(r, c) = m(r,c);
                }
            }
            return *this;
        }

        /// @brief
        ///    Matrix addition operator (self-modifying).
        ///
        /// @tparam OtherSP
        ///    Storage policy of other matrix.
        ///
        /// @param [in] m
        ///    Matrix whose data will be added to self.
        template <template<typename> class OtherSP>
        void operator+= (const FullMatrix<T, OtherSP, OrderingPolicy>& m)
        {
            ASSERT(numRows() == m.numRows() && numCols() == m.numCols());
            std::transform(data(), data() + this->size(),
                           m.data(), data(), std::plus<T>());
        }

        /// @brief
        ///    Multiply self by scalar.
        ///
        /// @param [in] scalar.
        ///    Scalar by which to multiply own data.
        void operator*= (const T& scalar)
        {
            std::transform(data(), data() + this->size(),
                           data(), boost::bind(std::multiplies<T>(), _1, scalar));
        }

        /// @brief
        ///    Matrix value type (i.e., element type).
        typedef T value_type;

        using StoragePolicy<T>::data;
        using OrderingPolicy::numRows;
        using OrderingPolicy::numCols;
        using OrderingPolicy::leadingDimension;

        /// @brief
        ///    Matrix read/write element access.
        ///
        /// @param [in] row
        ///    Row index of requested matrix element.
        ///
        /// @param [in] col
        ///    Column index of requested matrix element.
        ///
        /// @return
        ///    Matrix element at position (row,col).
        value_type&       operator()(int row, int col)
        {
            return this->operator[](this->idx(row, col));
        }

        /// @brief
        ///    Matrix read-only element access.
        ///
        /// @param [in] row
        ///    Row index of requested matrix element.
        ///
        /// @param [in] col
        ///    Column index of requested matrix element.
        ///
        /// @return
        ///    Matrix element at position (row,col).
        const value_type& operator()(int row, int col) const
        {
            return this->operator[](this->idx(row, col));
        }
    };



    // ----------------------------------------------------------------------
    // FullMatrix operations.
    //


    // Convenience typedefs.

    /// @brief
    ///    Convenience typedefs for C-ordered @code FullMatrix
    ///    @endcode types with 'Owning', 'Shared' and 'Immutable
    ///    Shared' matrix element storage semantics.
    typedef FullMatrix<double, OwnData,             COrdering>        OwnCMatrix;
    typedef FullMatrix<double, SharedData,          COrdering>        SharedCMatrix;
    typedef const FullMatrix<double, ImmutableSharedData, COrdering>  ImmutableCMatrix;


    /// @brief
    ///    Convenience typedefs for Fortran-ordered @code FullMatrix
    ///    @endcode types with 'Owning', 'Shared' and 'Immutable
    ///    Shared' matrix element storage semantics.
    typedef FullMatrix<double, OwnData,             FortranOrdering>       OwnFortranMatrix;
    typedef FullMatrix<double, SharedData,          FortranOrdering>       SharedFortranMatrix;
    typedef const FullMatrix<double, ImmutableSharedData, FortranOrdering> ImmutableFortranMatrix;


    /// @brief
    ///    Zero-fill a @code FullMatrix @endcode.
    ///
    /// @tparam Matrix
    ///    Matrix type.
    ///
    /// @param A
    ///    Specific matrix which will be zero-filled upon return.
    template<class Matrix>
    void zero(Matrix& A)
    {
        std::fill_n(A.data(), A.numRows() * A.numCols(),
                    typename Matrix::value_type(0.0));
    }

    /// @brief
    ///    Make an identity @code FullMatrix @endcode.
    ///
    /// @tparam Matrix
    ///    Matrix type.
    ///
    /// @param A
    ///    Specific matrix which will be zero-filled upon return.
    template<class Matrix>
    void eye(Matrix& A)
    {
        zero(A);
        for (int i = 0; i < std::min(A.numRows(), A.numCols()); ++i) {
            A(i, i) = 1.0;
        }
    }

    /// @brief
    ///    Compute matrix trace (i.e., sum(diag(A))).
    ///
    /// @tparam Matrix
    ///    Matrix type.
    ///
    /// @param [in] A
    ///    Matrix for which to compute the trace.
    ///
    /// @return
    ///    Trace of matrix.
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


    /// @brief
    ///    Matrix applied to a vector.
    ///
    /// @tparam Matrix
    ///    Matrix type.
    ///
    /// @tparam rows.
    ///    Number of matrix rows.
    ///
    /// @param [in] A
    ///    Matrix.
    ///
    /// @param [in] x
    ///    Vector
    ///
    /// @return
    ///    @f$ Ax @f$.
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

    /// @brief Compute C = AB. C must not overlap with A or B.
    /// @tparam Matrix1 a matrix type.
    /// @tparam Matrix2 a matrix type.
    /// @tparam MutableMatrix a matrix type with write access.
    /// @param[in] A left matrix of product.
    /// @param[in] B right matrix of product.
    /// @param[out] C resulting product matrix, it must already have the right size.
    template<class Matrix1, class Matrix2, class MutableMatrix>
    void prod(const Matrix1& A, const Matrix2& B, MutableMatrix& C)
    {
        int result_rows = A.numRows();
        int result_cols = B.numCols();
        int inner_dim = A.numCols();
        ASSERT (inner_dim == B.numRows());
        ASSERT(C.numRows() == result_rows);
        ASSERT(C.numCols() == result_cols);

        for (int c = 0; c < result_cols; ++c) {
            for (int r = 0; r < result_rows; ++r) {
                C(r,c) = 0.0;
                for (int i = 0; i < inner_dim; ++i) {
                    C(r,c) += A(r,i)*B(i,c);
                }
            }
        }
    }


    /// @brief
    ///    Construct orthonormal basis for matrix range (i.e., column
    ///    space).  Based on a QR factorization of the matrix.
    ///
    /// @tparam T
    ///    Matrix element type.
    ///
    /// @tparam StoragePolicy
    ///    Matrix storage policy.
    ///
    /// @param A
    ///    Matrix.  Will be overwritten by an orthogonal matrix, @f$ Q
    ///    @f$ whose columns represent an orthonormal basis for
    ///    range(A).
    ///
    /// @return
    ///    Zero for success, non-zero if an error occurred.
    template<typename T, template<typename> class StoragePolicy>
    int orthogonalizeColumns(FullMatrix<T,StoragePolicy,FortranOrdering>& A)
    {
        typedef typename FullMatrix<T,StoragePolicy,FortranOrdering>::value_type value_type;

        static std::vector<value_type> tau;
        static std::vector<value_type> work;

        if (int(tau .size()) <      A.numCols()) tau .resize(     A.numCols());
        if (int(work.size()) < 64 * A.numRows()) work.resize(64 * A.numRows());  // 64 from ILAENV

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


    /// @brief Matrix inversion, @f$A \leftarrow A^{-1} @f$.
    ///
    /// @tparam T
    ///    Matrix element type.
    ///
    /// @tparam StoragePolicy
    ///    Matrix storage policy.
    ///
    /// @tparam OrderingPolicy
    ///    Matrix ordering policy.
    ///
    /// @param A
    ///    Matrix.  Contains the inverse upon return from @code
    ///    invert() @endcode.
    ///
    /// @returns
    ///    The union of the LAPACK xGETRF and xGETRI 'INFO' subroutine
    ///    return flags.
    template<typename T, template<typename> class StoragePolicy, class OrderingPolicy>
    int invert(FullMatrix<T,StoragePolicy,OrderingPolicy>& A)
    {
        typedef typename FullMatrix<T,StoragePolicy,OrderingPolicy>::value_type value_type;

        ASSERT (A.numRows() == A.numCols());

        std::vector<int> ipiv(A.numRows());
        int info = 0;

        // Correct both for COrdering and FortranOrdering (inv(A)' == inv(A')).
        Dune::BLAS_LAPACK::GETRF(A.numRows(), A.numCols(), A.data(),
                                 A.leadingDimension(), &ipiv[0], info);

        if (info == 0) {
            std::vector<value_type> work(A.numRows());

            Dune::BLAS_LAPACK::GETRI(A.numRows(), A.data(), A.leadingDimension(),
                                     &ipiv[0], &work[0], int(work.size()), info);
        }

        return info;
    }


    /// @brief
    ///    Symmetric, rank @f$ k @f$ update of symmetric matrix.
    ///    Specifically, @f$ C \leftarrow a_1 AA^{\mathsf{T}} + a_2 C
    ///    @f$.
    ///
    /// @tparam T
    ///    Matrix element type.  Assumed to be an arithmetic type.
    ///    Typically @code T @endcode is an alias for @code double
    ///    @endcode.
    ///
    /// @tparam StoragePolicy.
    ///    Matrix storage policy.
    ///
    /// @param [in] a1
    ///    Scalar coefficient @f$ a_1 @f$.
    ///
    /// @param [in] A
    ///    Matrix @f$ A @f$.
    ///
    /// @param [in] a2
    ///    Scalar coefficient @f$ a_2 @f$.
    ///
    /// @param C
    ///    Matrix @f$ C @f$.
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
    /// @brief
    /// @todo Doc me!
    /// @tparam
    /// @param
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


    /// @brief
    ///    GEneral Matrix-Vector product (GAXPY operation).
    ///    Specifically, @f$ y \leftarrow a_1 Ax + a_2 y @f$.
    ///
    /// @tparam T
    ///    Matrix (and vector) element type.  Assumed to be an
    ///    arithmetic type and, typically, @code T @endcode is an
    ///    alias for @code double @endcode.
    ///
    /// @tparam SP
    ///    Matrix storage policy.
    ///
    /// @param [in] a1
    ///    Scalar coefficient @f$ a_1 @f$.
    ///
    /// @param [in] A
    ///    Matrix @f$ A @f$.
    ///
    /// @param [in] x
    ///    Vector @f$ x @f$.
    ///
    /// @param [in] a2
    ///    Scalar coefficient @f$ a_2 @f$.
    ///
    /// @param y
    ///    Vector @f$ y @f$.
    template<typename                 T ,
             template<typename> class SP>
    void vecMulAdd_N(const T&                                a1,
                     const FullMatrix<T,SP,FortranOrdering>& A ,
                     const std::vector<T>&                   x ,
                     const T&                                a2,
                     std::vector<T>&                         y)
    {
        ASSERT(A.numRows() == y.size());
        ASSERT(A.numCols() == x.size());

        Dune::BLAS_LAPACK::GEMV("No Transpose",
                                A.numRows(), A.numCols(),
                                a1, A.data(), A.leadingDimension(),
                                &x[0], 1, a2, &y[0], 1);
    }


    /// @brief
    ///    GEneral Matrix-Vector product (GAXPY operation).
    ///    Specifically, @f$ y \leftarrow a_1 Ax + a_2 y @f$.
    ///
    /// @tparam T
    ///    Matrix (and vector) element type.  Assumed to be an
    ///    arithmetic type and, typically, @code T @endcode is an
    ///    alias for @code double @endcode.
    ///
    /// @tparam SP
    ///    Matrix storage policy.
    ///
    /// @param [in] a1
    ///    Scalar coefficient @f$ a_1 @f$.
    ///
    /// @param [in] A
    ///    Matrix @f$ A @f$.
    ///
    /// @param [in] x
    ///    Vector @f$ x @f$.
    ///
    /// @param [in] a2
    ///    Scalar coefficient @f$ a_2 @f$.
    ///
    /// @param y
    ///    Vector @f$ y @f$.
    template<typename                 T ,
             template<typename> class SP>
    void vecMulAdd_N(const T&                                a1,
                     const FullMatrix<T,SP,FortranOrdering>& A ,
                     const T*                                x ,
                     const T&                                a2,
                     T*                                      y)
    {
        Dune::BLAS_LAPACK::GEMV("No Transpose",
                                A.numRows(), A.numCols(),
                                a1, A.data(), A.leadingDimension(),
                                x, 1, a2, y, 1);
    }


    /// @brief
    ///    GEneral Matrix-Vector product (GAXPY operation).
    ///    Specifically, @f$ y \leftarrow a_1 A^{\mathsf{T}}x + a_2 y
    ///    @f$.
    ///
    /// @tparam T
    ///    Matrix (and vector) element type.  Assumed to be an
    ///    arithmetic type and, typically, @code T @endcode is an
    ///    alias for @code double @endcode.
    ///
    /// @tparam SP
    ///    Matrix storage policy.
    ///
    /// @param [in] a1
    ///    Scalar coefficient @f$ a_1 @f$.
    ///
    /// @param [in] A
    ///    Matrix @f$ A @f$.
    ///
    /// @param [in] x
    ///    Vector @f$ x @f$.
    ///
    /// @param [in] a2
    ///    Scalar coefficient @f$ a_2 @f$.
    ///
    /// @param y
    ///    Vector @f$ y @f$.
    template<typename                 T ,
             template<typename> class SP>
    void vecMulAdd_T(const T&                                a1,
                     const FullMatrix<T,SP,FortranOrdering>& A ,
                     const std::vector<T>&                   x ,
                     const T&                                a2,
                     std::vector<T>&                         y)
    {
        ASSERT (A.numCols() == y.size());
        ASSERT (A.numRows() == x.size());

        Dune::BLAS_LAPACK::GEMV("Transpose",
                                A.numRows(), A.numCols(),
                                a1, A.data(), A.leadingDimension(),
                                &x[0], 1, a2, &y[0], 1);
    }


    /// @brief
    ///    GEneral Matrix-Vector product (GAXPY operation).
    ///    Specifically, @f$ y \leftarrow a_1 A^{\mathsf{T}}x + a_2 y
    ///    @f$.
    ///
    /// @tparam T
    ///    Matrix (and vector) element type.  Assumed to be an
    ///    arithmetic type and, typically, @code T @endcode is an
    ///    alias for @code double @endcode.
    ///
    /// @tparam SP
    ///    Matrix storage policy.
    ///
    /// @param [in] a1
    ///    Scalar coefficient @f$ a_1 @f$.
    ///
    /// @param [in] A
    ///    Matrix @f$ A @f$.
    ///
    /// @param [in] x
    ///    Vector @f$ x @f$.
    ///
    /// @param [in] a2
    ///    Scalar coefficient @f$ a_2 @f$.
    ///
    /// @param y
    ///    Vector @f$ y @f$.
    template<typename                 T ,
             template<typename> class SP>
    void vecMulAdd_T(const T&                                a1,
                     const FullMatrix<T,SP,FortranOrdering>& A ,
                     const T*                                x ,
                     const T&                                a2,
                     T*                                      y)
    {
        Dune::BLAS_LAPACK::GEMV("Transpose",
                                A.numRows(), A.numCols(),
                                a1, A.data(), A.leadingDimension(),
                                x, 1, a2, y, 1);
    }


    /// @brief
    ///    GEneral Matrix-Vector product (GAXPY operation).
    ///    Specifically, @f$ y \leftarrow a_1 Ax + a_2 y @f$.
    ///    Overload for C-ordered FullMatrix type.
    ///
    /// @tparam T
    ///    Matrix (and vector) element type.  Assumed to be an
    ///    arithmetic type and, typically, @code T @endcode is an
    ///    alias for @code double @endcode.
    ///
    /// @tparam SP
    ///    Matrix storage policy.
    ///
    /// @param [in] a1
    ///    Scalar coefficient @f$ a_1 @f$.
    ///
    /// @param [in] A
    ///    Matrix @f$ A @f$.
    ///
    /// @param [in] x
    ///    Vector @f$ x @f$.
    ///
    /// @param [in] a2
    ///    Scalar coefficient @f$ a_2 @f$.
    ///
    /// @param y
    ///    Vector @f$ y @f$.
    template<typename                 T ,
             template<typename> class SP>
    void vecMulAdd_N(const T&                          a1,
                     const FullMatrix<T,SP,COrdering>& A ,
                     const T*                          x ,
                     const T&                          a2,
                     T*                                y)
    {
        typedef FullMatrix<T, ImmutableSharedData, FortranOrdering> FMAT;

        const FMAT At(A.numCols(), A.numRows(), A.data());

        vecMulAdd_T(a1, At, x, a2, y);
    }


    /// @brief
    ///    GEneral Matrix-Matrix product update of other matrix.
    ///    Specificlly, @f$ C \leftarrow a_1AB + a_2C @f$.
    ///
    /// @tparam T
    ///    Matrix element type.  Assumed to be an arithmetic type and,
    ///    typically, @code T @endocde is an alias for @code double
    ///    @endcode.
    ///
    /// @tparam SP1
    ///    Storage policy of matrix @f$ A @f$.
    ///
    /// @tparam SP2
    ///    Storage policy of matrix @f$ B @f$.
    ///
    /// @tparam SP3
    ///    Storage policy of matrix @f$ C @f$.
    ///
    /// @param [in] a1
    ///    Scalar coefficient @f$ a_1 @f$.
    ///
    /// @param [in] A
    ///    Matrix @f$ A @f$.
    ///
    /// @param [in] B
    ///    Matrix @f$ B @f$.
    ///
    /// @param [in] a2
    ///    Scalar coefficient @f$ a_2 @f$.
    ///
    /// @param C
    ///    Matrix @f$ C @f$.
    template<typename                 T  ,
             template<typename> class SP1,
             template<typename> class SP2,
             template<typename> class SP3>
    void matMulAdd_NN(const T&                                 a1,
                      const FullMatrix<T,SP1,FortranOrdering>& A ,
                      const FullMatrix<T,SP2,FortranOrdering>& B ,
                      const T&                                 a2,
                      FullMatrix<T,SP3,FortranOrdering>&       C)
    {
        ASSERT(A.numRows() == C.numRows());
        ASSERT(A.numCols() == B.numRows());
        ASSERT(B.numCols() == C.numCols());

        int m = A.numRows();  // Number of *rows* in A
        int n = B.numCols();  // Number of *cols* in B
        int k = A.numCols();  // Number of *cols* in A (== numer of *rows* in B)

        Dune::BLAS_LAPACK::GEMM("No Transpose", "No Transpose", m, n, k,
                                a1, A.data(), A.leadingDimension(),
                                    B.data(), B.leadingDimension(),
                                a2, C.data(), C.leadingDimension());
    }


    /// @brief
    ///    GEneral Matrix-Matrix product update of other matrix.
    ///    Specificlly, @f$ C \leftarrow a_1AB^{\mathsf{T}} + a_2C
    ///    @f$.
    ///
    /// @tparam T
    ///    Matrix element type.  Assumed to be an arithmetic type and,
    ///    typically, @code T @endocde is an alias for @code double
    ///    @endcode.
    ///
    /// @tparam SP1
    ///    Storage policy of matrix @f$ A @f$.
    ///
    /// @tparam SP2
    ///    Storage policy of matrix @f$ B @f$.
    ///
    /// @tparam SP3
    ///    Storage policy of matrix @f$ C @f$.
    ///
    /// @param [in] a1
    ///    Scalar coefficient @f$ a_1 @f$.
    ///
    /// @param [in] A
    ///    Matrix @f$ A @f$.
    ///
    /// @param [in] B
    ///    Matrix @f$ B @f$.
    ///
    /// @param [in] a2
    ///    Scalar coefficient @f$ a_2 @f$.
    ///
    /// @param C
    ///    Matrix @f$ C @f$.
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

        int m = A.numRows();  // Number of *rows* in A
        int n = B.numRows();  // Number of *cols* in B'
        int k = A.numCols();  // Number of *cols* in A (== numer of *rows* in B')

        Dune::BLAS_LAPACK::GEMM("No Transpose", "Transpose", m, n, k,
                                a1, A.data(), A.leadingDimension(),
                                    B.data(), B.leadingDimension(),
                                a2, C.data(), C.leadingDimension());
    }


    /// @brief
    ///    GEneral Matrix-Matrix product update of other matrix.
    ///    Specificlly, @f$ C \leftarrow a_1A^{\mathsf{T}}B + a_2C
    ///    @f$.
    ///
    /// @tparam T
    ///    Matrix element type.  Assumed to be an arithmetic type and,
    ///    typically, @code T @endocde is an alias for @code double
    ///    @endcode.
    ///
    /// @tparam SP1
    ///    Storage policy of matrix @f$ A @f$.
    ///
    /// @tparam SP2
    ///    Storage policy of matrix @f$ B @f$.
    ///
    /// @tparam SP3
    ///    Storage policy of matrix @f$ C @f$.
    ///
    /// @param [in] a1
    ///    Scalar coefficient @f$ a_1 @f$.
    ///
    /// @param [in] A
    ///    Matrix @f$ A @f$.
    ///
    /// @param [in] B
    ///    Matrix @f$ B @f$.
    ///
    /// @param [in] a2
    ///    Scalar coefficient @f$ a_2 @f$.
    ///
    /// @param C
    ///    Matrix @f$ C @f$.
    template<typename                 T  ,
             template<typename> class SP1,
             template<typename> class SP2,
             template<typename> class SP3>
    void matMulAdd_TN(const T&                                 a1,
                      const FullMatrix<T,SP1,FortranOrdering>& A ,
                      const FullMatrix<T,SP2,FortranOrdering>& B ,
                      const T&                                 a2,
                      FullMatrix<T,SP3,FortranOrdering>&       C)
    {
        ASSERT (A.numCols() == C.numRows());
        ASSERT (A.numRows() == B.numRows());
        ASSERT (B.numCols() == C.numCols());

        int m = A.numCols();  // Number of *rows* in A'
        int n = B.numCols();  // Number of *cols* in B
        int k = A.numRows();  // Number of *cols* in A' (== numer of *rows* in B)

        Dune::BLAS_LAPACK::GEMM("Transpose", "No Transpose", m, n, k,
                                a1, A.data(), A.leadingDimension(),
                                    B.data(), B.leadingDimension(),
                                a2, C.data(), C.leadingDimension());
    }


    /// @brief
    ///    GEneral Matrix-Matrix product update of other matrix.
    ///    Specificlly, @f$ C \leftarrow a_1AB + a_2C @f$.  Overload
    ///    for C-ordered matrix @f$ B @f$.
    ///
    /// @tparam T
    ///    Matrix element type.  Assumed to be an arithmetic type and,
    ///    typically, @code T @endocde is an alias for @code double
    ///    @endcode.
    ///
    /// @tparam SP1
    ///    Storage policy of matrix @f$ A @f$.
    ///
    /// @tparam SP2
    ///    Storage policy of matrix @f$ B @f$.
    ///
    /// @tparam SP3
    ///    Storage policy of matrix @f$ C @f$.
    ///
    /// @param [in] a1
    ///    Scalar coefficient @f$ a_1 @f$.
    ///
    /// @param [in] A
    ///    Matrix @f$ A @f$.
    ///
    /// @param [in] B
    ///    Matrix @f$ B @f$.
    ///
    /// @param [in] a2
    ///    Scalar coefficient @f$ a_2 @f$.
    ///
    /// @param C
    ///    Matrix @f$ C @f$.
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


    /// @brief
    ///    Stream output operator for @code FullMatrix @endcode types.
    ///
    /// @tparam charT
    ///    Output stream character type.
    ///
    /// @tparam traits
    ///    Output stream character traits.
    ///
    /// @tparam T
    ///    Matrix element type.
    ///
    /// @tparam SP
    ///    Matrix storage policy.
    ///
    /// @tparam OP
    ///    Matrix ordering policy.
    ///
    /// @param os
    ///    Output stream.
    ///
    /// @param [in] A
    ///    Matrix.
    ///
    /// @return
    ///    Output stream (for output chaining).
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
