//===========================================================================
//
// File: blas_lapack.hpp
//
// Created: Sun Jun 21 18:56:51 2009
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

#ifndef OPENRS_BLAS_LAPACK_HEADER
#define OPENRS_BLAS_LAPACK_HEADER

#include <dune/common/ErrorMacros.hpp>
#include <dune/solvers/common/fortran.hpp>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef DGEMV
#undef DGEMV
#endif
#define  DGEMV F77_NAME(dgemv,DGEMV)

    /// @brief
    ///    Double precision general matrix-vector vector update.
    ///    Specifically, @f$ y \leftarrow a_1 \mathrm{op}(A) x + a_2 y
    ///    @f$.
    ///
    /// @param trans
    ///    Specifcation of @f$ \mathrm{op} @f$.  @code trans = 'N'
    ///    @endcode yields @f$ \mathrm{op}(A) = A @f$ while @code
    ///    trans = 'T' @endcode yields @f$ \mathrm{op}(A) =
    ///    A^{\mathsf{T}} @f$.
    ///
    /// @param [in] m
    ///    Number of matrix rows (and number of rows for which new
    ///    data will be assigned to the result vector @f$ y @f$.)
    ///
    /// @param [in] n
    ///    Number of matrix columns (and number of rows for which data
    ///    will be retrieved from the input vector @f$ x @f$).
    ///
    /// @param [in] a1
    ///    Scalar coefficient @f$ a_1 @f$.
    ///
    /// @param [in] A
    ///    Pointer to first data element of matrix @f$ A @f$.
    ///
    /// @param [in] ldA
    ///    Leading dimension of storage array containing the matrix
    ///    @f$ A @f$.
    ///
    /// @param [in] x
    ///    Input vector @f$ x @f$.
    ///
    /// @param [in] incX
    ///    Data element stride for input vector @f$ x @f$.
    ///
    /// @param [in] a2
    ///    Scalar coefficient @f$ a_2 @f$.
    ///
    /// @param y
    ///    Result vector @f$ y @f$.
    ///
    /// @param [in] incY
    ///    Data element stride for result vector @f$ y @f$.
    void DGEMV(F77_CHARACTER_TYPE,
               const int*    m   , const int*    n,
               const double* a1  , const double* A, const int* ldA ,
                                   const double* x, const int* incX,
               const double* a2  ,       double* y, const int* incY);


#ifdef DGEMM
#undef DGEMM
#endif
#define  DGEMM F77_NAME(dgemm,DGEMM)

    // C <- a1*op(A)*op(B) + a2*C  where op(X) \in {X, X.', X'}
    /// @brief
    ///    Double precision general matrix-matrix product matrix
    ///    update.  Specifically, @f$ C \leftarrow a_1 \mathrm{op}(A)
    ///    \mathrm{op}(B) + a_2 C @f$.
    ///
    /// @param [in] transA
    ///    Specifcation of @f$ \mathrm{op} @f$.  @code transA = 'N'
    ///    @endcode yields @f$ \mathrm{op}(A) = A @f$ while @code
    ///    transA = 'T' @endcode yields @f$ \mathrm{op}(A) =
    ///    A^{\mathsf{T}} @f$.
    ///
    /// @param [in] transB
    ///    Specifcation of @f$ \mathrm{op} @f$.  @code transB = 'N'
    ///    @endcode yields @f$ \mathrm{op}(B) = A @f$ while @code
    ///    transB = 'T' @endcode yields @f$ \mathrm{op}(B) =
    ///    B^{\mathsf{T}} @f$.
    ///
    /// @param [in] m
    ///    Number of rows of matrix @f$ \mathrm{op}(A) @f$ and of
    ///    matrix @f$ C @f$.
    ///
    /// @param [in] n
    ///    Number of columns of matrix @f$ \mathrm{op}(B) @f$ and of
    ///    matrix @f$ C @f$.
    ///
    /// @param [in] k
    ///    Number of columns of matrix @f$ \mathrm{op}(A) @f$ and
    ///    number of rows of matrix @f$ \mathrm{op}(B) @f$.
    ///
    /// @param [in] a1
    ///    Scalar coefficient @f$ a_1 @f$.
    ///
    /// @param [in] A
    ///    Pointer to first data element of matrix @f$ A @f$.
    ///
    /// @param [in] ldA
    ///    Leading dimension of storage array containing the matrix
    ///    @f$ A @f$.
    ///
    /// @param [in] B
    ///    Pointer to first data element of matrix @f$ B @f$.
    ///
    /// @param [in] ldB
    ///    Leading dimension of storage array containing the matrix
    ///    @f$ B @f$.
    ///
    /// @param [in] a2
    ///    Scalar coefficient @f$ a_2 @f$.
    ///
    /// @param C
    ///    Pointer to first data element of matrix @f$ C @f$.
    ///
    /// @param [in] ldC
    ///    Leading dimension of storage array containing the matrix
    ///    @f$ C @f$.
    void DGEMM(F77_CHARACTER_TYPE, F77_CHARACTER_TYPE,
               const int*    m   , const int*    n   , const int* k  ,
               const double* a1  , const double* A   , const int* ldA,
                                   const double* B   , const int* ldB,
               const double* a2  ,       double* C   , const int* ldC);


#ifdef DSYRK
#undef DSYRK
#endif
#define  DSYRK F77_NAME(dsyrk,DSYRK)

    // C <- a1*A*A' + a2*C   *or*   C <- a1*A'*A + a2*C
    void DSYRK(F77_CHARACTER_TYPE, F77_CHARACTER_TYPE,
               const int*    n   , const int*    k   ,
               const double* a1  , const double* A   , const int* ldA,
               const double* a2  ,       double* C   , const int* ldC);


#ifdef DTRMM
#undef DTRMM
#endif
#define  DTRMM F77_NAME(dtrmm,DTRMM)

    // B <- a*op(A)*B  *or*  B <- a*B*op(A)  where op(X) \in {X, X.', X'}
    void DTRMM(F77_CHARACTER_TYPE, F77_CHARACTER_TYPE,
               F77_CHARACTER_TYPE, F77_CHARACTER_TYPE,
               const int*    m   , const int* n      ,
               const double* a   ,
               const double* A   , const int* ldA    ,
                     double* B   , const int* ldB);


#ifdef DGEQRF
#undef DGEQRF
#endif
#define  DGEQRF F77_NAME(dgeqrf,DGEQRF)

    void DGEQRF(const int*    m    , const int*    n   ,
                      double* A    , const int*    ld  ,
                      double* tau  ,       double* work,
                const int*    lwork,       int*    info);


#ifdef DORGQR
#undef DORGQR
#endif
#define  DORGQR F77_NAME(dorgqr,DORGQR)

    void DORGQR(const int*    m   , const int* n    , const int*    k  ,
                      double* A   , const int* ld   , const double* tau,
                      double* work, const int* lwork,       int*    info);

#ifdef DGETRF
#undef DGETRF
#endif
#define  DGETRF F77_NAME(dgetrf,DGETRF)

    void DGETRF(const int*    m   , const int* n ,
                      double* A   , const int* ld,
                      int*    ipiv,       int* info);

#ifdef DGETRI
#undef DGETRI
#endif
#define  DGETRI F77_NAME(dgetri,DGETRI)

    void DGETRI(const int*    n   ,
                      double* A   , const int* ld,
                const int*    ipiv,
                      double* work,       int* lwork, int* info);

#ifdef __cplusplus
}
#endif

namespace Dune {
    namespace BLAS_LAPACK {
        //--------------------------------------------------------------------------
        /// @brief GEneral Matrix Vector product (Level 2 BLAS).
        ///
        /// @tparam T Element type of matrix/vector.
       template<typename T>
        void GEMV(const char* transA,
                  const int   m     , const int   n,
                  const T&    a1    , const T*    A, const int ldA ,
                                      const T*    x, const int incX,
                  const T&    a2    ,       T*    y, const int incY);

        /// @brief GEneral Matrix Vector product specialization for double.
        template<>
        void GEMV<double>(const char*   transA,
                          const int     m     , const int     n,
                          const double& a1    , const double* A, const int ldA,
                                                const double* x, const int incX,
                          const double& a2    ,       double* y, const int incY)
        {
            ASSERT((transA[0] == 'N') || (transA[0] == 'T'));

            DGEMV(F77_CHARACTER(transA[0]),
                  &m, &n, &a1, A, &ldA, x, &incX, &a2, y, &incY);
        }


        //--------------------------------------------------------------------------
        /// @brief GEneral Matrix Matrix product (Level 3 BLAS).
        ///
        /// @tparam T Element type of matrix.
        template<typename T>
        void GEMM(const char* transA, const char* transB,
                  const int   m     , const int   n     , const int k  ,
                  const T&    a1    , const T*    A     , const int ldA,
                                      const T*    B     , const int ldB,
                  const T&    a2    ,       T*    C     , const int ldC);

        /// @brief GEneral Matrix Matrix product specialization for double.
        template<>
        void GEMM<double>(const char*   transA, const char*   transB,
                          const int     m     , const int     n     , const int k  ,
                          const double& a1    , const double* A     , const int ldA,
                                                const double* B     , const int ldB,
                          const double& a2    ,       double* C     , const int ldC)
        {
            ASSERT((transA[0] == 'N') || (transA[0] == 'T'));
            ASSERT((transB[0] == 'N') || (transB[0] == 'T'));

            DGEMM(F77_CHARACTER(transA[0]), F77_CHARACTER(transB[0]),
                  &m, &n, &k, &a1, A, &ldA, B, &ldB, &a2, C, &ldC);
        }


        //--------------------------------------------------------------------------
        /// @brief SYmmetric Rank K update of symmetric matrix (Level 3 BLAS)
        ///
        /// @tparam T Matrix element type.
        template<typename T>
        void SYRK(const char* uplo, const char* trans,
                  const int   n   , const int   k    ,
                  const T&    a1  , const T*    A    , const int ldA,
                  const T&    a2  ,       T*    C    , const int ldC);

        /// @brief SYmmetric Rank K update of symmetric matrix specialization for double.
        template<>
        void SYRK<double>(const char*   uplo, const char*   trans,
                          const int     n   , const int     k    ,
                          const double& a1  , const double* A    , const int ldA,
                          const double& a2  ,       double* C    , const int ldC)
        {
            ASSERT((uplo[0]  == 'U') || (uplo[0]  == 'L'));
            ASSERT((trans[0] == 'N') || (trans[0] == 'T'));

            DSYRK(F77_CHARACTER(uplo[0]), F77_CHARACTER(trans[0]),
                  &n, &k, &a1, A, &ldA, &a2, C, &ldC);
        }



        //--------------------------------------------------------------------------
        /// @brief TRiangular Matrix Matrix product (Level 2 BLAS)
        ///
        /// @tparam T Matrix element type.
        template<typename T>
        void TRMM(const char* side  , const char* uplo,
                  const char* transA, const char* diag,
                  const int   m     , const int   n   , const T& a,
                  const T*    A     , const int   ldA ,
                        T*    B     , const int   ldB);

        /// @brief TRiangular Matrix Matrix product specialization for double.
        template<>
        void TRMM<double>(const char*   side  , const char* uplo,
                          const char*   transA, const char* diag,
                          const int     m     , const int   n   , const double& a,
                          const double* A     , const int   ldA ,
                                double* B     , const int   ldB)
        {
            ASSERT((side[0]   == 'L') || (side[0]   == 'R'));
            ASSERT((uplo[0]   == 'U') || (uplo[0]   == 'L'));
            ASSERT((transA[0] == 'N') || (transA[0] == 'T'));
            ASSERT((diag[0]   == 'N') || (diag[0]   == 'U'));

            DTRMM(F77_CHARACTER(side[0])  , F77_CHARACTER(uplo[0]),
                  F77_CHARACTER(transA[0]), F77_CHARACTER(diag[0]),
                  &m, &n, &a, A, &ldA, B, &ldB);
        }


        //--------------------------------------------------------------------------
        /// @brief GEneral matrix QR Factorization (LAPACK)
        ///
        /// @tparam T Matrix element type.
        template<typename T>
        void GEQRF(const int m    , const int  n   ,
                         T*  A    , const int  ld  ,
                         T*  tau  ,       T*   work,
                   const int lwork,       int& info);

        /// @brief GEneral matrix QR Factorization specialization for double.
        template<>
        void GEQRF<double>(const int     m    , const int     n   ,
                                 double* A    , const int     ld  ,
                                 double* tau  ,       double* work,
                           const int     lwork,       int&    info)
        {
            DGEQRF(&m, &n, A, &ld, tau, work, &lwork, &info);
        }


        //--------------------------------------------------------------------------
        /// @brief ORthogonal matrix Generator from QR factorization (LAPACK).
        ///
        /// @tparam T Matrix element type.
        template<typename T>
        void ORGQR(const int m   , const int n    , const int  k  ,
                         T*  A   , const int ld   , const T*   tau,
                         T*  work, const int lwork,       int& info);

        /// @brief
        ///    ORthogonal matrix Generator from QR factorization
        ///    specialization for double.
        template<>
        void ORGQR<double>(const int     m   , const int n    , const int     k  ,
                                 double* A   , const int ld   , const double* tau,
                                 double* work, const int lwork,       int&    info)
        {
            DORGQR(&m, &n, &k, A, &ld, tau, work, &lwork, &info);
        }


        //--------------------------------------------------------------------------
        /// @brief GEneral matrix TRiangular Factorization (LAPACK).
        ///
        /// @tparam T Matrix element type.
        template<typename T>
        void GETRF(const int m, const int n, T* A,
                   const int ld, int* ipiv, int& info);

        /// @brief
        ///    GEneral matrix TRiangular Factorization specialization for double.
        template<>
        void GETRF<double>(const int m, const int n , double* A,
                           const int ld, int* ipiv, int& info)
        {
            DGETRF(&m, &n, A, &ld, ipiv, &info);
        }



        //--------------------------------------------------------------------------
        /// @brief GEneral matrix TRiangular Inversion (LAPACK).
        ///
        /// @tparam T Matrix element type.
        template<typename T>
        void GETRI(const int  n   , T* A   , const int ld,
                   const int* ipiv, T* work, int lwork, int& info);

        /// @brief GEneral matrix TRiangular Inversion specialization for double.
        template<>
        void GETRI(const int  n   , double* A   , const int ld,
                   const int* ipiv, double* work, int lwork, int& info)
        {
            DGETRI(&n, A, &ld, ipiv, work, &lwork, &info);
        }
    } // namespace BLAS_LAPACK
} // namespace Dune

#endif // OPENRS_BLAS_LAPACK_HEADER
