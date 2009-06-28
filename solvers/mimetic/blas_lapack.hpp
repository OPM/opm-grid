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

#include <dune/grid/common/ErrorMacros.hpp>
#include <dune/solvers/mimetic/fortran.hpp>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef DGEMM
#undef DGEMM
#endif
#define  DGEMM F77_NAME(dgemm,DGEMM)

    // C <- a1*op(A)*op(B) + a2*C  where op(X) \in {X, X.', X'}
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

#ifdef __cplusplus
}
#endif

namespace Dune {
    namespace BLAS_LAPACK {
        //--------------------------------------------------------------------------
        template<typename T>
        void GEMM(const char* transA, const char* transB,
                  const int   m     , const int   n     , const int k  ,
                  const T&    a1    , const T*    A     , const int ldA,
                                      const T*    B     , const int ldB,
                  const T&    a2    ,       T*    C     , const int ldC);

        template<>
        void GEMM<double>(const char*   transA, const char*   transB,
                          const int     m     , const int     n     , const int k  ,
                          const double& a1    , const T*      A     , const int ldA,
                                                const double* B     , const int ldB,
                          const double& a2    ,       double* C     , const int ldC)
        {
            ASSERT((transA[0] == 'N') || (transA[0] == 'T'));
            ASSERT((transB[0] == 'N') || (transB[0] == 'T'));

            DGEMM(F77_CHARACTER(transA[0]), F77_CHARACTER(transB[0]),
                  &m, &n, &k, &a1, A, &ldA, B, &ldB, &a2, C, &ldC);
        }


        //--------------------------------------------------------------------------
        template<typename T>
        void SYRK(const char* uplo, const char* trans,
                  const int   n   , const int   k    ,
                  const T&    a1  , const T*    A    , const int ldA,
                  const T&    a2  ,       T*    C    , const int ldC);

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
        template<typename T>
        void TRMM(const char* side  , const char* uplo,
                  const char* transA, const char* diag,
                  const int   m     , const int   n   , const T& a,
                  const T*    A     , const int   ldA ,
                        T*    B     , const int   ldB);

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
        template<typename T>
        void GEQRF(const int m    , const int  n   ,
                         T*  A    , const int  ld  ,
                         T*  tau  ,       T*   work,
                   const int lwork,       int& info);

        template<>
        void GEQRF<double>(const int     m    , const int     n   ,
                                 double* A    , const int     ld  ,
                                 double* tau  ,       double* work,
                           const int     lwork,       int&    info)
        {
            DGEQRF(&m, &n, A, &ld, tau, work, &lwork, &info);
        }


        //--------------------------------------------------------------------------
        template<typename T>
        void ORGQR(const int m   , const int n    , const int  k  ,
                         T*  A   , const int ld   , const T*   tau,
                         T*  work, const int lwork,       int& info);

        template<>
        void ORGQR<double>(const int     m   , const int n    , const int     k  ,
                                 double* A   , const int ld   , const double* tau,
                                 double* work, const int lwork,       int&    info)
        {
            DORGQR(&m, &n, &k, A, &ld, tau, work, &lwork, &info);
        }
    } // namespace BLAS_LAPACK
} // namespace Dune

#endif // OPENRS_BLAS_LAPACK_HEADER
