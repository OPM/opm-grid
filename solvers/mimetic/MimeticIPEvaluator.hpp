//===========================================================================
//
// File: MimeticIPEvaluator.hpp
//
// Created: Thu Jun 25 13:09:23 2009
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

#ifndef OPENRS_MIMETICIPEVALUATOR_HEADER
#define OPENRS_MIMETICIPEVALUATOR_HEADER

#include <algorithm>   // for std::fill_n()
#include <vector>

#include <dune/solvers/mimetic/fortran.hpp>
#include <dune/solvers/mimetic/blas_lapack.hpp>
#include <dune/solvers/mimetic/FortranMatrix.hpp>

namespace Dune {
    template<class CellIter, bool computeInverseIP> class MimeticIPEvaluator;

    template<class CellIter>
    class MimeticIPEvaluator<CellIter,true> {
    public:
        typedef typename CellIter::Scalar Scalar;

        MimeticIPEvaluator(const int max_nf)
            : max_nf_(max_nf),
              work_  (64 * max_nf),  // 64 from ILAENV
              tau_   (3),            // Magic (only for <= 3D)
              fa_    (max_nf * max_nf),
              t1_    (max_nf * max_nf),
              t2_    (max_nf * max_nf)
        {}

        void evaluate(const CellIter& c, const std::vector<T>& perm,
                      FortranMatrix<Scalar>& Binv)
        {
            typedef typename CellIter::FaceIterator FI;
            typedef typename CellIter::Vector       CV;
            typedef typename FI      ::Vector       FV;

            const int nf = Binv.numRows();
            const CV  cc = c.centroid();
            const int nd = FV::size;    ASSERT(perm.size() == nd*nd);

            FortranMatrix<Scalar,false> T1(nf, nf, &t1_[0]);
            FortranMatrix<Scalar,false> T2(nf, nf, &t2_[0]);
            FortranMatrix<Scalar,false> FA(nf, nf, &fa_[0]);

            // Clear input matrix of any residual data.
            std::fill_n(Binv.data(), nf * nf, Scalar(0.0));

            int i = 0;
            for (FI f = c.facebegin(); f != c.faceend(); ++f, ++i) {
                Binv(i,i) = Scalar(1.0);
                FA(i,i)   = f.area();

                FV fc = f.centroid();  fc -= cc;  fc *= FA(i,i);
                FV fn = f.normal  ();             fn *= FA(i,i);

                for (int j = 0; j < cc.size(); ++j) {
                    T1(i,j) = fn[j];
                    T2(i,j) = fc[j];
                }
            }

            // T2 <- orth(T2)
            if (orthogonalizeColumns(T2, tau, work) != 0)
            {
                ASSERT (false);
            }

            // Binv <- Binv - T2*T2' == I - Q*Q'
            symmetricUpdate(-1, T2, 1, Binv);

            // Binv <- diag(A) * Binv * diag(A)
            symmetricUpdate(FA, Binv);

            // T2 <- N*K -- Assumes K (i.e., perm) is stored in C order
            // (i.e., transposed from the point of view of a FortranMatrix<T>).
            //
            FortranMatrix<Scalar,false> K(nd, nd, &perm[0]);
            MxM<Scalar,false,true>(Scalar(1.0), T1, K, Scalar(0.0), T2);

            // Binv <- (T2*N' + Binv) / vol(c)
            //      == (N*K*N' + t*(diag(A) * (I - Q*Q') * diag(A))) / vol(c)
            //
            // where t = 6/d * TRACE(K) (== 2*TRACE(K) for 3D).
            //
            Scalar t = 0.0;
            for (int j = 0; j < nd; ++j) t += K(j,j);

            MxM<Scalar,false,true>(Scalar(1.0)     /       c.volume() , T2, T1,
                                   Scalar(6.0) * t / (nd * c.volume()), Binv  );
        }

    private:
        int                 max_nf_      ;
        std::vector<Scalar> work_        ;
        std::vector<Scalar> tau_         ;
        std::vector<Scalar> fa_, t1_, t2_;
    };

} // namespace Dune

#endif // OPENRS_MIMETICIPEVALUATOR_HEADER
