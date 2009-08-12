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

#include <vector>

#include <dune/grid/common/ErrorMacros.hpp>

#include <dune/solvers/mimetic/fortran.hpp>
#include <dune/solvers/mimetic/blas_lapack.hpp>
#include <dune/solvers/mimetic/Matrix.hpp>

namespace Dune {
    template<class CellIter, int dim, bool computeInverseIP> class MimeticIPEvaluator;

    template<class CellIter, int dim>
    class MimeticIPEvaluator<CellIter,dim,true> {
    public:
        typedef typename CellIter::Scalar Scalar;

        MimeticIPEvaluator(const int max_nf)
            : max_nf_(max_nf),
              fa_    (max_nf * max_nf),
              t1_    (max_nf * dim),
              t2_    (max_nf * dim)
        {}


        template<class PermTensor, template<typename> class SP>
        void evaluate(const CellIter&                        c,
                      const PermTensor&                      K,
                      FullMatrix<Scalar,SP,FortranOrdering>& Binv)
        {
            typedef typename CellIter::FaceIterator FI;
            typedef typename CellIter::Vector       CV;
            typedef typename FI      ::Vector       FV;

            const int nf = Binv.numRows();

            ASSERT(Binv.numRows()  <= max_nf_);
            ASSERT(Binv.numRows()  == Binv.numCols());
            ASSERT(FV::size        == dim);
            ASSERT(int(t1_.size()) >= nf * dim);
            ASSERT(int(t2_.size()) >= nf * dim);
            ASSERT(int(fa_.size()) >= nf * nf);

            SharedFortranMatrix T1(nf, dim, &t1_[0]);
            SharedFortranMatrix T2(nf, dim, &t2_[0]);
            SharedFortranMatrix fa(nf, nf , &fa_[0]);

            // Clear matrices of any residual data.
            zero(Binv);  zero(T1);  zero(T2);  zero(fa);

            // Setup: Binv <- I, T1 <- N, T2 <- C
            const CV cc = c->centroid();
            int i = 0;
            for (FI f = c->facebegin(); f != c->faceend(); ++f, ++i) {
                Binv(i,i) = Scalar(1.0);
                fa(i,i)   = f->area();

                FV fc = f->centroid();  fc -= cc;  fc *= fa(i,i);
                FV fn = f->normal  ();             fn *= fa(i,i);

                for (int j = 0; j < dim; ++j) {
                    T1(i,j) = fn[j];
                    T2(i,j) = fc[j];
                }
            }
            ASSERT(i == nf);

            // T2 <- orth(T2)
            if (orthogonalizeColumns(T2) != 0) {
                ASSERT (false);
            }

            // Binv <- Binv - T2*T2' == I - Q*Q'
            symmetricUpdate(Scalar(-1.0), T2, Scalar(1.0), Binv);

            // Binv <- diag(A) * Binv * diag(A)
            symmetricUpdate(fa, Binv);

            // T2 <- N*K
            matMulAdd_NN(Scalar(1.0), T1, K, Scalar(0.0), T2);

            // Binv <- (T2*N' + Binv) / vol(c)
            //      == (N*K*N' + t*(diag(A) * (I - Q*Q') * diag(A))) / vol(c)
            //
            // where t = 6/d * TRACE(K) (== 2*TRACE(K) for 3D).
            //
            Scalar t = 0.0; //trace(K);
            for (int i = 0; i < dim; ++i) t += K(i,i);
            matMulAdd_NT(Scalar(1.0)     /        c->volume() , T2, T1,
                         Scalar(6.0) * t / (dim * c->volume()), Binv  );
        }

        template<class Point, class Vector>
        void gravityTerm(const CellIter& c,
                         const Point&    grav,
                         const Scalar    omega,
                         Vector&         gterm) const
        {
            typedef typename CellIter::FaceIterator FI;
            typedef FieldVector<double,dim> Point;

            ASSERT (gterm.size() <= max_nf_);

            const Point cc = c->centroid();
            int i = 0;
            for (FI f = c->facebegin(); f != c->faceend(); ++f, ++i) {
                Point fc = f->centroid();
                fc -= cc;
                gterm[i] = omega * (fc * grav);
            }
        }

    private:
        int                 max_nf_      ;
        std::vector<Scalar> fa_, t1_, t2_;
    };

} // namespace Dune

#endif // OPENRS_MIMETICIPEVALUATOR_HEADER
