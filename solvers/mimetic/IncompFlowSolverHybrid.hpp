//===========================================================================
//
// File: IncompFlowSolverHybrid.hpp
//
// Created: Tue Jun 30 10:25:40 2009
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

#ifndef OPENRS_INCOMPFLOWSOLVERHYBRID_HEADER
#define OPENRS_INCOMPFLOWSOLVERHYBRID_HEADER

#include <algorithm>
#include <functional>
#include <map>
#include <numeric>
#include <ostream>
#include <utility>
#include <vector>

#include <boost/bind.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/operators.hh>

#include <dune/grid/common/ErrorMacros.hpp>
#include <dune/grid/common/SparseTable.hpp>

namespace Dune {
    namespace {
        template<class GI>
        bool topologyIsSane(const GI& g)
        {
            typedef typename GI::CellIterator CI;
            typedef typename CI::FaceIterator FI;

            bool sane = g.numberOfCells() >= 0;

            for (CI c = g.cellbegin(); sane && c != g.cellend(); ++c) {
                std::vector<int> n;

                for (FI f = c->facebegin(); f != c->faceend(); ++f) {
                    if (!f->boundary()) {
                        n.push_back(f->neighbourCellIndex());
                    }
                }
                std::sort(n.begin(), n.end());

                sane = std::unique(n.begin(), n.end()) == n.end();
            }

            return sane;
        }


        template<typename T>
        class axpby : public std::binary_function<T,T,T> {
        public:
            axpby(const T& a, const T& b) : a_(a), b_(b) {}
            T operator()(const T& x, const T& y)
            {
                return a_*x + b_*y;
            }
        private:
            T a_, b_;
        };
    }


    template<class GridInterface, class InnerProduct>
    class IncompFlowSolverHybrid {
        typedef typename GridInterface::Scalar Scalar;
    public:
        void init(const GridInterface& g)
        {
            ASSERT (topologyIsSane(g));

            max_ncf_                = -1;
            num_internal_faces_     =  0;
            total_num_faces_        =  0;
            matrix_structure_valid_ = false;

            std::vector<int>(g.numberOfCells(), -1).swap(cellno_);
            cellFaces_.clear();
            F_        .clear();
            std::vector<Scalar>(g.numberOfCells()).swap(L_);
            Binv_     .clear();

            if (g.numberOfCells() > 0) {
                buildGridTopology(g);
                buildSystemStructure();
            }
        }


        template<class ReservoirInterface>
        void assembleStatic(const GridInterface&      g,
                            const ReservoirInterface& r)
        {
            ASSERT(matrix_structure_valid_);

            typedef typename GridInterface     ::CellIterator CI;
            typedef typename ReservoirInterface::PermTensor   PermTensor;

            InnerProduct ip(max_ncf_);
            int i = 0;
            for (CI c = g.cellbegin(); c != g.cellend(); ++c, ++i) {
                const int nf = cellFaces_[i].size();

                SharedFortranMatrix Binv(nf, nf, &Binv_[i][0]);

                ip.evaluate(c, r.permeability(c->index()), Binv);
            }
        }


        template<class ReservoirInterface, class BCInterface>
        void solve(const GridInterface&       g  ,
                   const ReservoirInterface&  r  ,
                   const BCInterface&         bc ,
                   const std::vector<double>& src,
                   const std::vector<double>& sat)
        {
            assembleDynamic(g, r, bc, src, sat);
#if 0
            solveLinSys();
            computePressureAndFluxes();
#endif
        }


        template<typename charT, class traits>
        void printStats(std::basic_ostream<charT,traits>& os)
        {
            os << "IncompFlowSolverHybrid<>:\n"
               << "\tMaximum number of cell faces = " << max_ncf_ << '\n'
               << "\tNumber of internal faces     = " << num_internal_faces_ << '\n'
               << "\tTotal number of faces        = " << total_num_faces_ << '\n';

            os << "cell index map = [";
            std::copy(cellno_.begin(), cellno_.end(),
                      std::ostream_iterator<int>(os, " "));
            os << "\b]\n";

            os << "cell faces     =\n";
            for (int i = 0; i < cellFaces_.size(); ++i)
            {
                os << "\t[" << i << "] -> [";
                std::copy(cellFaces_[i].begin(), cellFaces_[i].end(),
                          std::ostream_iterator<int>(os, ","));
                os << "\b]\n";
            }
        }

    private:
        typedef FieldMatrix<Scalar, 1, 1> BlockType;

        int                   max_ncf_;
        int                   num_internal_faces_;
        int                   total_num_faces_;

        bool                  matrix_structure_valid_;

        std::vector<int>      cellno_;
        SparseTable<int>      cellFaces_;

        SparseTable<Scalar>   F_;
        std::vector<Scalar>   L_;
        SparseTable<Scalar>   Binv_, f_;

        std::vector<Scalar>   g_;

        BCRSMatrix<BlockType> S_;
        std::vector<Scalar>   rhs_;



        // ----------------------------------------------------------------
        void buildGridTopology(const GridInterface& g)
        // ----------------------------------------------------------------
        {
            typedef typename GridInterface::CellIterator CI;
            typedef typename CI           ::FaceIterator FI;

            const int nc = g.numberOfCells();
            std::vector<int> fpos           ;   fpos  .reserve(nc + 1);
            std::vector<int> num_cf         ;   num_cf.reserve(nc);
            std::vector<int> faces          ;

            // First pass: enumerate internal faces.
            int cellno = 0; fpos.push_back(0);
            for (CI c = g.cellbegin(); c != g.cellend(); ++c, ++cellno) {
                const int c0 = c->index();
                ASSERT((0 <= c0) && (c0 < nc) && (cellno_[c0] == -1));

                cellno_[c0] = cellno;

                num_cf.push_back(0);
                int& ncf = num_cf.back();

                for (FI f = c->facebegin(); f != c-> faceend(); ++f) {
                    if (!f->boundary()) {
                        const int c1 = f->neighbourCellIndex();
                        ASSERT((0 <= c1) && (c1 < nc) && (c1 != c0));

                        if (cellno_[c1] == -1) {
                            // Previously undiscovered internal face.
                            faces.push_back(c1);
                        }
                    }
                    ++ncf;
                }

                fpos.push_back(int(faces.size()));
                max_ncf_ = std::max(max_ncf_, ncf);
            }
            ASSERT(cellno == nc);

            total_num_faces_ = num_internal_faces_ = int(faces.size());

            // Second pass: build cell-to-face mapping, including boundary.
            typedef std::vector<int>::iterator VII;
            for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
                const int c0 = c->index();

                ASSERT((0 <= c0         ) && (c0 < nc         ) &&
                       (0 <= cellno_[c0]) && (cellno_[c0] < nc));

                const int ncf = num_cf[cellno_[c0]];
                std::vector<int>    l2g; l2g.reserve(ncf);
                std::vector<Scalar> Binv_alloc(ncf * ncf);
                std::vector<Scalar> F_alloc(ncf);

                for (FI f = c->facebegin(); f != c->faceend(); ++f) {
                    if (f->boundary()) {
                        // External, not counted before.  Add new face...
                        l2g.push_back(total_num_faces_++);
                    } else {
                        // Internal face.  Need to determine during
                        // traversal of which cell we discovered this
                        // face first, and extract the face number
                        // from the 'faces' table range of that cell.

                        // Note: std::find() below is potentially
                        // *VERY* expensive (e.g., large number of
                        // seeks in moderately sized data in case of
                        // faulted cells).

                        const int c1 = f->neighbourCellIndex();
                        ASSERT((0 <= c1         ) && (c1 < nc         ) &&
                               (0 <= cellno_[c1]) && (cellno_[c1] < nc));

                        int t = c0, seek = c1;
                        if (cellno_[seek] < cellno_[t])
                            std::swap(t, seek);

                        int s = fpos[cellno_[t]], e = fpos[cellno_[t] + 1];

                        VII p = std::find(faces.begin() + s, faces.begin() + e, seek);
                        ASSERT(p != faces.begin() + e);

                        l2g.push_back(s + (p - (faces.begin() + s)));
                    }
                }
                ASSERT(int(l2g.size()) == ncf);

                cellFaces_.appendRow(l2g       .begin(), l2g       .end());
                F_        .appendRow(F_alloc   .begin(), F_alloc   .end());
                f_        .appendRow(F_alloc   .begin(), F_alloc   .end());
                Binv_     .appendRow(Binv_alloc.begin(), Binv_alloc.end());
            }
        }


        // ----------------------------------------------------------------
        void buildSystemStructure()
        // ----------------------------------------------------------------
        {
            ASSERT (!cellFaces_.empty());

            typedef SparseTable<int>::row_type::iterator fi;

            // Clear any residual data, prepare for assembling structure.
            S_.setSize(total_num_faces_, total_num_faces_);
            S_.setBuildMode(BCRSMatrix<BlockType>::random);

            // Compute row sizes
            for (int f = 0; f < total_num_faces_; ++f) {
                S_.setrowsize(f, 1);
            }

            for (int c = 0; c < cellFaces_.size(); ++c) {
                const int nf = cellFaces_[c].size();
                fi fb = cellFaces_[c].begin(), fe = cellFaces_[c].end();

                for (fi f = fb; f != fe; ++f) {
                    S_.incrementrowsize(*f, nf - 1);
                }
            }
            S_.endrowsizes();

            // Compute actual connections (the non-zero structure).
            for (int c = 0; c < cellFaces_.size(); ++c) {
                fi fb = cellFaces_[c].begin(), fe = cellFaces_[c].end();

                for (fi i = fb; i != fe; ++i) {
                    for (fi j = fb; j != fe; ++j) {
                        S_.addindex(*i, *j);
                    }
                }
            }
            S_.endindices();
            std::vector<Scalar>(total_num_faces_).swap(rhs_);
            std::vector<Scalar>(cellFaces_.size()).swap(g_);

            matrix_structure_valid_ = true;
        }


        template<class ReservoirInterface, class BCInterface>
        void assembleDynamic(const GridInterface&       g  ,
                             const ReservoirInterface&  r  ,
                             const BCInterface&         bc ,
                             const std::vector<double>& src,
                             const std::vector<double>& sat)
        {
            typedef typename GridInterface::CellIterator CI;

            std::vector<double> mob(ReservoirInterface::NumberOfPhases);
            std::vector<double> rho(ReservoirInterface::NumberOfPhases);

            std::vector<Scalar> data_store(max_ncf_ * max_ncf_);
            std::vector<Scalar> e  (max_ncf_);
            std::vector<Scalar> rhs(max_ncf_);

            std::vector<unsigned char> dirichlet_faces(max_ncf_);

            // Clear residual data
            S_ = 0;
            std::fill(g_  .begin(), g_  .end(), Scalar(0.0));
            std::fill(rhs_.begin(), rhs_.end(), Scalar(0.0));

            std::fill(e   .begin(), e   .end(), Scalar(1.0));

            // Assemble dynamic contributions for each cell
            for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
                const int ci = c->index();
                const int c0 = cellno_[ci];
                const int nf = cellFaces_[c0].size();

                r.phaseMobility(ci, sat[ci], mob);
                r.phaseDensity (ci,          rho);

                const double totmob = std::accumulate   (mob.begin(), mob.end(), 0.0);
                const double omega  = std::inner_product(rho.begin(), rho.end(),
                                                         mob.begin(), 0.0) / totmob;

                SharedFortranMatrix    S  (nf, nf, &data_store[0]);
                ImmutableFortranMatrix one(nf, 1 , &e[0]);

                setExternalContrib(c, c0, bc, src[ci], rhs, dirichlet_faces);

                buildCellContrib(c0, totmob, omega, one, S, rhs);

                addCellContrib(S, rhs, dirichlet_faces,
                               cellFaces_[c0].begin(),
                               cellFaces_[c0].end());
            }
        }


        void buildCellContrib(const int c, const Scalar totmob, const Scalar omega,
                              const ImmutableFortranMatrix& one,
                              SharedFortranMatrix& S, std::vector<Scalar>& rhs)
        {
            std::transform(Binv_[c].begin(), Binv_[c].end(), S.data(),
                           boost::bind(std::multiplies<Scalar>(), _1, totmob));

            // Ft <- B^{-t} * ones([size(S,2),1])
            SharedFortranMatrix Ft(S.numRows(), 1, &F_[c][0]);
            matMulAdd_TN(Scalar(1.0), S, one, Scalar(0.0), Ft);

            L_[c]  = std::accumulate   (Ft.data(), Ft.data() + Ft.numRows(), 0.0);
            g_[c] -= std::inner_product(Ft.data(), Ft.data() + Ft.numRows(),
                                        f_[c].begin(), Scalar(0.0));

            // rhs <- B^{-1}*f - r (==B^{-1}f + E\pi - f)
            vecMulAdd_N(Scalar(1.0), S, &f_[c][0], -Scalar(1.0), &rhs[0]);

            // rhs <- rhs + g_[c]/L_[c]*F
            std::transform(rhs.begin(), rhs.end(), Ft.data(), rhs.begin(),
                           axpby<Scalar>(Scalar(1.0), Scalar(g_[c] / L_[c])));

            // S <- S - F'*F/L_c
            symmetricUpdate(-1.0/L_[c], Ft, 1.0, S);
        }


        template<class BCInterface>
        void setExternalContrib(const typename GridInterface::CellIterator c,
                                const int c0, const BCInterface& bc,
                                const double src,
                                std::vector<unsigned char> dF,
                                std::vector<Scalar>& rhs)
        {
            typedef typename GridInterface::CellIterator::FaceIterator FI;
            typedef std::vector<unsigned char>::value_type uchar;

            std::fill(f_[c0].begin(), f_[c0].end(), Scalar(0.0));
            std::fill(rhs   .begin(), rhs   .end(), Scalar(0.0));
            std::fill(dF    .begin(), dF    .end(), uchar(0));

            g_[c0] = src;

            int k = 0;
            for (FI f = c->facebegin(); f != c->faceend(); ++f, ++k) {
                if (f->boundary()) {
                    const int bid = f->boundaryId();

                    if (bc[bid].isDirichlet()) {
                        f_[c0][k] = bc[bid].pressure();
                        dF    [k] = uchar(1);
                        rhs   [k] = -f_[c0][k];
                    } else {
                        ASSERT (bc[bid].isNeumann());
                        rhs[k] = bc[bid].outflux();
                    }
                }
            }
        }


        template<class L2GIterator>
        void addCellContrib(const SharedFortranMatrix&        S  ,
                            const std::vector<Scalar>&        rhs,
                            const std::vector<unsigned char>& dF ,
                            L2GIterator&                      b  ,
                            L2GIterator&                      e  )
        {
            int r = 0;
            for (L2GIterator i = b; i != e; ++i, ++r) {

                int c = 0;
                for (L2GIterator j = b; j != e; ++j, ++c) {
                    S_[*i][*j] += S(r,c);
                }

                // Handle Dirichlet (prescribed pressure) conditions.
                if (dF[r]) S_[*i][*i] += 1.0;

                rhs_[*i] += rhs[r];
            }
        }
    };
} // namespace Dune

#endif // OPENRS_INCOMPFLOWSOLVERHYBRID_HEADER
