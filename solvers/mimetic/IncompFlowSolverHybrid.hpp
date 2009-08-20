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

#include "config.h"

#include <algorithm>
#include <functional>
#include <map>
#include <numeric>
#include <ostream>
#include <utility>
#include <vector>

#include <tr1/unordered_map>

#include <boost/bind.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/ErrorMacros.hpp>
#include <dune/common/SparseTable.hpp>

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/io.hh>

#include <dune/istl/overlappingschwarz.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <dune/solvers/common/BoundaryConditions.hpp>

namespace Dune {
    namespace {
        /// @brief
        /// @todo Doc me!
        /// @tparam
        /// @param
        /// @return
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


        /// @brief
        /// @todo Doc me!
        /// @tparam
        template<typename T>
        class axpby : public std::binary_function<T,T,T> {
        public:
            /// @brief
            /// @todo Doc me!
            /// @param
            axpby(const T& a, const T& b) : a_(a), b_(b) {}

            /// @brief
            /// @todoc me!
            /// @return
            T operator()(const T& x, const T& y)
            {
                return a_*x + b_*y;
            }
        private:
            T a_, b_;
        };
    }


    /// @brief
    /// @todo Doc me!
    /// @tparam
    template<class GridInterface,
             class ReservoirInterface,
             class BCInterface,
             class InnerProduct>
    class IncompFlowSolverHybrid {
        typedef typename GridInterface::Scalar Scalar;
        enum FaceType { Internal, Dirichlet, Neumann, Periodic };

        class FlowSolution {
        public:
            /// @brief
            /// @todo Doc me!
            typedef typename GridInterface::Scalar       Scalar;
            typedef typename GridInterface::CellIterator CI;
            typedef typename CI           ::FaceIterator FI;

            friend class IncompFlowSolverHybrid;

            /// @brief
            /// @todo Doc me!
            /// @param
            /// @return
            Scalar pressure(const CI& c) const
            {
                return pressure_[cellno_[c->index()]];
            }
            /// @brief
            /// @todo Doc me!
            /// @param
            /// @return
            Scalar outflux (const FI& f) const
            {
                return outflux_[cellno_[f->cellIndex()]][f->localIndex()];
            }
        private:
            std::vector< int  > cellno_;
            SparseTable< int  > cellFaces_;
            std::vector<Scalar> pressure_;
            SparseTable<Scalar> outflux_;
        };

    public:
        /// @brief
        /// @todo Doc me!
        /// @param
        void init(const GridInterface&      g,
                  const ReservoirInterface& r,
                  const BCInterface&        bc)
        {
            clear();

            if (g.numberOfCells() > 0) {
                enumerateDof(g, bc);
                allocateConnections(bc);
                setConnections(bc);
                computeInnerProducts(r);
            }
        }

        /// @brief
        /// @todo Doc me!
        /// @param
        void clear()
        {
            pgrid_                  =  0;
            max_ncf_                = -1;
            num_internal_faces_     =  0;
            total_num_faces_        =  0;
            matrix_structure_valid_ = false;
            do_regularization_      = true; // Assume pure Neumann by default.

            bdry_id_map_.clear();

            std::vector<Scalar>().swap(L_);
            std::vector<Scalar>().swap(g_);
            Binv_.clear();   F_.clear();   f_.clear();

            std::vector<int>().swap(flowSolution_.cellno_);
            flowSolution_.cellFaces_.clear();

            cleared_state_ = true;
        }

        /// @brief
        /// @todo Doc me!
        /// @param
        void enumerateDof(const GridInterface& g, const BCInterface& bc)
        {
            ASSERT2 (cleared_state_,
                     "You must call clear() prior to enumerateDof()");
            ASSERT  (topologyIsSane(g));

            enumerateGridDof(g);
            enumerateBCDof(g, bc);

            pgrid_ = &g;
            cleared_state_ = false;
        }

        /// @brief
        /// @todo Doc me!
        /// @param
        void allocateConnections(const BCInterface& bc)
        {
            ASSERT2 (!cleared_state_,
                     "You must call enumerateDof() prior "
                     "to allocateConnections()");
            ASSERT  (!matrix_structure_valid_);

            // Clear any residual data, prepare for assembling structure.
            S_.setSize(total_num_faces_, total_num_faces_);
            S_.setBuildMode(BCRSMatrix<MatrixBlockType>::random);

            // Compute row sizes
            for (int f = 0; f < total_num_faces_; ++f) {
                S_.setrowsize(f, 1);
            }

            allocateGridConnections();
            allocateBCConnections(bc);

            S_.endrowsizes();

            rhs_ .resize(total_num_faces_);
            soln_.resize(total_num_faces_);
        }

        /// @brief
        /// @todo Doc me!
        /// @param
        void setConnections(const BCInterface& bc)
        {
            setGridConnections();
            setBCConnections(bc);

            S_.endindices();

            const int nc = pgrid_->numberOfCells();
            std::vector<Scalar>(nc).swap(flowSolution_.pressure_);
            std::vector<Scalar>(nc).swap(g_);
            std::vector<Scalar>(nc).swap(L_);

            matrix_structure_valid_ = true;
        }

        /// @brief
        /// @todo Doc me!
        /// @param
        void computeInnerProducts(const ReservoirInterface& r)
        {
            ASSERT2 (matrix_structure_valid_,
                     "You must call connectionsCompleted() prior "
                     "to computeInnerProducts()");

            typedef typename GridInterface     ::CellIterator CI;
            typedef typename ReservoirInterface::PermTensor   PermTensor;

            InnerProduct ip(max_ncf_);
            int i = 0;
            const SparseTable<int>& cellFaces = flowSolution_.cellFaces_;
            for (CI c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c, ++i) {
                const int nf = cellFaces[i].size();

                SharedFortranMatrix Binv(nf, nf, &Binv_[i][0]);

                ip.evaluate(c, r.permeability(c->index()), Binv);
            }
        }


        /// @brief
        /// @todo Doc me!
        /// @param
        void solve(const ReservoirInterface&  r  ,
                   const std::vector<double>& sat,
                   const BCInterface&         bc ,
                   const std::vector<double>& src,
                   const typename GridInterface::CellIterator::Vector& grav)
        {
            assembleDynamic(r, sat, bc, src, grav);
            // printSystem("linsys_mimetic");
#if 0
            solveLinearSystem();
#else
            solveLinearSystemAMG();
#endif
            computePressureAndFluxes(r, sat);
        }


        /// @brief
        /// @todo Doc me!
        typedef const FlowSolution& SolutionType;

        /// @brief
        /// @todo Doc me!
        /// @return
        SolutionType getSolution()
        {
            return flowSolution_;
        }


        /// @brief
        /// @todo Doc me!
        /// @tparam
        /// @param
        template<typename charT, class traits>
        void printStats(std::basic_ostream<charT,traits>& os)
        {
            os << "IncompFlowSolverHybrid<>:\n"
               << "\tMaximum number of cell faces = " << max_ncf_ << '\n'
               << "\tNumber of internal faces     = " << num_internal_faces_ << '\n'
               << "\tTotal number of faces        = " << total_num_faces_ << '\n';

            const std::vector<int>& cell = flowSolution_.cellno_;
            os << "cell index map = [";
            std::copy(cell.begin(), cell.end(),
                      std::ostream_iterator<int>(os, " "));
            os << "\b]\n";

            const SparseTable<int>& cf = flowSolution_.cellFaces_;
            os << "cell faces     =\n";
            for (int i = 0; i < cf.size(); ++i)
            {
                os << "\t[" << i << "] -> [";
                std::copy(cf[i].begin(), cf[i].end(),
                          std::ostream_iterator<int>(os, ","));
                os << "\b]\n";
            }
        }

        /// @brief
        /// @todo Doc me!
        /// @tparam
        /// @param
        template<class charT, class traits>
        void printIP(std::basic_ostream<charT,traits>& os)
        {
            const SparseTable<int>& cf = flowSolution_.cellFaces_;
            // Loop grid whilst building (and outputing) the inverse IP matrix.
            for (int c = 0; c != cf.size(); ++c) {
                const int nf = cf[c].size();
                ImmutableFortranMatrix Binv(nf, nf, &Binv_[c][0]);

                os << c << " -> Binv = [\n" << Binv << "]\n";
            }
        }

        /// @brief
        /// @todo Doc me!
        /// @param
        void printSystem(const std::string& prefix)
        {
            writeMatrixToMatlab(S_, prefix + "-mat.dat");

            std::string rhsfile(prefix + "-rhs.dat");
            std::ofstream rhs(rhsfile.c_str());
            std::copy(rhs_.begin(), rhs_.end(),
                      std::ostream_iterator<VectorBlockType>(rhs, "\n"));
        }

    private:
        typedef std::pair<int,int>                 DofID;
        typedef std::tr1::unordered_map<int,DofID> BdryIdMapType;
        typedef BdryIdMapType::const_iterator      BdryIdMapIterator;

        const GridInterface* pgrid_;
        BdryIdMapType        bdry_id_map_;

        // ----------------------------------------------------------------
        bool cleared_state_;
        int  max_ncf_;
        int  num_internal_faces_;
        int  total_num_faces_;

        // ----------------------------------------------------------------
        std::vector<Scalar> L_, g_;
        SparseTable<Scalar> Binv_, F_, f_;

        // ----------------------------------------------------------------
        // Actual, assembled system of linear equations
        typedef FieldVector<Scalar, 1   > VectorBlockType;
        typedef FieldMatrix<Scalar, 1, 1> MatrixBlockType;

        BCRSMatrix <MatrixBlockType>      S_;    // System matrix
        BlockVector<VectorBlockType>      rhs_;  // System RHS
        BlockVector<VectorBlockType>      soln_; // System solution (contact pressure)
        bool                              matrix_structure_valid_;
        bool                              do_regularization_;

        // ----------------------------------------------------------------
        // Physical quantities (derived)
        FlowSolution flowSolution_;


        // ----------------------------------------------------------------
        void enumerateGridDof(const GridInterface& g)
        // ----------------------------------------------------------------
        {
            typedef typename GridInterface::CellIterator CI;
            typedef typename CI           ::FaceIterator FI;

            const int nc = g.numberOfCells();
            std::vector<int> fpos           ;   fpos  .reserve(nc + 1);
            std::vector<int> num_cf         ;   num_cf.reserve(nc);
            std::vector<int> faces          ;

            // Allocate cell structures.
            std::vector<int>(nc, -1).swap(flowSolution_.cellno_);

            std::vector<int>& cell = flowSolution_.cellno_;

            // First pass: enumerate internal faces.
            int cellno = 0; fpos.push_back(0);
            for (CI c = g.cellbegin(); c != g.cellend(); ++c, ++cellno) {
                const int c0 = c->index();
                ASSERT((0 <= c0) && (c0 < nc) && (cell[c0] == -1));

                cell[c0] = cellno;

                num_cf.push_back(0);
                int& ncf = num_cf.back();

                for (FI f = c->facebegin(); f != c-> faceend(); ++f) {
                    if (!f->boundary()) {
                        const int c1 = f->neighbourCellIndex();
                        ASSERT((0 <= c1) && (c1 < nc) && (c1 != c0));

                        if (cell[c1] == -1) {
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

            SparseTable<int>& cf = flowSolution_.cellFaces_;

            // Second pass: build cell-to-face mapping, including boundary.
            typedef std::vector<int>::iterator VII;
            for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
                const int c0 = c->index();

                ASSERT ((0 <=      c0 ) && (     c0  < nc) &&
                        (0 <= cell[c0]) && (cell[c0] < nc));

                const int ncf = num_cf[cell[c0]];
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
                        ASSERT ((0 <=      c1 ) && (     c1  < nc) &&
                                (0 <= cell[c1]) && (cell[c1] < nc));

                        int t = c0, seek = c1;
                        if (cell[seek] < cell[t])
                            std::swap(t, seek);

                        int s = fpos[cell[t]], e = fpos[cell[t] + 1];

                        VII p = std::find(faces.begin() + s, faces.begin() + e, seek);
                        ASSERT(p != faces.begin() + e);

                        l2g.push_back(s + (p - (faces.begin() + s)));
                    }
                }
                ASSERT (int(l2g.size()) == ncf);

                cf   .appendRow(l2g       .begin(), l2g       .end());
                F_   .appendRow(F_alloc   .begin(), F_alloc   .end());
                f_   .appendRow(F_alloc   .begin(), F_alloc   .end());
                Binv_.appendRow(Binv_alloc.begin(), Binv_alloc.end());

                flowSolution_.outflux_
                    .appendRow (F_alloc   .begin(), F_alloc   .end());
            }
        }


        // ----------------------------------------------------------------
        void enumerateBCDof(const GridInterface& g, const BCInterface& bc)
        // ----------------------------------------------------------------
        {
            typedef typename GridInterface::CellIterator CI;
            typedef typename CI           ::FaceIterator FI;

            const std::vector<int>& cell = flowSolution_.cellno_;

            bdry_id_map_.clear();
            for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
                for (FI f = c->facebegin(); f != c->faceend(); ++f) {
                    if (f->boundary()) {
                        const int bid = f->boundaryId();
                        if (bc.flowCond(bid).isPeriodic()) {
                            DofID dof(cell[c->index()], f->localIndex());
                            bdry_id_map_.insert(std::make_pair(bid, dof));
                        }
                    }
                }
            }
        }



        // ----------------------------------------------------------------
        void allocateGridConnections()
        // ----------------------------------------------------------------
        {
            const   SparseTable<int>& cf = flowSolution_.cellFaces_;
            typedef SparseTable<int>::row_type::const_iterator fi;

            for (int c = 0; c < cf.size(); ++c) {
                const int nf = cf[c].size();
                fi fb = cf[c].begin(), fe = cf[c].end();

                for (fi f = fb; f != fe; ++f) {
                    S_.incrementrowsize(*f, nf - 1);
                }
            }
        }


        // ----------------------------------------------------------------
        void allocateBCConnections(const BCInterface& bc)
        // ----------------------------------------------------------------
        {
            // Include system connections for periodic boundary
            // conditions, if any.  We make an arbitrary choice in
            // that the face/degree-of-freedom with the minimum
            // numerical id of the two periodic partners represents
            // the coupling.  Suppose <i_p> is this minimum dof-id.
            // We then need to introduce a *symmetric* coupling to
            // <i_p> to each of the dof's of the cell *NOT* containing
            // <i_p>.  This choice is implemented in the following
            // loop by introducing couplings only when dof1 (self) is
            // less than dof2 (other).
            //
            // See also: setBCConnections() and addCellContrib().
            //
            typedef typename GridInterface::CellIterator CI;
            typedef typename CI           ::FaceIterator FI;

            const std::vector<int>& cell = flowSolution_.cellno_;
            const SparseTable<int>& cf   = flowSolution_.cellFaces_;

            if (!bdry_id_map_.empty()) {
                // At least one periodic BC.  Allocate corresponding
                // connections.
                for (CI c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c) {
                    for (FI f = c->facebegin(); f != c->faceend(); ++f) {
                        if (f->boundary()) {
                            const int bid = f->boundaryId();
                            if (bc.flowCond(bid).isPeriodic()) {
                                // dof-id of self
                                const int dof1 = cf[cell[c->index()]][f->localIndex()];

                                // dof-id of other
                                BdryIdMapIterator j =
                                    bdry_id_map_.find(bc.getPeriodicPartner(bid));
                                ASSERT (j != bdry_id_map_.end());
                                const int c2   = j->second.first;
                                const int dof2 = cf[c2][j->second.second];

                                if (dof1 < dof2) {
                                    // Allocate storage for the actual
                                    // couplings.
                                    //
                                    // Note: The equivalence
                                    // dof1==dof2 is handled
                                    // automatically in this loop.
                                    //
                                    const int ndof = cf.rowSize(c2);
                                    S_.incrementrowsize(dof1, ndof);          // self->other
                                    for (int dof = 0; dof < ndof; ++dof) {
                                        S_.incrementrowsize(cf[c2][dof], 1);  // other->self
                                    }
                               }
                            }
                        }
                    }
                }
            }
        }



        // ----------------------------------------------------------------
        void setGridConnections()
        // ----------------------------------------------------------------
        {
            const   SparseTable<int>& cf = flowSolution_.cellFaces_;
            typedef SparseTable<int>::row_type::const_iterator fi;

            // Compute actual connections (the non-zero structure).
            for (int c = 0; c < cf.size(); ++c) {
                fi fb = cf[c].begin(), fe = cf[c].end();

                for (fi i = fb; i != fe; ++i) {
                    for (fi j = fb; j != fe; ++j) {
                        S_.addindex(*i, *j);
                    }
                }
            }
        }


        // ----------------------------------------------------------------
        void setBCConnections(const BCInterface& bc)
        // ----------------------------------------------------------------
        {
            // Include system connections for periodic boundary
            // conditions, if any.  We make an arbitrary choice in
            // that the face/degree-of-freedom with the minimum
            // numerical id of the two periodic partners represents
            // the coupling.  Suppose <i_p> is this minimum dof-id.
            // We then need to introduce a *symmetric* coupling to
            // <i_p> to each of the dof's of the cell *NOT* containing
            // <i_p>.  This choice is implemented in the following
            // loop by introducing couplings only when dof1 (self) is
            // less than dof2 (other).
            //
            // See also: allocateBCConnections() and addCellContrib().
            //
            typedef typename GridInterface::CellIterator CI;
            typedef typename CI           ::FaceIterator FI;

            const std::vector<int>& cell = flowSolution_.cellno_;
            const SparseTable<int>& cf   = flowSolution_.cellFaces_;

            if (!bdry_id_map_.empty()) {
                // At least one periodic BC.  Assign periodic
                // connections.
                for (CI c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c) {
                    for (FI f = c->facebegin(); f != c->faceend(); ++f) {
                        if (f->boundary()) {
                            const int bid = f->boundaryId();
                            if (bc.flowCond(bid).isPeriodic()) {
                                // dof-id of self
                                const int dof1 = cf[cell[c->index()]][f->localIndex()];

                                // dof-id of other
                                BdryIdMapIterator j =
                                    bdry_id_map_.find(bc.getPeriodicPartner(bid));
                                ASSERT (j != bdry_id_map_.end());
                                const int c2   = j->second.first;
                                const int dof2 = cf[c2][j->second.second];

                                if (dof1 < dof2) {
                                    // Assign actual couplings.
                                    //
                                    // Note: The equivalence
                                    // dof1==dof2 is handled
                                    // automatically in this loop.
                                    //
                                    const int ndof = cf.rowSize(c2);
                                    for (int dof = 0; dof < ndof; ++dof) {
                                        S_.addindex(dof1, cf[c2][dof]);  // self->other
                                        S_.addindex(cf[c2][dof], dof1);  // other->self
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }



        // ----------------------------------------------------------------
        void assembleDynamic(const ReservoirInterface&  r  ,
                             const std::vector<double>& sat,
                             const BCInterface&         bc ,
                             const std::vector<double>& src,
                             const typename GridInterface::CellIterator::Vector& grav)
        // ----------------------------------------------------------------
        {
            typedef typename GridInterface::CellIterator CI;

            const std::vector<int>& cell = flowSolution_.cellno_;
            const SparseTable<int>& cf   = flowSolution_.cellFaces_;

            std::vector<double> mob(ReservoirInterface::NumberOfPhases);
            std::vector<double> rho(ReservoirInterface::NumberOfPhases);

            std::vector<Scalar> data_store(max_ncf_ * max_ncf_);
            std::vector<Scalar> e  (max_ncf_);
            std::vector<Scalar> rhs(max_ncf_);

            std::vector<FaceType> facetype(max_ncf_);
            std::vector<Scalar>   condval (max_ncf_);
            std::vector<int>      ppartner(max_ncf_);

            // Clear residual data
            S_   = 0.0;
            rhs_ = 0.0;

            std::fill(g_.begin(), g_.end(), Scalar(0.0));
            std::fill(L_.begin(), L_.end(), Scalar(0.0));
            std::fill(e .begin(), e .end(), Scalar(1.0));

            // We will have to regularize resulting system if there
            // are no prescribed pressures (i.e., Dirichlet BC's).
            do_regularization_ = true;

            InnerProduct ip(max_ncf_);

            // Assemble dynamic contributions for each cell
            for (CI c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c) {
                const int ci = c->index();
                const int c0 = cell[ci];            ASSERT (c0 < cf.size());
                const int nf = cf[c0].size();

                r.phaseMobility(ci, sat[ci], mob);
                r.phaseDensity (ci,          rho);

                const double totmob = std::accumulate   (mob.begin(), mob.end(), 0.0);
                const double omega  = std::inner_product(rho.begin(), rho.end(),
                                                         mob.begin(), 0.0) / totmob;

                SharedFortranMatrix    S  (nf, nf, &data_store[0]);
                ImmutableFortranMatrix one(nf, 1 , &e[0]);

                typename SparseTable<double>::mutable_row_type gterm = f_[c0];
                std::fill(gterm.begin(), gterm.end(), Scalar(0.0));
                ip.gravityTerm(c, grav, omega, gterm);

                setExternalContrib(c, c0, bc, src[ci], rhs,
                                   facetype, condval, ppartner);

                buildCellContrib(c0, totmob, omega, one, S, rhs);

                addCellContrib(S, rhs, facetype, condval, ppartner, cf[c0]);
            }
        }



        // ----------------------------------------------------------------
        void solveLinearSystem()
        // ----------------------------------------------------------------
        {
            // Adapted from DuMux...
            Scalar residTol = 1.0e-12;

            typedef BCRSMatrix <MatrixBlockType>        Matrix;
            typedef BlockVector<VectorBlockType>        Vector;
            typedef MatrixAdapter<Matrix,Vector,Vector> Adapter;

            // Regularize the matrix (only for pure Neumann problems...)
            if (do_regularization_) {
                S_[0][0] *= 2;
            }
            Adapter opS(S_);

            // initialize the preconditioner
            Dune::SeqILU0<Matrix,Vector,Vector> precond(S_, 1.0);

            // invert the linear equation system
            Dune::BiCGSTABSolver<Vector> linsolve(opS, precond, residTol, 500, 1);

            Dune::InverseOperatorResult result;
            soln_ = 0.0;
            linsolve.apply(soln_, rhs_, result);
        }



        // ----------------------------------------------------------------
        void solveLinearSystemAMG()
        // ----------------------------------------------------------------
        {
            // Adapted from upscaling.cc by Arne Rekdal, 2009
            Scalar residTol = 1.0e-8;

            // Representation types for linear system.
            typedef BCRSMatrix <MatrixBlockType>        Matrix;
            typedef BlockVector<VectorBlockType>        Vector;
            typedef MatrixAdapter<Matrix,Vector,Vector> Operator;

            // AMG specific types.
#define FIRST_DIAGONAL 1
#define SYMMETRIC 1

#if FIRST_DIAGONAL
            typedef Amg::FirstDiagonal CouplingMetric;
#else
            typedef Amg::RowSum        CouplingMetric;
#endif

#if SYMMETRIC
            typedef Amg::SymmetricCriterion<Matrix,CouplingMetric>   CriterionBase;
#else
            typedef Amg::UnSymmetricCriterion<Matrix,CouplingMetric> CriterionBase;
#endif

            typedef SeqILU0<Matrix,Vector,Vector>        Smoother;
            typedef Amg::CoarsenCriterion<CriterionBase> Criterion;
            typedef Amg::AMG<Operator,Vector,Smoother>   Precond;

            // Regularize the matrix (only for pure Neumann problems...)
            if (do_regularization_) {
                S_[0][0] *= 2;
            }
            Operator opS(S_);

            // initialize the preconditioner
            double relax = 1;
            typename Precond::SmootherArgs smootherArgs;
            smootherArgs.relaxationFactor = relax;

            Criterion criterion;
            Precond precond(opS, criterion, smootherArgs);

            // invert the linear equation system
            int verbose = 1;
            CGSolver<Vector> linsolve(opS, precond, residTol, S_.N(), verbose);

            InverseOperatorResult result;
            soln_ = 0.0;
            linsolve.apply(soln_, rhs_, result);
        }



        // ----------------------------------------------------------------
        void computePressureAndFluxes(const ReservoirInterface&  r  ,
                                      const std::vector<double>& sat)
        // ----------------------------------------------------------------
        {
            typedef typename GridInterface::CellIterator CI;

            const std::vector<int>& cell = flowSolution_.cellno_;
            const SparseTable<int>& cf   = flowSolution_.cellFaces_;

            std::vector<Scalar>& p = flowSolution_.pressure_;
            SparseTable<Scalar>& v = flowSolution_.outflux_;

            std::vector<double> mob(ReservoirInterface::NumberOfPhases);
            std::vector<double> pi (max_ncf_);

            // Assemble dynamic contributions for each cell
            for (CI c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c) {
                const int ci = c->index();
                const int c0 = cell[ci];
                const int nf = cf[c0].size();

                // Extract contact pressures for cell 'c'.
                for (int i = 0; i < nf; ++i) {
                    pi[i] = soln_[cf[c0][i]];
                }

                // Compute cell pressure in cell 'c'.
                p[c0] = (g_[c0] +
                         std::inner_product(F_[c0].begin(), F_[c0].end(),
                                            pi.begin(), 0.0)) / L_[c0];

                // Compute cell (out) fluxes for cell 'c'.
                // 1) Form system right hand side, r = f + Cp - D\pi
                std::transform(f_[c0].begin(), f_[c0].end(), pi.begin(),
                               pi    .begin(), //_1 + p[c0] - _2);
                               boost::bind(std::minus<Scalar>(),
                                           boost::bind(std::plus<Scalar>(),
                                                       _1,
                                                       p[c0]),
                                           _2));

                // 2) Solve system Bv = r
                r.phaseMobility(ci, sat[ci], mob);
                const double totmob = std::accumulate(mob.begin(), mob.end(), 0.0);

                ImmutableFortranMatrix Binv(nf, nf, &Binv_[c0][0]);
                vecMulAdd_N(totmob, Binv, &pi[0], Scalar(0.0), &v[c0][0]);
            }
        }




        // ----------------------------------------------------------------
        void setExternalContrib(const typename GridInterface::CellIterator c,
                                const int c0, const BCInterface& bc,
                                const double           src,
                                std::vector<Scalar>&   rhs,
                                std::vector<FaceType>& facetype,
                                std::vector<double>&   condval,
                                std::vector<int>&      ppartner)
        // ----------------------------------------------------------------
        {
            typedef typename GridInterface::CellIterator::FaceIterator FI;

            const std::vector<int>& cell = flowSolution_.cellno_;
            const SparseTable<int>& cf   = flowSolution_.cellFaces_;

            std::fill(rhs     .begin(), rhs     .end(), Scalar(0.0));
            std::fill(facetype.begin(), facetype.end(), Internal   );
            std::fill(condval .begin(), condval .end(), Scalar(0.0));
            std::fill(ppartner.begin(), ppartner.end(), -1         );

            g_[c0] = src;

            int k = 0;
            for (FI f = c->facebegin(); f != c->faceend(); ++f, ++k) {
                if (f->boundary()) {
                    const int bid = f->boundaryId();
		    const FlowBC& bcond = bc.flowCond(bid);
                    if (bcond.isDirichlet()) {
                        facetype[k]        = Dirichlet;
                        condval[k]         = bcond.pressure();
                        do_regularization_ = false;
                    } else if (bcond.isPeriodic()) {
                        BdryIdMapIterator j =
                            bdry_id_map_.find(bc.getPeriodicPartner(bid));
                        ASSERT (j != bdry_id_map_.end());

                        facetype[k] = Periodic;
                        condval[k]  = bcond.pressureDifference();
                        ppartner[k] = cf[j->second.first][j->second.second];
                    } else {
                        ASSERT (bcond.isNeumann());
                        facetype[k] = Neumann;
                        rhs[k]      = bcond.outflux();
                    }
                }
            }
        }




        // ----------------------------------------------------------------
        void buildCellContrib(const int c, const Scalar totmob, const Scalar omega,
                              const ImmutableFortranMatrix& one,
                              SharedFortranMatrix& S, std::vector<Scalar>& rhs)
        // ----------------------------------------------------------------
        {
            std::transform(Binv_[c].begin(), Binv_[c].end(), S.data(),
                           boost::bind(std::multiplies<Scalar>(), _1, totmob));

            // Ft <- B^{-t} * ones([size(S,2),1])
            SharedFortranMatrix Ft(S.numRows(), 1, &F_[c][0]);
            matMulAdd_TN(Scalar(1.0), S, one, Scalar(0.0), Ft);

            L_[c]  = std::accumulate   (Ft.data(), Ft.data() + Ft.numRows(), 0.0);
            g_[c] -= std::inner_product(Ft.data(), Ft.data() + Ft.numRows(),
                                        f_[c].begin(), Scalar(0.0));

            // rhs <- B^{-1}*f - r (==B^{-1}f + E\pi - h)
            vecMulAdd_N(Scalar(1.0), S, &f_[c][0], -Scalar(1.0), &rhs[0]);

            // rhs <- rhs + g_[c]/L_[c]*F
            std::transform(rhs.begin(), rhs.end(), Ft.data(), rhs.begin(),
                           axpby<Scalar>(Scalar(1.0), Scalar(g_[c] / L_[c])));

            // S <- S - F'*F/L_c
            symmetricUpdate(-1.0/L_[c], Ft, 1.0, S);
        }



        // ----------------------------------------------------------------
        /// \param l2g local-to-global face map.
        template<class L2G>
        void addCellContrib(const SharedFortranMatrix&   S       ,
                            const std::vector<Scalar>&   rhs     ,
                            const std::vector<FaceType>& facetype,
                            const std::vector<Scalar>&   condval ,
                            const std::vector<int>&      ppartner,
                            const L2G&                   l2g)
        // ----------------------------------------------------------------
        {
            typedef typename L2G::const_iterator it;

            int r = 0;
            for (it i = l2g.begin(); i != l2g.end(); ++i, ++r) {
                // Indirection for periodic BC handling.
                int ii = *i;

                switch (facetype[r]) {
                case Dirichlet:
                    // Pressure is already known.  Assemble trivial
                    // equation of the form: a*x = a*p where 'p' is
                    // the known pressure value (i.e., condval[r]).
                    //
                    S_[*i][*i] = S(r,r);
                    rhs_[*i]   = S(r,r) * condval[r];
                    break;
                case Periodic:
                    // Periodic boundary condition.  Contact pressures
                    // linked by relations of the form
                    //
                    //      a*(x0 - x1) = a*pd
                    //
                    // where 'pd' is the known pressure difference
                    // x0-x1 (condval[r]).  Preserve matrix symmetry
                    // by assembling both of the equations
                    //
                    //      a*(x0-x1) = a*  pd, and  (1)
                    //      a*(x1-x0) = a*(-pd)      (2)
                    //
                    ASSERT ((0 <= ppartner[r]) && (ppartner[r] < rhs_.size()));
                    ASSERT (ii != ppartner[r]);

                    {
                        const double a = S(r,r), b = a * condval[r];

                        // Equation (1)
                        S_  [         ii][         ii] += a;
                        S_  [         ii][ppartner[r]] -= a;
                        rhs_[         ii]              += b;

                        // Equation (2)
                        S_  [ppartner[r]][         ii] -= a;
                        S_  [ppartner[r]][ppartner[r]] += a;
                        rhs_[ppartner[r]]              -= b;
                    }

                    ii = std::min(ii, ppartner[r]);

                    // INTENTIONAL FALL-THROUGH!
                    // IOW: Don't insert <break;> here!
                    //
                default:
                    int c = 0;
                    for (it j = l2g.begin(); j != l2g.end(); ++j, ++c) {
                        // Indirection for periodic BC handling.
                        int jj = *j;

                        if (facetype[c] == Dirichlet) {
                            rhs_[ii] -= S(r,c) * condval[c];
                            continue;
                        }
                        if (facetype[c] == Periodic) {
                            ASSERT ((0 <= ppartner[c]) && (ppartner[c] < rhs_.size()));
                            ASSERT (jj != ppartner[c]);
                            rhs_[ii] -= S(r,c) * condval[c];
                            jj = std::min(jj, ppartner[c]);
                        }
                        S_[ii][jj] += S(r,c);
                    }
                    break;
                }
                rhs_[ii] += rhs[r];
            }
        }
    };
} // namespace Dune

#endif // OPENRS_INCOMPFLOWSOLVERHYBRID_HEADER
