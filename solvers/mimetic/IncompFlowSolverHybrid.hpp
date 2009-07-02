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
#include <map>
#include <ostream>
#include <utility>
#include <vector>

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
            Binv_     .clear();

            if (g.numberOfCells() > 0) {
                buildGridTopology(g);
                buildMatrixStructure();
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
                const int nf = cellFaces_[c].size();

                SharedFortranMatrix Binv(nf, nf, &Binv_[i][0]);

                ip.evaluate(c, r.permeability(c->index()), Binv);
            }
        }


        void printStats(std::ostream& os)
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
        std::vector<int>      cellno_;
        SparseTable<int>      cellFaces_;
        SparseTable<Scalar>   Binv_;
        BCRSMatrix<BlockType> S_;
        bool                  matrix_structure_valid_;



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
                Binv_     .appendRow(Binv_alloc.begin(), Binv_alloc.end());
            }
        }


        // ----------------------------------------------------------------
        void buildMatrixStructure()
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

            matrix_structure_valid_ = true;
        }
    };
} // namespace Dune

#endif // OPENRS_INCOMPFLOWSOLVERHYBRID_HEADER
