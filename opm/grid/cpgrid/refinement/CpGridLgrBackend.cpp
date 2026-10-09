/*
  Copyright 2026 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/refinement/conforming/ConformingBlockBuilder.hpp>
#include <opm/grid/cpgrid/refinement/GridStateWriter.hpp>
#include <opm/grid/cpgrid/refinement/RetainedCornerPointInput.hpp>
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Dune
{

void CpGrid::setLgrBackend(Opm::Refinement::Backend backend)
{
    if (maxLevel() > 0 && backend != lgr_backend_) {
        OPM_THROW(std::logic_error, "The LGR backend cannot change once the grid is refined.");
    }
    lgr_backend_ = backend;
    // Only the Conforming builder resamples the corner-point input.
    current_data_->front()->retain_cp_input_ = (backend == Opm::Refinement::Backend::Conforming);
}

void CpGrid::throwIfConforming_(const std::string& what) const
{
    if (lgr_backend_ == Opm::Refinement::Backend::Conforming) {
        OPM_THROW(std::logic_error, what + " is not available with the Conforming LGR backend; "
                  "refine with addLgrsUpdateLeafView() instead.");
    }
}

void CpGrid::globalRefineConforming_(int refCount)
{
    if (refCount < 0) {
        OPM_THROW(std::logic_error, "Invalid argument. Provide a nonnegative integer for global refinement.");
    }
    if (refCount == 0) {
        return;
    }
    // One builder pass yields one level; globalRefine(n) promises n nested levels.
    if (refCount > 1 || maxLevel() > 0) {
        OPM_THROW(std::logic_error, "The Conforming LGR backend refines an unrefined grid once; "
                  "use globalRefine(1) or addLgrsUpdateLeafView() with the factor you need.");
    }
    const auto dims = logicalCartesianSize();
    addLgrsUpdateLeafView({{2, 2, 2}}, {{0, 0, 0}}, {{dims[0], dims[1], dims[2]}}, {"GR1"});
}

// Same-level neighbours share a whole face of the cell, so they take its corner average like two
// coarse cells; an LGR boundary face takes its own points, or its centroid when it is not planar.
Dune::FieldVector<double,3> CpGrid::faceCenterEclConforming_(int cell_index, int face,
                                                             const Dune::cpgrid::Intersection& intersection) const
{
    static const int faceVxMap[6][4] = { {0, 2, 4, 6}, {1, 3, 5, 7}, {0, 1, 4, 5},
                                         {2, 3, 6, 7}, {0, 1, 2, 3}, {4, 5, 6, 7} };
    const bool sameLevelNeighbours = intersection.neighbor() &&
        (intersection.inside().level() == intersection.outside().level());
    const bool coarseOnBoundary = intersection.boundary() && !intersection.neighbor() &&
        (intersection.inside().level() == 0);
    if (!sameLevelNeighbours && !coarseOnBoundary) {
        const auto& fp = current_data_->back()->face_to_point_[intersection.id()];
        if (fp.size() == 4) {
            std::array<Dune::FieldVector<double,3>,4> v;
            int k = 0;
            for (auto it = fp.begin(); it != fp.end(); ++it) {
                v[k++] = vertexPosition(*it);
            }
            const auto e1 = v[1] - v[0];
            const auto e2 = v[2] - v[0];
            const Dune::FieldVector<double,3> nrm{ e1[1]*e2[2] - e1[2]*e2[1],
                                                   e1[2]*e2[0] - e1[0]*e2[2],
                                                   e1[0]*e2[1] - e1[1]*e2[0] };
            const double nn = nrm.two_norm();
            const double dev = (nn > 0.0) ? std::abs((v[3] - v[0]) * nrm) / nn : 0.0;
            if (dev <= 1e-9 * std::max(e1.two_norm() + e2.two_norm(), 1.0)) {
                auto center = v[0];
                center += v[1];
                center += v[2];
                center += v[3];
                center /= 4.0;
                return center;
            }
        }
        return intersection.geometry().center();
    }
    Dune::FieldVector<double,3> center(0.0);
    for (int i = 0; i < 4; ++i) {
        center += vertexPosition(current_data_->back()->cell_to_point_[cell_index][faceVxMap[face][i]]);
    }
    center /= 4.0;
    return center;
}

void CpGrid::setRefinementBuilder(std::shared_ptr<Opm::Refinement::Builder> builder)
{
    refinement_builder_ = std::move(builder);
}

// The input is retained on rank 0's undistributed level zero; every rank's builder resamples it.
// Entered on all ranks, whether or not the current view already holds it.
void CpGrid::broadcastRetainedInput_()
{
    auto input = std::make_shared<Opm::Refinement::RetainedCornerPointInput>();
    int present = 0;
    if (comm().rank() == 0 && !data_.empty()) {
        if (const auto source = Opm::Refinement::GridStateWriter::retainedCornerPointInput(*data_[0])) {
            *input = *source;
            present = 1;
        }
    }
    comm().broadcast(&present, 1, 0);
    if (!present) {
        return;
    }
    comm().broadcast(input->dims.data(), 3, 0);
    int edgeConformal = input->edgeConformal ? 1 : 0;
    comm().broadcast(&edgeConformal, 1, 0);
    input->edgeConformal = (edgeConformal != 0);
    std::array<int,3> sizes{static_cast<int>(input->coord.size()),
                            static_cast<int>(input->zcorn.size()),
                            static_cast<int>(input->actnum.size())};
    comm().broadcast(sizes.data(), 3, 0);
    input->coord.resize(sizes[0]);
    input->zcorn.resize(sizes[1]);
    input->actnum.resize(sizes[2]);
    if (sizes[0] > 0) comm().broadcast(input->coord.data(), sizes[0], 0);
    if (sizes[1] > 0) comm().broadcast(input->zcorn.data(), sizes[1], 0);
    if (sizes[2] > 0) comm().broadcast(input->actnum.data(), sizes[2], 0);
    Opm::Refinement::GridStateWriter::setRetainedCornerPointInput(*current_data_->front(), input);
}

void CpGrid::addLgrsUpdateLeafView(std::vector<Opm::Refinement::BlockRefinement> requests)
{
    Opm::Refinement::validateBlockRefinements(requests);

    if (lgr_backend_ == Opm::Refinement::Backend::Trilinear) {
        std::vector<std::array<int,3>> cellsPerDim, startIJK, endIJK;
        std::vector<std::string> names, parents;
        for (const auto& request : requests) {
            const bool graded = std::any_of(request.subdivision.begin(), request.subdivision.end(),
                                            [](const auto& sub) { return !sub.empty(); });
            if (graded || !request.minpvRemoved.empty() || request.pillarsFromBoxLayer) {
                OPM_THROW(std::invalid_argument,
                          "Refinement '" + request.name + "' is graded, removes cells by a block "
                          "MINPV or takes its pillars from the box layer; only the Conforming LGR "
                          "backend refines that.");
            }
            cellsPerDim.push_back(request.cellsPerDim);
            startIJK.push_back(request.startIJK);
            endIJK.push_back(request.endIJK);
            names.push_back(request.name);
            parents.push_back(request.parentGridName);
        }
        addLgrsUpdateLeafView(cellsPerDim, startIJK, endIJK, names, parents);
        return;
    }

    if (comm().size() > 1) {
        broadcastRetainedInput_();
    }
    auto builder = refinement_builder_;
    if (!builder) {
        if (const auto input = Opm::Refinement::GridStateWriter::retainedCornerPointInput(*current_data_->front())) {
            builder = std::make_shared<Opm::Refinement::ConformingBlockBuilder>(
                input->dims, input->coord, input->zcorn, input->actnum, input->edgeConformal);
        }
    }
    if (!builder) {
        OPM_THROW(std::logic_error, "The Conforming LGR backend has no refinement builder and no "
                  "retained corner-point input; select the backend before processEclipseFormat().");
    }
    const int preBuildMaxLevel = maxLevel();
    builder->build(*this, requests);

    for (std::size_t box = 0; box < requests.size(); ++box) {
        lgr_names_[requests[box].name] = preBuildMaxLevel + static_cast<int>(box) + 1;
    }
    if (global_id_set_ptr_) {
        auto& data = currentData();
        for (std::size_t gridIdx = preBuildMaxLevel + 1; gridIdx < data.size(); ++gridIdx) {
            global_id_set_ptr_->insertIdSet(*data[gridIdx]);
        }
    }
}

} // namespace Dune
