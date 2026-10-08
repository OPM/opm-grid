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
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <algorithm>
#include <array>
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
}

void CpGrid::setRefinementBuilder(std::shared_ptr<Opm::Refinement::Builder> builder)
{
    refinement_builder_ = std::move(builder);
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

    if (!refinement_builder_) {
        OPM_THROW(std::logic_error, "The Conforming LGR backend has no refinement builder.");
    }
    const int preBuildMaxLevel = maxLevel();
    refinement_builder_->build(*this, requests);

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
