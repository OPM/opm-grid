/*
  Copyright 2026 SINTEF Digital.

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

#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{

bool boxesOverlap(const Opm::Refinement::BlockRefinement& a,
                  const Opm::Refinement::BlockRefinement& b)
{
    for (int c = 0; c < 3; ++c) {
        if (a.endIJK[c] <= b.startIJK[c] || b.endIJK[c] <= a.startIJK[c]) {
            return false;
        }
    }
    return true;
}

} // anonymous namespace

namespace Opm
{
namespace Refinement
{

AxisSubdivision axisSubdivision(const BlockRefinement& request, const int dim)
{
    if (! request.subdivision[dim].empty()) {
        return request.subdivision[dim];
    }

    const int nparent = request.endIJK[dim] - request.startIJK[dim];
    const int factor = request.cellsPerDim[dim];

    AxisSubdivision uniform;
    uniform.parentOffset.reserve(nparent * factor);
    uniform.fracLo.reserve(nparent * factor);
    uniform.fracHi.reserve(nparent * factor);

    for (int parent = 0; parent < nparent; ++parent) {
        for (int sub = 0; sub < factor; ++sub) {
            uniform.parentOffset.push_back(parent);
            uniform.fracLo.push_back(static_cast<double>(sub) / factor);
            uniform.fracHi.push_back(static_cast<double>(sub + 1) / factor);
        }
    }

    return uniform;
}

AxisPositions axisPositions(const AxisSubdivision& sub)
{
    const auto n = sub.size();

    AxisPositions pos;
    pos.subIndex.resize(n);
    pos.parentCount.resize(n);

    // parentOffset is non-decreasing, so each parent's columns are one run.
    std::size_t column = 0;
    while (column < n) {
        const int parent = sub.parentOffset[column];

        std::size_t end = column;
        while ((end < n) && (sub.parentOffset[end] == parent)) {
            ++end;
        }

        pos.firstColumn.resize(parent + 1, static_cast<int>(column));
        pos.firstColumn[parent] = static_cast<int>(column);

        const int count = static_cast<int>(end - column);
        for (auto c = column; c < end; ++c) {
            pos.subIndex[c] = static_cast<int>(c - column);
            pos.parentCount[c] = count;
        }

        column = end;
    }

    return pos;
}

std::array<int,3> refinedDims(const BlockRefinement& request)
{
    std::array<int,3> dims{};
    for (int d = 0; d < 3; ++d) {
        dims[d] = request.subdivision[d].empty()
            ? (request.endIJK[d] - request.startIJK[d]) * request.cellsPerDim[d]
            : static_cast<int>(request.subdivision[d].size());
    }
    return dims;
}

std::vector<BlockRefinement>
blockRefinements(const std::vector<std::array<int,3>>& cellsPerDim,
                 const std::vector<std::array<int,3>>& startIJK,
                 const std::vector<std::array<int,3>>& endIJK,
                 const std::vector<std::string>& names,
                 const std::vector<std::string>& parentGridNames)
{
    const std::size_t numBoxes = startIJK.size();
    if (cellsPerDim.size() != numBoxes || endIJK.size() != numBoxes || names.size() != numBoxes
        || (!parentGridNames.empty() && parentGridNames.size() != numBoxes)) {
        throw std::invalid_argument("Sizes of provided refinement vectors must match.");
    }
    std::vector<BlockRefinement> requests(numBoxes);
    for (std::size_t box = 0; box < numBoxes; ++box) {
        requests[box].name = names[box];
        if (!parentGridNames.empty()) {
            requests[box].parentGridName = parentGridNames[box];
        }
        requests[box].cellsPerDim = cellsPerDim[box];
        requests[box].startIJK = startIJK[box];
        requests[box].endIJK = endIJK[box];
    }
    return requests;
}

void validateBlockRefinements(const std::vector<BlockRefinement>& requests)
{
    for (const auto& req : requests) {
        // "GLOBAL" is the reserved name of level zero in CpGrid::lgr_names_
        // (and the parentGridName convention for top-level boxes); a request
        // with that name would overwrite the level-zero map entry.
        if (req.name == "GLOBAL") {
            throw std::invalid_argument("Refinement request name 'GLOBAL' is reserved "
                                        "for the level-zero grid; choose another name.");
        }
        for (int c = 0; c < 3; ++c) {
            if (req.startIJK[c] < 0 || req.startIJK[c] >= req.endIJK[c]) {
                throw std::invalid_argument("Invalid IJK box in refinement '" + req.name
                                            + "': end I/J/K must be larger than start I/J/K (and start non-negative).");
            }
            if (req.subdivision[c].empty() && req.cellsPerDim[c] < 1) {
                throw std::invalid_argument("Invalid subdivisions in refinement '" + req.name
                                            + "': NX/NY/NZ must be positive.");
            }
            if (!req.subdivision[c].empty()
                && (req.subdivision[c].fracLo.size() != req.subdivision[c].size()
                    || req.subdivision[c].fracHi.size() != req.subdivision[c].size()))
            {
                throw std::invalid_argument("Refinement '" + req.name
                                            + "' has a graded subdivision whose extent tables "
                                              "do not match its column count.");
            }
        }
    }
    // Pairwise disjointness of boxes refining the same parent grid. Note:
    // upstream never validated this and overlapping CARFIN boxes corrupt
    // the result silently (review Part III.3) - here it is a hard error.
    for (std::size_t i = 0; i < requests.size(); ++i) {
        for (std::size_t j = i + 1; j < requests.size(); ++j) {
            if (requests[i].parentGridName != requests[j].parentGridName) {
                continue;
            }
            if (boxesOverlap(requests[i], requests[j])) {
                throw std::invalid_argument("Refinement boxes '" + requests[i].name + "' and '"
                                            + requests[j].name + "' overlap. CARFIN boxes must be disjoint.");
            }
        }
    }
}

} // namespace Refinement
} // namespace Opm
