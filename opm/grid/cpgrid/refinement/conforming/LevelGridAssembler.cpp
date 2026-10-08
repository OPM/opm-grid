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

#include <opm/grid/cpgrid/refinement/conforming/LevelGridAssembler.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/refinement/conforming/GrdeclRefinement.hpp>
#include <opm/grid/cpgrid/refinement/GridStateWriter.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <array>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Opm
{
namespace Refinement
{

namespace {

/// The parent's pinch connections inside the box, carried to the children: each
/// joins the bottom children of the upper cell to the top children of the lower.
std::set<std::pair<int,int>>
boxPinchPairs(const Dune::cpgrid::CpGridData& parentGrid,
              const std::array<int,3>& parentDims,
              const BlockRefinement& request,
              const std::array<int,3>& refinedDims)
{
    auto& parent = const_cast<Dune::cpgrid::CpGridData&>(parentGrid);
    const auto& faceToCell = GridStateWriter::faceToCell(parent);
    const auto& faceTag = GridStateWriter::faceTag(parent);
    const auto& globalCell = parentGrid.globalCell();

    const std::array<AxisSubdivision,3> subs = { axisSubdivision(request, 0),
                                                axisSubdivision(request, 1),
                                                axisSubdivision(request, 2) };
    const auto refinedOf = [&subs](int dim, int offset) {
        std::vector<int> refined;
        for (std::size_t r = 0; r < subs[dim].size(); ++r) {
            if (subs[dim].parentOffset[r] == offset) {
                refined.push_back(static_cast<int>(r));
            }
        }
        return refined;
    };
    const auto ijkOf = [&parentDims](int cart) {
        return std::array<int,3>{ cart % parentDims[0],
                                  (cart / parentDims[0]) % parentDims[1],
                                  cart / (parentDims[0] * parentDims[1]) };
    };
    const auto inBox = [&request](const std::array<int,3>& ijk) {
        for (int d = 0; d < 3; ++d) {
            if (ijk[d] < request.startIJK[d] || ijk[d] >= request.endIJK[d]) {
                return false;
            }
        }
        return true;
    };

    std::set<std::pair<int,int>> pairs;
    const int numFaces = faceToCell.size();
    for (int face = 0; face < numFaces; ++face) {
        const auto row = faceToCell[Dune::cpgrid::EntityRep<1>(face, true)];
        // A distributed grid marks an off-rank neighbour with an out-of-range index.
        if ((row.size() != 2) || (faceTag.get(face) != K_FACE)
            || (row[0].index() >= static_cast<int>(globalCell.size()))
            || (row[1].index() >= static_cast<int>(globalCell.size()))) {
            continue;
        }
        auto upper = ijkOf(globalCell[row[0].index()]);
        auto lower = ijkOf(globalCell[row[1].index()]);
        if (upper[2] > lower[2]) {
            std::swap(upper, lower);
        }
        if ((upper[0] != lower[0]) || (upper[1] != lower[1]) || (lower[2] - upper[2] < 2)
            || !inBox(upper) || !inBox(lower)) {
            continue;
        }
        const int kUpper = refinedOf(2, upper[2] - request.startIJK[2]).back();
        const int kLower = refinedOf(2, lower[2] - request.startIJK[2]).front();
        for (const int ir : refinedOf(0, upper[0] - request.startIJK[0])) {
            for (const int jr : refinedOf(1, upper[1] - request.startIJK[1])) {
                const int column = ir + refinedDims[0] * jr;
                const int layer = refinedDims[0] * refinedDims[1];
                pairs.emplace(column + layer * kUpper, column + layer * kLower);
            }
        }
    }
    return pairs;
}

} // Anonymous namespace

std::shared_ptr<Dune::cpgrid::CpGridData>
assembleBlockLevelGrid(const Dune::cpgrid::CpGridData& parentGrid,
                       const std::array<int,3>& parentDims,
                       const double* coord,
                       const double* zcorn,
                       const int* actnum,
                       const BlockRefinement& request,
                       int levelIndex,
                       int parentLevel,
                       std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& levelStorage,
                       Dune::MPIHelper::MPICommunicator comm,
                       RefinedBlockGrdecl* outRefined)
{
    // Stages 2+3: resample the block geometry onto sub-pillars.
    const RefinedBlockGrdecl refined = refineBlock(parentDims, coord, zcorn, actnum, request);

    // Stage 4: process the refined description with the corner-point
    // preprocessor. Faults/pinch-outs inside the block are matched here.
    // Refined level grids are rank-local (rank-interior LGR model): use a
    // self-communicator so none of their operations are collective over the
    // distributed communicator (avoids deadlock when only the owning rank
    // refines). The distributed leaf keeps the real communicator.
    (void)comm;
    auto level = std::make_shared<Dune::cpgrid::CpGridData>(
        Dune::MPIHelper::getLocalCommunicator(), levelStorage);
    grdecl raw;
    raw.dims[0] = refined.dims[0];
    raw.dims[1] = refined.dims[1];
    raw.dims[2] = refined.dims[2];
    raw.coord = refined.coord.data();
    raw.zcorn = refined.zcorn.data();
    raw.actnum = refined.actnum.data();

    // Index 0 holds pinch connections, as for level zero.
    std::array<std::set<std::pair<int,int>>, 2> nnc;
    nnc[0] = boxPinchPairs(parentGrid, parentDims, request, refined.dims);
    level->processEclipseFormat(raw,
                                nullptr,
                                nnc,
                                /* remove_ij_boundary = */ false,
                                /* turn_normals = */ false,
                                /* pinchActive = */ !nnc[0].empty(),
                                /* tolerance_unique_points = */ 0.0,
                                /* edge_conformal = */ false);

    // Stage 6a: parent relations. The processed grid's global_cell_ holds
    // the refined-local Cartesian index (the CARFIN convention for level
    // grids), which encodes the parent and the position inside it.
    const std::array<AxisSubdivision,3> subs = { axisSubdivision(request, 0),
                                                axisSubdivision(request, 1),
                                                axisSubdivision(request, 2) };
    const std::array<AxisPositions,3> pos = { axisPositions(subs[0]),
                                              axisPositions(subs[1]),
                                              axisPositions(subs[2]) };
    const auto& [nx, ny, nz] = parentDims;

    // Invert the parent grid's cartesian -> compressed mapping. For a
    // top-level box the parent is level zero; for a nested box it is the LGR
    // level grid being refined, whose globalCell() holds the refined-local
    // Cartesian index (the same CARFIN convention), so the same inversion
    // works in the parent's local index space.
    std::vector<int> parentCompressed(static_cast<std::size_t>(nx)*ny*nz, -1);
    const auto& parentGlobalCell = parentGrid.globalCell();
    for (std::size_t c = 0; c < parentGlobalCell.size(); ++c) {
        parentCompressed[parentGlobalCell[c]] = static_cast<int>(c);
    }

    const int numRefinedCells = level->size(0);
    std::vector<std::array<int,2>> childToParent(numRefinedCells);
    std::vector<int> idxInParent(numRefinedCells);

    for (int cell = 0; cell < numRefinedCells; ++cell) {
        const int refinedCart = level->globalCell()[cell];
        const int ir = refinedCart % refined.dims[0];
        const int jr = (refinedCart / refined.dims[0]) % refined.dims[1];
        const int kr = refinedCart / (refined.dims[0]*refined.dims[1]);

        const int ci = request.startIJK[0] + subs[0].parentOffset[ir];
        const int cj = request.startIJK[1] + subs[1].parentOffset[jr];
        const int ck = request.startIJK[2] + subs[2].parentOffset[kr];
        const int parentCart = ci + nx*cj + nx*ny*ck;
        const int parentIdx = parentCompressed[parentCart];
        if (parentIdx < 0) {
            // Children inherit ACTNUM, so an active child of an inactive
            // parent indicates an internal inconsistency.
            throw std::logic_error("Refined cell in '" + request.name
                                   + "' has no active parent cell.");
        }
        // Position within the parent, whose extent varies from parent to parent
        // once the box is graded, so the strides are that parent's own counts.
        const int cx = pos[0].parentCount[ir];
        const int cy = pos[1].parentCount[jr];
        childToParent[cell] = {parentLevel, parentIdx};
        idxInParent[cell] = pos[0].subIndex[ir]
                          + pos[1].subIndex[jr]*cx
                          + pos[2].subIndex[kr]*cx*cy;
    }

    GridStateWriter::setLevel(*level, levelIndex);
    GridStateWriter::setCellsPerDim(*level, request.cellsPerDim);
    // Only a graded box stores its tables; a uniform level keeps to
    // cells_per_dim, which everything reading it already understands.
    if (std::ranges::any_of(request.subdivision, [](const auto& s) { return !s.empty(); })) {
        GridStateWriter::setSubdivision(*level, subs);
    }
    GridStateWriter::setParentRelations(*level, std::move(childToParent), std::move(idxInParent));

    // Expose the resampled description so a nested child can be refined from
    // this level's geometry (its parent) rather than the global grid.
    if (outRefined != nullptr) {
        *outRefined = refined;
    }

    return level;
}

std::shared_ptr<Dune::cpgrid::CpGridData>
assembleEmptyLevelGrid(const BlockRefinement& request,
                       int levelIndex,
                       std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& levelStorage,
                       Dune::MPIHelper::MPICommunicator comm)
{
    (void)comm;
    auto level = std::make_shared<Dune::cpgrid::CpGridData>(
        Dune::MPIHelper::getLocalCommunicator(), levelStorage);
    GridStateWriter::setLevel(*level, levelIndex);
    GridStateWriter::setCellsPerDim(*level, request.cellsPerDim);
    if (std::ranges::any_of(request.subdivision, [](const auto& s) { return !s.empty(); })) {
        GridStateWriter::setSubdivision(*level, { axisSubdivision(request, 0),
                                                  axisSubdivision(request, 1),
                                                  axisSubdivision(request, 2) });
    }
    GridStateWriter::setLogicalCartesianSize(*level, Opm::Refinement::refinedDims(request));
    GridStateWriter::setGlobalCell(*level, {});
    GridStateWriter::setIndexSet(*level, 0, 0);
    GridStateWriter::setParentRelations(*level, {}, {});
    return level;
}

} // namespace Refinement
} // namespace Opm
