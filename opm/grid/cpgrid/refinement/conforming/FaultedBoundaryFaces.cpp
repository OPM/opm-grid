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

#include <opm/grid/cpgrid/refinement/conforming/FaultedBoundaryFaces.hpp>

#include <opm/grid/cpgrid/refinement/conforming/GrdeclRefinement.hpp>
#include <opm/grid/cpgrid/refinement/RefinementRequest.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <algorithm>
#include <array>

namespace Opm
{
namespace Refinement
{

std::vector<BoundaryConnection>
faultedBoundaryConnections(const std::array<int,3>& parentDims,
                           const double* coord,
                           const double* zcorn,
                           const int* actnum,
                           const BlockRefinement& box,
                           int axis,
                           int side,
                           bool edgeConformal,
                           const std::array<std::array<int,2>,3>* across)
{
    std::vector<BoundaryConnection> result;

    // Mini region: the box, plus a one-cell coarse "shell" on the (axis, side)
    // side, plus a one-cell halo in the two perpendicular axes so a fault that
    // throws the boundary onto a coarse cell outside the box footprint is still
    // captured. Refining box+shell+halo and processing reproduces the box's
    // boundary faces against the coarse neighbour(s), split at the fault.
    const std::array<int,3>& boxStartIJK = box.startIJK;
    const std::array<int,3>& boxEndIJK   = box.endIJK;

    std::array<int,3> miniStart = boxStartIJK;
    std::array<int,3> miniEnd   = boxEndIJK;
    if (side < 0) {
        miniStart[axis] -= 1;
        if (miniStart[axis] < 0) {
            return result; // shell off-grid: the box boundary is the domain edge
        }
    }
    else {
        miniEnd[axis] += 1;
        if (miniEnd[axis] > parentDims[axis]) {
            return result;
        }
    }
    for (int p = 0; p < 3; ++p) {
        if (p == axis) continue;
        miniStart[p] = std::max(0, miniStart[p] - 1);
        miniEnd[p]   = std::min(parentDims[p], miniEnd[p] + 1);
        if (across) {
            miniStart[p] = std::min(miniStart[p], (*across)[p][0]);
            miniEnd[p]   = std::max(miniEnd[p], (*across)[p][1] + 1);
        }
    }

    // The mini block keeps the box's own subdivision on the box's parent cells
    // and refines the surrounding shell uniformly, as before -- the shell's
    // resolution is what lets a box<->box faulted interface identify the
    // neighbour's child.
    BlockRefinement req;
    req.name = "FAULTBND";
    req.cellsPerDim = box.cellsPerDim;
    req.startIJK = miniStart;
    req.endIJK = miniEnd;
    for (int d = 0; d < 3; ++d) {
        const auto boxSub = axisSubdivision(box, d);
        const auto boxPos = axisPositions(boxSub);

        AxisSubdivision mini;
        for (int p = miniStart[d]; p < miniEnd[d]; ++p) {
            const int parent = p - miniStart[d];
            const bool inside = (p >= boxStartIJK[d]) && (p < boxEndIJK[d]);
            const int first = inside ? boxPos.firstColumn[p - boxStartIJK[d]] : 0;
            const int count = inside ? boxPos.parentCount[first] : box.cellsPerDim[d];

            for (int sub = 0; sub < count; ++sub) {
                mini.parentOffset.push_back(parent);
                mini.fracLo.push_back(inside ? boxSub.fracLo[first + sub]
                                             : static_cast<double>(sub) / count);
                mini.fracHi.push_back(inside ? boxSub.fracHi[first + sub]
                                             : static_cast<double>(sub + 1) / count);
            }
        }
        req.subdivision[d] = std::move(mini);
    }

    const auto miniSub = std::array { axisSubdivision(req, 0),
                                      axisSubdivision(req, 1),
                                      axisSubdivision(req, 2) };
    const auto miniPos = std::array { axisPositions(miniSub[0]),
                                      axisPositions(miniSub[1]),
                                      axisPositions(miniSub[2]) };

    const RefinedBlockGrdecl refined = refineBlock(parentDims, coord, zcorn, actnum, req);
    const std::array<int,3> mdims = refined.dims;

    grdecl raw;
    raw.dims[0] = mdims[0]; raw.dims[1] = mdims[1]; raw.dims[2] = mdims[2];
    raw.coord = refined.coord.data();
    raw.zcorn = refined.zcorn.data();
    raw.actnum = refined.actnum.empty() ? nullptr : refined.actnum.data();

    struct processed_grid out;
    if (!process_grdecl(/*pinchActive=*/1, edgeConformal ? 1 : 0, 1e-6, &raw, nullptr, &out)) {
        return result;
    }

    const enum face_tag axisFaceTag = (axis == 0) ? I_FACE
                                    : (axis == 1) ? J_FACE
                                                  : K_FACE;
    // Offset (in refined cells) from the mini origin to the box origin, and the
    // box's refined extent. The box's boundary layer in `axis` is at frame 0
    // (side<0) or its last layer (side>0).
    std::array<int,3> boxOffset{};
    std::array<int,3> boxRefDims{};
    for (int d = 0; d < 3; ++d) {
        boxOffset[d]  = miniPos[d].firstColumn[boxStartIJK[d] - miniStart[d]];
        boxRefDims[d] = static_cast<int>(axisSubdivision(box, d).size());
    }
    const int boundaryLayer = (side < 0) ? 0 : boxRefDims[axis] - 1;

    const auto latticeOf = [&](int activeCell) {
        const int g = out.local_cell_index[activeCell];   // active -> logical (mini) index
        return std::array<int,3>{ g % mdims[0],
                                  (g / mdims[0]) % mdims[1],
                                  g / (mdims[0] * mdims[1]) };
    };
    const auto parentOfMini = [&](const std::array<int,3>& L) {
        return std::array<int,3>{ miniStart[0] + miniSub[0].parentOffset[L[0]],
                                  miniStart[1] + miniSub[1].parentOffset[L[1]],
                                  miniStart[2] + miniSub[2].parentOffset[L[2]] };
    };
    const auto inBox = [&](const std::array<int,3>& p) {
        return p[0] >= boxStartIJK[0] && p[0] < boxEndIJK[0]
            && p[1] >= boxStartIJK[1] && p[1] < boxEndIJK[1]
            && p[2] >= boxStartIJK[2] && p[2] < boxEndIJK[2];
    };

    for (unsigned face = 0; face < out.number_of_faces; ++face) {
        if (out.face_tag[face] != axisFaceTag) {
            continue;
        }
        // For an axis-face the normal points from neighbor 0 to neighbor 1, i.e.
        // increasing `axis`. The box's boundary face toward `side` therefore has
        // the box cell on the high-axis slot for side<0, the low-axis slot for
        // side>0; the other slot is the coarse neighbour (or -1 at the domain).
        const int boxSlot   = (side < 0) ? out.face_neighbors[2*face + 1] : out.face_neighbors[2*face];
        const int otherSlot = (side < 0) ? out.face_neighbors[2*face]     : out.face_neighbors[2*face + 1];
        if (boxSlot < 0) {
            continue; // no box cell on this face
        }
        const std::array<int,3> boxL = latticeOf(boxSlot);
        std::array<int,3> boxFrame{ boxL[0] - boxOffset[0],
                                    boxL[1] - boxOffset[1],
                                    boxL[2] - boxOffset[2] };
        if (boxFrame[0] < 0 || boxFrame[0] >= boxRefDims[0]
            || boxFrame[1] < 0 || boxFrame[1] >= boxRefDims[1]
            || boxFrame[2] < 0 || boxFrame[2] >= boxRefDims[2]) {
            continue; // boxSlot is not a box cell
        }
        if (boxFrame[axis] != boundaryLayer) {
            continue; // not the box's (axis, side) boundary layer
        }

        int coarseCart = -1; // domain boundary part of the box boundary face
        std::array<int,3> coarseSub{};
        if (otherSlot >= 0) {
            const std::array<int,3> oL = latticeOf(otherSlot);
            const std::array<int,3> op = parentOfMini(oL);
            if (inBox(op)) {
                continue; // an internal box face, not a boundary face
            }
            coarseCart = op[0] + parentDims[0]*op[1] + parentDims[0]*parentDims[1]*op[2];
            // Sub-position of the neighbour cell within its parent column (this
            // box's refined frame); identifies the neighbour's refined child when
            // it is itself refined.
            for (int d = 0; d < 3; ++d) {
                coarseSub[d] = miniPos[d].subIndex[oL[d]];
            }
        }

        BoundaryConnection bc;
        bc.boxCell = boxFrame;
        bc.coarseNeighborCart = coarseCart;
        bc.coarseNeighborSub = coarseSub;
        for (unsigned n = out.face_node_ptr[face]; n < out.face_node_ptr[face + 1]; ++n) {
            const int node = out.face_nodes[n];
            bc.faceNodes.push_back({ out.node_coordinates[3*node],
                                     out.node_coordinates[3*node + 1],
                                     out.node_coordinates[3*node + 2] });
        }
        result.push_back(std::move(bc));
    }

    free_processed_grid(&out);
    return result;
}

} // namespace Refinement
} // namespace Opm
