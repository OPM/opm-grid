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

#include <opm/grid/cpgrid/refinement/conforming/LeafGridAssembler.hpp>

#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/cpgrid/EntityRep.hpp>
#include <opm/grid/cpgrid/Geometry.hpp>
#include <opm/grid/cpgrid/refinement/conforming/FaultedBoundaryFaces.hpp>
#include <opm/grid/cpgrid/refinement/conforming/RefinedIds.hpp>
#include <opm/grid/cpgrid/refinement/GridStateWriter.hpp>

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <cstdio>
#include <limits>

namespace
{

using Dune::cpgrid::CpGridData;
using Dune::cpgrid::EntityRep;
using Opm::Refinement::BlockRefinement;
using Opm::Refinement::GridStateWriter;

// Numerical aquifer cells are never refined (refused in a box): each is a
// level-zero cell of the leaf.
std::vector<int> leafAquiferCells(const CpGridData& level0,
                                  const std::vector<std::array<int,2>>& leafToLevel)
{
    const auto& aquifer0 = level0.sortedNumAquiferCells();
    std::vector<int> leaf;
    if (aquifer0.empty()) {
        return leaf;
    }
    for (int cell = 0; cell < static_cast<int>(leafToLevel.size()); ++cell) {
        if ((leafToLevel[cell][0] == 0) &&
            std::binary_search(aquifer0.begin(), aquifer0.end(), leafToLevel[cell][1])) {
            leaf.push_back(cell);
        }
    }
    return leaf;
}

struct SourceRef
{
    int grid;  // 0 = level zero, b+1 = level grid of box b
    int index; // cell or face index in that grid
};

// processEclipseFormat reverses face nodes for left-handed or upward-ZCORN input, which
// swaps the halves of cell_to_point along J. Returns the XOR mask from logical corner
// (bit0 = I, bit1 = J, bit2 = K) to cell_to_point slot, detected from logical neighbours.
int cornerSlotMask(CpGridData& grid)
{
    const auto& c2p = GridStateWriter::cellToPoint(grid);
    const auto& gc = grid.globalCell();
    const auto dims = grid.logicalCartesianSize();
    std::map<int,int> cellOfCart;
    for (int c = 0; c < static_cast<int>(gc.size()); ++c) {
        cellOfCart.emplace(gc[c], c);
    }
    int mask = 0;
    for (int axis = 0; axis < 2; ++axis) {
        const int bit = 1 << axis;
        const int stride = axis == 0 ? 1 : dims[0];
        int agree = 0;
        int flipped = 0;
        for (const auto& [cart, cell] : cellOfCart) {
            const int ijk = axis == 0 ? cart % dims[0] : (cart / dims[0]) % dims[1];
            if (ijk + 1 >= dims[axis]) {
                continue;
            }
            const auto nb = cellOfCart.find(cart + stride);
            if (nb == cellOfCart.end()) {
                continue;
            }
            const auto& nbPoints = c2p[nb->second];
            const auto shared = [&](int side) {
                for (int corner = 0; corner < 8; ++corner) {
                    if (((corner & bit) != 0) == (side == 1) &&
                        std::find(nbPoints.begin(), nbPoints.end(), c2p[cell][corner]) == nbPoints.end()) {
                        return false;
                    }
                }
                return true;
            };
            const bool high = shared(1);
            const bool low = shared(0);
            agree += high && !low;
            flipped += low && !high;
            if (agree + flipped >= 16) {
                break;
            }
        }
        if (flipped > agree) {
            mask |= bit;
        }
    }
    return mask;
}

#if HAVE_MPI
// On a distributed leaf, ids through the local id set follow local numbering and
// differ between ranks, so the leaf gets a mapping. Cells and points come from
// assignRefinedIds; coarse faces keep level zero's ids and new faces, which have
// no meaning across ranks in CpGrid, get per-rank ranges above every id in use.
void setDistributedLeafIds(CpGridData& leaf, const CpGridData& level0,
                           std::vector<int> cellIds,
                           const std::vector<SourceRef>& leafFaces,
                           int numLeafFaces, std::vector<int> pointIds, int idFloor,
                           const Dune::Communication<Dune::MPIHelper::MPICommunicator>& cc)
{
    const auto& faceIds0 = level0.globalIdSet().template getMapping<1>();
    int maxUsed = idFloor - 1;
    for (const int id : faceIds0) {
        maxUsed = std::max(maxUsed, id);
    }
    maxUsed = cc.max(maxUsed);

    std::vector<int> faceIds(numLeafFaces, -1);
    int numNew = 0;
    for (int f = 0; f < numLeafFaces; ++f) {
        const bool coarse = f < static_cast<int>(leafFaces.size()) && leafFaces[f].grid == 0;
        if (coarse) {
            faceIds[f] = faceIds0[leafFaces[f].index];
        } else {
            ++numNew;
        }
    }
    std::vector<int> counts(cc.size(), 0);
    cc.allgather(&numNew, 1, counts.data());
    int next = maxUsed + 1;
    for (int r = 0; r < cc.rank(); ++r) {
        next += counts[r];
    }
    for (auto& id : faceIds) {
        if (id < 0) {
            id = next++;
        }
    }
    GridStateWriter::setGlobalIdMapping(leaf, std::move(cellIds), std::move(faceIds),
                                        std::move(pointIds));
}
#endif

int axisOf(face_tag tag)
{
    switch (tag) {
    case I_FACE: return 0;
    case J_FACE: return 1;
    case K_FACE: return 2;
    default:
        throw std::logic_error("Unexpected face tag in refined level grid.");
    }
}

face_tag faceTagOf(int axis)
{
    return axis == 0 ? I_FACE : axis == 1 ? J_FACE : K_FACE;
}

// In a distributed grid, face_to_cell uses this sentinel for a neighbor cell
// that lives on another rank (not present locally).
constexpr int kRemoteCell = std::numeric_limits<int>::max();

// outsideNeighborOf returns this when the parent's boundary on a side is
// crossed by a fault (split or partial face); the boundary is then rebuilt
// from the corner-point processor instead of the simple mosaic.
constexpr int kFaultedSide = -2;

// Area-weighted normal (Newell), centroid and area of a planar-ish polygon.
struct PolyGeom { Dune::FieldVector<double,3> center; Dune::FieldVector<double,3> normal; double area; };
PolyGeom polygonGeometry(const std::vector<std::array<double,3>>& pts)
{
    PolyGeom g;
    g.center = 0.0;
    g.normal = 0.0;
    const int n = static_cast<int>(pts.size());
    for (int i = 0; i < n; ++i) {
        const auto& a = pts[i];
        const auto& b = pts[(i + 1) % n];
        g.normal[0] += (a[1] - b[1]) * (a[2] + b[2]);
        g.normal[1] += (a[2] - b[2]) * (a[0] + b[0]);
        g.normal[2] += (a[0] - b[0]) * (a[1] + b[1]);
        for (int d = 0; d < 3; ++d) g.center[d] += a[d];
    }
    for (int d = 0; d < 3; ++d) g.center[d] /= (n > 0 ? n : 1);
    g.area = 0.5 * g.normal.two_norm();
    if (g.area > 0.0) g.normal /= g.normal.two_norm();  // unit
    return g;
}

} // anonymous namespace

namespace Opm
{
namespace Refinement
{

std::shared_ptr<CpGridData>
assembleLeafGrid(std::vector<std::shared_ptr<CpGridData>>& storage,
                 const std::vector<BlockRefinement>& requests,
                 const std::array<int,3>& parentDims,
                 const double* coord,
                 const double* zcorn,
                 const int* actnum,
                 Dune::MPIHelper::MPICommunicator comm)
{
    CpGridData& level0 = *storage[0];
    const int numBoxes = static_cast<int>(requests.size());
    const auto& dims0 = level0.logicalCartesianSize();
    const int numCells0 = level0.size(0);
    const int numCorners0 = level0.size(3);

    // ----- Basic level-zero lookups ------------------------------------
    std::vector<int> compressed0(static_cast<std::size_t>(dims0[0])*dims0[1]*dims0[2], -1);
    for (int c = 0; c < numCells0; ++c) {
        compressed0[level0.globalCell()[c]] = c;
    }
    const auto boxOfIJK = [&](int i, int j, int k) {
        for (int b = 0; b < numBoxes; ++b) {
            const auto& req = requests[b];
            if (i >= req.startIJK[0] && i < req.endIJK[0]
                && j >= req.startIJK[1] && j < req.endIJK[1]
                && k >= req.startIJK[2] && k < req.endIJK[2]) {
                return b;
            }
        }
        return -1;
    };
    // Box membership of each level-zero cell (-1 if unrefined).
    std::vector<int> boxOfCell(numCells0, -1);
    for (int c = 0; c < numCells0; ++c) {
        const int cart = level0.globalCell()[c];
        const int i = cart % dims0[0];
        const int j = (cart / dims0[0]) % dims0[1];
        const int k = cart / (dims0[0]*dims0[1]);
        boxOfCell[c] = boxOfIJK(i, j, k);
    }

    auto& cellToPoint0 = GridStateWriter::cellToPoint(level0);
    auto& cellToFace0 = GridStateWriter::cellToFace(level0);
    auto& faceToCell0 = GridStateWriter::faceToCell(level0);
    auto& faceToPoint0 = GridStateWriter::faceToPoint(level0);
    auto& faceTag0 = GridStateWriter::faceTag(level0);
    auto& faceNormals0 = GridStateWriter::faceNormals(level0);
    const int numFaces0 = faceToCell0.size();

    // ----- Per-box data -------------------------------------------------
    struct BoxData
    {
        CpGridData* level;
        std::array<int,3> factors;
        std::array<AxisSubdivision,3> subs;
        std::array<AxisPositions,3> pos;
        bool graded{false};
        std::array<int,3> refinedDims;
        // children of each level-zero parent, ordered by idxInParent
        std::map<int, std::vector<int>> childrenOfParent;
        // level corner -> level-zero corner (or -1)
        std::vector<int> cornerEquiv;
        // level corner -> leaf pool index
        std::vector<int> cornerToLeaf;
        // level cell/face -> leaf index
        std::vector<int> cellToLeaf;
        std::vector<int> faceToLeaf;
    };
    std::vector<BoxData> boxes(numBoxes);
    const int mask0 = cornerSlotMask(level0);

    for (int b = 0; b < numBoxes; ++b) {
        BoxData& box = boxes[b];
        box.level = storage[b + 1].get();
        box.factors = requests[b].cellsPerDim;
        for (int d = 0; d < 3; ++d) {
            box.subs[d] = axisSubdivision(requests[b], d);
            box.pos[d] = axisPositions(box.subs[d]);
            box.graded = box.graded || !requests[b].subdivision[d].empty();
        }
        box.refinedDims = box.level->logicalCartesianSize();

        const auto& childToParent = GridStateWriter::childToParent(*box.level);
        const auto& idxInParent = GridStateWriter::idxInParent(*box.level);
        std::map<int, std::vector<std::pair<int,int>>> tmp; // parent -> (idxInParent, child)
        for (int cell = 0; cell < box.level->size(0); ++cell) {
            tmp[childToParent[cell][1]].emplace_back(idxInParent[cell], cell);
        }
        for (auto& [parent, kids] : tmp) {
            std::sort(kids.begin(), kids.end());
            auto& target = box.childrenOfParent[parent];
            target.reserve(kids.size());
            for (const auto& [pos, cell] : kids) {
                target.push_back(cell);
            }
        }

        // Corner identification with level zero, through the parent cells
        // (per-cell identification stays consistent across faults inside
        // the block, where lattice positions are not unique).
        // A refined corner coincides with a parent corner when it sits on the
        // parent's boundary in all three directions. Taken from the cell's
        // refined lattice position rather than from idxInParent, whose strides
        // are that parent's own counts once the box is graded.
        box.cornerEquiv.assign(box.level->size(3), -1);
        auto& cellToPointL = GridStateWriter::cellToPoint(*box.level);
        const auto& rd = box.refinedDims;
        const int maskL = cornerSlotMask(*box.level);
        for (int cell = 0; cell < box.level->size(0); ++cell) {
            const int parent = childToParent[cell][1];
            const int refinedCart = box.level->globalCell()[cell];
            const std::array<int,3> lattice = { refinedCart % rd[0],
                                                (refinedCart / rd[0]) % rd[1],
                                                refinedCart / (rd[0]*rd[1]) };
            for (int corner = 0; corner < 8; ++corner) {
                const std::array<int,3> d = { corner & 1, (corner >> 1) & 1, (corner >> 2) & 1 };
                std::array<int,3> parentCornerIJK{};
                bool onParentCorner = true;
                for (int c = 0; c < 3 && onParentCorner; ++c) {
                    const int edge = box.pos[c].subIndex[lattice[c]] + d[c];
                    if (edge == 0) {
                        parentCornerIJK[c] = 0;
                    }
                    else if (edge == box.pos[c].parentCount[lattice[c]]) {
                        parentCornerIJK[c] = 1;
                    }
                    else {
                        onParentCorner = false;
                    }
                }
                if (onParentCorner) {
                    const int parentCorner = cellToPoint0[parent][(parentCornerIJK[0]
                                                                  + 2*parentCornerIJK[1]
                                                                  + 4*parentCornerIJK[2]) ^ mask0];
                    box.cornerEquiv[cellToPointL[cell][corner ^ maskL]] = parentCorner;
                }
            }
        }
    }

    // ----- Corner pool (D1) ----------------------------------------------
    // Pool = all level-zero corners (indices preserved) + new refined corners.
    auto leaf = std::make_shared<CpGridData>(comm, storage);
    auto& leafGeometry = GridStateWriter::geometry(*leaf);
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& leafCorners =
        leafGeometry.geomVector(std::integral_constant<int,3>());
    const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& corners0 =
        GridStateWriter::geometry(level0).geomVector(std::integral_constant<int,3>());
    leafCorners.assign(corners0.begin(), corners0.end());

    std::vector<std::array<int,2>> leafCornerHistory(numCorners0);
    for (int corner = 0; corner < numCorners0; ++corner) {
        leafCornerHistory[corner] = {0, corner};
    }

    // Refined corners shared between touching boxes (edge/corner-sharing)
    // coincide bitwise: both boxes resample the same parent description with
    // equal subdivisions, so the arithmetic is identical. Exact-coordinate
    // matching is therefore correct (not a fragile floating-point heuristic)
    // and merges those corners into one pool entry. Disjoint boxes only
    // coincide on shared boundaries, so a global map over refined corners is
    // safe.
    std::map<std::array<double,3>, int> refinedCornerPool;
    const auto coordKey = [](const Dune::cpgrid::Geometry<0,3>& corner) {
        const auto& c = corner.center();
        return std::array<double,3>{ c[0], c[1], c[2] };
    };

    for (int b = 0; b < numBoxes; ++b) {
        BoxData& box = boxes[b];
        const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& cornersL =
            GridStateWriter::geometry(*box.level).geomVector(std::integral_constant<int,3>());
        box.cornerToLeaf.assign(box.level->size(3), -1);
        std::vector<std::array<int,2>> levelCornerHistory(box.level->size(3), {-1, -1});
        for (int corner = 0; corner < box.level->size(3); ++corner) {
            if (box.cornerEquiv[corner] >= 0) {
                box.cornerToLeaf[corner] = box.cornerEquiv[corner];
                levelCornerHistory[corner] = {0, box.cornerEquiv[corner]};
            }
            else {
                const auto key = coordKey(cornersL[corner]);
                const auto existing = refinedCornerPool.find(key);
                if (existing != refinedCornerPool.end()) {
                    // Shared with an already-processed box: reuse, and keep
                    // the corner's birth identity from that earlier box.
                    box.cornerToLeaf[corner] = existing->second;
                    levelCornerHistory[corner] = leafCornerHistory[existing->second];
                }
                else {
                    const int leafCorner = static_cast<int>(leafCorners.size());
                    box.cornerToLeaf[corner] = leafCorner;
                    leafCorners.push_back(cornersL[corner]);
                    leafCornerHistory.push_back({b + 1, corner});
                    refinedCornerPool.emplace(key, leafCorner);
                }
            }
        }
        GridStateWriter::setCornerHistory(*box.level, std::move(levelCornerHistory));
    }
    // May grow below if faulted box boundaries introduce new split vertices.
    int numLeafCorners = static_cast<int>(leafCorners.size());

    // ----- Leaf cells -----------------------------------------------------
    // Level-zero order, parents replaced by their children (idxInParent order).
    std::vector<SourceRef> leafCells;
    std::vector<int> leafIdxOfCell0(numCells0, -1);
    std::vector<std::vector<int>> leafIdxOfLevelCell(numBoxes);
    for (int b = 0; b < numBoxes; ++b) {
        leafIdxOfLevelCell[b].assign(boxes[b].level->size(0), -1);
    }
    std::vector<std::tuple<int, std::vector<int>>> parentToChildren(numCells0,
                                                                    std::make_tuple(-1, std::vector<int>{}));
    for (int c = 0; c < numCells0; ++c) {
        const int b = boxOfCell[c];
        if (b < 0) {
            leafIdxOfCell0[c] = static_cast<int>(leafCells.size());
            leafCells.push_back({0, c});
        }
        else {
            auto childrenIt = boxes[b].childrenOfParent.find(c);
            if (childrenIt == boxes[b].childrenOfParent.end()) {
                continue; // active parent without active children cannot happen; inactive handled by absence
            }
            for (const int child : childrenIt->second) {
                leafIdxOfLevelCell[b][child] = static_cast<int>(leafCells.size());
                leafCells.push_back({b + 1, child});
            }
            parentToChildren[c] = std::make_tuple(b + 1, childrenIt->second);
        }
    }
    const int numLeafCells = static_cast<int>(leafCells.size());

    // ----- Identify replaced level-zero faces and mosaic pairings ---------
    // A level-zero face is dropped if any adjacent cell is a refined parent.
    // For parent faces toward an unrefined neighbor, the refined boundary
    // faces of the box take over, paired with that neighbor.
    std::vector<bool> faceDropped(numFaces0, false);
    for (int face = 0; face < numFaces0; ++face) {
        const auto row = faceToCell0[EntityRep<1>(face, true)];
        for (int e = 0; e < row.size(); ++e) {
            const int cell = row[e].index();
            // A remote (off-rank) neighbor is never a local box cell.
            if (cell != kRemoteCell && boxOfCell[cell] >= 0) {
                faceDropped[face] = true;
                break;
            }
        }
    }

    // Conformity bookkeeping: (parent, axis, side) -> verified once.
    // The parent's boundary in that direction must consist of exactly one
    // level-zero face that is the parent's *full* logical face (its four
    // lattice corners); otherwise the block boundary is faulted/degenerate
    // and the refined mosaic cannot replace it (unsupported, throw).
    Dune::cpgrid::EntityVariableBase<enum face_tag>& tags0 = faceTag0;
    std::map<std::tuple<int,int,int>, int> verifiedNeighbor;
    const auto outsideNeighborOf = [&](int b, int parent, int axis, int side) {
        (void)b;
        const auto key = std::make_tuple(parent, axis, side);
        auto it = verifiedNeighbor.find(key);
        if (it != verifiedNeighbor.end()) {
            return it->second;
        }

        // All faces of the parent in the (axis, side) direction.
        std::vector<int> directionFaces;
        const auto row = cellToFace0[EntityRep<0>(parent, true)];
        for (int e = 0; e < row.size(); ++e) {
            const int face = row[e].index();
            const enum face_tag tag = tags0.get(face);
            if (tag != NNC_FACE && axisOf(tag) == axis
                && row[e].orientation() == (side > 0)) {
                directionFaces.push_back(face);
            }
        }

        // Below a pinch-out the parent keeps its own top face, now without a
        // neighbour, and gains the face that connects it to the cell above.
        if ((axis == 2) && (directionFaces.size() == 2)) {
            int open = 0;
            int pinched = -1;
            const auto ijk = [&dims0](int cart) {
                return std::array<int,3>{ cart % dims0[0], (cart / dims0[0]) % dims0[1],
                                          cart / (dims0[0] * dims0[1]) };
            };
            const auto p = ijk(level0.globalCell()[parent]);
            for (const int face : directionFaces) {
                const auto cells = faceToCell0[EntityRep<1>(face, true)];
                if (cells.size() == 1) {
                    ++open;
                    continue;
                }
                for (int q = 0; q < cells.size(); ++q) {
                    const int cell = cells[q].index();
                    if ((cell == parent) || (cell == kRemoteCell)) {
                        continue;
                    }
                    const auto c = ijk(level0.globalCell()[cell]);
                    if ((c[0] == p[0]) && (c[1] == p[1]) && (std::abs(c[2] - p[2]) > 1)) {
                        pinched = cell;
                    }
                }
            }
            if ((open == 1) && (pinched >= 0)) {
                verifiedNeighbor[key] = pinched;
                return pinched;
            }
        }

        int neighbor = -1;
        if (directionFaces.size() > 1) {
            neighbor = kFaultedSide;            // fault-split face
        }
        else if (directionFaces.size() == 1) {
            const int face = directionFaces[0];
            // The face must span the parent's full logical face: its corner
            // set must equal the parent's four lattice corners on that side.
            const int fixedBit = (side > 0) ? 1 : 0;
            std::set<int> expected;
            for (int corner = 0; corner < 8; ++corner) {
                const int bit = (axis == 0) ? (corner & 1)
                              : (axis == 1) ? ((corner >> 1) & 1)
                                            : ((corner >> 2) & 1);
                if (bit == fixedBit) {
                    expected.insert(cellToPoint0[parent][corner ^ mask0]);
                }
            }
            std::set<int> actual;
            auto facePoints = faceToPoint0[face];
            for (auto pt = facePoints.begin(); pt != facePoints.end(); ++pt) {
                actual.insert(*pt);
            }
            if (actual != expected) {
                neighbor = kFaultedSide;        // partial (faulted/degenerate) face
            }
            else {
                const auto cells = faceToCell0[EntityRep<1>(face, true)];
                for (int q = 0; q < cells.size(); ++q) {
                    const int cell = cells[q].index();
                    if (cell != parent && cell != kRemoteCell) {
                        neighbor = cell;
                    }
                }
            }
        }
        // No direction faces: domain boundary, inactive, or off-rank neighbor.

        verifiedNeighbor[key] = neighbor;
        return neighbor;
    };

    // ----- Leaf faces ------------------------------------------------------
    std::vector<SourceRef> leafFaces;
    // outside coarse leaf cell for mosaic faces, aligned with leafFaces (-1 otherwise)
    std::vector<int> mosaicOutside;
    std::vector<int> leafIdxOfFace0(numFaces0, -1);
    for (int face = 0; face < numFaces0; ++face) {
        if (!faceDropped[face]) {
            leafIdxOfFace0[face] = static_cast<int>(leafFaces.size());
            leafFaces.push_back({0, face});
            mosaicOutside.push_back(-1);
        }
    }
    // Face-sharing LGRs: a boundary face of box A toward another refined
    // box B pairs with B's matching boundary face. Both have the same
    // merged corner set (corners already merged), so we key by it. The box
    // seen first owns the leaf face; the partner maps onto it and patches
    // in its cell as the second side (the same mechanism as a coarse
    // mosaic neighbor, only the outside cell is refined).
    std::map<std::vector<int>, int> sharedBoundaryFaceLeaf; // corner set -> owner leaf face

    // Per box: refined (box-global) Cartesian index -> level cell index. Used by
    // the compatible sub-face mosaic below to find, for a finer box's boundary
    // sub-face, the coarser box's child cell it abuts (LGR_GAPS A2).
    std::vector<std::map<int,int>> refinedCartToCell(numBoxes);
    for (int b = 0; b < numBoxes; ++b) {
        for (int c = 0; c < boxes[b].level->size(0); ++c) {
            refinedCartToCell[b][boxes[b].level->globalCell()[c]] = c;
        }
    }

    // Pre-pass: which (box, axis, side) boundaries are faulted (some parent on
    // that side has a split/partial level-zero face). Decided before the face
    // loop so that ALL box boundary faces on a faulted side are suppressed
    // consistently -- including parents whose throw clears the neighbour
    // entirely (those alone look like a plain domain boundary). The whole side
    // is then rebuilt from the corner-point processor below.
    std::set<std::tuple<int,int,int>> faultedSides;
    // Parent IJK range of the cells really across each box side, from the level-zero
    // faces: a throw taller than the box puts them outside the box's own layers.
    std::map<std::tuple<int,int,int>, std::array<std::array<int,2>,3>> acrossRange;
    const auto extendAcross = [&](int b, int parent, int axis, int side) {
        const auto row = cellToFace0[EntityRep<0>(parent, true)];
        for (int e = 0; e < row.size(); ++e) {
            const int face = row[e].index();
            const enum face_tag tag = tags0.get(face);
            if (tag == NNC_FACE || axisOf(tag) != axis || row[e].orientation() != (side > 0)) {
                continue;
            }
            const auto cells = faceToCell0[EntityRep<1>(face, true)];
            for (int q = 0; q < cells.size(); ++q) {
                const int cell = cells[q].index();
                if (cell == parent || cell == kRemoteCell) {
                    continue;
                }
                const int cart = level0.globalCell()[cell];
                const std::array<int,3> nijk = { cart % dims0[0],
                                                 (cart / dims0[0]) % dims0[1],
                                                 cart / (dims0[0]*dims0[1]) };
                const auto key = std::make_tuple(b, axis, side);
                auto [it, fresh] = acrossRange.try_emplace(key);
                for (int d = 0; d < 3; ++d) {
                    it->second[d][0] = fresh ? nijk[d] : std::min(it->second[d][0], nijk[d]);
                    it->second[d][1] = fresh ? nijk[d] : std::max(it->second[d][1], nijk[d]);
                }
            }
        }
    };
    for (int c = 0; c < numCells0; ++c) {
        const int b = boxOfCell[c];
        if (b < 0) {
            continue;
        }
        const int cart = level0.globalCell()[c];
        const std::array<int,3> ijk = { cart % dims0[0],
                                        (cart / dims0[0]) % dims0[1],
                                        cart / (dims0[0]*dims0[1]) };
        for (int axis = 0; axis < 3; ++axis) {
            for (int side = -1; side <= 1; side += 2) {
                const bool onBoundary = (side < 0)
                    ? (ijk[axis] == requests[b].startIJK[axis])
                    : (ijk[axis] == requests[b].endIJK[axis] - 1);
                if (!onBoundary) {
                    continue;
                }
                if (outsideNeighborOf(b, c, axis, side) == kFaultedSide) {
                    faultedSides.emplace(b, axis, side);
                }
                extendAcross(b, c, axis, side);
            }
        }
    }

    for (int b = 0; b < numBoxes; ++b) {
        BoxData& box = boxes[b];
        auto& faceToCellL = GridStateWriter::faceToCell(*box.level);
        auto& faceTagL = GridStateWriter::faceTag(*box.level);
        auto& faceToPointL = GridStateWriter::faceToPoint(*box.level);
        Dune::cpgrid::EntityVariableBase<enum face_tag>& tagsL = faceTagL;
        const auto& childToParent = GridStateWriter::childToParent(*box.level);
        box.faceToLeaf.assign(faceToCellL.size(), -1);

        const auto mergedCornerSet = [&](int face) {
            std::vector<int> set;
            auto pts = faceToPointL[face];
            for (auto p = pts.begin(); p != pts.end(); ++p) {
                set.push_back(box.cornerToLeaf[*p]);
            }
            std::sort(set.begin(), set.end());
            return set;
        };

        for (int face = 0; face < faceToCellL.size(); ++face) {
            const auto row = faceToCellL[EntityRep<1>(face, true)];
            int outside = -1;
            if (row.size() == 1) {
                // Boundary face of the refined block: internal hole, domain
                // boundary, block boundary toward a coarse neighbor, or
                // toward another refined box (face-sharing).
                const int cell = row[0].index();
                const bool normalOut = row[0].orientation();
                const int axis = axisOf(tagsL.get(face));
                const int side = normalOut ? 1 : -1;
                const int cart = box.level->globalCell()[cell];
                std::array<int,3> lattice = { cart % box.refinedDims[0],
                                              (cart / box.refinedDims[0]) % box.refinedDims[1],
                                              cart / (box.refinedDims[0]*box.refinedDims[1]) };
                lattice[axis] += side;
                if (lattice[axis] < 0 || lattice[axis] >= box.refinedDims[axis]) {
                    // Faulted side: suppress this box face; the whole side is
                    // rebuilt below from the corner-point processor.
                    if (faultedSides.count(std::make_tuple(b, axis, side))) {
                        continue;
                    }
                    // Block boundary: find the neighbor parent (or none).
                    const int parent = childToParent[cell][1];
                    const int neighbor = outsideNeighborOf(b, parent, axis, side);
                    if (neighbor >= 0 && boxOfCell[neighbor] < 0) {
                        // Toward a coarse (unrefined) neighbor: mosaic.
                        outside = leafIdxOfCell0[neighbor];
                    }
                    else if (neighbor >= 0) {
                        // Toward another refined box. The two boxes share this
                        // 2-D face; its in-face directions are the two axes
                        // other than `axis` (the perpendicular touch direction).
                        const int nbBox = boxOfCell[neighbor];
                        if (box.graded || boxes[nbBox].graded) {
                            // The conformity rule below compares one subdivision
                            // factor per in-face direction; a graded box has a
                            // column-by-column subdivision instead.
                            throw std::invalid_argument(
                                "Refinement boxes touch and at least one is graded "
                                "(N*FIN/H*FIN). Touching graded boxes are not supported.");
                        }
                        const auto& facA = box.factors;
                        const auto& facB = boxes[nbBox].factors;
                        const int u = (axis == 0) ? 1 : 0;
                        const int v = (axis == 2) ? 1 : 2;
                        if (facA[u] == facB[u] && facA[v] == facB[v]) {
                            // Matching in-face subdivisions: 1:1 face pairing by
                            // corner set (both sides have identical sub-faces).
                            const auto key = mergedCornerSet(face);
                            const auto it = sharedBoundaryFaceLeaf.find(key);
                            const int thisCellLeaf = leafIdxOfLevelCell[b][cell];
                            if (it != sharedBoundaryFaceLeaf.end()) {
                                // Partner: map onto the owner's face, supply our
                                // cell as its second side, do not emit a face.
                                box.faceToLeaf[face] = it->second;
                                mosaicOutside[it->second] = thisCellLeaf;
                                continue;
                            }
                            // Owner: emit now, partner patches its cell later.
                            sharedBoundaryFaceLeaf.emplace(key, static_cast<int>(leafFaces.size()));
                        }
                        else {
                            // Compatible different subdivisions (guard-enforced:
                            // one box is uniformly finer). The finer side emits
                            // every sub-face and references the coarser side's
                            // covering cell as the outside; the coarser side
                            // suppresses its interface faces, so its interface
                            // cells become >6-face hexes tiled by the finer
                            // sub-faces (LGR_GAPS A2 sub-face mosaic).
                            const bool thisFiner = (facA[u] >= facB[u] && facA[v] >= facB[v]);
                            if (!thisFiner) {
                                continue;       // coarser side: the finer references us
                            }
                            // Sub-position of this finer cell within its parent;
                            // map it to the coarser box's covering child cell.
                            std::array<int,3> subB{};
                            subB[axis] = (side > 0) ? 0 : (facB[axis] - 1);
                            for (const int c : {u, v}) {
                                const int subA = lattice[c] % facA[c];
                                subB[c] = subA * facB[c] / facA[c];
                            }
                            const int ncart = level0.globalCell()[neighbor];
                            const std::array<int,3> nijk = { ncart % dims0[0],
                                                             (ncart / dims0[0]) % dims0[1],
                                                             ncart / (dims0[0]*dims0[1]) };
                            const auto& startB = requests[nbBox].startIJK;
                            const auto& rdB = boxes[nbBox].refinedDims;
                            std::array<int,3> latB{};
                            for (int c = 0; c < 3; ++c) {
                                latB[c] = (nijk[c] - startB[c]) * facB[c] + subB[c];
                            }
                            const int cartB = latB[0] + rdB[0]*latB[1] + rdB[0]*rdB[1]*latB[2];
                            const auto cit = refinedCartToCell[nbBox].find(cartB);
                            if (cit != refinedCartToCell[nbBox].end()) {
                                outside = leafIdxOfLevelCell[nbBox][cit->second];
                            }
                        }
                    }
                }
                // else: internal hole boundary -> stays a boundary face.
            }
            box.faceToLeaf[face] = static_cast<int>(leafFaces.size());
            leafFaces.push_back({b + 1, face});
            mosaicOutside.push_back(outside);
        }
    }

    // Every face-sharing owner face must have been claimed by its partner box.
    // An unclaimed one means the shared-face corners did not merge exactly in
    // the corner pool, and the box<->box interface would silently seal.
    for (const auto& entry : sharedBoundaryFaceLeaf) {
        const int leafFace = entry.second;
        if (mosaicOutside[leafFace] < 0) {
            const int ownerBox = leafFaces[leafFace].grid - 1;
            throw std::logic_error("Internal error: refined box '"
                                   + requests[ownerBox].name
                                   + "' has an interface face toward a face-sharing "
                                     "refined box with no matching partner face. "
                                     "The shared boundary corners failed to merge "
                                     "exactly; the interface would be sealed.");
        }
    }
    const int numSourceFaces = static_cast<int>(leafFaces.size());

    // ----- Faulted box boundaries: synthetic split faces -------------------
    // Phase 1 (docs/PLAN.md Track 1): for each faulted (box, axis, side), the
    // corner-point processor (faultedBoundaryConnections) gives the split
    // connections between the refined boundary cells and the neighbour(s). The
    // neighbour is a coarse cell, or -- when it is itself refined -- the
    // neighbour box's child cell across the fault (box<->box faulted interface,
    // LGR_GAPS A5). New split vertices join the leaf corner pool; each connection
    // becomes a leaf face with real polygon geometry. (Reference-matching transmissibility
    // geometry at these faces is a separate, later layer.)
    struct SyntheticFace
    {
        int boxCell;
        int neighborCell;                        // coarse, or a refined child (A5)
        std::vector<int> points;                 // leaf corner indices
        Dune::FieldVector<double,3> normal;      // unit, along +axis (CpGrid convention)
        Dune::FieldVector<double,3> center;
        double area;
        enum face_tag tag;
        int side;                                // box boundary side: -1 low, +1 high
        int box;
    };
    std::vector<SyntheticFace> syntheticFaces;
    std::vector<std::array<int,2>> splitKeys;    // {box, ordinal} per leaf corner past the history
    if (!faultedSides.empty()) {
        // coord -> leaf corner index, over the current pool (exact match: the
        // box-boundary corners are produced by identical resampling arithmetic).
        std::map<std::array<double,3>, int> cornerByCoord;
        for (int c = 0; c < static_cast<int>(leafCorners.size()); ++c) {
            cornerByCoord[coordKey(leafCorners[c])] = c;
        }

        // Map a coordinate to a leaf-pool corner (dedup; append if new).
        const auto poolCorner = [&](const std::array<double,3>& nd) {
            const std::array<double,3> key3{ nd[0], nd[1], nd[2] };
            auto pit = cornerByCoord.find(key3);
            if (pit != cornerByCoord.end()) {
                return pit->second;
            }
            const int idx = static_cast<int>(leafCorners.size());
            leafCorners.push_back(Dune::cpgrid::Geometry<0,3>(
                Dune::FieldVector<double,3>{ nd[0], nd[1], nd[2] }));
            cornerByCoord.emplace(key3, idx);
            return idx;
        };
        // Emit one synthetic face from a polygon given by node coordinates.
        // The normal points along +axis, never box -> neighbour: CpGrid derives a
        // face's min/max side (and hence Intersection::indexInInside) from the
        // stored normal together with each cell's orientation flag, so a normal
        // pointing down-axis makes both neighbours report the opposite face of
        // themselves. The side is kept for the face_to_cell ordering below.
        const auto emitPolygon = [&](int box, int boxLeaf, int neighborLeaf, int axis, int side,
                                     const std::vector<std::array<double,3>>& nodes) {
            std::vector<int> pts;
            pts.reserve(nodes.size());
            for (const auto& nd : nodes) {
                pts.push_back(poolCorner(nd));
            }
            PolyGeom pg = polygonGeometry(nodes);
            if (pg.normal[axis] < 0.0) {
                pg.normal *= -1.0;
                std::reverse(pts.begin(), pts.end());
            }
            syntheticFaces.push_back(SyntheticFace{
                boxLeaf, neighborLeaf, std::move(pts), pg.normal, pg.center,
                pg.area, faceTagOf(axis), side, box});
        };
        // Emit a box cell's WHOLE (axis, side) boundary face from its level grid.
        // Used to merge an over-refined shell's split pieces (unequal factors
        // across a fault) that all abut the SAME coarser neighbour child back into
        // one face, so the cell pair shares exactly one face.
        const auto emitWholeBoxFace = [&](int bIdx, int boxCellLevel, int boxLeaf,
                                          int neighborLeaf, int axis, int side) {
            auto& cf = GridStateWriter::cellToFace(*boxes[bIdx].level);
            auto& f2p = GridStateWriter::faceToPoint(*boxes[bIdx].level);
            Dune::cpgrid::EntityVariableBase<enum face_tag>& tags =
                GridStateWriter::faceTag(*boxes[bIdx].level);
            const Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& fnorm =
                GridStateWriter::faceNormals(*boxes[bIdx].level);
            const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<2,3>>& fgeom =
                GridStateWriter::geometry(*boxes[bIdx].level)
                  .geomVector(std::integral_constant<int,1>());
            const auto row = cf[EntityRep<0>(boxCellLevel, true)];
            for (int e = 0; e < row.size(); ++e) {
                const int f = row[e].index();
                if (axisOf(tags.get(f)) != axis || row[e].orientation() != (side > 0)) {
                    continue;
                }
                std::vector<int> pts;
                auto fp = f2p[f];
                for (int n = 0; n < fp.size(); ++n) {
                    pts.push_back(boxes[bIdx].cornerToLeaf[fp[n]]);
                }
                Dune::FieldVector<double,3> normal = fnorm[f];
                if (normal[axis] < 0.0) {
                    normal *= -1.0;
                    std::reverse(pts.begin(), pts.end());
                }
                syntheticFaces.push_back(SyntheticFace{
                    boxLeaf, neighborLeaf, std::move(pts), normal,
                    fgeom[f].center(), fgeom[f].volume(), faceTagOf(axis), side, bIdx});
                return;
            }
        };

        // Box<->box faulted connections, grouped by (box cell, neighbour child).
        // The emitting box's shell is refined to ITS factors; when it is finer
        // than the neighbour, one box-cell face can be split into several pieces
        // that all abut the SAME coarser neighbour child -- those are merged into
        // the box cell's whole face below. Coarse-neighbour and domain
        // connections are emitted directly.
        struct BoxBoxPiece {
            int boxLeaf, neighborLeaf, boxIdx, boxCellLevel, axis, side;
            std::vector<std::array<double,3>> nodes;
        };
        std::map<std::pair<int,int>, std::vector<BoxBoxPiece>> boxBoxGroups;

        for (const auto& [b, axis, side] : faultedSides) {
            const auto rangeIt = acrossRange.find(std::make_tuple(b, axis, side));
            const auto conns = faultedBoundaryConnections(
                parentDims, coord, zcorn, actnum, requests[b],
                axis, side, /*edgeConformal=*/true,
                (rangeIt != acrossRange.end()) ? &rangeIt->second : nullptr);
            const auto& rd = boxes[b].refinedDims;
            for (const auto& conn : conns) {
                const int refinedCart = conn.boxCell[0]
                                      + rd[0]*conn.boxCell[1]
                                      + rd[0]*rd[1]*conn.boxCell[2];
                auto cit = refinedCartToCell[b].find(refinedCart);
                if (cit == refinedCartToCell[b].end()) {
                    continue; // refined boundary cell inactive
                }
                const int boxLeaf = leafIdxOfLevelCell[b][cit->second];
                if (boxLeaf < 0) {
                    continue;
                }
                // coarseNeighborCart == -1: this part of the box boundary faces
                // the domain (the fault scarp) -> a boundary face (no outside).
                int neighborLeaf = -1;   // cell opposite the box across the fault
                bool boxBox = false;
                if (conn.coarseNeighborCart >= 0) {
                    const int l0 = compressed0[conn.coarseNeighborCart];
                    if (l0 < 0) {
                        continue; // neighbour inactive
                    }
                    neighborLeaf = leafIdxOfCell0[l0];
                    if (neighborLeaf < 0) {
                        // Neighbour is refined away: a box<->box faulted interface
                        // (LGR_GAPS A5). Connect to the neighbour box's child cell
                        // this sub-face abuts. Emit each interface face once: the
                        // FINER box owns it (equal in-face factors -> the lower box
                        // index), so the coarser side's faces are taken over by the
                        // finer side's sub-faces (a >6-face hex, like the A2 mosaic)
                        // -- now with the fault throw baked into the polygon. The
                        // mini-grid in faultedBoundaryConnections refines the shell
                        // to the emitting box's factors, so emitting from the finer
                        // side makes the sub-faces match the finer resolution; each
                        // maps to the coarser neighbour child that contains it.
                        const int nbBox = boxOfCell[l0];
                        if (nbBox < 0) {
                            continue;
                        }
                        const int u = (axis == 0) ? 1 : 0;
                        const int v = (axis == 2) ? 1 : 2;
                        if (boxes[b].graded || boxes[nbBox].graded) {
                            // Pairing the two sides' sub-faces below compares one
                            // subdivision factor per in-face direction, which a
                            // graded box does not have.
                            throw std::invalid_argument(
                                "Refinement '" + requests[b].name + "' meets '"
                                + requests[nbBox].name + "' across a fault, and at least "
                                "one of them is graded (N*FIN/H*FIN). Box-to-box "
                                "interfaces between graded boxes are not supported.");
                        }
                        const auto& facThis = boxes[b].factors;
                        const auto& facN    = boxes[nbBox].factors;
                        // Boxes that meet only through the throw miss the builder's
                        // IJK-contact guard; without nesting both sides would skip.
                        const auto nests = [u, v](const auto& fine, const auto& coarse) {
                            return fine[u] % coarse[u] == 0 && fine[v] % coarse[v] == 0;
                        };
                        if (!nests(facThis, facN) && !nests(facN, facThis)) {
                            throw std::invalid_argument(
                                "Refinement '" + requests[b].name + "' meets '"
                                + requests[nbBox].name + "' across a fault with in-face "
                                "subdivisions that do not nest. Use equal subdivisions, or "
                                "make one a multiple of the other.");
                        }
                        if (facThis[u] == facN[u] && facThis[v] == facN[v]) {
                            if (b >= nbBox) {
                                continue;        // equal: lower box index owns
                            }
                        }
                        else if (!(facThis[u] >= facN[u] && facThis[v] >= facN[v])) {
                            continue;            // coarser side: the finer references us
                        }
                        // (Nesting is checked above, so facThis is a multiple of facN
                        // in the in-face directions.)
                        const std::array<int,3> pijk = {
                            conn.coarseNeighborCart % dims0[0],
                            (conn.coarseNeighborCart / dims0[0]) % dims0[1],
                            conn.coarseNeighborCart / (dims0[0]*dims0[1]) };
                        const auto& startN = requests[nbBox].startIJK;
                        const auto& rdN    = boxes[nbBox].refinedDims;
                        // Neighbour child sub-position in the NEIGHBOUR's frame: the
                        // in-face components scale from the emitting box's sub-position
                        // by the factor ratio (finer -> coarser); the touch axis is the
                        // neighbour's boundary layer facing box `b`.
                        std::array<int,3> sN{};
                        sN[axis] = (side > 0) ? 0 : (facN[axis] - 1);
                        sN[u] = conn.coarseNeighborSub[u] * facN[u] / facThis[u];
                        sN[v] = conn.coarseNeighborSub[v] * facN[v] / facThis[v];
                        std::array<int,3> latN{};
                        for (int d = 0; d < 3; ++d) {
                            latN[d] = (pijk[d] - startN[d]) * facN[d] + sN[d];
                        }
                        const int cartN = latN[0] + rdN[0]*latN[1] + rdN[0]*rdN[1]*latN[2];
                        const auto nit = refinedCartToCell[nbBox].find(cartN);
                        if (nit == refinedCartToCell[nbBox].end()) {
                            continue; // neighbour child inactive
                        }
                        neighborLeaf = leafIdxOfLevelCell[nbBox][nit->second];
                        if (neighborLeaf < 0) {
                            continue;
                        }
                        boxBox = true;
                    }
                }

                if (boxBox) {
                    boxBoxGroups[{boxLeaf, neighborLeaf}].push_back(
                        BoxBoxPiece{boxLeaf, neighborLeaf, b, cit->second, axis, side,
                                    conn.faceNodes});
                }
                else {
                    emitPolygon(b, boxLeaf, neighborLeaf, axis, side, conn.faceNodes);
                }
            }
        }

        // Emit the box<->box interface: one face per (box cell, neighbour child).
        // A single piece is emitted as-is; multiple pieces (the over-refined shell
        // split one face) are merged into the box cell's whole boundary face.
        for (auto& [keyPair, pieces] : boxBoxGroups) {
            if (pieces.size() == 1) {
                const auto& p = pieces.front();
                emitPolygon(p.boxIdx, p.boxLeaf, p.neighborLeaf, p.axis, p.side, p.nodes);
            }
            else {
                const auto& p = pieces.front();
                emitWholeBoxFace(p.boxIdx, p.boxCellLevel, p.boxLeaf, p.neighborLeaf,
                                 p.axis, p.side);
            }
        }
        numLeafCorners = static_cast<int>(leafCorners.size());

        // Split corners belong to the lowest box using them and are ordered by coordinate
        // there, so their ids do not depend on the order the faces were emitted in.
        const int firstSplit = static_cast<int>(leafCornerHistory.size());
        std::vector<int> owner(numLeafCorners - firstSplit, numBoxes);
        for (const auto& sf : syntheticFaces) {
            for (const int p : sf.points) {
                if (p >= firstSplit) {
                    owner[p - firstSplit] = std::min(owner[p - firstSplit], sf.box);
                }
            }
        }
        splitKeys.assign(owner.size(), {-1, -1});
        for (int b = 0; b < numBoxes; ++b) {
            std::vector<int> mine;
            for (std::size_t s = 0; s < owner.size(); ++s) {
                if (owner[s] == b) {
                    mine.push_back(static_cast<int>(s));
                }
            }
            std::ranges::sort(mine, [&](int x, int y) {
                return coordKey(leafCorners[firstSplit + x]) < coordKey(leafCorners[firstSplit + y]);
            });
            for (std::size_t n = 0; n < mine.size(); ++n) {
                splitKeys[mine[n]] = {b, static_cast<int>(n)};
            }
        }
    }

    const int numLeafFaces = numSourceFaces + static_cast<int>(syntheticFaces.size());

    // ----- Leaf topology ---------------------------------------------------
    auto& leafFaceToCell = GridStateWriter::faceToCell(*leaf);
    auto& leafFaceToPoint = GridStateWriter::faceToPoint(*leaf);
    auto& leafCellToPoint = GridStateWriter::cellToPoint(*leaf);
    auto& leafCellToFace = GridStateWriter::cellToFace(*leaf);
    auto& leafFaceTag = GridStateWriter::faceTag(*leaf);
    auto& leafFaceNormals = GridStateWriter::faceNormals(*leaf);

    Dune::cpgrid::EntityVariableBase<enum face_tag>& leafTags = leafFaceTag;
    Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& leafNormals = leafFaceNormals;
    leafTags.resize(numLeafFaces);
    leafNormals.resize(numLeafFaces);

    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<2,3>>& leafFaceGeom =
        leafGeometry.geomVector(std::integral_constant<int,1>());
    leafFaceGeom.resize(numLeafFaces);

    std::vector<EntityRep<0>> rowBuffer;
    std::vector<int> pointBuffer;
    for (int face = 0; face < numSourceFaces; ++face) {
        const SourceRef src = leafFaces[face];
        CpGridData& srcGrid = (src.grid == 0) ? level0 : *boxes[src.grid - 1].level;
        auto& srcFaceToCell = (src.grid == 0) ? faceToCell0 : GridStateWriter::faceToCell(srcGrid);
        auto& srcFaceToPoint = (src.grid == 0) ? faceToPoint0 : GridStateWriter::faceToPoint(srcGrid);
        const Dune::cpgrid::EntityVariableBase<enum face_tag>& srcTags =
            (src.grid == 0) ? static_cast<const Dune::cpgrid::EntityVariableBase<enum face_tag>&>(faceTag0)
                            : static_cast<const Dune::cpgrid::EntityVariableBase<enum face_tag>&>(GridStateWriter::faceTag(srcGrid));
        const Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& srcNormals =
            (src.grid == 0) ? static_cast<const Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>&>(faceNormals0)
                            : static_cast<const Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>&>(GridStateWriter::faceNormals(srcGrid));
        const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<2,3>>& srcFaceGeom =
            GridStateWriter::geometry(srcGrid).geomVector(std::integral_constant<int,1>());

        // face_to_cell row (cells mapped to leaf indices; mosaic outside added)
        rowBuffer.clear();
        const auto srcRow = srcFaceToCell[EntityRep<1>(src.index, true)];
        for (int e = 0; e < srcRow.size(); ++e) {
            const int srcCell = srcRow[e].index();
            // Preserve the off-rank sentinel: a face to a remote cell stays
            // a remote-neighbor face on the distributed leaf.
            const int leafCell = (srcCell == kRemoteCell)
                ? kRemoteCell
                : (src.grid == 0)
                    ? leafIdxOfCell0[srcCell]
                    : leafIdxOfLevelCell[src.grid - 1][srcCell];
            rowBuffer.push_back(EntityRep<0>(leafCell, srcRow[e].orientation()));
        }
        if (mosaicOutside[face] >= 0) {
            // The refined cell entry already exists; add the coarse neighbor
            // with opposite orientation, keeping normal-side first.
            const bool refinedNormalOut = srcRow[0].orientation();
            const EntityRep<0> outsideEntry(mosaicOutside[face], !refinedNormalOut);
            if (refinedNormalOut) {
                rowBuffer.push_back(outsideEntry);
            }
            else {
                rowBuffer.insert(rowBuffer.begin(), outsideEntry);
            }
        }
        leafFaceToCell.appendRow(rowBuffer.begin(), rowBuffer.end());

        // face_to_point row (pool indices)
        pointBuffer.clear();
        auto srcPoints = srcFaceToPoint[src.index];
        for (int e = 0; e < srcPoints.size(); ++e) {
            pointBuffer.push_back((src.grid == 0) ? srcPoints[e]
                                                  : boxes[src.grid - 1].cornerToLeaf[srcPoints[e]]);
        }
        leafFaceToPoint.appendRow(pointBuffer.begin(), pointBuffer.end());

        leafTags[face] = srcTags[src.index];
        leafNormals[face] = srcNormals[src.index];
        leafFaceGeom[face] = srcFaceGeom[src.index];
    }

    // Synthetic faulted-boundary faces: real polygon geometry, the cell on the
    // face's negative side first with orientation true (the mosaic convention,
    // and what makeInverseRelation would produce). The normal runs along +axis,
    // so that cell is the box cell on the box's high side and the neighbour
    // (coarse cell, or a refined child for a box<->box faulted interface) on its
    // low side -- not the box cell either way.
    for (int s = 0; s < static_cast<int>(syntheticFaces.size()); ++s) {
        const int face = numSourceFaces + s;
        const SyntheticFace& sf = syntheticFaces[s];
        std::vector<EntityRep<0>> cells;
        const bool boxOnNegativeSide = (sf.side > 0);
        if (boxOnNegativeSide) {
            cells.emplace_back(sf.boxCell, true);
            if (sf.neighborCell >= 0) {
                cells.emplace_back(sf.neighborCell, false);  // interior split face
            }
        }
        else {
            if (sf.neighborCell >= 0) {
                cells.emplace_back(sf.neighborCell, true);   // interior split face
            }
            cells.emplace_back(sf.boxCell, false);
        }
        // A missing neighbour means the box boundary faces the domain (fault
        // scarp) -> boundary face.
        leafFaceToCell.appendRow(cells.begin(), cells.end());
        leafFaceToPoint.appendRow(sf.points.begin(), sf.points.end());
        leafTags[face] = sf.tag;
        leafNormals[face] = sf.normal;
        leafFaceGeom[face] = Dune::cpgrid::Geometry<2,3>(sf.center, sf.area);
    }

    // cell_to_face_ as the inverse of face_to_cell. Built manually rather
    // than via makeInverseRelation because face_to_cell may carry the
    // off-rank sentinel (distributed grids), which the generic inverse — it
    // sizes the table from the maximum cell index — cannot handle. The
    // orientation convention matches makeInverseRelation: a cell keeps the
    // face's orientation when it is on the face's normal side, else flips it.
    {
        std::vector<std::vector<EntityRep<1>>> cellFaces(numLeafCells);
        for (int face = 0; face < numLeafFaces; ++face) {
            const auto row = leafFaceToCell[EntityRep<1>(face, true)];
            for (int e = 0; e < row.size(); ++e) {
                const int cell = row[e].index();
                if (cell == kRemoteCell) {
                    continue; // remote neighbor: no local cell owns this side
                }
                const EntityRep<1> faceEnt(face, true);
                cellFaces[cell].push_back(row[e].orientation() ? faceEnt : faceEnt.opposite());
            }
        }
        for (int cell = 0; cell < numLeafCells; ++cell) {
            leafCellToFace.appendRow(cellFaces[cell].begin(), cellFaces[cell].end());
        }
    }

    // ----- Leaf cell data ----------------------------------------------------
    leafCellToPoint.resize(numLeafCells);
    std::vector<int> leafGlobalCell(numLeafCells);
    std::vector<std::array<int,2>> leafChildToParent(numLeafCells, std::array<int,2>{-1, -1});
    std::vector<int> leafIdxInParent(numLeafCells, -1);
    std::vector<std::array<int,2>> leafToLevel(numLeafCells);

    for (int cell = 0; cell < numLeafCells; ++cell) {
        const SourceRef src = leafCells[cell];
        if (src.grid == 0) {
            leafCellToPoint[cell] = cellToPoint0[src.index];
            leafGlobalCell[cell] = level0.globalCell()[src.index];
            leafToLevel[cell] = {0, src.index};
        }
        else {
            const BoxData& box = boxes[src.grid - 1];
            const auto& srcPoints = GridStateWriter::cellToPoint(*box.level)[src.index];
            for (int corner = 0; corner < 8; ++corner) {
                leafCellToPoint[cell][corner] = box.cornerToLeaf[srcPoints[corner]];
            }
            const int parent = GridStateWriter::childToParent(*box.level)[src.index][1];
            leafGlobalCell[cell] = level0.globalCell()[parent];
            leafChildToParent[cell] = {0, parent};
            leafIdxInParent[cell] = GridStateWriter::idxInParent(*box.level)[src.index];
            leafToLevel[cell] = {src.grid, src.index};
        }
    }

    // Cell geometries: centers/volumes from the source grids, corner list
    // pointing into the leaf pool via the leaf's cell_to_point_.
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<3,3>>& leafCellGeom =
        leafGeometry.geomVector(std::integral_constant<int,0>());
    leafCellGeom.resize(numLeafCells);
    auto leafCornersPtr = leafGeometry.pointGeometryPtr();
    for (int cell = 0; cell < numLeafCells; ++cell) {
        const SourceRef src = leafCells[cell];
        CpGridData& srcGrid = (src.grid == 0) ? level0 : *boxes[src.grid - 1].level;
        const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<3,3>>& srcCellGeom =
            GridStateWriter::geometry(srcGrid).geomVector(std::integral_constant<int,0>());
        const auto& geom = srcCellGeom[src.index];
        leafCellGeom[cell] = Dune::cpgrid::Geometry<3,3>(geom.center(), geom.volume(),
                                                         leafCornersPtr,
                                                         leafCellToPoint[cell].data());
    }

    // ----- Parallel leaf structures (distributed grids only) -------------------
    // Build the leaf's parallel cell index set, remote indices, partition types,
    // and communication interfaces so that overlap leaf cells take part in
    // collective communication (e.g. opm-models' dofTotalVolume border sync,
    // which sums interior volumes and then communicates them to the overlap).
    // Coarse leaf cells inherit their level-zero global index and attribute;
    // refined cells are rank-interior (owner only) and get assignRefinedIds' ids,
    // which are above every coarse global id.
    const auto ids = assignRefinedIds(storage, numBoxes, leafToLevel, leafCornerHistory,
                                      splitKeys, comm);
#if HAVE_MPI
    Dune::Communication<Dune::MPIHelper::MPICommunicator> cc(comm);
    if (cc.size() > 1) {
        using ParallelIndexSet = CpGridData::ParallelIndexSet;
        using AttributeSet = CpGridData::AttributeSet;

        std::vector<int> leafCellGlobalId(numLeafCells, -1);
        std::vector<AttributeSet> leafCellAttr(numLeafCells, AttributeSet::owner);

        for (const auto& entry : level0.cellIndexSet()) {
            const int l0local = entry.local().local();
            const int leafC = leafIdxOfCell0[l0local];   // -1 if refined away
            if (leafC >= 0) {
                leafCellGlobalId[leafC] = entry.global();
                leafCellAttr[leafC] = entry.local().attribute();
            }
        }

        for (int cell = 0; cell < numLeafCells; ++cell) {
            if (leafCells[cell].grid != 0) {
                leafCellGlobalId[cell] = ids.cells[cell];
                leafCellAttr[cell] = AttributeSet::owner;
            }
        }

        auto& cis = leaf->cellIndexSet();
        cis.beginResize();
        for (int cell = 0; cell < numLeafCells; ++cell) {
            cis.add(leafCellGlobalId[cell],
                    ParallelIndexSet::LocalIndex(cell, leafCellAttr[cell], true));
        }
        cis.endResize();

        leaf->cellRemoteIndices().template rebuild<false>();
        leaf->computeCellPartitionType();
        leaf->computePointPartitionType();
        leaf->computeCommunicationInterfaces(numLeafCorners);
        setDistributedLeafIds(*leaf, level0, leafCellGlobalId, leafFaces, numLeafFaces,
                              ids.points, ids.end, cc);
    }
    else
#endif
    {
        std::vector<int> faceIds(numLeafFaces);
        std::iota(faceIds.begin(), faceIds.end(), numLeafCells);
        GridStateWriter::setGlobalIdMapping(*leaf, ids.cells, std::move(faceIds), ids.points);
    }

    // ----- Leaf metadata -------------------------------------------------------
    GridStateWriter::setLogicalCartesianSize(*leaf, dims0);
    GridStateWriter::setGlobalCell(*leaf, std::move(leafGlobalCell));
    GridStateWriter::setIndexSet(*leaf, numLeafCells, numLeafCorners);
    GridStateWriter::setParentRelations(*leaf, std::move(leafChildToParent), std::move(leafIdxInParent));
    GridStateWriter::setAquiferCells(*leaf, leafAquiferCells(level0, leafToLevel));
    GridStateWriter::setLeafToLevel(*leaf, std::move(leafToLevel));
    GridStateWriter::setCornerHistory(*leaf, std::move(leafCornerHistory));
    GridStateWriter::setRefinementMaxLevel(*leaf, numBoxes);
    for (int b = 0; b < numBoxes; ++b) {
        GridStateWriter::setRefinementMaxLevel(*boxes[b].level, numBoxes);
    }

    // Level zero learns about its children.
    GridStateWriter::setParentToChildren(level0, std::move(parentToChildren));

    return leaf;
}

std::shared_ptr<CpGridData>
assembleNestedLeafGrid(std::vector<std::shared_ptr<CpGridData>>& storage,
                       const std::vector<BlockRefinement>& requests,
                       const std::vector<LevelGeom>& levelGeom,
                       Dune::MPIHelper::MPICommunicator comm)
{
    const int numBoxes = static_cast<int>(requests.size());
    CpGridData& level0 = *storage[0];
    const auto& dims0 = level0.logicalCartesianSize();
    const int numCells0 = level0.size(0);
    const int numCorners0 = level0.size(3);

    // ----- Box tree (parentGridName -> level index b+1) --------------------
    std::map<std::string,int> nameToBox;
    for (int b = 0; b < numBoxes; ++b) {
        nameToBox[requests[b].name] = b;
    }
    std::vector<int> parentBoxOf(numBoxes, -1);   // -1 == parent is GLOBAL (level 0)
    std::vector<int> parentLevelOf(numBoxes, 0);  // grid level of the parent (0 == level0)
    for (int b = 0; b < numBoxes; ++b) {
        if (requests[b].parentGridName != "GLOBAL") {
            parentBoxOf[b] = nameToBox.at(requests[b].parentGridName);
            parentLevelOf[b] = parentBoxOf[b] + 1;
        }
    }

    // ----- Per-box data ----------------------------------------------------
    struct BoxData
    {
        CpGridData* level;
        std::array<int,3> factors;
        std::array<AxisPositions,3> pos;
        std::array<int,3> refinedDims;
        std::map<int, std::vector<int>> childrenOfParent;  // parent cell (in parent grid) -> children
        std::vector<int> cornerEquiv;                      // level corner -> parent leaf corner (or -1)
        std::vector<int> cornerToLeaf;                     // level corner -> leaf pool index
    };
    std::vector<BoxData> boxes(numBoxes);
    // childBoxOf[L][cell] = the box (index) refining that cell of grid level L, or -1.
    std::vector<std::vector<int>> childBoxOf(numBoxes + 1);
    childBoxOf[0].assign(numCells0, -1);
    for (int b = 0; b < numBoxes; ++b) {
        childBoxOf[b + 1].assign(storage[b + 1]->size(0), -1);
    }
    for (int b = 0; b < numBoxes; ++b) {
        BoxData& box = boxes[b];
        box.level = storage[b + 1].get();
        box.factors = requests[b].cellsPerDim;
        for (int d = 0; d < 3; ++d) {
            box.pos[d] = axisPositions(axisSubdivision(requests[b], d));
        }
        box.refinedDims = box.level->logicalCartesianSize();
        const auto& childToParent = GridStateWriter::childToParent(*box.level);
        const auto& idxInParent = GridStateWriter::idxInParent(*box.level);
        std::map<int, std::vector<std::pair<int,int>>> tmp;
        for (int cell = 0; cell < box.level->size(0); ++cell) {
            tmp[childToParent[cell][1]].emplace_back(idxInParent[cell], cell);
        }
        for (auto& [parent, kids] : tmp) {
            std::sort(kids.begin(), kids.end());
            auto& target = box.childrenOfParent[parent];
            target.reserve(kids.size());
            for (const auto& [pos, cell] : kids) {
                target.push_back(cell);
            }
            childBoxOf[parentLevelOf[b]][parent] = b;   // parent cell (in its grid) is refined by b
        }
    }

    // ----- Containment validation (fully-contained, one nesting level) -----
    for (int b = 0; b < numBoxes; ++b) {
        if (parentBoxOf[b] < 0) {
            continue;
        }
        // Rank-interior ("refine-after") model: a box that is not present on this
        // rank is an empty placeholder whose parent LGR level grid is also empty,
        // so levelGeom carries no dims to validate against (validating here would
        // throw spuriously on the ranks that do not own the box, deadlocking
        // against the ranks that do). The box is validated on the rank that owns
        // it. A contained child is always co-located with its parent, so the
        // child's own emptiness is a sufficient test.
        if (boxes[b].level->size(0) == 0) {
            continue;
        }
        if (parentBoxOf[parentBoxOf[b]] >= 0) {
            throw std::logic_error("Nested LGR deeper than one level is not implemented "
                "yet (docs/NESTED_LGR_PLAN.md Phase C).");
        }
        const auto& req = requests[b];
        const auto& pdims = levelGeom[parentLevelOf[b]].dims;
        const auto ijkRange = [&]() {
            std::string s;
            for (int d = 0; d < 3; ++d) {
                s += (d ? " x [" : "[") + std::to_string(req.startIJK[d]) + ","
                   + std::to_string(req.endIJK[d]) + ")";
            }
            return s;
        };
        const std::string parentExtent = std::to_string(pdims[0]) + "x"
            + std::to_string(pdims[1]) + "x" + std::to_string(pdims[2]);
        for (int d = 0; d < 3; ++d) {
            if (req.startIJK[d] < 0 || req.endIJK[d] > pdims[d]
                || req.startIJK[d] >= req.endIJK[d]) {
                throw std::logic_error(
                    "Nested refinement box '" + req.name + "' extends outside its parent LGR '"
                    + req.parentGridName + "' (parent has " + parentExtent + " refined cells): "
                    "its parent-local IJK range " + ijkRange()
                    + " must be non-empty and lie within the parent.");
            }
            if (req.startIJK[d] < 1 || req.endIJK[d] > pdims[d] - 1) {
                throw std::logic_error(
                    "Nested refinement box '" + req.name + "' touches the boundary of its "
                    "parent LGR '" + req.parentGridName + "' (parent " + parentExtent
                    + "; box parent-local IJK " + ijkRange() + "). Only a child fully "
                    "contained in the interior of its parent LGR is supported yet "
                    "(docs/NESTED_LGR_PLAN.md Phase C): every direction must satisfy "
                    "1 <= startIJK and endIJK <= parentDim-1.");
            }
        }
    }

    // Parent grid of box b, and the leaf-pool index of one of its corners.
    const auto parentGridOf = [&](int b) -> CpGridData& {
        return parentLevelOf[b] == 0 ? level0 : *boxes[parentBoxOf[b]].level;
    };
    const auto parentCornerLeaf = [&](int b, int parentCornerIdx) -> int {
        // level-0 corners occupy leaf-pool indices [0, numCorners0) unchanged.
        return parentLevelOf[b] == 0 ? parentCornerIdx
                                     : boxes[parentBoxOf[b]].cornerToLeaf[parentCornerIdx];
    };

    // ----- Corner pool (level0 corners + refined corners), parent-before-child
    auto leaf = std::make_shared<CpGridData>(comm, storage);
    auto& leafGeometry = GridStateWriter::geometry(*leaf);
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& leafCorners =
        leafGeometry.geomVector(std::integral_constant<int,3>());
    const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& corners0 =
        GridStateWriter::geometry(level0).geomVector(std::integral_constant<int,3>());
    leafCorners.assign(corners0.begin(), corners0.end());
    std::vector<std::array<int,2>> leafCornerHistory(numCorners0);
    for (int corner = 0; corner < numCorners0; ++corner) {
        leafCornerHistory[corner] = {0, corner};
    }
    std::map<std::array<double,3>, int> refinedCornerPool;
    const auto coordKey = [](const Dune::cpgrid::Geometry<0,3>& corner) {
        const auto& c = corner.center();
        return std::array<double,3>{ c[0], c[1], c[2] };
    };

    for (int b = 0; b < numBoxes; ++b) {   // request order is parent-before-child
        BoxData& box = boxes[b];
        CpGridData& pgrid = parentGridOf(b);
        auto& parentCellToPoint = GridStateWriter::cellToPoint(pgrid);
        const auto& childToParent = GridStateWriter::childToParent(*box.level);
        auto& cellToPointL = GridStateWriter::cellToPoint(*box.level);

        // Per-cell corner identification against the PARENT grid (level0 or a
        // box): a refined corner coincides with a parent corner when it sits on
        // the parent's boundary in all three directions. Read from the cell's
        // refined lattice position, whose relation to the parent holds whether
        // or not the box is graded.
        const auto& rdb = box.refinedDims;
        const int maskP = cornerSlotMask(pgrid);
        const int maskL = cornerSlotMask(*box.level);
        box.cornerEquiv.assign(box.level->size(3), -1);
        for (int cell = 0; cell < box.level->size(0); ++cell) {
            const int parent = childToParent[cell][1];
            const int refinedCart = box.level->globalCell()[cell];
            const std::array<int,3> lattice = { refinedCart % rdb[0],
                                                (refinedCart / rdb[0]) % rdb[1],
                                                refinedCart / (rdb[0]*rdb[1]) };
            for (int corner = 0; corner < 8; ++corner) {
                const std::array<int,3> d = { corner & 1, (corner >> 1) & 1, (corner >> 2) & 1 };
                std::array<int,3> parentCornerIJK{};
                bool onParentCorner = true;
                for (int c = 0; c < 3 && onParentCorner; ++c) {
                    const int edge = box.pos[c].subIndex[lattice[c]] + d[c];
                    if (edge == 0) {
                        parentCornerIJK[c] = 0;
                    }
                    else if (edge == box.pos[c].parentCount[lattice[c]]) {
                        parentCornerIJK[c] = 1;
                    }
                    else {
                        onParentCorner = false;
                    }
                }
                if (onParentCorner) {
                    const int parentCorner = parentCellToPoint[parent][(parentCornerIJK[0]
                                                                       + 2*parentCornerIJK[1]
                                                                       + 4*parentCornerIJK[2]) ^ maskP];
                    box.cornerEquiv[cellToPointL[cell][corner ^ maskL]] = parentCornerLeaf(b, parentCorner);
                }
            }
        }

        // Pool: equivalent corners reuse the parent leaf corner; the rest are
        // appended (exact-coordinate dedup with already-appended refined corners).
        const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& cornersL =
            GridStateWriter::geometry(*box.level).geomVector(std::integral_constant<int,3>());
        box.cornerToLeaf.assign(box.level->size(3), -1);
        std::vector<std::array<int,2>> levelCornerHistory(box.level->size(3), {-1, -1});
        for (int corner = 0; corner < box.level->size(3); ++corner) {
            if (box.cornerEquiv[corner] >= 0) {
                box.cornerToLeaf[corner] = box.cornerEquiv[corner];
                levelCornerHistory[corner] = leafCornerHistory[box.cornerEquiv[corner]];
            }
            else {
                const auto key = coordKey(cornersL[corner]);
                const auto existing = refinedCornerPool.find(key);
                if (existing != refinedCornerPool.end()) {
                    box.cornerToLeaf[corner] = existing->second;
                    levelCornerHistory[corner] = leafCornerHistory[existing->second];
                }
                else {
                    const int leafCorner = static_cast<int>(leafCorners.size());
                    box.cornerToLeaf[corner] = leafCorner;
                    leafCorners.push_back(cornersL[corner]);
                    leafCornerHistory.push_back({b + 1, corner});
                    refinedCornerPool.emplace(key, leafCorner);
                }
            }
        }
        GridStateWriter::setCornerHistory(*box.level, std::move(levelCornerHistory));
    }
    const int numLeafCorners = static_cast<int>(leafCorners.size());

    // ----- Leaf cells: tree walk, parents replaced by children -------------
    std::vector<SourceRef> leafCells;
    std::vector<std::vector<int>> leafIdxOf(numBoxes + 1);
    leafIdxOf[0].assign(numCells0, -1);
    for (int b = 0; b < numBoxes; ++b) {
        leafIdxOf[b + 1].assign(boxes[b].level->size(0), -1);
    }
    std::vector<std::tuple<int, std::vector<int>>> parentToChildren0(numCells0,
                                                                    std::make_tuple(-1, std::vector<int>{}));
    std::function<void(int,int)> emitCell = [&](int level, int cell) {
        const int b = childBoxOf[level][cell];
        if (b < 0) {
            leafIdxOf[level][cell] = static_cast<int>(leafCells.size());
            leafCells.push_back({level, cell});
            return;
        }
        const auto it = boxes[b].childrenOfParent.find(cell);
        if (it == boxes[b].childrenOfParent.end()) {
            return;  // active parent without active children cannot happen
        }
        if (level == 0) {
            parentToChildren0[cell] = std::make_tuple(b + 1, it->second);
        }
        for (const int child : it->second) {
            emitCell(b + 1, child);
        }
    };
    for (int c = 0; c < numCells0; ++c) {
        emitCell(0, c);
    }
    const int numLeafCells = static_cast<int>(leafCells.size());

    // ----- Faces -----------------------------------------------------------
    // A face of grid level L is dropped if a side cell is refined by a child
    // box (the child's mosaic/interior faces take over). Applies to level0 and
    // to any parent box level. Per-level dropped flags.
    std::vector<std::vector<char>> faceDropped(numBoxes + 1);
    {
        auto& f2c0 = GridStateWriter::faceToCell(level0);
        faceDropped[0].assign(f2c0.size(), 0);
        for (int face = 0; face < f2c0.size(); ++face) {
            const auto row = f2c0[EntityRep<1>(face, true)];
            for (int e = 0; e < row.size(); ++e) {
                const int cell = row[e].index();
                if (cell != kRemoteCell && childBoxOf[0][cell] >= 0) {
                    faceDropped[0][face] = 1; break;
                }
            }
        }
        for (int b = 0; b < numBoxes; ++b) {
            auto& f2cL = GridStateWriter::faceToCell(*boxes[b].level);
            faceDropped[b + 1].assign(f2cL.size(), 0);
            for (int face = 0; face < f2cL.size(); ++face) {
                const auto row = f2cL[EntityRep<1>(face, true)];
                for (int e = 0; e < row.size(); ++e) {
                    const int cell = row[e].index();
                    if (cell != kRemoteCell && childBoxOf[b + 1][cell] >= 0) {
                        faceDropped[b + 1][face] = 1; break;
                    }
                }
            }
        }
    }

    // Lattice neighbour in a Cartesian parent grid (contained case has no faults
    // at the nested boundary, so the simple lattice neighbour is correct).
    std::vector<std::vector<int>> cartToCell(numBoxes + 1);
    {
        cartToCell[0].assign(static_cast<std::size_t>(dims0[0])*dims0[1]*dims0[2], -1);
        for (int c = 0; c < numCells0; ++c) cartToCell[0][level0.globalCell()[c]] = c;
        for (int b = 0; b < numBoxes; ++b) {
            const auto& rd = boxes[b].refinedDims;
            cartToCell[b + 1].assign(static_cast<std::size_t>(rd[0])*rd[1]*rd[2], -1);
            for (int c = 0; c < boxes[b].level->size(0); ++c)
                cartToCell[b + 1][boxes[b].level->globalCell()[c]] = c;
        }
    }
    const auto latticeNeighbour = [&](int gridLevel, int cell, int axis, int side) -> int {
        const auto& dims = (gridLevel == 0) ? dims0 : boxes[gridLevel - 1].refinedDims;
        CpGridData& g = (gridLevel == 0) ? level0 : *boxes[gridLevel - 1].level;
        const int cart = g.globalCell()[cell];
        std::array<int,3> ijk = { cart % dims[0], (cart / dims[0]) % dims[1],
                                  cart / (dims[0]*dims[1]) };
        ijk[axis] += side;
        if (ijk[axis] < 0 || ijk[axis] >= dims[axis]) return -1;          // domain boundary
        const int ncart = ijk[0] + dims[0]*(ijk[1] + dims[1]*ijk[2]);
        return cartToCell[gridLevel][ncart];                              // -1 if inactive
    };

    // Leaf faces: level0 non-dropped, then each box's faces (interior + mosaic).
    std::vector<SourceRef> leafFaces;
    std::vector<int> mosaicOutside;                 // outside leaf cell for a box boundary face
    {
        auto& f2c0 = GridStateWriter::faceToCell(level0);
        for (int face = 0; face < f2c0.size(); ++face) {
            if (!faceDropped[0][face]) {
                leafFaces.push_back({0, face});
                mosaicOutside.push_back(-1);
            }
        }
    }
    for (int b = 0; b < numBoxes; ++b) {
        BoxData& box = boxes[b];
        auto& faceToCellL = GridStateWriter::faceToCell(*box.level);
        Dune::cpgrid::EntityVariableBase<enum face_tag>& tagsL = GridStateWriter::faceTag(*box.level);
        const auto& childToParent = GridStateWriter::childToParent(*box.level);
        for (int face = 0; face < faceToCellL.size(); ++face) {
            if (faceDropped[b + 1][face]) {
                continue;                            // toward a deeper-refined cell
            }
            const auto row = faceToCellL[EntityRep<1>(face, true)];
            int outside = -1;
            if (row.size() == 1) {
                const int cell = row[0].index();
                const bool normalOut = row[0].orientation();
                const int axis = axisOf(tagsL.get(face));
                const int side = normalOut ? 1 : -1;
                const int cart = box.level->globalCell()[cell];
                std::array<int,3> lattice = { cart % box.refinedDims[0],
                                              (cart / box.refinedDims[0]) % box.refinedDims[1],
                                              cart / (box.refinedDims[0]*box.refinedDims[1]) };
                lattice[axis] += side;
                if (lattice[axis] < 0 || lattice[axis] >= box.refinedDims[axis]) {
                    // Block boundary: stitch against the parent-grid neighbour.
                    const int parent = childToParent[cell][1];
                    const int neighbour = latticeNeighbour(parentLevelOf[b], parent, axis, side);
                    if (neighbour >= 0 && childBoxOf[parentLevelOf[b]][neighbour] < 0) {
                        outside = leafIdxOf[parentLevelOf[b]][neighbour];
                    }
                    // else: domain boundary / deeper-refined neighbour (cannot
                    // happen for a contained child) -> stays a boundary face.
                }
                // else: internal hole boundary -> stays a boundary face.
            }
            leafFaces.push_back({b + 1, face});
            mosaicOutside.push_back(outside);
        }
    }
    const int numLeafFaces = static_cast<int>(leafFaces.size());

    // ----- Leaf topology + geometry write-out ------------------------------
    auto& leafFaceToCell = GridStateWriter::faceToCell(*leaf);
    auto& leafFaceToPoint = GridStateWriter::faceToPoint(*leaf);
    auto& leafCellToPoint = GridStateWriter::cellToPoint(*leaf);
    auto& leafCellToFace = GridStateWriter::cellToFace(*leaf);
    Dune::cpgrid::EntityVariableBase<enum face_tag>& leafTags = GridStateWriter::faceTag(*leaf);
    Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& leafNormals = GridStateWriter::faceNormals(*leaf);
    leafTags.resize(numLeafFaces);
    leafNormals.resize(numLeafFaces);
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<2,3>>& leafFaceGeom =
        leafGeometry.geomVector(std::integral_constant<int,1>());
    leafFaceGeom.resize(numLeafFaces);

    std::vector<EntityRep<0>> rowBuffer;
    std::vector<int> pointBuffer;
    for (int face = 0; face < numLeafFaces; ++face) {
        const SourceRef src = leafFaces[face];
        CpGridData& srcGrid = (src.grid == 0) ? level0 : *boxes[src.grid - 1].level;
        auto& srcFaceToCell = GridStateWriter::faceToCell(srcGrid);
        auto& srcFaceToPoint = GridStateWriter::faceToPoint(srcGrid);
        const Dune::cpgrid::EntityVariableBase<enum face_tag>& srcTags = GridStateWriter::faceTag(srcGrid);
        const Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& srcNormals = GridStateWriter::faceNormals(srcGrid);
        const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<2,3>>& srcFaceGeom =
            GridStateWriter::geometry(srcGrid).geomVector(std::integral_constant<int,1>());

        rowBuffer.clear();
        const auto srcRow = srcFaceToCell[EntityRep<1>(src.index, true)];
        for (int e = 0; e < srcRow.size(); ++e) {
            const int srcCell = srcRow[e].index();
            const int leafCell = (srcCell == kRemoteCell) ? kRemoteCell
                                                          : leafIdxOf[src.grid][srcCell];
            rowBuffer.push_back(EntityRep<0>(leafCell, srcRow[e].orientation()));
        }
        if (mosaicOutside[face] >= 0) {
            const bool refinedNormalOut = srcRow[0].orientation();
            const EntityRep<0> outsideEntry(mosaicOutside[face], !refinedNormalOut);
            if (refinedNormalOut) rowBuffer.push_back(outsideEntry);
            else rowBuffer.insert(rowBuffer.begin(), outsideEntry);
        }
        leafFaceToCell.appendRow(rowBuffer.begin(), rowBuffer.end());

        pointBuffer.clear();
        auto srcPoints = srcFaceToPoint[src.index];
        for (int e = 0; e < srcPoints.size(); ++e) {
            pointBuffer.push_back((src.grid == 0) ? srcPoints[e]
                                                  : boxes[src.grid - 1].cornerToLeaf[srcPoints[e]]);
        }
        leafFaceToPoint.appendRow(pointBuffer.begin(), pointBuffer.end());
        leafTags[face] = srcTags[src.index];
        leafNormals[face] = srcNormals[src.index];
        leafFaceGeom[face] = srcFaceGeom[src.index];
    }

    // cell_to_face_ = inverse of face_to_cell_.
    {
        std::vector<std::vector<EntityRep<1>>> cellFaces(numLeafCells);
        for (int face = 0; face < numLeafFaces; ++face) {
            const auto row = leafFaceToCell[EntityRep<1>(face, true)];
            for (int e = 0; e < row.size(); ++e) {
                const int cell = row[e].index();
                if (cell == kRemoteCell) continue;
                const EntityRep<1> faceEnt(face, true);
                cellFaces[cell].push_back(row[e].orientation() ? faceEnt : faceEnt.opposite());
            }
        }
        for (int cell = 0; cell < numLeafCells; ++cell) {
            leafCellToFace.appendRow(cellFaces[cell].begin(), cellFaces[cell].end());
        }
    }

    // ----- Leaf cell data (multi-level father chain) -----------------------
    leafCellToPoint.resize(numLeafCells);
    std::vector<int> leafGlobalCell(numLeafCells);
    std::vector<std::array<int,2>> leafChildToParent(numLeafCells, std::array<int,2>{-1, -1});
    std::vector<int> leafIdxInParent(numLeafCells, -1);
    std::vector<std::array<int,2>> leafToLevel(numLeafCells);

    // Cartesian index of a refined cell's GLOBAL (level-0) ancestor.
    const auto globalAncestorCart = [&](int gridLevel, int cell) -> int {
        int L = gridLevel, c = cell;
        while (L > 0) {
            c = GridStateWriter::childToParent(*boxes[L - 1].level)[c][1];
            L = parentLevelOf[L - 1];
        }
        return level0.globalCell()[c];
    };

    for (int cell = 0; cell < numLeafCells; ++cell) {
        const SourceRef src = leafCells[cell];
        if (src.grid == 0) {
            leafCellToPoint[cell] = GridStateWriter::cellToPoint(level0)[src.index];
            leafGlobalCell[cell] = level0.globalCell()[src.index];
            leafToLevel[cell] = {0, src.index};
        }
        else {
            const BoxData& box = boxes[src.grid - 1];
            const auto& srcPoints = GridStateWriter::cellToPoint(*box.level)[src.index];
            for (int corner = 0; corner < 8; ++corner) {
                leafCellToPoint[cell][corner] = box.cornerToLeaf[srcPoints[corner]];
            }
            leafChildToParent[cell] = { parentLevelOf[src.grid - 1],
                                        GridStateWriter::childToParent(*box.level)[src.index][1] };
            leafIdxInParent[cell] = GridStateWriter::idxInParent(*box.level)[src.index];
            leafGlobalCell[cell] = globalAncestorCart(src.grid, src.index);
            leafToLevel[cell] = {src.grid, src.index};
        }
    }

    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<3,3>>& leafCellGeom =
        leafGeometry.geomVector(std::integral_constant<int,0>());
    leafCellGeom.resize(numLeafCells);
    auto leafCornersPtr = leafGeometry.pointGeometryPtr();
    for (int cell = 0; cell < numLeafCells; ++cell) {
        const SourceRef src = leafCells[cell];
        CpGridData& srcGrid = (src.grid == 0) ? level0 : *boxes[src.grid - 1].level;
        const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<3,3>>& srcCellGeom =
            GridStateWriter::geometry(srcGrid).geomVector(std::integral_constant<int,0>());
        const auto& geom = srcCellGeom[src.index];
        leafCellGeom[cell] = Dune::cpgrid::Geometry<3,3>(geom.center(), geom.volume(),
                                                         leafCornersPtr,
                                                         leafCellToPoint[cell].data());
    }

    // Build the leaf's parallel cell index set, remote indices, partition types
    // and communication interfaces in the rank-interior ("refine-after") model,
    // so the distributed leaf view knows the cells it owns (coarse + all nested
    // refined cells on this rank) and overlap cells take part in collective
    // communication. This mirrors the single-level assembleLeafGrid branch:
    // coarse leaf cells inherit their level-zero global index/attribute; every
    // refined cell (any non-level-0 source, i.e. all nested levels) is rank-
    // interior (owner only) and gets its assignRefinedIds id. In refine-before-redistribute the leaf is
    // assembled with a self-communicator (cc.size() == 1), so this is skipped and
    // the distribution happens later in scatterGrid.
    const auto ids = assignRefinedIds(storage, numBoxes, leafToLevel, leafCornerHistory,
                                      {}, comm);
#if HAVE_MPI
    Dune::Communication<Dune::MPIHelper::MPICommunicator> cc(comm);
    if (cc.size() > 1) {
        using ParallelIndexSet = CpGridData::ParallelIndexSet;
        using AttributeSet = CpGridData::AttributeSet;

        std::vector<int> leafCellGlobalId(numLeafCells, -1);
        std::vector<AttributeSet> leafCellAttr(numLeafCells, AttributeSet::owner);

        for (const auto& entry : level0.cellIndexSet()) {
            const int l0local = entry.local().local();
            const int leafC = leafIdxOf[0][l0local];   // -1 if refined away
            if (leafC >= 0) {
                leafCellGlobalId[leafC] = entry.global();
                leafCellAttr[leafC] = entry.local().attribute();
            }
        }

        for (int cell = 0; cell < numLeafCells; ++cell) {
            if (leafCells[cell].grid != 0) {
                leafCellGlobalId[cell] = ids.cells[cell];
                leafCellAttr[cell] = AttributeSet::owner;
            }
        }

        auto& cis = leaf->cellIndexSet();
        cis.beginResize();
        for (int cell = 0; cell < numLeafCells; ++cell) {
            cis.add(leafCellGlobalId[cell],
                    ParallelIndexSet::LocalIndex(cell, leafCellAttr[cell], true));
        }
        cis.endResize();

        leaf->cellRemoteIndices().template rebuild<false>();
        leaf->computeCellPartitionType();
        leaf->computePointPartitionType();
        leaf->computeCommunicationInterfaces(numLeafCorners);
        setDistributedLeafIds(*leaf, level0, leafCellGlobalId, leafFaces, numLeafFaces,
                              ids.points, ids.end, cc);
    }
    else
#endif
    {
        std::vector<int> faceIds(numLeafFaces);
        std::iota(faceIds.begin(), faceIds.end(), numLeafCells);
        GridStateWriter::setGlobalIdMapping(*leaf, ids.cells, std::move(faceIds), ids.points);
    }

    // ----- Leaf metadata ---------------------------------------------------
    GridStateWriter::setLogicalCartesianSize(*leaf, dims0);
    GridStateWriter::setGlobalCell(*leaf, std::move(leafGlobalCell));
    GridStateWriter::setIndexSet(*leaf, numLeafCells, numLeafCorners);
    GridStateWriter::setParentRelations(*leaf, std::move(leafChildToParent), std::move(leafIdxInParent));
    GridStateWriter::setAquiferCells(*leaf, leafAquiferCells(level0, leafToLevel));
    GridStateWriter::setLeafToLevel(*leaf, std::move(leafToLevel));
    GridStateWriter::setCornerHistory(*leaf, std::move(leafCornerHistory));
    GridStateWriter::setRefinementMaxLevel(*leaf, numBoxes);
    for (int b = 0; b < numBoxes; ++b) {
        GridStateWriter::setRefinementMaxLevel(*boxes[b].level, numBoxes);
    }
    // Each parent grid learns about its children (level0 -> top boxes;
    // a parent box -> its nested children).
    GridStateWriter::setParentToChildren(level0, std::move(parentToChildren0));
    for (int b = 0; b < numBoxes; ++b) {
        const int pl = parentLevelOf[b];
        if (pl == 0) {
            continue;  // already recorded on level0 above
        }
        CpGridData& pgrid = *boxes[parentBoxOf[b]].level;
        std::vector<std::tuple<int,std::vector<int>>> p2c(pgrid.size(0),
                                                          std::make_tuple(-1, std::vector<int>{}));
        for (const auto& [parent, kids] : boxes[b].childrenOfParent) {
            p2c[parent] = std::make_tuple(b + 1, kids);
        }
        GridStateWriter::setParentToChildren(pgrid, std::move(p2c));
    }

    return leaf;
}

} // namespace Refinement
} // namespace Opm
