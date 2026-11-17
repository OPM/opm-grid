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

#include <opm/grid/cpgrid/refinement/conforming/LeafGridAssembler.hpp>

#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/cpgrid/EntityRep.hpp>
#include <opm/grid/cpgrid/Geometry.hpp>
#include <opm/grid/cpgrid/LgrHelpers.hpp>
#include <opm/grid/cpgrid/refinement/conforming/FaultedBoundaryFaces.hpp>
#include <opm/grid/cpgrid/refinement/conforming/GrdeclRefinement.hpp>
#include <opm/grid/cpgrid/refinement/GridStateWriter.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <dune/common/dotproduct.hh>
#include <dune/common/fmatrixev.hh>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <map>
#include <optional>
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

// Terms. A box is one refinement request: the cells of its parent grid with
// startIJK <= (i,j,k) < endIJK. These are its host cells, or parents, and each is split
// into child cells. The parent grid is level zero, or, for a nested box, the level grid
// of another box, its parent LGR. The level grid of a box holds its children. The leaf
// grid holds the cells of level zero that no box refines, and the children in place of
// their hosts. A box side is the boundary of the box in one direction (axis) and sense
// (side, -1 or +1). A fault offsets the cells on its two sides; the offset is the throw.
// A pinch-out joins the cells above and below collapsed layers. A side is faulted where
// the face of a host on it is not the host's whole face toward one cell: several faces,
// or a face that covers part of the host's side only (a fault, or a pinch-out). The
// corner-point processing of a faulted side cuts it into pieces, each between one box
// cell and one cell across, or the domain. On a side that is not faulted, the cell
// across a host gets the faces of the host's children on the side in place of its face
// toward the host: a mosaic.

// The leaf cells that are numerical aquifer cells. No box holds one, so each is a
// level-zero cell of the leaf: opm-common's EclipseState refuses a numerical aquifer
// cell inside a refinement box, and refuseDeckConnectionsInsideBoxes, which the
// simulator calls before refining, refuses an aquifer connection onto a box cell.
std::vector<int> leafAquiferCells(const CpGridData& level0,
                                  const std::vector<std::array<int,2>>& leafToLevel)
{
    // Refinement keeps the order of the level-zero cells (children replace their parent in
    // place), so the sorted aquifer cells are found in one forward pass.
    const auto& aquifer0 = level0.sortedNumAquiferCells();
    std::vector<int> leaf;
    leaf.reserve(aquifer0.size());
    auto it = leafToLevel.begin();
    for (const int aquiferCell : aquifer0) {
        it = std::find(it, leafToLevel.end(), std::array<int,2>{0, aquiferCell});
        if (it == leafToLevel.end()) {
            OPM_THROW(std::logic_error, "Numerical aquifer cell " + std::to_string(aquiferCell)
                      + " of level zero is not a cell of the leaf grid.");
        }
        leaf.push_back(static_cast<int>(it - leafToLevel.begin()));
    }
    return leaf;
}

using SourceRef = Opm::Refinement::LeafGridAssembler::SourceRef;

// A grid that *either* is left-handed *or* has its K axis reversed
// (ZCORN decreasing with K) will have a different ordering of its corners
// than a right-handed grid with a standard K axis. We will call such a grid
// "flipped". (If both, then it cancels out and the ordering is the normal one.)
// For normal right-handed grids with a standard K axis, the corner ordering is
// the lexicographic order of the I, J, and K indices, I running fastest.
// For a flipped grid, the corner ordering can be found by XORing the number of
// the regular corner ordering with a fixed mask that flips the J axis bit.
//
// Note that this function does not look at geometry, yet handedness or K axis reversal
// comes from the embedding of the grid in R^3, not from the topology by itself. This is
// possible because the processing of the cornerpoint data into a (level) grid takes this
// into account in order to make the corners appear in the correct order (to e.g. ensure
// positive volumes from straight-forward volume calculations) and to make face-vertices
// ordered so that a right-hand rule is satisfied for the face normals.
//
// The mask produced by this function will be used as follows:
//     base_corneridx = I + 2*J + 4 * K;
//     possibly_corrected_corneridx = base_corneridx ^ mask;
// The mask is 0 for an unflipped grid and 2 for a flipped grid. All three axes
// are checked: I and K must retain their logical order, and J must have the same
// ordering throughout the grid.
int cornerSlotMask(const CpGridData& grid)
{
    const auto& c2p = grid.cellToPoint();
    int mask = 0;
    // An empty grid has no corners to translate (e.g. on an empty MPI rank).
    if (c2p.empty()) {
        return mask;
    }
    // The positive-side face is identified by its axis tag and positive cell incidence,
    // independently of corner slots. Its four corners must occupy slots whose
    // axis bit is either all set (logical order) or all clear (flipped order).
    // One usable cell determines the bit: preprocessing applies one permutation
    // to the whole grid. Check all other usable cells for contradictions.
    constexpr std::array<face_tag,3> axisTags = {I_FACE, J_FACE, K_FACE};
    for (int axis = 0; axis < 3; ++axis) {
        const int bit = 1 << axis;
        const enum face_tag tag = axisTags[axis];
        std::optional<bool> axisFlipped;
        for (int cell = 0; cell < static_cast<int>(c2p.size()); ++cell) {
            const auto row = grid.cellToFace(cell);
            int highFace = -1;
            int highFaces = 0;
            for (int e = 0; e < row.size(); ++e) {
                if ((grid.faceTag(row[e].index()) == tag) && row[e].orientation()) {
                    highFace = row[e].index();
                    ++highFaces;
                }
            }
            if (highFaces != 1) {
                continue;   // A split face cannot establish the corner permutation.
            }
            const auto points = grid.faceToPoint(highFace);
            if (points.size() != 4) {
                continue;
            }
            int withBit = 0;
            int withoutBit = 0;
            for (const int point : points) {
                const auto slot = std::find(c2p[cell].begin(), c2p[cell].end(), point);
                if (slot == c2p[cell].end()) {
                    break;
                }
                // The distance is this point's slot number (0..7). Test its
                // I, J, or K bit and count whether that bit is set or clear.
                // On a positive-side face, four set bits mean logical order;
                // four clear bits mean the slot order is flipped on this axis.
                ((std::distance(c2p[cell].begin(), slot) & bit) != 0) ? ++withBit : ++withoutBit;
            }
            // Partial or degenerate faces need not contain the four logical
            // side corners, so they cannot establish the permutation either.
            if ((withBit != 4) && (withoutBit != 4)) {
                continue;
            }
            const bool flipped = (withoutBit == 4);
            if (axisFlipped.has_value() && (*axisFlipped != flipped)) {
                OPM_THROW(std::logic_error, "Inconsistent cell corner ordering along axis "
                          + std::to_string(axis) + " at cell " + std::to_string(cell) + ".");
            }
            axisFlipped = flipped;
        }
        if (!axisFlipped.has_value()) {
            OPM_THROW(std::logic_error, "Cannot determine cell corner ordering along axis "
                      + std::to_string(axis) + ": no usable whole face.");
        }
        if (*axisFlipped) {
            if (axis != 1) {
                OPM_THROW(std::logic_error, "Unexpected cell corner flip along axis "
                          + std::to_string(axis) + ": only the J axis may be flipped.");
            }
            mask |= bit;
        }
    }
    return mask;
}

#if HAVE_MPI
// Sets the global ids of a distributed leaf: @p cellIds for the cells, and ids for the
// faces and points. The local id set numbers entities by their local index, which
// differs between ranks. A face or point of level zero keeps its level-zero global id.
// Every other face and point is new: each rank counts its new ones and gives them the
// next block of ids above every id in use, after the blocks of the ranks below it. No
// new entity gets two ids: a box and the cells next to it lie on one rank, since
// classifyBoxPresence refuses a box that any rank holds as overlap, so every refined
// face and point exists on exactly one rank.
void setDistributedLeafIds(CpGridData& leaf, const CpGridData& level0,
                           std::vector<int> cellIds,
                           const std::vector<SourceRef>& leafFaces,
                           int numLeafFaces, int numLeafCorners,
                           const Dune::Communication<Dune::MPIHelper::MPICommunicator>& cc)
{
    const auto& ids0 = level0.globalIdSet();
    const auto& faceIds0 = ids0.template getMapping<1>();
    const auto& pointIds0 = ids0.template getMapping<3>();
    const int numCorners0 = level0.size(3);

    // The largest id in use: level zero's (over cells and points; every face id lies
    // below every point id), and the leaf's cell ids.
    int maxUsed = static_cast<int>(ids0.getMaxGlobalId());
    if (!cellIds.empty()) {
        maxUsed = std::max(maxUsed, std::ranges::max(cellIds));
    }
    maxUsed = cc.max(maxUsed);

    std::vector<int> faceIds(numLeafFaces, -1);
    std::vector<int> pointIds(numLeafCorners, -1);
    int numNew = 0;
    // numLeafFaces counts the faces of leafFaces, taken from level zero and the level
    // grids, and after them the faces the leaf assembly builds on box sides: the split
    // faces where a fault crosses a side, and the faces where a box's MINPV removes a
    // child on a side. These built faces are new and get new ids.
    for (int f = 0; f < numLeafFaces; ++f) {
        const bool coarse = f < static_cast<int>(leafFaces.size()) && leafFaces[f].grid == 0;
        if (coarse) {
            faceIds[f] = faceIds0[leafFaces[f].index];
        } else {
            ++numNew;
        }
    }
    // The leaf stores level zero's corners first, at the same indices, so p < numCorners0
    // is a level-zero corner.
    for (int p = 0; p < numLeafCorners; ++p) {
        if (p < numCorners0) {
            pointIds[p] = pointIds0[p];
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
    for (auto& id : pointIds) {
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
        OPM_THROW(std::logic_error, "Unexpected face tag in refined level grid.");
    }
}

// The faces of a cell on one side along an axis. orientation() is true for a face toward
// +axis; side is -1 or +1, so it is compared as side > 0 (a bool compared with -1 never
// matches).
std::vector<int> facesOnSide(const CpGridData& grid, int cell, int axis, int side)
{
    std::vector<int> faces;
    const auto row = grid.cellToFace(cell);
    for (int e = 0; e < row.size(); ++e) {
        const enum face_tag tag = grid.faceTag(row[e].index());
        if ((tag != NNC_FACE) && (axisOf(tag) == axis) && (row[e].orientation() == (side > 0))) {
            faces.push_back(row[e].index());
        }
    }
    return faces;
}

// The four corners of a cell on one side (side -1 or +1) along an axis; mask is the
// grid's cornerSlotMask.
std::set<int> cornersOnSide(const CpGridData& grid, int cell, int mask, int axis, int side)
{
    std::set<int> corners;
    for (int corner = 0; corner < 8; ++corner) {
        if ((((corner >> axis) & 1) == 1) == (side > 0)) {
            corners.insert(grid.cellToPoint(cell)[corner ^ mask]);
        }
    }
    return corners;
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

// Area-weighted normal (Newell), centroid and area of a polygon that is close to planar.
struct PolyGeom { Dune::FieldVector<double,3> center; Dune::FieldVector<double,3> normal; double area; };
PolyGeom polygonGeometry(const std::vector<Dune::FieldVector<double,3>>& pts)
{
    PolyGeom g;
    g.center = 0.0;
    g.normal = 0.0;
    const int n = static_cast<int>(pts.size());
    Dune::FieldVector<double,3> mean(0.0);
    for (int i = 0; i < n; ++i) {
        const auto& a = pts[i];
        const auto& b = pts[(i + 1) % n];
        g.normal[0] += (a[1] - b[1]) * (a[2] + b[2]);
        g.normal[1] += (a[2] - b[2]) * (a[0] + b[0]);
        g.normal[2] += (a[0] - b[0]) * (a[1] + b[1]);
        mean += a;
    }
    mean /= (n > 0 ? n : 1);
    g.area = 0.5 * g.normal.two_norm();
    if (g.area > 0.0) g.normal /= g.normal.two_norm();  // unit
    // Centroid: the triangles from the vertex mean, weighted by their area along the normal.
    double weight = 0.0;
    for (int i = 0; i < n; ++i) {
        const auto& a = pts[i];
        const auto& b = pts[(i + 1) % n];
        const auto e1 = a - mean;
        const auto e2 = b - mean;
        const double w = 0.5 * Dune::dot(Dune::FMatrixHelp::Impl::crossProduct(e1, e2), g.normal);
        auto centre = mean;
        centre += a;
        centre += b;
        centre *= w / 3.0;
        g.center += centre;
        weight += w;
    }
    if (weight > 0.0) {
        g.center /= weight;
    }
    else {
        g.center = mean;
    }
    return g;
}

// The outline of polygons that tile one planar region, all wound the same way: edges two
// polygons share cancel, and what remains must close into a single loop.
std::optional<std::vector<int>> mergedOutline(const std::vector<std::vector<int>>& polygons)
{
    std::map<std::pair<int,int>, int> edges;
    for (const auto& polygon : polygons) {
        const int n = static_cast<int>(polygon.size());
        for (int i = 0; i < n; ++i) {
            const int a = polygon[i];
            const int b = polygon[(i + 1) % n];
            // Two consecutive corners of a piece can be the same leaf corner: where two
            // lines of the faulted side cross very close to a pillar, the computed crossing
            // point can round onto that pillar's corner. Such an edge has zero length and
            // is not part of the outline.
            if (a == b) {
                continue;
            }
            const auto reverse = edges.find({b, a});
            if (reverse != edges.end()) {
                if (--reverse->second == 0) {
                    edges.erase(reverse);
                }
            }
            else {
                ++edges[{a, b}];
            }
        }
    }
    std::map<int,int> next;
    for (const auto& [edge, count] : edges) {
        if ((count != 1) || !next.emplace(edge.first, edge.second).second) {
            return std::nullopt;
        }
    }
    if (next.empty()) {
        return std::nullopt;
    }
    std::vector<int> outline;
    int node = next.begin()->first;
    do {
        outline.push_back(node);
        const auto it = next.find(node);
        if ((it == next.end()) || (outline.size() > next.size())) {
            return std::nullopt;
        }
        node = it->second;
    } while (node != outline.front());
    if (outline.size() != next.size()) {
        return std::nullopt;
    }
    return outline;
}

// Points with an index each, found by their coordinates: a point within tol of a stored
// point is that point. Two computations of one point, such as a corner of a box and the
// same corner computed by the box across a throw, can differ by round-off. tol is 1e-9
// times the extent of the given corners (their largest span along x, y or z), at least
// 1e-9: far below the size of a cell, and above the round-off of the coordinates. The
// points are kept in cubic buckets of size tol, so a point within tol of x lies in the
// bucket of x or in one of its 26 neighbours.
class PointLookup
{
public:
    explicit PointLookup(const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& corners)
    {
        double extent = 1.0;
        if (!corners.empty()) {
            auto lo = corners[0].center();
            auto hi = lo;
            for (const auto& corner : corners) {
                for (int d = 0; d < 3; ++d) {
                    lo[d] = std::min(lo[d], corner.center()[d]);
                    hi[d] = std::max(hi[d], corner.center()[d]);
                }
            }
            extent = std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2], extent});
        }
        tol_ = 1e-9 * extent;
    }

    // The index of a stored point within tol of x, or -1.
    int find(const Dune::FieldVector<double,3>& x) const
    {
        const auto bucket = bucketOf(x);
        for (long long di = -1; di <= 1; ++di) {
            for (long long dj = -1; dj <= 1; ++dj) {
                for (long long dk = -1; dk <= 1; ++dk) {
                    const auto it = buckets_.find({bucket[0] + di, bucket[1] + dj, bucket[2] + dk});
                    if (it == buckets_.end()) {
                        continue;
                    }
                    for (const auto& [point, index] : it->second) {
                        if ((point - x).two_norm() <= tol_) {
                            return index;
                        }
                    }
                }
            }
        }
        return -1;
    }

    void add(const Dune::FieldVector<double,3>& x, int index)
    {
        buckets_[bucketOf(x)].emplace_back(x, index);
    }

private:
    std::array<long long,3> bucketOf(const Dune::FieldVector<double,3>& x) const
    {
        return { static_cast<long long>(std::floor(x[0] / tol_)),
                 static_cast<long long>(std::floor(x[1] / tol_)),
                 static_cast<long long>(std::floor(x[2] / tol_)) };
    }

    double tol_;
    std::map<std::array<long long,3>, std::vector<std::pair<Dune::FieldVector<double,3>, int>>> buckets_;
};

} // anonymous namespace

namespace Opm
{
namespace Refinement
{

LeafGridAssembler::LeafGridAssembler(const std::vector<std::shared_ptr<CpGridData>>& storage,
                                     const std::vector<BlockRefinement>& requests,
                                     const std::array<int,3>& parentDims,
                                     const double* coord,
                                     const double* zcorn,
                                     const int* actnum,
                                     Dune::MPIHelper::MPICommunicator comm)
    : storage_(storage)
    , requests_(requests)
    , parentDims_(parentDims)
    , coord_(coord)
    , zcorn_(zcorn)
    , actnum_(actnum)
    , comm_(comm)
    , level0_(*storage[0])
    , numBoxes_(static_cast<int>(requests.size()))
    , dims0_(level0_.logicalCartesianSize())
    , numCells0_(level0_.size(0))
    , numCorners0_(level0_.size(3))
    , parentBoxOf_(computeParentBoxOf())
    , compressed0_(computeCompressed0())
    , boxOfCell_(computeBoxOfCell())
    , cellToPoint0_(level0_.cellToPoint())
    , faceToCell0_(GridStateWriter::faceToCell(level0_))
    , faceToPoint0_(GridStateWriter::faceToPoint(level0_))
    , faceTag0_(GridStateWriter::faceTag(level0_))
    , faceNormals0_(GridStateWriter::faceNormals(level0_))
    , numFaces0_(faceToCell0_.size())
{
}

std::shared_ptr<CpGridData> LeafGridAssembler::build()
{
    // The stages fill the members once; a second run would add to them.
    assert(!built_);
    built_ = true;

    // The stages before buildParallelStructures do no collective operation, and they
    // hold every refusal of the assembly.
    const auto localStages = [this] {
        collectBoxData();
        buildCornerPool();
        buildLeafCells();
        refuseThrowsLeavingBoxes();
        findDroppedFaces();
        findFaultedSides();
        buildLeafFaces();
        buildSyntheticFaces();
        buildLeafTopology();
        computeCellData();
    };
    const Dune::Communication<Dune::MPIHelper::MPICommunicator> cc(comm_);
    if (cc.size() > 1) {
        // In a distributed run a refusal fires only on the rank that holds the box. The
        // other ranks would go on into the collective operations of
        // buildParallelStructures and wait there forever. So an error on this rank is
        // caught, one reduction tells every rank the lowest rank that failed, and then
        // every rank throws. The lowest rank that failed throws with its own message and
        // logs it; the other ranks throw without logging.
        bool failed = false;
        std::string message;
        try {
            localStages();
        } catch (const std::exception& e) {
            failed = true;
            message = e.what();
        }
        // The lowest rank that failed, or the number of ranks if none did.
        const int failedRank = cc.min(failed ? cc.rank() : cc.size());
        if (failedRank < cc.size()) {
            if (cc.rank() == failedRank) {
                OPM_THROW(std::runtime_error, "Parallel LGR leaf assembly failed on rank "
                          + std::to_string(cc.rank()) + ": " + message);
            }
            OPM_THROW_NOLOG(std::runtime_error, "Parallel LGR leaf assembly aborted on rank "
                            + std::to_string(cc.rank()) + " because rank "
                            + std::to_string(failedRank) + " failed (see that rank's message).");
        }
    } else {
        localStages();
    }
    buildParallelStructures();
    recordLeafMetadata();

    return leaf_;
}

std::vector<int> LeafGridAssembler::computeParentBoxOf() const
{
    // ----- Box tree ---------------------------------------------------------
    // A nested box refines cells of another box's level grid (its parent LGR),
    // which comes before it in the requests.
    std::map<std::string,int> boxOfName;
    std::vector<int> parentBoxOf(numBoxes_, -1);
    for (int b = 0; b < numBoxes_; ++b) {
        boxOfName[requests_[b].name] = b;
        if (requests_[b].parentGridName != "GLOBAL") {
            parentBoxOf[b] = boxOfName.at(requests_[b].parentGridName);
        }
    }
    return parentBoxOf;
}

bool LeafGridAssembler::nested(int b) const
{
    return parentBoxOf_[b] >= 0;
}

std::vector<int> LeafGridAssembler::computeCompressed0() const
{
    // ----- Basic level-zero lookups ------------------------------------
    std::vector<int> compressed0(static_cast<std::size_t>(dims0_[0])*dims0_[1]*dims0_[2], -1);
    for (int c = 0; c < numCells0_; ++c) {
        compressed0[level0_.globalCell()[c]] = c;
    }
    return compressed0;
}

std::vector<int> LeafGridAssembler::computeBoxOfCell() const
{
    const auto boxOfIJK = [&](const std::array<int,3>& ijk) {
        for (int b = 0; b < numBoxes_; ++b) {
            if (!nested(b) && contains(requests_[b], ijk)) {
                return b;
            }
        }
        return -1;
    };
    // Box membership of each level-zero cell (-1 if unrefined).
    std::vector<int> boxOfCell(numCells0_, -1);
    for (int c = 0; c < numCells0_; ++c) {
        std::array<int,3> ijk;
        level0_.getIJK(c, ijk);
        boxOfCell[c] = boxOfIJK(ijk);
    }
    return boxOfCell;
}

void LeafGridAssembler::collectBoxData()
{
    // ----- Per-box data -------------------------------------------------
    boxes_.resize(numBoxes_);
    mask0_ = cornerSlotMask(level0_);

    for (int b = 0; b < numBoxes_; ++b) {
        BoxData& box = boxes_[b];
        box.level = storage_[b + 1].get();
        box.factors = requests_[b].cellsPerDim;
        for (int d = 0; d < 3; ++d) {
            box.subs[d] = axisSubdivision(requests_[b], d);
            box.pos[d] = axisPositions(box.subs[d]);
            box.graded = box.graded || !requests_[b].subdivision[d].empty();
        }
        box.refinedDims = box.level->logicalCartesianSize();

        const auto& childToParent = GridStateWriter::childToParent(*box.level);
        const auto& cartesianIndexInParentCell = GridStateWriter::cartesianIndexInParentCell(*box.level);
        for (int cell = 0; cell < box.level->size(0); ++cell) {
            box.childrenOfParent[childToParent[cell][1]].push_back(cell);
        }
        for (auto& [parent, kids] : box.childrenOfParent) {
            std::sort(kids.begin(), kids.end(), [&](int i, int j) {
                return cartesianIndexInParentCell[i] < cartesianIndexInParentCell[j];
            });
        }

        // The corners of the level grid that are corners of the parent grid, found
        // through each child's host: across a fault inside the box, one lattice
        // position has several corners, so the position alone does not name one. A
        // child's corner is a corner of its host when it lies on the host's boundary
        // in all three directions. The child's place in its host comes from its
        // lattice position, not from cartesianIndexInParentCell, whose strides are the
        // host's own numbers of children once the box is graded.
        box.cornerEquiv.assign(box.level->size(3), -1);
        const auto& cellToPointL = box.level->cellToPoint();
        const int maskL = cornerSlotMask(*box.level);
        box.cornerMask = maskL;
        const auto& cellToPointP = nested(b) ? boxes_[parentBoxOf_[b]].level->cellToPoint()
                                             : cellToPoint0_;
        const int maskP = nested(b) ? boxes_[parentBoxOf_[b]].cornerMask : mask0_;
        for (int cell = 0; cell < box.level->size(0); ++cell) {
            const int parent = childToParent[cell][1];
            std::array<int,3> lattice;
            box.level->getIJK(cell, lattice);
            for (int corner = 0; corner < 8; ++corner) {
                // The corner's lattice offset from the cell's lowest corner: 0 or 1 along
                // I, J and K (bits 0, 1 and 2 of corner).
                const std::array<int,3> latticeOffset = { corner & 1, (corner >> 1) & 1, (corner >> 2) & 1 };
                std::array<int,3> parentCornerIJK{};
                bool onParentCorner = true;
                for (int axis = 0; axis < 3 && onParentCorner; ++axis) {
                    const int edge = box.pos[axis].subIndex[lattice[axis]] + latticeOffset[axis];
                    if (edge == 0) {
                        parentCornerIJK[axis] = 0;
                    }
                    else if (edge == box.pos[axis].parentCount[lattice[axis]]) {
                        parentCornerIJK[axis] = 1;
                    }
                    else {
                        onParentCorner = false;
                    }
                }
                if (onParentCorner) {
                    const int parentCorner = cellToPointP[parent][(parentCornerIJK[0]
                                                                  + 2*parentCornerIJK[1]
                                                                  + 4*parentCornerIJK[2]) ^ maskP];
                    box.cornerEquiv[cellToPointL[cell][corner ^ maskL]] = parentCorner;
                }
            }
        }
    }

    // A nested box refines a top-level LGR and stays at least one cell inside it.
    // checkConformingNesting(), which CpGrid::addLgrsUpdateLeafView() calls before the level
    // grids are built, refuses other requests; the assembly below relies on it.
#ifndef NDEBUG
    for (int b = 0; b < numBoxes_; ++b) {
        if (nested(b)) {
            assert(!nested(parentBoxOf_[b]));
            for (int d = 0; d < 3; ++d) {
                assert((requests_[b].startIJK[d] >= 1)
                       && (requests_[b].endIJK[d] <= boxes_[parentBoxOf_[b]].refinedDims[d] - 1));
            }
        }
    }
#endif
    // nestedBoxOf_[pb][cell]: for a cell of the level grid of the parent box pb, the index
    // of the nested (child) box that refines it, or -1. A cell whose children the nested
    // box's MINPV removes all has that box's index too: it is refined into nothing.
    nestedBoxOf_.resize(numBoxes_);
    for (int b = 0; b < numBoxes_; ++b) {
        nestedBoxOf_[b].assign(boxes_[b].level->size(0), -1);
    }
    for (int b = 0; b < numBoxes_; ++b) {
        if (nested(b)) {
            const CpGridData& parentLgr = *boxes_[parentBoxOf_[b]].level;
            for (int cell = 0; cell < parentLgr.size(0); ++cell) {
                std::array<int,3> ijk;
                parentLgr.getIJK(cell, ijk);
                if (contains(requests_[b], ijk)) {
                    nestedBoxOf_[parentBoxOf_[b]][cell] = b;
                }
            }
        }
    }
}

void LeafGridAssembler::buildCornerPool()
{
    // ----- Corner pool ---------------------------------------------------
    // Pool = all level-zero corners (indices preserved) + new refined corners.
    leaf_ = std::make_shared<CpGridData>(comm_, storage_);
    auto& leafGeometry = GridStateWriter::mutableGeometry(*leaf_);
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& leafCorners =
        leafGeometry.geomVector(std::integral_constant<int,3>());
    const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& corners0 =
        GridStateWriter::geometry(level0_).geomVector(std::integral_constant<int,3>());
    leafCorners.assign(corners0.begin(), corners0.end());

    leafCornerOrigin_.resize(numCorners0_);
    for (int corner = 0; corner < numCorners0_; ++corner) {
        leafCornerOrigin_[corner] = {0, corner};
    }

    // The refined corners in the leaf, by their coordinates. A corner on the face or
    // edge two touching boxes share has equal coordinates in both boxes: both refine the
    // same parent description with uniform subdivisions, which put the corner at the
    // same fraction of the parent cell (graded boxes that touch are refused). The lookup
    // makes it one leaf corner. Boxes do not overlap, so refined corners of two boxes
    // are equal only on what the boxes share. A refined corner can also lie on a
    // level-zero point of its host's faces that is not a corner of the host, such as a
    // corner of the cell across a fault: it is then that point. So the lookup starts with
    // the points of the faces of every host cell.
    PointLookup cornerPool(corners0);
    for (int c = 0; c < numCells0_; ++c) {
        if (boxOfCell_[c] < 0) {
            continue;
        }
        const auto row = level0_.cellToFace(c);
        for (int e = 0; e < row.size(); ++e) {
            for (const int p : faceToPoint0_[row[e].index()]) {
                if (cornerPool.find(corners0[p].center()) < 0) {
                    cornerPool.add(corners0[p].center(), p);
                }
            }
        }
    }

    for (int b = 0; b < numBoxes_; ++b) {
        BoxData& box = boxes_[b];
        const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& cornersL =
            GridStateWriter::geometry(*box.level).geomVector(std::integral_constant<int,3>());
        box.cornerToLeaf.assign(box.level->size(3), -1);
        std::vector<std::array<int,2>> levelCornerOrigin(box.level->size(3), {-1, -1});
        for (int corner = 0; corner < box.level->size(3); ++corner) {
            if (box.cornerEquiv[corner] >= 0) {
                // Level-zero corners keep their index in the leaf; a nested box's
                // parent corners are already in the pool.
                const int leafCorner = nested(b)
                    ? boxes_[parentBoxOf_[b]].cornerToLeaf[box.cornerEquiv[corner]]
                    : box.cornerEquiv[corner];
                box.cornerToLeaf[corner] = leafCorner;
                levelCornerOrigin[corner] = leafCornerOrigin_[leafCorner];
            }
            else {
                const auto x = cornersL[corner].center();
                const int existing = cornerPool.find(x);
                if (existing >= 0) {
                    // A point of a host's faces, or shared with a box processed before:
                    // reuse the corner, with its origin (the grid it comes from).
                    box.cornerToLeaf[corner] = existing;
                    levelCornerOrigin[corner] = leafCornerOrigin_[existing];
                }
                else {
                    const int leafCorner = static_cast<int>(leafCorners.size());
                    box.cornerToLeaf[corner] = leafCorner;
                    leafCorners.push_back(cornersL[corner]);
                    leafCornerOrigin_.push_back({b + 1, corner});
                    cornerPool.add(x, leafCorner);
                }
            }
        }
        GridStateWriter::setCornerOrigin(*box.level, std::move(levelCornerOrigin));
    }
    // May grow below if faulted box boundaries introduce new split vertices.
    numLeafCorners_ = static_cast<int>(leafCorners.size());
}

void LeafGridAssembler::buildLeafCells()
{
    // ----- Leaf cells -----------------------------------------------------
    // Level-zero order, parents replaced by their children (cartesianIndexInParentCell order).
    leafIdxOfCell0_.assign(numCells0_, -1);
    leafIdxOfLevelCell_.resize(numBoxes_);
    for (int b = 0; b < numBoxes_; ++b) {
        leafIdxOfLevelCell_[b].assign(boxes()[b].level->size(0), -1);
    }
    level0ToChildren_.assign(numCells0_, std::make_tuple(-1, std::vector<int>{}));
    for (int c = 0; c < numCells0_; ++c) {
        const int b = boxOfCell_[c];
        if (b < 0) {
            leafIdxOfCell0_[c] = static_cast<int>(leafCells_.size());
            leafCells_.push_back({0, c});
        }
        else {
            auto childrenIt = boxes()[b].childrenOfParent.find(c);
            if (childrenIt == boxes()[b].childrenOfParent.end()) {
                // The box's MINPV removed every child: refined, into nothing.
                level0ToChildren_[c] = std::make_tuple(b + 1, std::vector<int>{});
                continue;
            }
            for (const int child : childrenIt->second) {
                const int nb = nestedBoxOf()[b][child];
                if (nb < 0) {
                    leafIdxOfLevelCell_[b][child] = static_cast<int>(leafCells_.size());
                    leafCells_.push_back({b + 1, child});
                    continue;
                }
                // Refined again by a nested box: its children in its place, none if the
                // nested box's MINPV removed them all. The conforming refinement supports
                // one level of nesting, so these children are not refined again. With
                // more levels, this loop needs to be generalised.
                const auto grandchildren = boxes()[nb].childrenOfParent.find(child);
                if (grandchildren == boxes()[nb].childrenOfParent.end()) {
                    continue;
                }
                for (const int grandchild : grandchildren->second) {
                    leafIdxOfLevelCell_[nb][grandchild] = static_cast<int>(leafCells_.size());
                    leafCells_.push_back({nb + 1, grandchild});
                }
            }
            level0ToChildren_[c] = std::make_tuple(b + 1, childrenIt->second);
        }
    }
    numLeafCells_ = static_cast<int>(leafCells_.size());
    // The parent LGRs of nested boxes learn about their children.
    for (int nb = 0; nb < numBoxes_; ++nb) {
        if (!nested(nb)) {
            continue;
        }
        const int pb = parentBoxOf_[nb];
        // Several nested boxes can refine one LGR: its table is made by the first of them,
        // and each fills in the cells it refines.
        auto [it, fresh] = nestedParentToChildren_.try_emplace(pb);
        if (fresh) {
            it->second.assign(boxes()[pb].level->size(0), std::make_tuple(-1, std::vector<int>{}));
        }
        // Each cell nb refines gets (nb + 1, its children): the level of nb, and the cell's
        // children there, none when nb's MINPV removed them all. The cells no nested box
        // refines keep (-1, {}).
        for (int parent = 0; parent < boxes()[pb].level->size(0); ++parent) {
            if (nestedBoxOf()[pb][parent] == nb) {
                const auto kids = boxes()[nb].childrenOfParent.find(parent);
                it->second[parent] = std::make_tuple(nb + 1, (kids != boxes()[nb].childrenOfParent.end())
                                                                 ? kids->second : std::vector<int>{});
            }
        }
    }
}

// The grid a box refines is level zero, or for a nested box its parent LGR.
int LeafGridAssembler::boxOfParentGridCell(int b, int cell) const
{
    return nested(b) ? nestedBoxOf_[parentBoxOf_[b]][cell] : boxOfCell_[cell];
}

int LeafGridAssembler::leafOfParentGridCell(int b, int cell) const
{
    return nested(b) ? leafIdxOfLevelCell_[parentBoxOf_[b]][cell] : leafIdxOfCell0_[cell];
}

std::array<int,3> LeafGridAssembler::parentGridIJK(int b, int cell) const
{
    const CpGridData& grid = nested(b) ? *boxes_[parentBoxOf_[b]].level : level0_;
    std::array<int,3> ijk;
    grid.getIJK(cell, ijk);
    return ijk;
}

void LeafGridAssembler::refuseThrowsLeavingBoxes() const
{
    // Refuse a connection from a cell inside the box to a cell outside it through a face
    // that is internal to the box's logical IJK extent, not on a box side. Faults and
    // pinch-outs create such connections when the box does not cover every layer they
    // connect. Connections across the box's logical boundary are handled by
    // buildLeafFaces() and buildSyntheticFaces(). Internal ones are not: refining the box
    // leaves an open face, and the leaf assembly removes the outside cell's original face
    // without supplying a replacement.
    const auto text = [](const std::array<int,3>& c) {
        return "(" + std::to_string(c[0] + 1) + "," + std::to_string(c[1] + 1) + ","
            + std::to_string(c[2] + 1) + ")";
    };
    for (int b = 0; b < numBoxes_; ++b) {
        const CpGridData& grid = nested(b) ? *boxes()[parentBoxOf_[b]].level : level0_;
        const Dune::cpgrid::EntityVariableBase<enum face_tag>& tags = GridStateWriter::faceTag(grid);
        for (int cell = 0; cell < grid.size(0); ++cell) {
            if (boxOfParentGridCell(b, cell) != b) {
                continue;
            }
            const auto ijk = parentGridIJK(b, cell);
            const auto row = GridStateWriter::cellToFace(grid)[EntityRep<0>(cell, true)];
            for (int e = 0; e < row.size(); ++e) {
                const int face = row[e].index();
                if (tags.get(face) == NNC_FACE) {
                    continue;
                }
                const int axis = axisOf(tags.get(face));
                // By construction, a cell's faces toward -I, -J and -K have orientation
                // false and those toward +I, +J and +K orientation true, whatever the grid's
                // handedness. If the processing changes this, it must change here as well.
                const bool onSide = row[e].orientation() ? (ijk[axis] == requests_[b].endIJK[axis] - 1)
                                                         : (ijk[axis] == requests_[b].startIJK[axis]);
                const auto cells = GridStateWriter::faceToCell(grid)[EntityRep<1>(face, true)];
                for (int q = 0; (q < cells.size()) && !onSide; ++q) {
                    const int other = cells[q].index();
                    if ((other == cell) || (other == kRemoteCell) || (boxOfParentGridCell(b, other) == b)) {
                        continue;
                    }
                    const auto outside = parentGridIJK(b, other);
                    OPM_THROW(std::invalid_argument, "Refinement box '" + requests_[b].name + "' has a "
                              "fault or a pinch-out inside it that joins its cell " + text(ijk)
                              + " to the cell " + text(outside) + " outside the box (I,J,K from 1 in '"
                              + requests_[b].parentGridName + "'). This is not supported: the box "
                              "must cover every layer the throw joins.");
                }
            }
        }
    }
}

// The cell across a nested box's side: across its parent cell's face there, which must
// be the parent's whole face (a fault or pinch-out on it is not supported), or -1.
// The check reads the parent LGR's level grid, because level zero cannot decide it
// exactly: a throw that lines up with the parent's child layers gives whole faces in
// the parent's level grid, and such a side is built correctly.
int LeafGridAssembler::nestedNeighborOf(int b, int parent, int axis, int side) const
{
    const CpGridData& grid = *boxes_[parentBoxOf_[b]].level;
    const auto faces = facesOnSide(grid, parent, axis, side);
    const auto expected = cornersOnSide(grid, parent, boxes_[parentBoxOf_[b]].cornerMask, axis, side);
    // One face whose points are the parent's four corners on this side is the parent's
    // whole face. Several faces (a fault), one face with other points (a partial face, on
    // a fault or at a pinch-out with a gap), or no face, are not supported.
    std::set<int> actual;
    if (faces.size() == 1) {
        const auto points = GridStateWriter::faceToPoint(grid)[faces[0]];
        actual.insert(points.begin(), points.end());
    }
    if (actual != expected) {
        OPM_THROW(std::logic_error, "Nested refinement '" + requests_[b].name + "' has a "
                  "side on a fault or pinch-out inside its parent LGR '"
                  + requests_[b].parentGridName + "'; that is not supported.");
    }
    const auto cells = GridStateWriter::faceToCell(grid)[EntityRep<1>(faces[0], true)];
    for (int q = 0; q < cells.size(); ++q) {
        if (cells[q].index() != parent) {
            return cells[q].index();
        }
    }
    return -1;
}

void LeafGridAssembler::findDroppedFaces()
{
    // ----- Identify replaced level-zero faces and mosaic pairings ---------
    // A level-zero face is dropped if any adjacent cell is a refined parent.
    // For parent faces toward an unrefined neighbor, the refined boundary
    // faces of the box take over, paired with that neighbor. A parent with no
    // child left leaves the face to the unrefined neighbor, as its boundary face.
    // Level zero only: buildLeafFaces() drops a parent LGR's faces toward the cells
    // its nested boxes refine, by the same rule.
    faceDropped_.assign(numFaces0_, false);
    for (int face = 0; face < numFaces0_; ++face) {
        const auto row = faceToCell0_[EntityRep<1>(face, true)];
        bool refinedSide = false;   // a parent with children
        bool emptySide = false;     // a parent with no child left
        bool coarseSide = false;
        for (int e = 0; e < row.size(); ++e) {
            const int cell = row[e].index();
            // A remote (off-rank) neighbor is never a local box cell.
            if (cell == kRemoteCell) {
                continue;
            }
            if (boxOfCell_[cell] < 0) {
                coarseSide = true;
            }
            else if (std::get<1>(level0ToChildren()[cell]).empty()) {
                emptySide = true;
            }
            else {
                refinedSide = true;
            }
        }
        faceDropped_[face] = refinedSide || (emptySide && !coarseSide);
    }
}

// outsideNeighborOf(parent, axis, side): the level-zero cell across a host cell's
// face on a box side, found once per (parent, axis, side). When the host has one face
// there with its four corners on that side, that is its whole face, and the cell
// across it is returned. Next to a pinch-out, the cell across the collapsed layers is
// returned. With no face there, or no cell across, -1 is returned. Otherwise the side
// is faulted: kFaultedSide is returned, and the faulted-side stage below rebuilds the
// side from the corner-point processing.
int LeafGridAssembler::outsideNeighborOf(int parent, int axis, int side) const
{
    const auto key = std::make_tuple(parent, axis, side);
    auto it = verifiedNeighbor_.find(key);
    if (it != verifiedNeighbor_.end()) {
        return it->second;
    }

    // All faces of the parent in the (axis, side) direction.
    const auto directionFaces = facesOnSide(level0_, parent, axis, side);

    // Next to a pinch-out the parent has two faces on this side: its own face, with
    // no cell across, and the face that joins it to the cell across the collapsed
    // layers.
    if ((axis == 2) && (directionFaces.size() == 2)) {
        int open = 0;
        int pinched = -1;
        std::array<int,3> p;
        level0_.getIJK(parent, p);
        for (const int face : directionFaces) {
            const auto cells = faceToCell0_[EntityRep<1>(face, true)];
            if (cells.size() == 1) {
                ++open;
                continue;
            }
            for (int q = 0; q < cells.size(); ++q) {
                const int cell = cells[q].index();
                if ((cell == parent) || (cell == kRemoteCell)) {
                    continue;
                }
                std::array<int,3> c;
                level0_.getIJK(cell, c);
                if ((c[0] == p[0]) && (c[1] == p[1]) && (std::abs(c[2] - p[2]) > 1)) {
                    pinched = cell;
                }
            }
        }
        if ((open == 1) && (pinched >= 0)) {
            verifiedNeighbor_[key] = pinched;
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
        const auto expected = cornersOnSide(level0_, parent, mask0_, axis, side);
        std::set<int> actual;
        auto facePoints = faceToPoint0_[face];
        for (auto pt = facePoints.begin(); pt != facePoints.end(); ++pt) {
            actual.insert(*pt);
        }
        if (actual != expected) {
            neighbor = kFaultedSide;        // partial (faulted/degenerate) face
        }
        else {
            const auto cells = faceToCell0_[EntityRep<1>(face, true)];
            for (int q = 0; q < cells.size(); ++q) {
                const int cell = cells[q].index();
                if (cell != parent && cell != kRemoteCell) {
                    neighbor = cell;
                }
            }
        }
    }
    // No direction faces: domain boundary, inactive, or off-rank neighbor.

    verifiedNeighbor_[key] = neighbor;
    return neighbor;
}

void LeafGridAssembler::findFaultedSides()
{
    // Identify top-level box sides that require reconstruction by the corner-point
    // processor rather than assembly from existing refined boundary faces. If any
    // host cell has a split or partial face on a side, mark the entire box side in
    // faultedSides_. This includes cells whose fault displacement leaves no neighbour:
    // viewed alone, their faces can look like ordinary domain-boundary faces.
    // buildLeafFaces() uses this set to exclude the existing refined boundary faces;
    // buildSyntheticFaces() reconstructs those sides, including their exposed parts.
    //
    // Record the bounding IJK range of cells across each side in acrossRange_.
    // buildSyntheticFaces() passes this range to boxBoundaryConnections() so the
    // corner-point processor includes the cells connected across the fault, or
    // across a throw of whole layers on a side with a child the box's MINPV removes.
    // The local boxesAcross map records neighbouring boxes and propagates marked
    // sides across box-to-box connections. Both sides then use the reconstruction
    // path; ownership rules in buildSyntheticFaces() prevent duplicate interface faces.
    //
    // Only top-level boxes are examined for faults here, using level-zero connectivity.
    // Nested box boundaries are checked separately by nestedNeighborOf(), called from
    // buildLeafFaces(); it rejects sides on faults or pinch-outs that do not form
    // a whole face of the parent LGR cell, rather than reconstructing those sides.
    // acrossRange_ covers nested boxes too, from their parent LGR's connectivity:
    // a throw of whole layers inside the parent LGR leaves the faces whole.

    // Boxes across each box side: a side shared with a faulted side is rebuilt with it.
    std::map<std::tuple<int,int,int>, std::set<int>> boxesAcross;
    const auto extendAcross = [&](int b, int parent, int axis, int side) {
        const CpGridData& grid = nested(b) ? *boxes()[parentBoxOf_[b]].level : level0_;
        for (const int face : facesOnSide(grid, parent, axis, side)) {
            const auto cells = GridStateWriter::faceToCell(grid)[EntityRep<1>(face, true)];
            for (int q = 0; q < cells.size(); ++q) {
                const int cell = cells[q].index();
                if (cell == parent || cell == kRemoteCell) {
                    continue;
                }
                const std::array<int,3> nijk = parentGridIJK(b, cell);
                const auto key = std::make_tuple(b, axis, side);
                const int nb = boxOfParentGridCell(b, cell);
                if ((nb >= 0) && (nb != b)) {
                    boxesAcross[key].insert(nb);
                }
                auto [it, fresh] = acrossRange_.try_emplace(key);
                for (int d = 0; d < 3; ++d) {
                    it->second.startIJK[d] = fresh ? nijk[d] : std::min(it->second.startIJK[d], nijk[d]);
                    it->second.endIJK[d] = fresh ? nijk[d] + 1 : std::max(it->second.endIJK[d], nijk[d] + 1);
                }
            }
        }
    };
    for (int c = 0; c < numCells0_; ++c) {
        const int b = boxOfCell_[c];
        if (b < 0) {
            continue;
        }
        std::array<int,3> ijk;
        level0_.getIJK(c, ijk);
        for (int axis = 0; axis < 3; ++axis) {
            for (int side = -1; side <= 1; side += 2) {
                const bool onBoundary = (side < 0)
                    ? (ijk[axis] == requests_[b].startIJK[axis])
                    : (ijk[axis] == requests_[b].endIJK[axis] - 1);
                if (!onBoundary) {
                    continue;
                }
                if (outsideNeighborOf(c, axis, side) == kFaultedSide) {
                    faultedSides_.emplace(b, axis, side);
                }
                extendAcross(b, c, axis, side);
            }
        }
    }
    for (int b = 0; b < numBoxes_; ++b) {
        if (!nested(b)) {
            continue;
        }
        const auto& inParent = nestedBoxOf()[parentBoxOf_[b]];
        for (int c = 0; c < static_cast<int>(inParent.size()); ++c) {
            if (inParent[c] != b) {
                continue;
            }
            const std::array<int,3> ijk = parentGridIJK(b, c);
            for (int axis = 0; axis < 3; ++axis) {
                for (int side = -1; side <= 1; side += 2) {
                    if (ijk[axis] == ((side < 0) ? requests_[b].startIJK[axis]
                                                 : requests_[b].endIJK[axis] - 1)) {
                        extendAcross(b, c, axis, side);
                    }
                }
            }
        }
    }
    // A side across a faulted side is rebuilt with it, so both sides of a box-to-box
    // interface take the same path. buildSyntheticFaces() emits each face there once: by
    // the finer box, or by the lower box index when the in-face factors are equal. Each
    // marked side is processed once; newly marked sides are appended to propagate the
    // marking.
    std::vector<std::tuple<int,int,int>> pending(faultedSides_.begin(), faultedSides_.end());
    for (std::size_t i = 0; i < pending.size(); ++i) {
        // Copy the entry: appending to pending can invalidate its references.
        const auto [b, axis, side] = pending[i];
        const auto across = boxesAcross.find(std::make_tuple(b, axis, side));
        if (across == boxesAcross.end()) {
            continue;
        }
        for (const int nb : across->second) {
            if (faultedSides_.emplace(nb, axis, -side).second) {
                pending.emplace_back(nb, axis, -side);
            }
        }
    }
}

int LeafGridAssembler::childAcross(int b, const std::array<int,3>& lattice, int axis, int side,
                                   int nb, int neighbor) const
{
    // The sub-position of b's child within its host, mapped to nb's frame: the in-face
    // components scale by the factor ratio (finer -> coarser); along the axis it is nb's
    // boundary layer facing b.
    const auto& facA = boxes()[b].factors;
    const auto& facB = boxes()[nb].factors;
    std::array<int,3> subB{};
    for (int c = 0; c < 3; ++c) {
        subB[c] = (c == axis) ? ((side > 0) ? 0 : (facB[axis] - 1))
                              : (lattice[c] % facA[c]) * facB[c] / facA[c];
    }
    const std::array<int,3> nijk = parentGridIJK(b, neighbor);
    const auto& startB = requests_[nb].startIJK;
    const auto& rdB = boxes()[nb].refinedDims;
    std::array<int,3> latB{};
    for (int c = 0; c < 3; ++c) {
        latB[c] = (nijk[c] - startB[c]) * facB[c] + subB[c];
    }
    const int cartB = latB[0] + rdB[0]*latB[1] + rdB[0]*rdB[1]*latB[2];
    const auto cit = refinedCartToCell()[nb].find(cartB);
    return (cit != refinedCartToCell()[nb].end()) ? cit->second : -1;
}

void LeafGridAssembler::buildLeafFaces()
{
    // ----- Leaf faces ------------------------------------------------------
    for (int face = 0; face < numFaces0_; ++face) {
        if (!faceDropped()[face]) {
            leafFaces_.push_back({0, face});
            mosaicOutside_.push_back(-1);
        }
    }
    // Face-sharing LGRs: a boundary face of box A toward another refined
    // box B pairs with B's matching boundary face. Both have the same leaf
    // corners (corners already merged); the key is them sorted, so the two
    // orderings match. The box seen first, the one with the lower index, owns
    // the leaf face; the partner maps onto it and patches in its cell as the
    // second side (the same mechanism as a coarse mosaic neighbor, only the
    // outside cell is refined).
    std::map<std::vector<int>, int> sharedBoundaryFaceLeaf; // sorted leaf corners -> owner leaf face

    // Per box: refined (box-global) Cartesian index -> level cell index. Used by
    // the compatible sub-face mosaic below to find, for a finer box's boundary
    // sub-face, the coarser box's child cell it abuts.
    refinedCartToCell_.resize(numBoxes_);
    for (int b = 0; b < numBoxes_; ++b) {
        for (int c = 0; c < boxes()[b].level->size(0); ++c) {
            refinedCartToCell_[b][boxes()[b].level->globalCell()[c]] = c;
        }
    }

    for (int b = 0; b < numBoxes_; ++b) {
        BoxData& box = boxes_[b];
        auto& faceToCellL = GridStateWriter::faceToCell(*box.level);
        auto& faceTagL = GridStateWriter::faceTag(*box.level);
        auto& faceToPointL = GridStateWriter::faceToPoint(*box.level);
        const Dune::cpgrid::EntityVariableBase<enum face_tag>& tagsL = faceTagL;
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
            // Toward a cell a nested box refines: that box's faces take over. A cell with
            // no child left leaves the face to the cell on the other side, as its
            // boundary face. The same rule as findDroppedFaces() for level zero.
            bool towardNested = false;
            bool towardEmpty = false;
            bool towardUnrefined = false;
            for (int e = 0; e < row.size(); ++e) {
                const int cell = row[e].index();
                const int nb = nestedBoxOf()[b][cell];
                if (nb < 0) {
                    towardUnrefined = true;
                }
                else if (boxes()[nb].childrenOfParent.count(cell) == 0) {
                    towardEmpty = true;
                }
                else {
                    towardNested = true;
                }
            }
            if (towardNested || (towardEmpty && !towardUnrefined)) {
                continue;
            }
            if (row.size() == 1) {
                // Boundary face of the refined block: internal hole, domain
                // boundary, block boundary toward a coarse neighbor, or
                // toward another refined box (face-sharing).
                const int cell = row[0].index();
                const bool normalOut = row[0].orientation();
                const int axis = axisOf(tagsL.get(face));
                const int side = normalOut ? 1 : -1;
                const std::array<int,3> lattice =
                    Opm::Lgr::getIJK(box.level->globalCell()[cell], box.refinedDims);
                if (side < 0 ? lattice[axis] == 0 : lattice[axis] == box.refinedDims[axis] - 1) {
                    // Faulted side: suppress this box face; the whole side is
                    // rebuilt below from the corner-point processor.
                    if (faultedSides().count(std::make_tuple(b, axis, side))) {
                        continue;
                    }
                    // Block boundary: find the neighbor parent (or none).
                    const int parent = childToParent[cell][1];
                    const int neighbor = nested(b) ? nestedNeighborOf(b, parent, axis, side)
                                                   : outsideNeighborOf(parent, axis, side);
                    if (neighbor >= 0 && boxOfParentGridCell(b, neighbor) < 0) {
                        // Toward a coarse (unrefined) neighbor: mosaic.
                        outside = leafOfParentGridCell(b, neighbor);
                    }
                    else if (neighbor >= 0) {
                        // Toward another refined box. The two boxes share this
                        // 2-D face; its in-face directions are the two axes
                        // other than `axis` (the perpendicular touch direction).
                        const int nbBox = boxOfParentGridCell(b, neighbor);
                        if (box.graded || boxes()[nbBox].graded) {
                            // checkTouchingBoxes() refuses, before the build, graded
                            // boxes whose index ranges share a face or an edge. Boxes
                            // can also meet across a pinch-out or a throw, where the
                            // cell across is not the index neighbour; those are
                            // refused here. The conformity rule below compares one
                            // subdivision factor per in-face direction, which a
                            // graded box does not have.
                            OPM_THROW(std::invalid_argument,
                                "Refinement boxes '" + requests_[b].name + "' and '"
                                + requests_[nbBox].name + "' touch and at least one is "
                                "graded. Touching graded boxes are not supported.");
                        }
                        const auto& facA = box.factors;
                        const auto& facB = boxes()[nbBox].factors;
                        const int u = (axis == 0) ? 1 : 0;
                        const int v = (axis == 2) ? 1 : 2;
                        if (facA[u] == facB[u] && facA[v] == facB[v]) {
                            // Matching in-face subdivisions: 1:1 face pairing by
                            // corner set (both sides have identical sub-faces). A
                            // partner the other box's MINPV removed leaves the face
                            // to this cell, as a boundary face.
                            if (childAcross(b, lattice, axis, side, nbBox, neighbor) >= 0) {
                                const auto key = mergedCornerSet(face);
                                const auto it = sharedBoundaryFaceLeaf.find(key);
                                const int thisCellLeaf = leafIdxOfLevelCell()[b][cell];
                                if (it != sharedBoundaryFaceLeaf.end()) {
                                    // Partner: map onto the owner's face, supply our
                                    // cell as its second side, do not emit a face.
                                    box.faceToLeaf[face] = it->second;
                                    mosaicOutside_[it->second] = thisCellLeaf;
                                    continue;
                                }
                                // Owner: emits the face; the partner adds its cell to it.
                                sharedBoundaryFaceLeaf.emplace(key, static_cast<int>(leafFaces_.size()));
                            }
                        }
                        else {
                            // Different factors, one box finer in both in-face
                            // directions (checkTouchingBoxes refuses other pairs).
                            // The finer side emits every sub-face, with the
                            // coarser box's child that covers it as the outside
                            // cell, or none where its MINPV removed that child. The
                            // coarser side emits no face there, so each of its
                            // children on the face has the finer sub-faces in place
                            // of its face: a mosaic, with more than six faces per
                            // cell. Where the finer box's MINPV removes a child,
                            // buildSyntheticFaces() gives the coarser child a
                            // boundary face in its place; where it removes every
                            // child of the host across, the coarser side keeps its
                            // faces, as boundary faces.
                            const bool thisFiner = (facA[u] >= facB[u] && facA[v] >= facB[v]);
                            if (thisFiner) {
                                const int across = childAcross(b, lattice, axis, side, nbBox, neighbor);
                                if (across >= 0) {
                                    outside = leafIdxOfLevelCell()[nbBox][across];
                                }
                            }
                            else if (boxes()[nbBox].childrenOfParent.count(neighbor) > 0) {
                                continue;       // coarser side: the finer references us
                            }
                        }
                    }
                }
                // else: internal hole boundary -> stays a boundary face.
            }
            box.faceToLeaf[face] = static_cast<int>(leafFaces_.size());
            leafFaces_.push_back({b + 1, face});
            mosaicOutside_.push_back(outside);
        }
    }

    // Every face-sharing owner face must have been claimed by its partner box: an
    // owner is recorded only where the partner box has the child across. An
    // unclaimed one means the shared-face corners did not merge exactly in the
    // corner pool, and the box<->box interface would silently seal.
    for (const auto& entry : sharedBoundaryFaceLeaf) {
        const int leafFace = entry.second;
        if (mosaicOutside_[leafFace] < 0) {
            const int ownerBox = leafFaces_[leafFace].grid - 1;
            OPM_THROW(std::logic_error, "Internal error: refined box '"
                                   + requests_[ownerBox].name
                                   + "' has an interface face toward a face-sharing "
                                     "refined box with no matching partner face. "
                                     "The shared boundary corners failed to merge "
                                     "exactly; the interface would be sealed.");
        }
    }
    numSourceFaces_ = static_cast<int>(leafFaces_.size());
}

void LeafGridAssembler::buildSyntheticFaces()
{
    // The cell at a Cartesian index of the grid box b refines, or -1 if that grid does not
    // have it.
    const auto parentCellOfCart = [&](int b, int cart) {
        if (!nested(b)) {
            return compressed0_[cart];
        }
        const auto& cells = refinedCartToCell()[parentBoxOf_[b]];
        const auto it = cells.find(cart);
        return (it != cells.end()) ? it->second : -1;
    };
    const auto parentGridOf = [&](int b) -> const CpGridData& {
        return nested(b) ? *boxes()[parentBoxOf_[b]].level : level0_;
    };
    // The host of the child at a lattice position of box b when the box's MINPV removed
    // that child and the host is in the grid, else -1.
    const auto hostOfRemovedChild = [&](int b, const std::array<int,3>& lattice) {
        const auto& removed = requests_[b].minpvRemoved;
        const auto& rd = boxes()[b].refinedDims;
        // minpvRemoved: 1 = removed, 0 = stays on the leaf; empty when the box set no MINPV.
        const auto child = lattice[0] + rd[0]*(lattice[1] + rd[1]*lattice[2]);
        if (removed.empty() || (removed[child] == 0)) {
            return -1;
        }
        const auto& dims = nested(b) ? boxes()[parentBoxOf_[b]].refinedDims : dims0_;
        int cart = 0;
        for (int d = 2; d >= 0; --d) {
            cart = cart*dims[d] + requests_[b].startIJK[d] + boxes()[b].subs[d].parentOffset[lattice[d]];
        }
        return parentCellOfCart(b, cart);
    };
    // Box sides the fault does not cross with a child the box's MINPV removed, of a host
    // that keeps other children: the cell across has lost its face toward the host, and
    // keeps a boundary face where the child would be.
    std::set<std::tuple<int,int,int>> removedChildSides;
    for (int b = 0; b < numBoxes_; ++b) {
        const auto& rd = boxes()[b].refinedDims;
        for (std::size_t child = 0; child < requests_[b].minpvRemoved.size(); ++child) {
            const std::array<int,3> lattice = Opm::Lgr::getIJK(static_cast<int>(child), rd);
            const int host = hostOfRemovedChild(b, lattice);
            if ((host < 0) || (boxes()[b].childrenOfParent.count(host) == 0)) {
                continue;
            }
            for (int axis = 0; axis < 3; ++axis) {
                for (const int side : {-1, 1}) {
                    const bool onSide = (lattice[axis] == ((side < 0) ? 0 : rd[axis] - 1));
                    if (onSide && !faultedSides().count(std::make_tuple(b, axis, side))) {
                        removedChildSides.emplace(b, axis, side);
                    }
                }
            }
        }
    }

    // ----- Faulted box boundaries: synthetic split faces -------------------
    // For each faulted (box, axis, side), the
    // corner-point processor (boxBoundaryConnections) gives the split
    // connections between the refined boundary cells and the neighbour(s). The
    // neighbour is a coarse cell, or -- when it is itself refined -- the
    // neighbour box's child cell across the fault (box<->box faulted interface).
    // New split vertices join the leaf corner pool; each connection becomes a
    // leaf face with real polygon geometry.
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& leafCorners =
        GridStateWriter::mutableGeometry(*leaf_).geomVector(std::integral_constant<int,3>());
    if (!faultedSides().empty() || !removedChildSides.empty()) {
        // The leaf corners by their coordinates. The processing of a side refines the
        // box's host cells with refineBlock, as the box's level grid does, so a corner of
        // a piece that is a corner of a box cell has that corner's coordinates, up to the
        // round-off where the box across a throw computes it.
        PointLookup cornerByCoord(
            GridStateWriter::geometry(level0_).geomVector(std::integral_constant<int,3>()));
        for (int c = 0; c < static_cast<int>(leafCorners.size()); ++c) {
            cornerByCoord.add(leafCorners[c].center(), c);
        }

        // Map a coordinate to a leaf-pool corner (dedup; append if new).
        const auto poolCorner = [&](const std::array<double,3>& nd) {
            const Dune::FieldVector<double,3> x{ nd[0], nd[1], nd[2] };
            const int existing = cornerByCoord.find(x);
            if (existing >= 0) {
                return existing;
            }
            const int idx = static_cast<int>(leafCorners.size());
            leafCorners.push_back(Dune::cpgrid::Geometry<0,3>(x));
            cornerByCoord.add(x, idx);
            return idx;
        };
        // The faceGroups map stores polygon pieces that have the same adjacent leaf cells
        // and logical side. The key is (cell, neighbour, axis, side), with neighbour == -1
        // for a boundary. The side is -1 or +1 along the given axis, relative to the first
        // cell, which may be a box child or a surviving neighbour of a removed child.
        //
        // Boundary processing can produce several pieces for the same key. These will be
        // merged into one face where possible. Disconnected boundary pieces remain
        // separate, but disconnected pieces for an internal connection are rejected.
        struct FaceGroup
        {
            int box;
            Dune::FieldVector<double,3> cellCentre;  // of the cell the pieces are a side of
            std::vector<std::vector<int>> pieces;   // leaf corner indices
        };
        std::map<std::tuple<int,int,int,int>, FaceGroup> faceGroups;
        const auto addPiece = [&](int b, const Dune::FieldVector<double,3>& cellCentre,
                                  int leafCell, int neighborLeaf, int axis, int side,
                                  const std::vector<std::array<double,3>>& nodes) {
            auto& group = faceGroups[std::make_tuple(leafCell, neighborLeaf, axis, side)];
            group.box = b;
            group.cellCentre = cellCentre;
            std::vector<int> pts;
            pts.reserve(nodes.size());
            for (const auto& nd : nodes) {
                pts.push_back(poolCorner(nd));
            }
            group.pieces.push_back(std::move(pts));
        };

        const auto centreOf = [](const CpGridData& grid, int cell) {
            const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<3,3>>& cells =
                GridStateWriter::geometry(grid).geomVector(std::integral_constant<int,0>());
            return cells[cell].center();
        };
        for (const auto& [b, axis, side] : faultedSides()) {
            const auto rangeIt = acrossRange().find(std::make_tuple(b, axis, side));
            const auto conns = boxBoundaryConnections(
                parentDims_, coord_, zcorn_, actnum_, requests_[b],
                axis, side, /*edgeConformal=*/true,
                (rangeIt != acrossRange().end()) ? &rangeIt->second : nullptr);
            const auto& rd = boxes()[b].refinedDims;
            const auto pieceArea = [](const BoundaryConnection& conn) {
                std::vector<Dune::FieldVector<double,3>> xs;
                for (const auto& nd : conn.faceNodes) {
                    xs.push_back({ nd[0], nd[1], nd[2] });
                }
                return polygonGeometry(xs).area;
            };
            std::map<std::array<int,3>, double> sideArea;
            for (const auto& conn : conns) {
                sideArea[conn.boxCell] += pieceArea(conn);
            }
            for (const auto& conn : conns) {
                // Where the sub-layers of the two sides line up across the throw, the same
                // depth comes out of two computations that can differ in the last bit, and
                // the processing cuts a piece of no area between them, which would join
                // cells that meet only along an edge. A piece below 1e-8 of its box cell's
                // side, which the cell's pieces tile, is dropped: such pieces have an area
                // of 1e-16 of the side or less, and the face geometry itself rounds at
                // about 1e-9 at field coordinates.
                if (pieceArea(conn) <= 1e-8 * sideArea[conn.boxCell]) {
                    continue;
                }
                const int refinedCart = conn.boxCell[0]
                                      + rd[0]*conn.boxCell[1]
                                      + rd[0]*rd[1]*conn.boxCell[2];
                const auto cit = refinedCartToCell()[b].find(refinedCart);
                // A child the box's MINPV removed leaves its piece to the cell across, as
                // a boundary face of that cell, where its face toward the host is gone.
                const bool childExists = cit != refinedCartToCell()[b].end();
                const int removedHost = childExists
                    ? -1
                    : hostOfRemovedChild(b, conn.boxCell);
                const int boxLeaf = childExists
                    ? leafIdxOfLevelCell()[b][cit->second]
                    : -1;
                if (removedHost < 0 && boxLeaf < 0) {
                    continue;
                }
                // coarseNeighborCart == -1: this piece of the side faces the domain,
                // where the throw leaves no cell across; it becomes a boundary face.
                int neighborLeaf = -1;   // cell opposite the box across the fault
                Dune::FieldVector<double,3> neighborCentre(0.0);
                // A neighbour that is not in the grid (inactive or removed) leaves the
                // piece facing the domain.
                const int l0 = (conn.coarseNeighborCart >= 0) ? compressed0_[conn.coarseNeighborCart] : -1;
                if (l0 >= 0) {
                    neighborLeaf = leafIdxOfCell0()[l0];
                    neighborCentre = centreOf(level0_, l0);
                    if ((neighborLeaf >= 0) && (removedHost >= 0)
                        && (boxes()[b].childrenOfParent.count(removedHost) == 0)) {
                        continue; // a host with no child left keeps its faces to coarse cells
                    }
                    if (neighborLeaf < 0) {
                        // The cell across is a host of another box: the piece joins
                        // this box cell to the child of the other box that it touches,
                        // across the fault. Each such face is emitted once, by the
                        // finer box, or by the lower box index when the in-face
                        // factors are equal. The coarser box's children on the side
                        // then have the finer box's pieces in place of their faces: a
                        // mosaic whose faces follow the throw. The processing of the
                        // side refines the cells across with the emitting box's
                        // factors, so each piece lies inside one face of a child of
                        // the coarser box.
                        const int nbBox = boxOfCell_[l0];
                        if (nbBox < 0) {
                            continue;
                        }
                        const int u = (axis == 0) ? 1 : 0;
                        const int v = (axis == 2) ? 1 : 2;
                        // checkTouchingBoxes() compares boxes by index and does not see
                        // boxes that meet only across the fault, so the two checks below
                        // are made here.
                        if (boxes()[b].graded || boxes()[nbBox].graded) {
                            // Pairing the two sides' sub-faces below compares one
                            // subdivision factor per in-face direction, which a
                            // graded box does not have.
                            OPM_THROW(std::invalid_argument,
                                "Refinement '" + requests_[b].name + "' meets '"
                                + requests_[nbBox].name + "' across a fault, and at least "
                                "one of them is graded. Box-to-box "
                                "interfaces between graded boxes are not supported.");
                        }
                        const auto& facThis = boxes()[b].factors;
                        const auto& facN    = boxes()[nbBox].factors;
                        // With factors that do not nest, both sides would skip the face
                        // below.
                        const auto nests = [u, v](const auto& fine, const auto& coarse) {
                            return fine[u] % coarse[u] == 0 && fine[v] % coarse[v] == 0;
                        };
                        if (!nests(facThis, facN) && !nests(facN, facThis)) {
                            OPM_THROW(std::invalid_argument,
                                "Refinement '" + requests_[b].name + "' meets '"
                                + requests_[nbBox].name + "' across a fault with in-face "
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
                        const std::array<int,3> pijk =
                            Opm::Lgr::getIJK(conn.coarseNeighborCart, dims0_);
                        const auto& startN = requests_[nbBox].startIJK;
                        const auto& rdN    = boxes()[nbBox].refinedDims;
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
                        // A child across that is not in the grid (its box's MINPV
                        // removed it) leaves the piece facing the domain.
                        const auto nit = refinedCartToCell()[nbBox].find(cartN);
                        if (nit != refinedCartToCell()[nbBox].end()) {
                            neighborLeaf = leafIdxOfLevelCell()[nbBox][nit->second];
                            neighborCentre = centreOf(*boxes()[nbBox].level, nit->second);
                            if (neighborLeaf < 0) {
                                continue;
                            }
                        }
                    }
                }

                if (removedHost >= 0) {
                    if (neighborLeaf >= 0) {
                        addPiece(b, neighborCentre, neighborLeaf, -1, axis, -side, conn.faceNodes);
                    }
                    continue;
                }
                addPiece(b, centreOf(*boxes()[b].level, cit->second), boxLeaf, neighborLeaf,
                         axis, side, conn.faceNodes);
            }
        }

        // A side the fault does not cross, with a child the box's MINPV removed: the
        // processing of the side, which keeps the child, gives the child's face toward the
        // cell across, which becomes that cell's boundary face. A nested box's side is
        // processed in its parent LGR's description, as its level grid was refined.
        std::map<int, RefinedBlockGrdecl> parentLgrGeometry;
        for (const auto& [b, axis, side] : removedChildSides) {
            if (nested(b) && !parentLgrGeometry.count(parentBoxOf_[b])) {
                parentLgrGeometry[parentBoxOf_[b]] = refineBlock(parentDims_, coord_, zcorn_, actnum_,
                                                                 requests_[parentBoxOf_[b]]);
            }
            // A throw of two or more whole layers leaves the side unfaulted, with the cells
            // across outside the box's own layers.
            const auto rangeIt = acrossRange().find(std::make_tuple(b, axis, side));
            const CellRange* cellsAcross = (rangeIt != acrossRange().end()) ? &rangeIt->second : nullptr;
            const auto conns = nested(b)
                ? boxBoundaryConnections(parentLgrGeometry[parentBoxOf_[b]].intDims(),
                                         parentLgrGeometry[parentBoxOf_[b]].coord.data(),
                                         parentLgrGeometry[parentBoxOf_[b]].zcorn.data(),
                                         parentLgrGeometry[parentBoxOf_[b]].actnum.data(),
                                         requests_[b], axis, side, /*edgeConformal=*/true, cellsAcross)
                : boxBoundaryConnections(parentDims_, coord_, zcorn_, actnum_, requests_[b],
                                         axis, side, /*edgeConformal=*/true, cellsAcross);
            const auto& rd = boxes()[b].refinedDims;
            for (const auto& conn : conns) {
                const int refinedCart = conn.boxCell[0] + rd[0]*(conn.boxCell[1] + rd[1]*conn.boxCell[2]);
                if (refinedCartToCell()[b].count(refinedCart) || (conn.coarseNeighborCart < 0)) {
                    continue;   // a child in the grid has its own face; none across
                }
                const int host = hostOfRemovedChild(b, conn.boxCell);
                if ((host < 0) || (boxes()[b].childrenOfParent.count(host) == 0)) {
                    continue;
                }
                const int across = parentCellOfCart(b, conn.coarseNeighborCart);
                if (across < 0) {
                    continue;
                }
                const int nb = boxOfParentGridCell(b, across);
                if (nb < 0) {
                    // Across: a cell of the grid the box refines that no box refines.
                    addPiece(b, centreOf(parentGridOf(b), across), leafOfParentGridCell(b, across), -1,
                             axis, -side, conn.faceNodes);
                    continue;
                }
                // Across: a host of another box. When that box is coarser along the side,
                // its children there have this box's faces in place of theirs
                // (buildLeafFaces()), so its child under the piece keeps it as a boundary
                // face. A box with equal or finer subdivisions keeps its own face there.
                const auto& fac = boxes()[b].factors;
                const auto& facN = boxes()[nb].factors;
                const int u = (axis == 0) ? 1 : 0;
                const int v = (axis == 2) ? 1 : 2;
                const bool coarserAcross = (fac[u] >= facN[u]) && (fac[v] >= facN[v])
                                         && ((fac[u] != facN[u]) || (fac[v] != facN[v]));
                const int child = coarserAcross ? childAcross(b, conn.boxCell, axis, side, nb, across) : -1;
                if (child >= 0) {
                    addPiece(b, centreOf(*boxes()[nb].level, child), leafIdxOfLevelCell()[nb][child], -1,
                             axis, -side, conn.faceNodes);
                }
            }
        }

        // One face per group. The normal points along +axis, never box -> neighbour:
        // CpGrid derives a face's min/max side (and hence Intersection::indexInInside)
        // from the stored normal together with each cell's orientation flag, so a
        // normal pointing down-axis makes both neighbours report the opposite face of
        // themselves. +axis is out of the box cell on its high side and into it on its
        // low side, whatever the grid's handedness. The side is kept for the
        // face_to_cell ordering below.
        const auto coordinates = [&](const std::vector<int>& pts) {
            std::vector<Dune::FieldVector<double,3>> xs;
            xs.reserve(pts.size());
            for (const int p : pts) {
                xs.push_back(leafCorners[p].center());
            }
            return xs;
        };
        for (auto& [key, group] : faceGroups) {
            const auto& [leafCell, neighborLeaf, axis, side] = key;
            const auto& cellCentre = group.cellCentre;
            for (auto& piece : group.pieces) {
                const PolyGeom pg = polygonGeometry(coordinates(piece));
                auto outward = pg.center;
                outward -= cellCentre;
                if ((pg.normal * outward) * side < 0.0) {
                    std::reverse(piece.begin(), piece.end());
                }
            }
            std::vector<std::vector<int>> faces;
            if (group.pieces.size() == 1) {
                faces = std::move(group.pieces);
            }
            else if (auto outline = mergedOutline(group.pieces)) {
                faces.push_back(std::move(*outline));
            }
            else if (neighborLeaf < 0) {
                faces = std::move(group.pieces);    // separate stretches facing the domain
            }
            else {
                OPM_THROW(std::logic_error, "The faulted side of refinement '"
                          + requests_[group.box].name + "' joins a refined cell to one "
                          "neighbour through separate pieces, which cannot form one face.");
            }
            for (auto& pts : faces) {
                const PolyGeom pg = polygonGeometry(coordinates(pts));
                syntheticFaces_.push_back(SyntheticFace{
                    leafCell, neighborLeaf, std::move(pts), pg.normal, pg.center,
                    pg.area, faceTagOf(axis), side});
            }
        }
        numLeafCorners_ = static_cast<int>(leafCorners.size());
    }

    numLeafFaces_ = numSourceFaces() + static_cast<int>(syntheticFaces_.size());
}

void LeafGridAssembler::buildLeafTopology()
{
    // ----- Leaf topology ---------------------------------------------------
    auto& leafFaceToCell = GridStateWriter::mutableFaceToCell(*leaf_);
    auto& leafFaceToPoint = GridStateWriter::mutableFaceToPoint(*leaf_);
    auto& leafCellToFace = GridStateWriter::mutableCellToFace(*leaf_);
    auto& leafFaceTag = GridStateWriter::mutableFaceTag(*leaf_);
    auto& leafFaceNormals = GridStateWriter::mutableFaceNormals(*leaf_);
    auto& leafGeometry = GridStateWriter::mutableGeometry(*leaf_);

    Dune::cpgrid::EntityVariableBase<enum face_tag>& leafTags = leafFaceTag;
    Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& leafNormals = leafFaceNormals;
    leafTags.resize(numLeafFaces());
    leafNormals.resize(numLeafFaces());

    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<2,3>>& leafFaceGeom =
        leafGeometry.geomVector(std::integral_constant<int,1>());
    leafFaceGeom.resize(numLeafFaces());

    std::vector<EntityRep<0>> rowBuffer;
    std::vector<int> pointBuffer;
    for (int face = 0; face < numSourceFaces(); ++face) {
        const SourceRef src = leafFaces()[face];
        CpGridData& srcGrid = (src.grid == 0) ? level0_ : *boxes()[src.grid - 1].level;
        auto& srcFaceToCell = (src.grid == 0) ? faceToCell0_ : GridStateWriter::faceToCell(srcGrid);
        auto& srcFaceToPoint = (src.grid == 0) ? faceToPoint0_ : GridStateWriter::faceToPoint(srcGrid);
        const Dune::cpgrid::EntityVariableBase<enum face_tag>& srcTags =
            (src.grid == 0) ? static_cast<const Dune::cpgrid::EntityVariableBase<enum face_tag>&>(faceTag0_)
                            : static_cast<const Dune::cpgrid::EntityVariableBase<enum face_tag>&>(GridStateWriter::faceTag(srcGrid));
        const Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& srcNormals =
            (src.grid == 0) ? static_cast<const Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>&>(faceNormals0_)
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
                    ? leafIdxOfCell0()[srcCell]
                    : leafIdxOfLevelCell()[src.grid - 1][srcCell];
            if (leafCell < 0) {
                continue;   // a parent with no child left: the face is the other cell's boundary
            }
            rowBuffer.push_back(EntityRep<0>(leafCell, srcRow[e].orientation()));
        }
        if (mosaicOutside()[face] >= 0) {
            // The refined cell entry already exists; add the cell across
            // with opposite orientation, keeping normal-side first.
            const bool refinedNormalOut = srcRow[0].orientation();
            const EntityRep<0> outsideEntry(mosaicOutside()[face], !refinedNormalOut);
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
        for (int e = 0; e < static_cast<int>(srcPoints.size()); ++e) {
            pointBuffer.push_back((src.grid == 0) ? srcPoints[e]
                                                  : boxes()[src.grid - 1].cornerToLeaf[srcPoints[e]]);
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
    for (int s = 0; s < static_cast<int>(syntheticFaces().size()); ++s) {
        const int face = numSourceFaces() + s;
        const SyntheticFace& sf = syntheticFaces()[s];
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
        // With no neighbour, the piece faces the domain: a boundary face.
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
        std::vector<std::vector<EntityRep<1>>> cellFaces(numLeafCells());
        for (int face = 0; face < numLeafFaces(); ++face) {
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
        for (int cell = 0; cell < numLeafCells(); ++cell) {
            leafCellToFace.appendRow(cellFaces[cell].begin(), cellFaces[cell].end());
        }
    }
}

void LeafGridAssembler::computeCellData()
{
    // ----- Leaf cell data ----------------------------------------------------
    auto& leafCellToPoint = GridStateWriter::mutableCellToPoint(*leaf_);
    leafCellToPoint.resize(numLeafCells());
    leafGlobalCell_.resize(numLeafCells());
    leafChildToParent_.assign(numLeafCells(), std::array<int,2>{-1, -1});
    leafCartesianIndexInParentCell_.assign(numLeafCells(), -1);
    leafToLevel_.resize(numLeafCells());

    for (int cell = 0; cell < numLeafCells(); ++cell) {
        const SourceRef src = leafCells()[cell];
        if (src.grid == 0) {
            // The leaf keeps level zero's corners first, at their own indices
            // (buildCornerPool), so a level-zero cell's corners need no mapping.
            leafCellToPoint[cell] = cellToPoint0_[src.index];
            leafGlobalCell_[cell] = level0_.globalCell()[src.index];
            leafToLevel_[cell] = {0, src.index};
        }
        else {
            const BoxData& box = boxes()[src.grid - 1];
            const auto& srcPoints = box.level->cellToPoint(src.index);
            for (int corner = 0; corner < 8; ++corner) {
                leafCellToPoint[cell][corner] = box.cornerToLeaf[srcPoints[corner]];
            }
            const int b = src.grid - 1;
            const int parent = GridStateWriter::childToParent(*box.level)[src.index][1];
            // A nested child's level-zero host is its parent LGR cell's.
            const int host = nested(b)
                ? GridStateWriter::childToParent(*boxes()[parentBoxOf_[b]].level)[parent][1]
                : parent;
            leafGlobalCell_[cell] = level0_.globalCell()[host];
            leafChildToParent_[cell] = {nested(b) ? parentBoxOf_[b] + 1 : 0, parent};
            leafCartesianIndexInParentCell_[cell] = GridStateWriter::cartesianIndexInParentCell(*box.level)[src.index];
            leafToLevel_[cell] = {src.grid, src.index};
        }
    }

    // Cell geometries: centers/volumes from the source grids, corner list
    // pointing into the leaf pool via the leaf's cell_to_point_.
    auto& leafGeometry = GridStateWriter::mutableGeometry(*leaf_);
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<3,3>>& leafCellGeom =
        leafGeometry.geomVector(std::integral_constant<int,0>());
    leafCellGeom.resize(numLeafCells());
    auto leafCornersPtr = leafGeometry.pointGeometryPtr();
    for (int cell = 0; cell < numLeafCells(); ++cell) {
        const SourceRef src = leafCells()[cell];
        CpGridData& srcGrid = (src.grid == 0) ? level0_ : *boxes()[src.grid - 1].level;
        const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<3,3>>& srcCellGeom =
            GridStateWriter::geometry(srcGrid).geomVector(std::integral_constant<int,0>());
        const auto& geom = srcCellGeom[src.index];
        leafCellGeom[cell] = Dune::cpgrid::Geometry<3,3>(geom.center(), geom.volume(),
                                                         leafCornersPtr,
                                                         leafCellToPoint[cell].data());
    }
}

void LeafGridAssembler::buildParallelStructures()
{
    // ----- Parallel leaf structures (distributed grids only) -------------------
    // Build the leaf's parallel cell index set, remote indices, partition types,
    // and communication interfaces so that overlap leaf cells take part in
    // collective communication (e.g. opm-models' dofTotalVolume border sync,
    // which sums interior volumes and then communicates them to the overlap).
    // Coarse leaf cells inherit their level-zero global index and attribute;
    // refined cells are rank-interior (owner only) and get fresh, globally
    // unique indices in a per-rank range above every coarse global id, so the
    // remote-index matching never pairs them with a cell on another rank.
#if HAVE_MPI
    Dune::Communication<Dune::MPIHelper::MPICommunicator> cc(comm_);
    if (cc.size() > 1) {
        using ParallelIndexSet = CpGridData::ParallelIndexSet;
        using AttributeSet = CpGridData::AttributeSet;

        std::vector<int> leafCellGlobalId(numLeafCells(), -1);
        std::vector<AttributeSet> leafCellAttr(numLeafCells(), AttributeSet::owner);

        int localMaxCoarseGlobal = 0;
        for (const auto& entry : level0_.cellIndexSet()) {
            const int l0local = entry.local();
            localMaxCoarseGlobal = std::max(localMaxCoarseGlobal, static_cast<int>(entry.global()));
            const int leafC = leafIdxOfCell0()[l0local];   // -1 if refined away
            if (leafC >= 0) {
                leafCellGlobalId[leafC] = entry.global();
                leafCellAttr[leafC] = entry.local().attribute();
            }
        }

        // Fresh global ids for refined (interior-only) cells in a per-rank
        // disjoint range above every coarse global id on any rank.
        int numRefinedLocal = 0;
        for (int cell = 0; cell < numLeafCells(); ++cell) {
            if (leafCells()[cell].grid != 0) ++numRefinedLocal;
        }
        const int globalMaxCoarse = cc.max(localMaxCoarseGlobal);
        std::vector<int> refinedCounts(cc.size(), 0);
        cc.allgather(&numRefinedLocal, 1, refinedCounts.data());
        int refinedNext = globalMaxCoarse + 1;
        for (int r = 0; r < cc.rank(); ++r) {
            refinedNext += refinedCounts[r];
        }
        for (int cell = 0; cell < numLeafCells(); ++cell) {
            if (leafCells()[cell].grid != 0) {
                leafCellGlobalId[cell] = refinedNext++;
                leafCellAttr[cell] = AttributeSet::owner;
            }
        }

        auto& cis = leaf_->cellIndexSet();
        cis.beginResize();
        for (int cell = 0; cell < numLeafCells(); ++cell) {
            cis.add(leafCellGlobalId[cell],
                    ParallelIndexSet::LocalIndex(cell, leafCellAttr[cell], true));
        }
        cis.endResize();

        leaf_->cellRemoteIndices().template rebuild<false>();
        leaf_->computeCellPartitionType();
        leaf_->computePointPartitionType();
        leaf_->computeCommunicationInterfaces(numLeafCorners());
        setDistributedLeafIds(*leaf_, level0_, leafCellGlobalId, leafFaces(), numLeafFaces(),
                              numLeafCorners(), cc);
    }
#endif
}

void LeafGridAssembler::recordLeafMetadata()
{
    // ----- Leaf metadata -------------------------------------------------------
    GridStateWriter::setLogicalCartesianSize(*leaf_, dims0_);
    GridStateWriter::setGlobalCell(*leaf_, std::move(leafGlobalCell_));
    GridStateWriter::setIndexSet(*leaf_, numLeafCells(), numLeafCorners());
    GridStateWriter::setParentRelations(*leaf_, std::move(leafChildToParent_),
                                        std::move(leafCartesianIndexInParentCell_));
    GridStateWriter::setAquiferCells(*leaf_, leafAquiferCells(level0_, leafToLevel_));
    // Level zero and each level grid learn the leaf index of their cells, -1 for a cell
    // the leaf does not have, as getLeafIdxFromLevelIdx() reads it.
    std::vector<std::vector<int>> levelToLeaf(numBoxes_ + 1);
    levelToLeaf[0].assign(numCells0_, -1);
    for (int b = 0; b < numBoxes_; ++b) {
        levelToLeaf[b + 1].assign(boxes()[b].level->size(0), -1);
    }
    for (int cell = 0; cell < numLeafCells(); ++cell) {
        const auto& [level, index] = leafToLevel_[cell];
        levelToLeaf[level][index] = cell;
    }
    GridStateWriter::setLevelToLeaf(level0_, std::move(levelToLeaf[0]));
    for (int b = 0; b < numBoxes_; ++b) {
        GridStateWriter::setLevelToLeaf(*boxes()[b].level, std::move(levelToLeaf[b + 1]));
    }
    GridStateWriter::setLeafToLevel(*leaf_, std::move(leafToLevel_));
    GridStateWriter::setCornerOrigin(*leaf_, std::move(leafCornerOrigin_));
    GridStateWriter::setRefinementMaxLevel(*leaf_, numBoxes_);
    for (int b = 0; b < numBoxes_; ++b) {
        GridStateWriter::setRefinementMaxLevel(*boxes()[b].level, numBoxes_);
    }

    // Level zero learns about its children, and each parent LGR about its nested ones.
    GridStateWriter::setParentToChildren(level0_, std::move(level0ToChildren_));
    for (auto& [pb, children] : nestedParentToChildren_) {
        GridStateWriter::setParentToChildren(*boxes()[pb].level, std::move(children));
    }
}

} // namespace Refinement
} // namespace Opm
