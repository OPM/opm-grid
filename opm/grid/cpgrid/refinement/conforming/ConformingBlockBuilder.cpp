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

#include <opm/grid/cpgrid/refinement/conforming/ConformingBlockBuilder.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/refinement/conforming/EdgeConformal.hpp>
#include <opm/grid/cpgrid/refinement/GridStateWriter.hpp>
#include <opm/grid/cpgrid/refinement/conforming/LeafGridAssembler.hpp>
#include <opm/grid/cpgrid/refinement/conforming/LevelGridAssembler.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <algorithm>
#include <array>
#include <cassert>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Opm
{
namespace Refinement
{

ConformingBlockBuilder::ConformingBlockBuilder(const std::array<int,3>& parentDims,
                                               std::vector<double> coord,
                                               std::vector<double> zcorn,
                                               std::vector<int> actnum,
                                               bool edgeConformal)
    : ConformingBlockBuilder {
          std::make_shared<const RetainedCornerPointInput>(parentDims, std::move(coord), std::move(zcorn),
                                                           std::move(actnum), edgeConformal)
      }
{
}

ConformingBlockBuilder::ConformingBlockBuilder(std::shared_ptr<const RetainedCornerPointInput> input)
    : input_(std::move(input))
{
}

namespace
{

// Terms. A box is one request: the cells of its parent grid with
// startIJK <= (i,j,k) < endIJK. Those cells are its host cells, and each host cell
// is split into child cells. The parent grid is level zero, or, for a nested box,
// the level grid of another box, its parent LGR. In a distributed grid a rank holds
// interior cells, which it owns, and overlap cells: copies of the cells next to its
// interior cells that other ranks own.

// The corner-point description of one level: index 0 is level zero, index b+1 the
// level grid of requests[b]. A nested box is refined from its parent LGR's
// description. The pointers refer to data held elsewhere.
struct LevelGeom
{
    std::array<int,3> dims{};
    const double* coord{nullptr};
    const double* zcorn{nullptr};
    const int* actnum{nullptr};
};

// Whether this rank refines a box (Owned), or holds none of its cells and builds an
// empty level grid for it (Absent), so that every rank has the same levels.
enum class BoxPresence { Owned, Absent };

// The presence of a top-level box on this rank of a distributed grid. Owned: the
// rank holds cells of the box, all of them interior. Absent: the rank holds no cell
// of the box. Otherwise the rank holds a cell of the box as an overlap cell, so the
// box and the cells next to it are not all owned by one rank: this throws.
BoxPresence classifyBoxPresence(const Dune::cpgrid::CpGridData& level0,
                                const Opm::Refinement::BlockRefinement& req)
{
    int interior = 0;
    int nonInterior = 0;
    for (int c = 0; c < level0.size(0); ++c) {
        std::array<int,3> ijk;
        level0.getIJK(c, ijk);
        if (contains(req, ijk)) {
            const Dune::cpgrid::Entity<0> e(level0, c, true);
            if (e.partitionType() == Dune::InteriorEntity) {
                ++interior;
            } else {
                ++nonInterior;
            }
        }
    }
    // Presence is decided by the partition type of the box cells this rank holds. It
    // is not decided by comparing their number with the number of cells of the box:
    // that number is not known here. The box's Cartesian size counts inactive cells,
    // and the retained ACTNUM does not mark the cells MINPV removes.
    if (interior == 0 && nonInterior == 0) {
        return BoxPresence::Absent;
    }
    if (interior > 0 && nonInterior == 0) {
        return BoxPresence::Owned;
    }
    OPM_THROW(std::logic_error, "Refinement box '" + req.name + "' is split across MPI ranks or "
                           "touches the overlap region. Only rank-interior LGRs are supported; "
                           "use CpGrid::setPartitionCellGroups to keep the box on one rank.");
}

// Number of cells of the box's level grid: the children of the active cells of the
// parent grid that are inside the refinement box, less the children that the box's
// own MINPV removes. Inactive cells, and cells the parent grid's MINPV removes, are
// not in the parent grid.
int childCellsLeft(const Dune::cpgrid::CpGridData& parent,
                   const Opm::Refinement::BlockRefinement& req)
{
    std::array<int,3> size{};
    for (int d = 0; d < 3; ++d) {
        size[d] = req.endIJK[d] - req.startIJK[d];
    }
    const auto hostIndex = [&size](int i, int j, int k) {
        return i + static_cast<std::size_t>(size[0])*(j + static_cast<std::size_t>(size[1])*k);
    };
    std::vector<bool> present(static_cast<std::size_t>(size[0])*size[1]*size[2], false);
    for (int c = 0; c < parent.size(0); ++c) {
        std::array<int,3> ijk;
        parent.getIJK(c, ijk);
        if (contains(req, ijk)) {
            present[hostIndex(ijk[0] - req.startIJK[0], ijk[1] - req.startIJK[1],
                              ijk[2] - req.startIJK[2])] = true;
        }
    }
    const std::array<AxisSubdivision,3> subs = { axisSubdivision(req, 0),
                                                axisSubdivision(req, 1),
                                                axisSubdivision(req, 2) };
    const auto dims = refinedDims(req);
    // A mask of the wrong size is refused when the box is refined.
    const bool masked = req.minpvRemoved.size()
                     == static_cast<std::size_t>(dims[0])*dims[1]*dims[2];
    int left = 0;
    std::size_t child = 0;
    for (int k = 0; k < dims[2]; ++k) {
        for (int j = 0; j < dims[1]; ++j) {
            for (int i = 0; i < dims[0]; ++i, ++child) {
                const bool removed = masked && (req.minpvRemoved[child] != 0);
                if (!removed && present[hostIndex(subs[0].parentOffset[i], subs[1].parentOffset[j],
                                                  subs[2].parentOffset[k])]) {
                    ++left;
                }
            }
        }
    }
    return left;
}

// Refuses two boxes of one parent grid whose subdivisions do not fit where the boxes
// touch. In each direction, the index ranges of two boxes either have a gap between
// them, overlap, or abut. Boxes with a gap in some direction share nothing. Otherwise
// they share a face, an edge or a corner of the parent grid, and the factors (the
// numbers of children per host cell) must fit in the directions where the ranges
// overlap:
// - on a shared face, the face's two directions: in each, one box's factor is a
//   multiple of the other's, and the same box has the larger factor, or an equal
//   one, in both;
// - along a shared edge, its direction: the factors are equal;
// - a shared corner is a corner of the parent grid, with no rule.
// The factors in a direction where the boxes abut do not matter: the boxes meet there
// on a boundary between host cells. A graded box may share only a corner.
void checkTouchingBoxes(const std::vector<Opm::Refinement::BlockRefinement>& requests)
{
    for (std::size_t i = 0; i < requests.size(); ++i) {
        for (std::size_t j = i + 1; j < requests.size(); ++j) {
            const auto& a = requests[i];
            const auto& b = requests[j];
            // Boxes of different parent grids have their (i,j,k) in different index
            // spaces. They do not touch: a nested box stays at least one cell inside its
            // parent LGR.
            if (a.parentGridName != b.parentGridName) {
                continue;
            }
            // A gap in some direction: the boxes share nothing.
            bool interacts = true;
            for (int c = 0; c < 3; ++c) {
                const bool gap = (a.endIJK[c] < b.startIJK[c]) || (b.endIJK[c] < a.startIJK[c]);
                interacts = interacts && !gap;
            }
            if (!interacts) {
                continue;
            }
            // The number of directions in which the index ranges overlap; in the others
            // the boxes abut. It is at most 2: boxes whose ranges overlap in all three
            // directions share cells, and validateBlockRefinements() refuses those before
            // this. Two overlaps: the boxes share a face. One overlap: an edge. None: a
            // corner.
            int overlapCount = 0;
            for (int c = 0; c < 3; ++c) {
                const bool overlap = (a.startIJK[c] < b.endIJK[c]) && (b.startIJK[c] < a.endIJK[c]);
                if (overlap) {
                    ++overlapCount;
                }
            }
            // A graded box is subdivided column by column, not by one factor per
            // direction, so on a face or an edge it shares, the two boxes need not have
            // the same points: the leaf would have a gap there, or points on one side
            // only. A shared corner is a corner of the parent grid.
            const auto graded = [](const Opm::Refinement::BlockRefinement& request) {
                return std::any_of(request.subdivision.begin(), request.subdivision.end(),
                                   [](const auto& sub) { return !sub.empty(); });
            };
            if ((overlapCount > 0) && (graded(a) || graded(b))) {
                OPM_THROW(std::invalid_argument,
                    "Refinement boxes '" + a.name + "' and '" + b.name + "' touch and at least "
                    "one is graded. Touching graded boxes are not supported; keep a graded box "
                    "apart from other boxes (meeting at a corner is fine), or refine both "
                    "uniformly.");
            }
            // Compare the factors in the directions where the ranges overlap.
            bool anyDiffer = false;
            bool anyIncompatible = false;
            int incompatDir = -1;
            // Whether a, or b, has a factor at least as large as the other's in every
            // such direction: that box is then the finer side of a shared face.
            bool aFinerOrEqualAll = true;
            bool bFinerOrEqualAll = true;
            for (int c = 0; c < 3; ++c) {
                const bool overlap = (a.startIJK[c] < b.endIJK[c]) && (b.startIJK[c] < a.endIJK[c]);
                if (!overlap || a.cellsPerDim[c] == b.cellsPerDim[c]) {
                    continue;
                }
                anyDiffer = true;
                const int na = a.cellsPerDim[c];
                const int nb = b.cellsPerDim[c];
                const int hi = (na > nb) ? na : nb;
                const int lo = (na > nb) ? nb : na;
                // Graded boxes that touch are refused above, and uniform factors are at least 1.
                assert(lo > 0);
                if (hi % lo != 0) {
                    anyIncompatible = true;
                    incompatDir = c;
                }
                if (na < nb) { aFinerOrEqualAll = false; }
                if (nb < na) { bFinerOrEqualAll = false; }
            }
            if (!anyDiffer) {
                continue;                       // equal in every overlap direction
            }
            if (anyIncompatible) {
                const int na = a.cellsPerDim[incompatDir];
                const int nb = b.cellsPerDim[incompatDir];
                OPM_THROW(std::logic_error,
                    "Refinement boxes '" + a.name + "' and '" + b.name +
                    "' meet with incompatible subdivisions in direction " +
                    std::to_string(incompatDir) + " (" + std::to_string(na) + " vs " +
                    std::to_string(nb) + "); neither is a multiple of the other, so the shared "
                    "interface cannot be made conformal. Use equal subdivisions in the shared "
                    "direction(s).");
            }
            // The factors differ, and in each direction one is a multiple of the other.
            // This is supported on a shared face where one box is the finer side in both
            // directions: a cell of the coarser box on that face gets, in place of its
            // face there, the sub-faces of the finer box that cover it, so it has more
            // than six faces (LeafGridAssembler builds them). Different factors along a
            // shared edge, and a face where each box is finer in a different direction,
            // are not supported.
            const bool cleanlyNested = aFinerOrEqualAll || bFinerOrEqualAll;
            if (overlapCount == 2 && cleanlyNested) {
                continue;                       // a shared face, one box finer in both
            }
            OPM_THROW(std::logic_error,
                "Refinement boxes '" + a.name + "' and '" + b.name +
                "' meet with different but compatible subdivisions, but in an unsupported "
                "configuration (" + (overlapCount == 2
                    ? std::string("mixed nesting -- each box is finer in a different in-face "
                                  "direction")
                    : std::string("edge/corner contact, which would need a 1-D mosaic")) +
                "). Only a face-sharing pair with one box uniformly finer is supported; that "
                "is not implemented for this case yet. Use equal subdivisions in the shared "
                "direction(s).");
        }
    }
}

} // anonymous namespace

void ConformingBlockBuilder::build(Dune::CpGrid& grid,
                                   const std::vector<BlockRefinement>& requests)
{
    if (grid.maxLevel() != 0) {
        OPM_THROW(std::logic_error, "ConformingBlockBuilder requires an unrefined starting grid.");
    }
    checkTouchingBoxes(requests);

    auto& storage = grid.currentData();
    // Distributed: loadBalance() has scattered the grid over the ranks, and each box
    // is refined by the one rank that holds it.
    const bool distributed = grid.comm().size() > 1 && grid.isDistributed();
    // Before load balancing, on more than one rank: the whole grid is on rank 0 and
    // the other ranks hold no cells. Rank 0 refines every box, the other ranks build
    // empty level grids, and the leaf is assembled with the self-communicator. The
    // refined grid is distributed afterwards.
    const bool refineBefore = grid.comm().size() > 1 && !grid.isDistributed();

    // The corner-point description of each level, filled in as the level grids are
    // built.
    const auto& dims = input_->dims;
    const int* actnum0 = input_->actnum.empty() ? nullptr : input_->actnum.data();
    std::vector<LevelGeom> levelGeom(requests.size() + 1);
    levelGeom[0] = LevelGeom{ dims, input_->coord.data(), input_->zcorn.data(), actnum0 };
    // Each refined level's own description, for the boxes nested in it.
    std::vector<RefinedBlockGrdecl> refinedGeom(requests.size());
    std::map<std::string,int> nameToLevel{ {"GLOBAL", 0} };
    // The presence of each box on this rank, by name: a nested box takes its parent
    // LGR's.
    std::map<std::string,BoxPresence> presenceOf;

    // Builds the level grids on this rank. In a distributed run this can throw on
    // some ranks only, for example where classifyBoxPresence refuses a box. The other
    // ranks would go on into the collective operations of the leaf assembly and wait
    // there forever, so the call below is guarded.
    const auto buildLevelGrids = [&]() {
        for (std::size_t b = 0; b < requests.size(); ++b) {
            const auto& req = requests[b];
            const auto itParent = nameToLevel.find(req.parentGridName);
            if (itParent == nameToLevel.end()) {
                OPM_THROW(std::logic_error, "Refinement box '" + req.name + "' names parent grid '"
                    + req.parentGridName + "', which has not been built yet. Parent LGRs "
                    "must be ordered before their children.");
            }
            const int parentLevel = itParent->second;
            const LevelGeom& pg = levelGeom[parentLevel];

            // Whether this rank refines the box. Serial: always. Distributed: the rank
            // that holds the box. Before load balancing: rank 0. A nested box's (i,j,k)
            // are in its parent LGR's index space, not level zero's. It lies inside its
            // parent LGR, so it takes the parent's presence, which is known: a parent
            // LGR comes before the boxes nested in it.
            const bool nested = (req.parentGridName != "GLOBAL");
            BoxPresence presence = BoxPresence::Owned;
            if (distributed) {
                presence = nested ? presenceOf.at(req.parentGridName)
                                  : classifyBoxPresence(*storage[0], req);
            } else if (refineBefore) {
                presence = nested ? presenceOf.at(req.parentGridName)
                                  : ((grid.getLevelZeroComm().rank() == 0) ? BoxPresence::Owned
                                                                           : BoxPresence::Absent);
            }
            // A box with no host cell left (all inactive, or removed by MINPV), or
            // whose own MINPV removes every child, has nothing to refine: it
            // becomes an empty LGR, as on a rank that does not hold it. The
            // corner-point processing below needs at least one cell.
            if ((presence == BoxPresence::Owned)
                && (childCellsLeft(*storage[parentLevel], req) == 0)) {
                presence = BoxPresence::Absent;
            }
            presenceOf[req.name] = presence;
            RefinedBlockGrdecl& childRefined = refinedGeom[b];
            auto level = (presence == BoxPresence::Owned)
                ? assembleBlockLevelGrid(*storage[parentLevel], pg.dims,
                                         pg.coord, pg.zcorn, pg.actnum,
                                         req, static_cast<int>(b) + 1, parentLevel,
                                         storage, &childRefined)
                : assembleEmptyLevelGrid(req, static_cast<int>(b) + 1, storage);
            storage.push_back(std::move(level));
            levelGeom[b + 1] = LevelGeom{ childRefined.intDims(),
                                          childRefined.coord.data(),
                                          childRefined.zcorn.data(),
                                          childRefined.actnum.empty() ? nullptr
                                                                      : childRefined.actnum.data() };
            nameToLevel[req.name] = static_cast<int>(b) + 1;
        }
    };

    if (distributed) {
        // An error on this rank is caught. One reduction tells every rank the lowest
        // rank that failed, and then every rank throws, so no rank waits in the
        // collective operations of the leaf assembly. The lowest rank that failed
        // throws with its own message and logs it; the other ranks throw without
        // logging.
        // TODO: an error thrown inside this phase is still logged on every rank
        // that hits it. Per-rank logging, as with opm-simulators' DeferredLogger,
        // would collect those messages.
        int localBuildError = 0;
        std::string localBuildMsg;
        try {
            buildLevelGrids();
        } catch (const std::exception& e) {
            localBuildError = 1;
            localBuildMsg = e.what();
        }
        const auto& cc = grid.getLevelZeroComm();
        // The lowest rank that failed, or the number of ranks if none did.
        const int failedRank = cc.min(localBuildError != 0 ? cc.rank() : cc.size());
        if (failedRank < cc.size()) {
            if (cc.rank() == failedRank) {
                OPM_THROW(std::runtime_error, "Parallel LGR build failed on rank "
                          + std::to_string(cc.rank()) + ": " + localBuildMsg);
            }
            OPM_THROW_NOLOG(std::runtime_error, "Parallel LGR build aborted on rank "
                            + std::to_string(cc.rank()) + " because rank "
                            + std::to_string(failedRank) + " failed (see that rank's message).");
        }
    } else {
        buildLevelGrids();
    }

    // A distributed refinement assembles its leaf with level zero's
    // communicator. A serial one, also the one before load balancing, uses a
    // self-communicator, so the assembly does no collective operations.
    std::shared_ptr<Dune::cpgrid::CpGridData> leaf =
        LeafGridAssembler(storage, requests, dims,
                          input_->coord.data(), input_->zcorn.data(), actnum0,
                          distributed
                              ? static_cast<Dune::MPIHelper::MPICommunicator>(grid.getLevelZeroComm())
                              : Dune::MPIHelper::getLocalCommunicator()).build();
    if (input_->edgeConformal) {
        edgeConformalizeLeaf(*leaf);
    }
    storage.push_back(std::move(leaf));

    if (refineBefore) {
        // The new level grids and the leaf have the self-communicator. They get level
        // zero's, so that the refined grid, which is load balanced next, reports every
        // rank: CpGrid::comm() returns the leaf's communicator.
        for (std::size_t i = 1; i < storage.size(); ++i) {
            GridStateWriter::setCommunicator(*storage[i], grid.getLevelZeroComm());
        }
    }
}

} // namespace Refinement
} // namespace Opm
