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

#include <opm/grid/cpgrid/refinement/conforming/ConformingBlockBuilder.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/refinement/conforming/EdgeConformal.hpp>
#include <opm/grid/cpgrid/refinement/GridStateWriter.hpp>
#include <opm/grid/cpgrid/refinement/conforming/LeafGridAssembler.hpp>
#include <opm/grid/cpgrid/refinement/conforming/LevelGridAssembler.hpp>

#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <array>
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
    : dims_(parentDims)
    , coord_(std::move(coord))
    , zcorn_(std::move(zcorn))
    , actnum_(std::move(actnum))
    , edgeConformal_(edgeConformal)
{
}

namespace
{

// Per-box presence on this rank (distributed runs). The rank-interior model
// (PLAN Track 1 step 6): a box is either fully owned here (all its cells
// local and interior) or absent (no cell local, not even as overlap).
// Anything in between means the box crosses a rank boundary / touches the
// overlap, which is not supported.
enum class BoxPresence { Owned, Absent };

BoxPresence classifyBox(const Dune::cpgrid::CpGridData& level0,
                        const std::array<int,3>& dims,
                        const Opm::Refinement::BlockRefinement& req)
{
    int interior = 0;
    int nonInterior = 0;
    const auto& globalCell = level0.globalCell();
    for (int c = 0; c < level0.size(0); ++c) {
        const int cart = globalCell[c];
        const int i = cart % dims[0];
        const int j = (cart / dims[0]) % dims[1];
        const int k = cart / (dims[0]*dims[1]);
        if (i >= req.startIJK[0] && i < req.endIJK[0]
            && j >= req.startIJK[1] && j < req.endIJK[1]
            && k >= req.startIJK[2] && k < req.endIJK[2]) {
            const Dune::cpgrid::Entity<0> e(level0, c, true);
            if (e.partitionType() == Dune::InteriorEntity) {
                ++interior;
            } else {
                ++nonInterior;
            }
        }
    }
    // Presence is judged by what this rank holds, not by counting the box: a
    // rank that owns every box cell it can see, and sees no box cell it does not
    // own, holds the whole box -- a split box puts the far side in this rank's
    // overlap, which is exactly the non-interior case. Counting instead needs a
    // total that is right for the distributed grid, and the two candidates are
    // both wrong: the Cartesian extent ignores inactive cells, and the retained
    // ACTNUM predates MINPV. On Norne those two errors are 623 and 97 cells.
    if (interior == 0 && nonInterior == 0) {
        return BoxPresence::Absent;
    }
    if (interior > 0 && nonInterior == 0) {
        return BoxPresence::Owned;
    }
    throw std::logic_error("Refinement box '" + req.name + "' is split across MPI ranks or "
                           "touches the overlap region. Only rank-interior LGRs are supported; "
                           "use CpGrid::setPartitionCellGroups to keep the box on one rank.");
}

// refine-before-redistribute classification: before load balancing the global
// grid lives on rank 0 (all cells) and is empty on the other ranks. Classify a
// box purely by cell *presence* (ignoring partition type, which is not yet set):
// the rank that holds all of the box's cells refines it; ranks with none carry
// an empty placeholder level grid. A partial presence (some but not all cells)
// is not expected here and is rejected.
BoxPresence boxCellPresence(const Dune::cpgrid::CpGridData& level0,
                            const std::array<int,3>& dims,
                            const Opm::Refinement::BlockRefinement& req)
{
    const int boxCells = (req.endIJK[0] - req.startIJK[0])
                       * (req.endIJK[1] - req.startIJK[1])
                       * (req.endIJK[2] - req.startIJK[2]);
    int present = 0;
    const auto& globalCell = level0.globalCell();
    for (int c = 0; c < level0.size(0); ++c) {
        const int cart = globalCell[c];
        const int i = cart % dims[0];
        const int j = (cart / dims[0]) % dims[1];
        const int k = cart / (dims[0]*dims[1]);
        if (i >= req.startIJK[0] && i < req.endIJK[0]
            && j >= req.startIJK[1] && j < req.endIJK[1]
            && k >= req.startIJK[2] && k < req.endIJK[2]) {
            ++present;
        }
    }
    if (present == 0) {
        return BoxPresence::Absent;
    }
    if (present == boxCells) {
        return BoxPresence::Owned;
    }
    throw std::logic_error("Refinement box '" + req.name + "' is only partially present on "
                           "this rank before load balancing; refine-before-redistribute "
                           "expects the full global grid on a single rank.");
}

} // anonymous namespace

void ConformingBlockBuilder::build(Dune::CpGrid& grid,
                                   const std::vector<BlockRefinement>& requests)
{
    if (grid.maxLevel() != 0) {
        throw std::logic_error("ConformingBlockBuilder requires an unrefined starting grid.");
    }
    // Conformity check for touching boxes. With disjoint boxes, a pair
    // either is fully separated, shares an edge/corner, or shares a 2D
    // face. The shared refined entities (corners along a shared edge,
    // corners+faces on a shared face) coincide bitwise only if the
    // subdivision factors match in the *shared* directions: the overlap
    // dimension(s). The touch dimension(s) are parent-cell boundaries where
    // corners are identified with level zero regardless. Subdivisions in a
    // dimension where the boxes neither touch nor overlap are unconstrained
    // (e.g. the out-of-face direction of two face-sharing boxes).
    for (std::size_t i = 0; i < requests.size(); ++i) {
        for (std::size_t j = i + 1; j < requests.size(); ++j) {
            const auto& a = requests[i];
            const auto& b = requests[j];
            // Only boxes refining the same parent grid live in a common
            // Cartesian index space; their IJK extents are comparable. Boxes
            // with different parents (e.g. a nested box vs a top-level box)
            // are indexed in different spaces and cannot touch by construction.
            if (a.parentGridName != b.parentGridName) {
                continue;
            }
            // The boxes share an entity (corner/edge/face) only if they
            // touch-or-overlap in every dimension; a gap in any dimension
            // means they are fully separated and impose no constraint.
            bool interacts = true;
            for (int c = 0; c < 3; ++c) {
                const bool gap = (a.endIJK[c] < b.startIJK[c]) || (b.endIJK[c] < a.startIJK[c]);
                interacts = interacts && !gap;
            }
            if (!interacts) {
                continue;
            }
            // Classify how the two boxes meet. They interact (no gap in any
            // direction), so each direction is either an *overlap* (their
            // projections overlap) or a *touch* (they abut). A pair that
            // overlaps in two directions and touches in the third shares a 2-D
            // face; the two overlap directions are the *in-face* directions and
            // the touch direction is perpendicular. Conformity constrains only
            // the in-face directions.
            int overlapCount = 0;
            for (int c = 0; c < 3; ++c) {
                const bool overlap = (a.startIJK[c] < b.endIJK[c]) && (b.startIJK[c] < a.endIJK[c]);
                if (overlap) {
                    ++overlapCount;
                }
            }
            // Scan the overlap directions where the subdivisions differ.
            bool anyDiffer = false;
            bool anyIncompatible = false;
            int incompatDir = -1;
            // Cleanly nested: the *same* box is finer-or-equal in every differing
            // overlap direction (so one box is uniformly the finer side and the
            // mosaic has a well-defined coarse side).
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
                if (lo <= 0 || hi % lo != 0) {
                    anyIncompatible = true;
                    incompatDir = c;
                }
                if (na < nb) { aFinerOrEqualAll = false; }
                if (nb < na) { bFinerOrEqualAll = false; }
            }
            if (!anyDiffer) {
                continue;                       // equal in every shared direction
            }
            if (anyIncompatible) {
                const int na = a.cellsPerDim[incompatDir];
                const int nb = b.cellsPerDim[incompatDir];
                throw std::logic_error(
                    "Refinement boxes '" + a.name + "' and '" + b.name +
                    "' meet with incompatible subdivisions in direction " +
                    std::to_string(incompatDir) + " (" + std::to_string(na) + " vs " +
                    std::to_string(nb) + "); neither is a multiple of the other, so the shared "
                    "interface cannot be made conformal. Use equal subdivisions in the shared "
                    "direction(s).");
            }
            // Compatible (every differing overlap factor is a multiple of the
            // other). Supported only as a 2-D sub-face mosaic on a shared face
            // (overlapCount == 2) with one box uniformly the finer side
            // (cleanly nested): the coarser side's interface cells become
            // >6-face hexes tiled by the finer side's sub-faces (LGR_GAPS A2,
            // assembled in assembleLeafGrid). Edge/corner contact (a 1-D mosaic)
            // and a mixed nesting (each box finer in a different direction) are
            // not implemented.
            const bool cleanlyNested = aFinerOrEqualAll || bFinerOrEqualAll;
            if (overlapCount == 2 && cleanlyNested) {
                continue;                       // supported compatible face mosaic
            }
            throw std::logic_error(
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

    auto& storage = grid.currentData();
    // The refinement is "distributed" (rank-interior: only the owning rank
    // refines each box) only once the grid has actually been scattered. Before
    // load balancing - the refine-before-redistribute path - the grid is still
    // the replicated global grid; refine it serially on every rank instead, so
    // the builder refines every cell (no rank-interior classifyBox) and the
    // leaf assembly does no collectives, producing the full refined grid a
    // serial run yields, ready to be distributed afterwards.
    const bool distributed = grid.comm().size() > 1 && grid.isDistributed();
    // Capture the grid's actual (world) communicator before any new level or
    // leaf grid is pushed - CpGrid::comm() reads the back() entry of the
    // storage, which is the coarse level-zero grid at this point and therefore
    // still carries the world communicator on every rank.
    const auto worldComm = grid.comm();
    // Use a self-communicator for serial refinement (refine-before-redistribute
    // and the self-comm output/reference grid) so no operation is collective
    // over the world communicator; the distributed simulation grid keeps the
    // world communicator.
    Dune::MPIHelper::MPICommunicator comm = Dune::MPIHelper::getLocalCommunicator();
    if (distributed) {
        comm = grid.comm();
    }
    // refine-before-redistribute: comm().size() > 1 but the grid is not yet
    // scattered. The global grid is on rank 0 and empty on the others, so we
    // classify boxes by cell presence (rank 0 refines, the rest get empty
    // placeholders) and refine serially with the self-communicator above.
    const bool refineBefore = grid.comm().size() > 1 && !grid.isDistributed();

    // Nested refinement (parent != GLOBAL): build the level grids over their
    // parent LGR and stitch them with assembleNestedLeafGrid. A fully contained
    // nested box rides on its parent rank in both parallel models (the rank-
    // interior "refine-after" path and refine-before-redistribute), so detect
    // nesting and route to the nested assembler; the containment constraints are
    // enforced inside assembleNestedLeafGrid.
    bool anyNested = false;
    for (const auto& req : requests) {
        if (req.parentGridName != "GLOBAL") {
            anyNested = true;
        }
    }

    // Parent corner-point description per built level (index 0 = GLOBAL base,
    // index b+1 = the level grid of request b). A nested box is refined from
    // its parent level's resampled geometry, not the global grid. LevelGeom is
    // declared in LeafGridAssembler.hpp so the nested leaf assembler can reuse
    // the per-level geometry.
    std::vector<LevelGeom> levelGeom(requests.size() + 1);
    levelGeom[0] = LevelGeom{ dims_, coord_, zcorn_, actnum_ };
    std::map<std::string,int> nameToLevel{ {"GLOBAL", 0} };
    // Per-box presence (refine-before-redistribute), keyed by box name, so a
    // nested box can inherit its parent box's presence (see below).
    std::map<std::string,BoxPresence> presenceOf;

    // Per-rank construction of the level grids. In a distributed run this can
    // throw on one rank (e.g. classifyBox rejecting a box that ends up touching
    // the overlap region) while other ranks succeed and march on into the leaf
    // assembly's collectives (cc.max / cc.allgather). Those would then block
    // forever waiting for the rank that already threw -- the classic asymmetric
    // exception deadlock. Run the loop under a collective-safe guard below.
    const auto buildLevelGrids = [&]() {
        for (std::size_t b = 0; b < requests.size(); ++b) {
            const auto& req = requests[b];
            const auto itParent = nameToLevel.find(req.parentGridName);
            if (itParent == nameToLevel.end()) {
                throw std::logic_error("Refinement box '" + req.name + "' names parent grid '"
                    + req.parentGridName + "', which has not been built yet. Parent LGRs "
                    "must be ordered before their children.");
            }
            const int parentLevel = itParent->second;
            const LevelGeom& pg = levelGeom[parentLevel];

            // In a distributed run, only the rank that owns the (rank-interior)
            // box refines it; other ranks carry an empty placeholder level grid.
            // Nested boxes in a distributed run are still serial-only (guarded
            // above), so a distributed run is always top-level here and classifies
            // against level zero. The refine-before-redistribute path puts the full
            // global grid on a single rank and classifies a top-level box by its
            // cell presence there. A NESTED box addresses cells in its parent LGR's
            // own (local) Cartesian space, so it cannot be classified against the
            // global grid; instead it inherits its parent box's presence - a child
            // fully contained in its parent lives on whatever rank refines the
            // parent (and the parent-before-child ordering guarantees the parent's
            // presence is already known). Plain serial stays Owned.
            const bool nested = (req.parentGridName != "GLOBAL");
            BoxPresence presence = BoxPresence::Owned;
            if (distributed) {
                // Rank-interior model: a nested box's I/J/K are parent-LGR-local, so
                // it cannot be classified against level zero. A contained child lives
                // on whichever rank owns its parent box, so inherit the parent's
                // classification (parent processed before child).
                presence = nested ? presenceOf.at(req.parentGridName)
                                  : classifyBox(*storage[0], dims_, req);
            } else if (refineBefore) {
                presence = nested ? presenceOf.at(req.parentGridName)
                                  : boxCellPresence(*storage[0], dims_, req);
            }
            presenceOf[req.name] = presence;
            RefinedBlockGrdecl childRefined;
            auto level = (presence == BoxPresence::Owned)
                ? assembleBlockLevelGrid(*storage[parentLevel], pg.dims,
                                         pg.coord.data(), pg.zcorn.data(),
                                         pg.actnum.empty() ? nullptr : pg.actnum.data(),
                                         req, static_cast<int>(b) + 1, parentLevel,
                                         storage, comm, &childRefined)
                : assembleEmptyLevelGrid(req, static_cast<int>(b) + 1, storage, comm);
            storage.push_back(std::move(level));
            levelGeom[b + 1] = LevelGeom{ childRefined.dims,
                                          std::move(childRefined.coord),
                                          std::move(childRefined.zcorn),
                                          std::move(childRefined.actnum) };
            nameToLevel[req.name] = static_cast<int>(b) + 1;
        }
    };

    if (distributed) {
        // Collective-safe: catch a per-rank failure, agree across ranks with a
        // single reduction, and throw symmetrically (with the offending rank's
        // message) so the run fails cleanly instead of deadlocking in the leaf
        // assembly's collectives below.
        int localBuildError = 0;
        std::string localBuildMsg;
        try {
            buildLevelGrids();
        } catch (const std::exception& e) {
            localBuildError = 1;
            localBuildMsg = e.what();
        }
        Dune::Communication<Dune::MPIHelper::MPICommunicator> cc(comm);
        if (cc.max(localBuildError) != 0) {
            throw std::runtime_error(localBuildError != 0
                ? ("Parallel LGR build failed on rank " + std::to_string(cc.rank())
                   + ": " + localBuildMsg)
                : ("Parallel LGR build aborted on rank " + std::to_string(cc.rank())
                   + " because another rank failed (see that rank's message)."));
        }
    } else {
        buildLevelGrids();
    }

    // The nested hierarchy's level grids are built. Nested leaf stitching is
    // isolated in assembleNestedLeafGrid so the single-level assembler below
    // (and therefore every GLOBAL-parent deck) stays byte-identical.
    std::shared_ptr<Dune::cpgrid::CpGridData> leaf = anyNested
        ? assembleNestedLeafGrid(storage, requests, levelGeom, comm)
        : assembleLeafGrid(storage, requests, dims_,
                           coord_.data(), zcorn_.data(),
                           actnum_.empty() ? nullptr : actnum_.data(),
                           comm);
    if (edgeConformal_) {
        edgeConformalizeLeaf(*leaf);
    }
    storage.push_back(std::move(leaf));

    if (refineBefore) {
        // The new level grids (indices 1..N) and the leaf (back) were assembled
        // with a self-communicator so the serial refinement does no collectives.
        // Restore the world communicator on them: CpGrid::comm() reads the leaf,
        // so without this the refined grid would report a size-1 communicator
        // and the subsequent load balancing would treat it as a serial grid and
        // skip distributing the refined leaf. Level zero (index 0) is the
        // pre-existing coarse grid and already carries the world communicator.
        for (std::size_t i = 1; i < storage.size(); ++i) {
            GridStateWriter::setCommunicator(*storage[i], worldComm);
        }
    }
}

} // namespace Refinement
} // namespace Opm
