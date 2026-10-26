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
#ifndef OPM_GRID_REFINEMENT_LEAFGRIDASSEMBLER_HEADER_INCLUDED
#define OPM_GRID_REFINEMENT_LEAFGRIDASSEMBLER_HEADER_INCLUDED

#include <opm/grid/cpgrid/refinement/RefinementRequest.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <memory>
#include <vector>

namespace Dune
{
namespace cpgrid
{
class CpGridData;
}
}

namespace Opm
{
namespace Refinement
{

/// Stages 5+6b of the build pipeline (docs/DESIGN-builder.md D4): merge
/// level zero and the refined level grids into the leaf view.
///
/// Leaf cells are the unrefined level-zero cells plus all refined cells
/// (children replace their parent at its position in the cell ordering).
/// Parent faces on block boundaries are replaced by the refined mosaic,
/// pairing refined cells with their level-zero neighbors. Corners follow
/// the shared-pool design (D1): the leaf corner vector starts with all
/// level-zero corners (indices preserved), refined corners that coincide
/// with parent-lattice corners are identified with the level-zero corner
/// through the parent cells (so faults inside blocks identify per cell,
/// not per lattice key), and the rest are appended.
///
/// Ids work through the existing delegation machinery: leaf_to_level_cells_
/// and corner_history_ are populated so IdSet routes every leaf entity to
/// the id of its birth-level entity.
///
/// Block-boundary parent faces crossed by a fault are rebuilt from the
/// corner-point processor (split faces against the coarse neighbour(s)); the
/// parent description is needed for that.
///
/// A nested box (one refining another box's level grid) replaces the cells
/// it refines by its own: its corners on its parents' corners are theirs, and
/// its sides face the parent LGR's cells. It must lie inside its parent LGR,
/// away from the parent's boundary, one level deep, with no fault or pinch-out
/// on its sides.
///
/// @param storage   The hierarchy [level0, level1..levelB]; level grids as
///                  produced by assembleBlockLevelGrid for requests[b].
/// @param requests  The block requests, one per level grid; a parent LGR comes
///                  before the boxes nested in it.
/// @param parentDims, coord, zcorn, actnum  The parent corner-point
///                  description (actnum may be null).
/// @param comm      Communicator for the leaf grid object.
/// @return the assembled leaf; the caller appends it to storage.
std::shared_ptr<Dune::cpgrid::CpGridData>
assembleLeafGrid(std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& storage,
                 const std::vector<BlockRefinement>& requests,
                 const std::array<int,3>& parentDims,
                 const double* coord,
                 const double* zcorn,
                 const int* actnum,
                 Dune::MPIHelper::MPICommunicator comm);

} // namespace Refinement
} // namespace Opm

#endif // OPM_GRID_REFINEMENT_LEAFGRIDASSEMBLER_HEADER_INCLUDED
