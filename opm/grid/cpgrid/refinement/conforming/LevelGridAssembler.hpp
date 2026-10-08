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
#ifndef OPM_GRID_REFINEMENT_LEVELGRIDASSEMBLER_HEADER_INCLUDED
#define OPM_GRID_REFINEMENT_LEVELGRIDASSEMBLER_HEADER_INCLUDED

#include <opm/grid/cpgrid/refinement/RefinementRequest.hpp>
#include <opm/grid/cpgrid/refinement/conforming/GrdeclRefinement.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <array>
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

/// Stage 4+6a of the build pipeline (docs/DESIGN-builder.md D4): turn one
/// block request into a fully processed refined *level grid*.
///
/// The block is resampled onto sub-pillars (refineBlock), run through the
/// corner-point preprocessor (so faults and pinch-outs inside the block are
/// matched by the same code that handles them on level zero), and equipped
/// with the parent relations the rest of CpGrid reads (father(),
/// geometryInFather(), CARFIN-local Cartesian indices).
///
/// The returned grid has level_ = levelIndex and its level_data_ptr_ bound
/// to levelStorage; the caller is responsible for placing it at position
/// levelIndex in that vector (and level zero at position 0).
///
/// @param parentGrid   The grid the parents live in (level zero for a
///                     top-level box, or an LGR level grid for a nested box).
/// @param parentDims   Cartesian dimensions of the parent description.
/// @param coord,zcorn  Parent COORD/ZCORN arrays.
/// @param actnum       Parent ACTNUM or nullptr (all active).
/// @param request      The block to refine (IJK in @p parentGrid's space).
/// @param levelIndex   Index this level grid will occupy in the hierarchy.
/// @param parentLevel  Level index of @p parentGrid (0 = GLOBAL); recorded in
///                     the child->parent relation so nested grids point at the
///                     LGR level they refine rather than always level zero.
/// @param levelStorage The hierarchy vector the grid will live in.
/// @param comm         Communicator for the new grid object.
/// @param outRefined   If non-null, receives the resampled refined corner-point
///                     description of this level (its own children, if any, are
///                     refined from it).
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
                       RefinedBlockGrdecl* outRefined = nullptr);

/// An empty (zero-cell) refined level grid for a box that has no cells on
/// this rank (distributed runs, rank-interior LGRs). All ranks must carry
/// the same number of level grids, so absent boxes still need a placeholder
/// with consistent level index, subdivision factors and logical Cartesian
/// size.
std::shared_ptr<Dune::cpgrid::CpGridData>
assembleEmptyLevelGrid(const BlockRefinement& request,
                       int levelIndex,
                       std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& levelStorage,
                       Dune::MPIHelper::MPICommunicator comm);

} // namespace Refinement
} // namespace Opm

#endif // OPM_GRID_REFINEMENT_LEVELGRIDASSEMBLER_HEADER_INCLUDED
