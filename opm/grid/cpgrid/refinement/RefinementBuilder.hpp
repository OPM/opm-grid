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
#ifndef OPM_GRID_REFINEMENT_BUILDER_HEADER_INCLUDED
#define OPM_GRID_REFINEMENT_BUILDER_HEADER_INCLUDED

#include <opm/grid/cpgrid/refinement/RefinementRequest.hpp>

#include <memory>
#include <vector>

namespace Dune
{
class CpGrid;
}

namespace Opm
{
namespace Refinement
{

/// Interface for static refinement backends ("the seam", docs/PLAN.md
/// Track 1 milestone 2; design rationale in docs/lgr_review.md §6, §IV.2).
///
/// CpGrid::addLgrsUpdateLeafView() validates the request and forwards it
/// here. A backend is free to refine in place or to rebuild the grid's
/// level hierarchy and leaf view from scratch (e.g. the planned
/// preprocessor-level builder, review §8.2-S2). Contract:
///  - on success, grid.maxLevel() reflects the new refined level grids and
///    the leaf view is consistent (parent/child maps, ids, index sets);
///  - on failure, throw; the grid must be left unchanged;
///  - requests have already passed validateBlockRefinements().
class Builder
{
public:
    virtual ~Builder() = default;

    virtual void build(Dune::CpGrid& grid,
                       const std::vector<BlockRefinement>& requests) = 0;
};

/// Which algorithm refines a CpGrid: upstream's cell-by-cell trilinear
/// refinement (the default), or a conforming Builder.
enum class Backend { Trilinear, Conforming };

} // namespace Refinement
} // namespace Opm

#endif // OPM_GRID_REFINEMENT_BUILDER_HEADER_INCLUDED
