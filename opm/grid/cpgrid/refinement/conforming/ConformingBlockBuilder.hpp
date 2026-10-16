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
#ifndef OPM_GRID_REFINEMENT_CONFORMINGBLOCKBUILDER_HEADER_INCLUDED
#define OPM_GRID_REFINEMENT_CONFORMINGBLOCKBUILDER_HEADER_INCLUDED

#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>

#include <array>
#include <vector>

namespace Opm
{
namespace Refinement
{

/// First refinement backend (docs/PLAN.md Track 1 milestone 3): conforming
/// block refinement built through the preprocessor pipeline
/// (refineBlock -> assembleBlockLevelGrid -> assembleLeafGrid).
///
/// The builder owns a copy of the parent corner-point description, since
/// the grid does not reliably retain COORD/ZCORN after construction.
///
/// Restrictions (throw): an unrefined starting grid; in parallel each box must be
/// rank-interior; boxes that touch, in IJK or across a fault, need in-face
/// subdivisions that nest; graded boxes may not touch another box.
class ConformingBlockBuilder : public Builder
{
public:
    /// @param edgeConformal When true, the refined leaf is made
    ///        edge-conformal (an edge-conformalization post-pass that
    ///        inserts the refinement's boundary nodes into the coarse faces
    ///        sharing those edges). Default false keeps the current,
    ///        face-conformal-only behaviour.
    ConformingBlockBuilder(const std::array<int,3>& parentDims,
                           std::vector<double> coord,
                           std::vector<double> zcorn,
                           std::vector<int> actnum,
                           bool edgeConformal = false);

    void build(Dune::CpGrid& grid,
               const std::vector<BlockRefinement>& requests) override;

private:
    std::array<int,3> dims_;
    std::vector<double> coord_;
    std::vector<double> zcorn_;
    std::vector<int> actnum_;
    bool edgeConformal_;
};

} // namespace Refinement
} // namespace Opm

#endif // OPM_GRID_REFINEMENT_CONFORMINGBLOCKBUILDER_HEADER_INCLUDED
