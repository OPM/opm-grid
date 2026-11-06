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
#ifndef OPM_GRID_REFINEMENT_CONFORMINGBLOCKBUILDER_HEADER_INCLUDED
#define OPM_GRID_REFINEMENT_CONFORMINGBLOCKBUILDER_HEADER_INCLUDED

#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>
#include <opm/grid/cpgrid/refinement/RetainedCornerPointInput.hpp>

#include <array>
#include <memory>
#include <vector>

namespace Opm
{
namespace Refinement
{

/// The Conforming refinement backend. A box is a block of cells of its parent
/// grid that one request refines. Each box is resampled from the parent grid's
/// corner-point description (refineBlock), processed into the box's level grid
/// (assembleBlockLevelGrid), and the level grids are joined with level zero into
/// the leaf grid (LeafGridAssembler).
///
/// A builder made from a grid shares the corner-point description that the grid
/// keeps of its level zero; it does not copy it. The grid drops its own pointer
/// once it is refined, so the description is freed with the builder.
///
/// Restrictions, each refused with an exception:
/// - The grid is not refined yet.
/// - In a distributed grid, each box and the cells next to it are owned by one
///   rank.
/// - Boxes that meet need subdivisions that fit together. Boxes meet when their
///   index ranges are adjacent, so that they share a face, an edge or a corner
///   of the parent grid (boxes never share cells), or when their cells meet
///   through a fault face although their index ranges are not adjacent. On a
///   shared face, in each of the face's two directions, one box's number of
///   children per parent cell is a whole multiple of the other's (equal, or 4
///   and 2), and the same box is the finer one in both directions: the finer
///   box's sub-faces then tile the coarser box's faces. Along a shared edge the
///   numbers are equal.
/// - A graded box, one that gives the extent of each child column instead of
///   one number of children per parent cell, meets another box only at a corner.
class ConformingBlockBuilder : public Builder
{
public:
    /// Refines the grid that @p parentDims, @p coord, @p zcorn and @p actnum
    /// describe; the builder keeps the description.
    ///
    /// @param edgeConformal When true, the leaf is made edge-conformal: the
    ///        points a box adds on the edges of its sides are inserted into the
    ///        coarse faces that share those edges (edgeConformalizeLeaf). When
    ///        false, the default, those faces keep their own corners only.
    ConformingBlockBuilder(const std::array<int,3>& parentDims,
                           std::vector<double> coord,
                           std::vector<double> zcorn,
                           std::vector<int> actnum,
                           bool edgeConformal = false);

    /// Shares @p input, the corner-point description that the grid keeps of its
    /// level zero, which the grid drops once it is refined. The pointer is taken
    /// by value: the caller's pointer stays set unless the caller moves it in.
    explicit ConformingBlockBuilder(std::shared_ptr<const RetainedCornerPointInput> input);

    void build(Dune::CpGrid& grid,
               const std::vector<BlockRefinement>& requests) override;

private:
    std::shared_ptr<const RetainedCornerPointInput> input_;
};

} // namespace Refinement
} // namespace Opm

#endif // OPM_GRID_REFINEMENT_CONFORMINGBLOCKBUILDER_HEADER_INCLUDED
