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
#ifndef OPM_GRID_REFINEMENT_REQUEST_HEADER_INCLUDED
#define OPM_GRID_REFINEMENT_REQUEST_HEADER_INCLUDED

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace Opm
{
namespace Refinement
{

/// Where each refined column of one direction sits: the parent cell of the box
/// it lies in (0-based within the box) and its normalised extent within that
/// cell. One entry per refined column.
///
/// This is what N*FIN/H*FIN mean. A uniform box leaves it empty and is
/// described by BlockRefinement::cellsPerDim alone; axisSubdivision() fills in
/// the uniform tables so the builder has one code path.
struct AxisSubdivision
{
    std::vector<int> parentOffset{};
    std::vector<double> fracLo{};
    std::vector<double> fracHi{};

    bool empty() const { return parentOffset.empty(); }
    std::size_t size() const { return parentOffset.size(); }
};

/// One block-shaped static refinement request (one CARFIN box).
///
/// IJK indices refer to the logical Cartesian space of the parent grid
/// ("GLOBAL" for level zero). The box is half-open in cell terms: the last
/// refined cell column is endIJK - 1, matching CpGrid::addLgrsUpdateLeafView.
struct BlockRefinement
{
    std::string name;
    std::string parentGridName{"GLOBAL"};
    std::array<int,3> cellsPerDim{};
    std::array<int,3> startIJK{};
    std::array<int,3> endIJK{};
    /// Graded subdivision, one entry per direction. Empty means uniform.
    std::array<AxisSubdivision,3> subdivision{};
    /// Refined cells the block's own MINPV removes, one entry per refined
    /// Cartesian cell (1 = removed). Empty when the block set no MINPV.
    /// Decided in opm-common, where the father's pore volume lives, so that the
    /// simulation grid and the LGR grids written to the EGRID agree on which
    /// refined cells exist.
    std::vector<int> minpvRemoved{};
    /// Refined pillars from the box's own first layer (as the reference) rather
    /// than the column's; boxes stacked on sheared pillars then do not conform.
    bool pillarsFromBoxLayer{false};
};

/// The request's subdivision in one direction, uniform tables filled in when
/// the request carries none.
AxisSubdivision axisSubdivision(const BlockRefinement& request, int dim);

/// Refined dimensions of the block: the subdivision sizes, which for a uniform
/// request are box size times cellsPerDim.
std::array<int,3> refinedDims(const BlockRefinement& request);

/// Where each refined column sits inside its parent cell. subIndex and
/// parentCount are indexed by refined column; firstColumn by parent cell.
struct AxisPositions
{
    std::vector<int> subIndex{};
    std::vector<int> parentCount{};
    std::vector<int> firstColumn{};
};

AxisPositions axisPositions(const AxisSubdivision& sub);

/// Upstream's vector form of a refinement as block requests; parent names
/// default to "GLOBAL". Throws std::invalid_argument when the sizes differ.
std::vector<BlockRefinement>
blockRefinements(const std::vector<std::array<int,3>>& cellsPerDim,
                 const std::vector<std::array<int,3>>& startIJK,
                 const std::vector<std::array<int,3>>& endIJK,
                 const std::vector<std::string>& names,
                 const std::vector<std::string>& parentGridNames);

/// Validate a set of block requests: matching sizes are the caller's
/// responsibility; this checks per-box consistency (start < end, positive
/// subdivisions) and pairwise disjointness of boxes sharing a parent grid.
/// Throws std::invalid_argument with a message naming the offending box(es).
void validateBlockRefinements(const std::vector<BlockRefinement>& requests);

} // namespace Refinement
} // namespace Opm

#endif // OPM_GRID_REFINEMENT_REQUEST_HEADER_INCLUDED
