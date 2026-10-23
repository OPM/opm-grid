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
#ifndef OPM_GRID_REFINEMENT_RETAINEDCORNERPOINTINPUT_HEADER_INCLUDED
#define OPM_GRID_REFINEMENT_RETAINEDCORNERPOINTINPUT_HEADER_INCLUDED

#include <array>
#include <set>
#include <utility>
#include <vector>

namespace Opm
{
namespace Refinement
{

/// The corner-point description a grid was built from, retained only when
/// the deck requests LGRs (docs/DESIGN-builder.md D4 input note): the
/// refinement builder needs the post-MINPV COORD/ZCORN that level zero was
/// actually constructed from, and the grid does not otherwise keep them.
struct RetainedCornerPointInput
{
    std::array<int,3> dims{};
    std::vector<double> coord;
    std::vector<double> zcorn;
    std::vector<int> actnum;
    /// True if level zero was built edge-conformal: the refined leaf is then
    /// made edge-conformal too (an edge-conformalization post-pass).
    bool edgeConformal{false};
    /// Level zero's pinch-out connections (Cartesian pairs) and whether PINCH
    /// was active, so a grid rebuilt from this description has them too.
    std::set<std::pair<int,int>> pinchNnc{};
    bool pinchActive{false};
};

} // namespace Refinement
} // namespace Opm

#endif // OPM_GRID_REFINEMENT_RETAINEDCORNERPOINTINPUT_HEADER_INCLUDED
