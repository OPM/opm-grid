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
#ifndef OPM_GRID_REFINEMENT_REFINEDDISTRIBUTION_HEADER_INCLUDED
#define OPM_GRID_REFINEMENT_REFINEDDISTRIBUTION_HEADER_INCLUDED

#include <opm/grid/cpgrid/CpGridDataTraits.hpp>

#include <vector>

namespace Dune
{
class CpGrid;
namespace cpgrid { class CpGridData; }
}

/// Refine-before-redistribute: a Conforming grid refined before loadBalance() is partitioned on
/// level zero and its refined leaf is distributed, so the partitioner balances the refined cells.
namespace Opm::Refinement
{

/// The rank of each leaf cell: that of its level-zero ancestor in level0Part.
std::vector<int> leafPartitionFromLevelZero(const Dune::CpGrid& grid,
                                            const std::vector<int>& level0Part);

/// Give the serially assembled leaf the identity cell index set the scatter needs, and log the
/// number of leaf cells per rank.
void prepareLeafForScatter(Dune::cpgrid::CpGridData& leaf,
                           const std::vector<int>& leafPart,
                           const Dune::cpgrid::CpGridDataTraits::Communication& comm);

} // namespace Opm::Refinement

#endif // OPM_GRID_REFINEMENT_REFINEDDISTRIBUTION_HEADER_INCLUDED
