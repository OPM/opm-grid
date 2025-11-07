/*
  Copyright 2025 Equinor ASA.

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
#include "config.h"

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/cpgrid/LgrHelpers.hpp>
#include <opm/grid/cpgrid/LgrOutputHelpers.hpp>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>
#include <opm/grid/cpgrid/ParentToChildCellToPointGlobalIdHandle.hpp>

#include <algorithm>    // for std::max
#include <array>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>  // for std::integral_constant
#include <unordered_map>
#include <utility> // for std::pair
#include <vector>

namespace Opm
{
namespace Lgr
{

std::vector<int> mapLevelIndicesToCartesianOutputOrder(const Dune::CpGrid& grid,
                                                       const Opm::LevelCartesianIndexMapper<Dune::CpGrid>& levelCartMapp,
                                                       int level)
{
    const int lgr_cells = grid.levelGridView(level).size(0);

    std::vector<std::pair<int,int>> sorted_levelIdxToLevelCartIdx;
    sorted_levelIdxToLevelCartIdx.reserve(lgr_cells);

    for (const auto& element : Dune::elements(grid.levelGridView(level))) {
        sorted_levelIdxToLevelCartIdx.push_back(std::make_pair(element.index(),
                                                               levelCartMapp.cartesianIndex(element.index(),level)));
    }

    std::sort(sorted_levelIdxToLevelCartIdx.begin(), sorted_levelIdxToLevelCartIdx.end(),
              [](std::pair<int,int>& a, std::pair<int,int>& b) {
                  return a.second < b.second;
              });

    std::vector<int> toOutput; // Consecutive numbers, from 0 to total elemts in LGR1 -1
    toOutput.reserve(lgr_cells);

    // Redefinition of level Cartesian indices (output-style)
    for (const auto& [sorted_elemIdx, sorted_cartIdx] : sorted_levelIdxToLevelCartIdx) {
        toOutput.push_back(sorted_elemIdx);
    }
    return toOutput;
}


} // namespace Lgr
} // namespace Opm
