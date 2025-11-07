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



std::vector<std::vector<int>> ifLeafGetLevelIndexToLeafIndex(const Dune::CpGrid& grid)
{
    int maxLevel = grid.maxLevel();

    int rubbish = -1; // value for vanished cells (i.e. parent cells, not appearing on the leaf)
    
    std::vector<std::vector<int>> levelIdx_to_leafIdx{};
    levelIdx_to_leafIdx.resize(maxLevel +1);
    
    for (int level = 0; level <= maxLevel; ++level) {
        levelIdx_to_leafIdx[level].resize(grid.levelGridView(level).size(0), rubbish);
    }
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        levelIdx_to_leafIdx[ element.level() ][ element.getLevelElem().index() ] = element.index();
    }
    return levelIdx_to_leafIdx;
}

std::vector<std::array<int,2>> ifVanishedGetDescendantsPresentOnLeaf(const Dune::cpgrid::Entity<0>& element,
                                                                     int maxLevel)
{
    if (element.isLeaf()) {
        throw;
    }
    
    std::vector<std::array<int,2>> refinedChildrenPresentOnLeaf{}; // {level, child level index}
   
         
    auto it = element.hbegin(maxLevel);
    const auto& endIt = element.hend(maxLevel);
    
    for (; it != endIt; ++it) {
        if (it->isLeaf())
            refinedChildrenPresentOnLeaf.emplace_back(std::array{it->level(), it-> index()});
    }
    return refinedChildrenPresentOnLeaf;
}

std::vector<std::vector<std::vector<int>>> getDescendantsPresentOnLeaf(const Dune::CpGrid& grid)
{
    int maxLevel = grid.maxLevel();
    std::vector<std::vector<std::vector<int>>> descendantsOnLeaf{};
    descendantsOnLeaf.resize(maxLevel+1);
    
    const auto levelIdx_to_leafIdx = ifLeafGetLevelIndexToLeafIndex(grid);

    for (int level = 0; level <= maxLevel; ++level) {
        descendantsOnLeaf[level].resize(grid.levelGridView(level).size(0));
        
        for (const auto& element : Dune::elements(grid.levelGridView(level))) {
            if (element.isLeaf()) {
                descendantsOnLeaf[level][element.index()].emplace_back(levelIdx_to_leafIdx[ level ][element.index()]);
            }
            else {
                for (const auto& level_levelIdx : ifVanishedGetDescendantsPresentOnLeaf(element, grid.maxLevel())) {
                    descendantsOnLeaf[level][element.index()].emplace_back(levelIdx_to_leafIdx[ level_levelIdx[0] ][level_levelIdx[1]]);
                }
            }
        }
    }
    return descendantsOnLeaf;
}




} // namespace Lgr
} // namespace Opm
