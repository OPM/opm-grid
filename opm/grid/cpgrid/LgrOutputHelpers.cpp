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

#include <opm/input/eclipse/Units/UnitSystem.hpp>
#include <opm/output/data/Cells.hpp>
#include <opm/output/data/Solution.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/LgrOutputHelpers.hpp>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>

#include <algorithm> // for std::sort
#include <utility>   // for std::pair
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

void extractSolutionLevelGrids(const Dune::CpGrid& grid,
                               const Opm::data::Solution& leafSolution,
                               std::vector<Opm::data::Solution>& levelSolutions)

{
    int maxLevel = grid.maxLevel();
    const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);

    // To store in the order expected by outout files (increasing level Cartesian indices)
    std::vector<std::vector<int>> toOutput_levels{};
    toOutput_levels.resize(maxLevel); // exclude level zero (does not need reordering)
    for (int level = 1; level <= maxLevel; ++level) {
        toOutput_levels[level-1] = mapLevelIndicesToCartesianOutputOrder(grid, levelCartMapp, level);
    }

    // To restrict/create the level cell data, based on the leaf cells and the hierarchy
    levelSolutions.resize(maxLevel+1);

    for (const auto& [name, leafCellData] : leafSolution)
    {
        leafCellData.visit([&grid,
                            &maxLevel,
                            &toOutput_levels,
                            &levelSolutions,
                            &name,
                            &leafCellData](const auto& leafVector) {
            using T = std::decay_t<decltype(leafVector)>;

            if constexpr (std::is_same_v<T, std::monostate>) {
                // do nothing
            }
            else {
                if (!leafVector.empty()) {
                    std::vector<T> levelVectors{};
                    levelVectors.resize(maxLevel+1);
                    using ScalarType = std::decay_t<decltype(leafVector[0])>;

                    for (int level = 0; level <= maxLevel; ++level) {
                        // Parent cells get rubbish value: std::numeric_limist<ScalarType>::max()
                        levelVectors[level].resize(grid.levelGridView(level).size(0),
                                                   std::numeric_limits<ScalarType>::max());
                    }
                    // For level cells that appear in the leaf, extract the data value from leafVector
                    // and assign it the the equivalent level cell.
                    for (const auto& element : Dune::elements(grid.leafGridView())) {
                        levelVectors[element.level()][element.getLevelElem().index()] = leafVector[element.index()];
                    }
                    // Use toOutput_levels to reorder in ascending level cartesian indices
                    for (int level = 1; level<=maxLevel; ++level) { // exclude level zero (does not need reordering)
                        levelVectors[level] = Opm::Lgr::reorderForOutput(levelVectors[level], toOutput_levels[level-1]);
                    }

                    for (int level = 0; level <= maxLevel; ++level) {
                        if constexpr (std::is_same_v<T, std::vector<double>>) {
                            levelSolutions[level].insert(name,
                                                         leafCellData.dim,  // Opm::UnitSystem::measure
                                                         std::move(levelVectors[level]),
                                                         leafCellData.target); // Opm::data::TargetType>
                        }
                        else if constexpr (std::is_same_v<T, std::vector<int>>) {
                            levelSolutions[level].insert(name,
                                                         std::move(levelVectors[level]),
                                                         leafCellData.target); // Opm::data::TargetType>
                        }
                    }
                }
            }
        });
    }
}

} // namespace Lgr
} // namespace Opm
