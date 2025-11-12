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
#include <opm/grid/cpgrid/LgrOutputHelpers.hpp>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>

#include <algorithm> // for std::sort
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

std::vector<Opm::data::Solution> extractSolutionLevelGrids(const Dune::CpGrid& grid,
                                                           const Opm::data::Solution& leafSolution)

{
    int maxLevel = grid.maxLevel();

    // To restrict/create the level cell data, based on the leaf cells and the hierarchy
    std::vector<Opm::data::Solution> levelSolutions{};
    levelSolutions.resize(maxLevel+1);

    for (const auto& [name, leafCellData] : leafSolution)
    {
        const auto& leafVector = leafCellData.template data<double>();
        const auto& measure =  leafCellData.dim;  // Opm::UnitSystem::measure;
        const auto& target = leafCellData.target; // Opm::data::TargetType>;

        if (leafVector.empty()) {
            continue;
        }

        std::vector<std::vector<double>> levelVectors{};
        levelVectors.resize(maxLevel+1);

        const auto rubbish = -1;
        for (int level = 0; level <= maxLevel; ++level) {
            levelVectors[level].resize(grid.levelGridView(level).size(0), rubbish);
        }

        // For level cells that appear in the leaf, extract the data value from leafVector
        // and assign it the the equivalent level cell.
        // Notice that cells that vanished (parent cells) get the rubbish value.
        // Store in the order expected by outout files (increasing level Cartesian indices)
        const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
        for (const auto& element : Dune::elements(grid.leafGridView())) {
            int levelCartIdx = levelCartMapp.cartesianIndex(element.getLevelElem().index(), element.level());
            levelVectors[element.level()][levelCartIdx] = leafVector[element.index()];
        }

        for (int level = 0; level <= maxLevel; ++level) {
            levelSolutions[level].insert(name,
                                         measure,
                                         std::move(levelVectors[level]),
                                         target);
        }
    }
    return levelSolutions;
}

std::vector<Opm::RestartValue> getRestartValueLevelGrids(const Dune::CpGrid& grid,
                                                         const Opm::RestartValue& leafRestartValue)
{
    int maxLevel = grid.maxLevel();
    std::vector<Opm::RestartValue> restartValue_levels{};
    restartValue_levels.resize(maxLevel+1); // level 0, 1, ..., max level

    const auto dataSolutionLevels = extractSolutionLevelGrids(grid, leafRestartValue.solution);

    for (int level = 0; level <= maxLevel; ++level) {
        restartValue_levels[level] = Opm::RestartValue(dataSolutionLevels[level],
                                                       leafRestartValue.wells,
                                                       leafRestartValue.grp_nwrk,
                                                       leafRestartValue.aquifer,
                                                       level);
    }

    for (const auto& [rst_key, leafVector] : leafRestartValue.extra) {

        std::vector<std::vector<double>> levelVectors{};
        levelVectors.resize(maxLevel+1);

        const auto rubbish = -1;
        for (int level = 0; level <= maxLevel; ++level) {
            levelVectors[level].resize(grid.levelGridView(level).size(0), rubbish);
        }

        // For level cells that appear in the leaf, extract the data value from leafVector
        // and assign it the the equivalent level cell.
        // Notice that cells that vanished (parent cells) get the rubbish value.
        // Store in the order expected by outout files (increasing level Cartesian indices)
        const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
        for (const auto& element : Dune::elements(grid.leafGridView())) {
            int levelCartIdx = levelCartMapp.cartesianIndex(element.getLevelElem().index(), element.level());
            levelVectors[element.level()][levelCartIdx] = leafVector[element.index()];
        }

        for (int level = 0; level <= maxLevel; ++level) {
            restartValue_levels[level].addExtra(rst_key.key, rst_key.dim, levelVectors[level]);
        }
    }

    return restartValue_levels;
}


} // namespace Lgr
} // namespace Opm
