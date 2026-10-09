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

#include <opm/common/ErrorMacros.hpp>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/LgrOutputHelpers.hpp>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>

#include <algorithm> // for std::sort
#include <utility>   // for std::pair
#include <string>
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

    std::ranges::sort(sorted_levelIdxToLevelCartIdx,
                      [](const auto& a, const auto& b)
                      { return a.second < b.second; });

    std::vector<int> toOutput; // Consecutive numbers, from 0 to total elemts in LGR1 -1
    toOutput.reserve(lgr_cells);

    // Redefinition of level Cartesian indices (output-style)
    for (const auto& [sorted_elemIdx, sorted_cartIdx] : sorted_levelIdxToLevelCartIdx) {
        toOutput.push_back(sorted_elemIdx);
    }
    return toOutput;
}

std::vector<std::unordered_map<int,int>> levelCartesianToLevelCompressedMaps(const Dune::CpGrid& grid,
                                                                             const Opm::LevelCartesianIndexMapper<Dune::CpGrid>& levelCartMapp)
{
    std::vector<std::unordered_map<int,int>> levelCartToLevelCompressed{};
    levelCartToLevelCompressed.resize(grid.maxLevel()+1);
    for (int level = 0; level <= grid.maxLevel(); ++level) {
        levelCartToLevelCompressed[level].reserve(grid.levelGridView(level).size(0));
        for (const auto& element : Dune::elements(grid.levelGridView(level))) {
            levelCartToLevelCompressed[level].emplace(levelCartMapp.cartesianIndex(element.index(), level), element.index());
        }
    }
    return levelCartToLevelCompressed;
}

void assembleSolutionFromLevelGrids(const Dune::CpGrid& grid,
                                    const std::vector<Opm::data::Solution>& levelSolutions,
                                    Opm::data::Solution& leafSolution)
{
    const int maxLevel = grid.maxLevel();
    if (static_cast<int>(levelSolutions.size()) != maxLevel + 1) {
        OPM_THROW(std::invalid_argument,
                  "Restart solution has " + std::to_string(levelSolutions.size())
                  + " level sections for a grid with " + std::to_string(maxLevel + 1)
                  + " levels.");
    }

    // The permutation the writer applied to every level above zero.
    const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
    std::vector<std::vector<int>> toOutput_refinedLevels(maxLevel);
    for (int level = 1; level <= maxLevel; ++level) {
        toOutput_refinedLevels[level-1] =
            mapLevelIndicesToCartesianOutputOrder(grid, levelCartMapp, level);
    }

    const auto leafSize = static_cast<std::size_t>(grid.leafGridView().size(0));

    for (const auto& [name, levelZeroData] : levelSolutions.front()) {
        // A key has to be on every level to be assembled; the writer emits them
        // together, so a missing one means the file was not written by this path.
        const bool everywhere =
            std::all_of(levelSolutions.begin() + 1, levelSolutions.end(),
                        [&name = name](const Opm::data::Solution& sol)
                        { return sol.has(name); });
        if (! everywhere) {
            continue;
        }

        levelZeroData.visit([&](const auto& levelZeroVector) {
            using T = std::decay_t<decltype(levelZeroVector)>;

            if constexpr (! std::is_same_v<T, std::monostate>) {
                // Undo the output ordering, then let each leaf cell take the
                // value its own level holds for it.
                std::vector<T> levelVectors(maxLevel + 1);
                levelVectors[0] = levelZeroVector;
                for (int level = 1; level <= maxLevel; ++level) {
                    levelVectors[level] =
                        reorderFromOutput(levelSolutions[level].data<typename T::value_type>(name),
                                          toOutput_refinedLevels[level-1]);
                }

                T leafVector(leafSize);
                for (const auto& element : Dune::elements(grid.leafGridView())) {
                    leafVector[element.index()] =
                        levelVectors[element.level()][element.getLevelElem().index()];
                }

                if constexpr (std::is_same_v<T, std::vector<double>>) {
                    leafSolution.insert(name, levelZeroData.dim,
                                        std::move(leafVector), levelZeroData.target);
                }
                else {
                    leafSolution.insert(name, std::move(leafVector), levelZeroData.target);
                }
            }
        });
    }
}

void extractSolutionLevelGrids(const Dune::CpGrid& grid,
                               const std::vector<std::vector<int>>& toOutput_refinedLevels,
                               const Opm::data::Solution& leafSolution,
                               std::vector<Opm::data::Solution>& levelSolutions)

{
    int maxLevel = grid.maxLevel();
    // To restrict/create the level cell data, based on the leaf cells and the hierarchy
    levelSolutions.resize(maxLevel+1);

    for (const auto& leaf : leafSolution)
    {
        const auto& name = leaf.first;
        const auto& leafCellData = leaf.second;
        leafCellData.visit([&grid,
                            &maxLevel,
                            &toOutput_refinedLevels,
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

                    populateDataVectorLevelGrids<ScalarType>(grid,
                                                             maxLevel,
                                                             leafVector,
                                                             toOutput_refinedLevels,
                                                             levelVectors);

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
