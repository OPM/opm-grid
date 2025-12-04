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

#ifndef OPM_GRID_CPGRID_LGROUTPUTHELPERS_HEADER_INCLUDED
#define OPM_GRID_CPGRID_LGROUTPUTHELPERS_HEADER_INCLUDED

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>
#include <opm/output/data/Cells.hpp>
#include <opm/output/data/Solution.hpp>
#include <opm/output/eclipse/RestartValue.hpp>

#include <cstddef>      // for std::size_t
#include <unordered_map>
#include <utility>      // for std::move
#include <type_traits>  // for std::is_same_v
#include <vector>

namespace Opm
{
namespace Lgr
{

/// @brief Builds a mapping from level element indices to the Cartesian ordering required by output files.
///
/// Output file solution data containers for LGRs expect elements to be ordered in strict Cartesian order,
/// with the i-direction varying fastest, followed by j, then k.
///
/// @param [in] grid          The CpGrid instance representing the simulation grid.
/// @param [in] levelCartMapp The LevelCartesianIndexMapper providing Cartesian index information across all levels.
/// @param [in] level         The grid level (integer) for which to build the index mapping.
/// @return A vector where each position corresponds to the output ordering, and the value is the element index
///         in that position. This ensures compatibility with the Cartesian ordering expected by output files.
std::vector<int> mapLevelIndicesToCartesianOutputOrder(const Dune::CpGrid& grid,
                                                       const Opm::LevelCartesianIndexMapper<Dune::CpGrid>& levelCartMapp,
                                                       int level);

/// @brief Reorder data from a simulation container into the order assumed by output for refined level grids.
///
/// @param [in] simulatorContainer  Container with simulation data ordered by compressed indices.
/// @param [in] toOutput            A vector where each position corresponds to the output ordering, and the
///                                 value is the element index in that position.
/// @return container with reordered data as expected by output, i.e. in strict Cartesian order, with the
///         i-direction varying fastest, followed by j, then k.
template <typename Container>
Container reorderForOutput(const Container& simulatorContainer,
                           const std::vector<int>& toOutput);

/// @brief Map level Cartesian index to level compressed index (active cell)
///
/// @param [in] grid
/// @param [in] levelCartesianIndexMapper
/// @return Vector of unordered maps, each entry holds the map for a level grid.
///         levelCartesianIndex (key) -> levelCompressedIndex (value).
std::vector<std::unordered_map<int,int>>
levelCartesianToLevelCompressedMaps(const Dune::CpGrid& grid,
                                    const Opm::LevelCartesianIndexMapper<Dune::CpGrid>& levelCartMapp);

/// @brief Extracts and organizes solution data for all grid refinement levels.
///
/// It derives these level-specific solutions from a given leaf-solution. For cells
/// that no longer exist in the leaf grid (i.e., parent cells that were refined away),
/// the default Scalar value (or a "rubbish std::numeric_limits<T>::max()") is assigned.
/// This behavior is temporary and will be replaced with proper restriction in a future
/// implementation.
///
/// Each resulting ScalarBuffer (std::vector<Scalar>) follows the ordering expected by
/// output files, i.e., increasing level Cartesian indices.
///
/// @param [in]       grid
/// @param [in]       toOutput_refinedLevels For level grids 1,2,..,maxLevel, a map to store
///                                          data in the order expected by outout files
///                                          (increasing level Cartesian indices)
/// @param [in]       leafSolution The complete solution defined on the leaf grid, containing
///                                one or more named data fields (e.g., pressure, saturation).
/// @return   A vector of Opm::data::Solution objects, one for each refinement level
///           (from level 0 to grid.maxLevel()), where each entry contains data reordered
///           according to increasing level Cartesian indices for output.
void extractSolutionLevelGrids(const Dune::CpGrid& grid,
                               const std::vector<std::vector<int>>& toOutput_refinedLevels,
                               const Opm::data::Solution& leafSolution,
                               std::vector<Opm::data::Solution>&);

/// @brief Constructs restart-value containers for all grid refinement levels.
///
/// The level-specific solution data are first derived from the leaf solution
/// using extractSolutionLevelGrids(...). Other data components (such as wells,
/// group/network values, and aquifers) are passed unchanged to each level.
///
/// @param [template] Grid The function has no effect for grids other than CpGrid.
/// @param [in]       grid
/// @param [in]       leafRestartValue
/// @param [out]      A vector of RestartValue objects, one for each refinement level
///                   (from level 0 to grid.maxLevel()).
template <typename Grid>
void extractRestartValueLevelGrids(const Grid& grid,
                                   const Opm::RestartValue& leafRestartValue,
                                   std::vector<Opm::RestartValue>& restartValue_levels);

template <typename Grid, typename TransmissibilityType, typename DirectVerticalNeighborsFunc>
void extractTransLevelGrids(const Grid& grid,
                            const TransmissibilityType& leafTrans,
                            std::vector<Opm::data::Solution>& outputTrans_levels,
                            const DirectVerticalNeighborsFunc& directVerticalNeighbors);

} // namespace Lgr
} // namespace Opm

template <typename Container>
Container Opm::Lgr::reorderForOutput(const Container& simulatorContainer,
                                     const std::vector<int>& toOutput)
{
    // Use toOutput to reorder simulatorContainer
    Container outputContainer;
    outputContainer.resize(toOutput.size());
    for (std::size_t i = 0; i < toOutput.size(); ++i) {
        outputContainer[i] = simulatorContainer[toOutput[i]];
    }
    return outputContainer;
}

template <typename Grid>
void Opm::Lgr::extractRestartValueLevelGrids(const Grid& grid,
                                             const Opm::RestartValue& leafRestartValue,
                                             std::vector<Opm::RestartValue>& restartValue_levels)
{
    if constexpr (std::is_same_v<Grid, Dune::CpGrid>) {

        int maxLevel = grid.maxLevel();
        restartValue_levels.resize(maxLevel+1); // level 0, 1, ..., max level
        
        // To store leafRestartValue.extra data in the order expected
        // by outout files (increasing level Cartesian indices)
        std::vector<std::vector<int>> toOutput_refinedLevels{};
        toOutput_refinedLevels.resize(maxLevel); // exclude level zero (does not need reordering)

        const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
        for (int level = 1; level <= maxLevel; ++level) { // exclude level zero (does not need reordering)
            toOutput_refinedLevels[level-1] = mapLevelIndicesToCartesianOutputOrder(grid, levelCartMapp, level);
        }

        std::vector<Opm::data::Solution> dataSolutionLevels{};
        extractSolutionLevelGrids(grid,
                                  toOutput_refinedLevels,
                                  leafRestartValue.solution,
                                  dataSolutionLevels);

        for (int level = 0; level <= maxLevel; ++level) { // exclude level zero (does not need reordering)
            restartValue_levels[level] = Opm::RestartValue(std::move(dataSolutionLevels[level]),
                                                           leafRestartValue.wells,
                                                           leafRestartValue.grp_nwrk,
                                                           leafRestartValue.aquifer,
                                                           level);
        }

        for (const auto& [rst_key, leafVector] : leafRestartValue.extra) {

            std::vector<std::vector<double>> levelVectors{};
            levelVectors.resize(maxLevel+1);

            // Parent cells do not appear on the leaf grid -> get rubbish value
            // std::numeric_limits<double>::max()
            for (int level = 0; level <= maxLevel; ++level) {
                levelVectors[level].resize(grid.levelGridView(level).size(0),
                                           std::numeric_limits<double>::max());
            }
            // For level cells that appear in the leaf, extract the data value from leafVector
            // and assign it to the equivalent level cell.
            for (const auto& element : Dune::elements(grid.leafGridView())) {
                levelVectors[element.level()][element.getLevelElem().index()] = leafVector[element.index()];
            }

            // Use toOutput_levels to reorder in ascending level cartesian indices
            for (int level = 1; level<=maxLevel; ++level) { // exclude level zero (does not need reordering)
                levelVectors[level] = Opm::Lgr::reorderForOutput(levelVectors[level], toOutput_refinedLevels[level-1]);
            }

            for (int level = 0; level <= maxLevel; ++level) {
                restartValue_levels[level].addExtra(rst_key.key, rst_key.dim, std::move(levelVectors[level]));
            }
        }
    }
}

template <typename Grid, typename TransmissibilityType, typename DirectVecticalNeighborsFunc>
void Opm::Lgr::extractTransLevelGrids(const Grid& grid,
                                      const TransmissibilityType& leafTrans,
                                      std::vector<Opm::data::Solution>& outputTrans_levels,
                                      const DirectVecticalNeighborsFunc& directVerticalNeighbors)

{
    if constexpr (std::is_same_v<Grid, Dune::CpGrid>) {

        int maxLevel = grid.maxLevel();
        outputTrans_levels.resize(maxLevel+1); // level 0, 1, ..., max level

        const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
        const auto levelCartToLevelCompressed = levelCartesianToLevelCompressedMaps(grid, levelCartMapp);

        // Extract trans values for level cells that appear in the leaf grid
        for (int level = 0; level <= maxLevel; ++level) {

            const auto& levelCartDims = levelCartMapp.cartesianDimensions(level);

            auto createCellData = [&levelCartDims]() {
                return Opm::data::CellData{
                    Opm::UnitSystem::measure::transmissibility,
                    std::vector<double>(levelCartDims[0] * levelCartDims[1] * levelCartDims[2], 0.0),
                    Opm::data::TargetType::INIT
                };
            };

            outputTrans_levels[level].clear();
            outputTrans_levels[level].emplace("TRANX", createCellData());
            outputTrans_levels[level].emplace("TRANY", createCellData());
            outputTrans_levels[level].emplace("TRANZ", createCellData());

            auto& tranx = outputTrans_levels[level].at("TRANX");
            auto& trany = outputTrans_levels[level].at("TRANY");
            auto& tranz = outputTrans_levels[level].at("TRANZ");

            for (const auto& element : Dune::elements(grid.levelGridView(level))) {
                if (!element.isLeaf()) // will be considered later
                    continue;

                for (const auto& intersection : Dune::intersections(grid.levelGridView(level), element)) {
                    if (!intersection.neighbor())
                        continue; // intersection is on the level-domain boundary
                    if (!intersection.outside().isLeaf())
                        continue; // now, we only care about pair of level cells that are leaf

                    assert(intersection.inside().level() == intersection.outside().level());

                    const unsigned levelIdxIn = intersection.inside().index();
                    const unsigned levelIdxOut= intersection.outside().index();

                    if (levelIdxIn > levelIdxOut)
                        continue; // we only need to handle each connection once.

                    const int levelCartIdxIn = levelCartMapp.cartesianIndex( levelIdxIn, level);
                    const int levelCartIdxOut = levelCartMapp.cartesianIndex( levelIdxOut, level);

                    int minLevelCartIdx = std::min(levelCartIdxIn, levelCartIdxOut);
                    int maxLevelCartIdx = std::max(levelCartIdxIn, levelCartIdxOut);

                    int leafIdxIn = grid.currentData()[level]->getLeafIdxFromLevelIdx(levelIdxIn);
                    int leafIdxOut = grid.currentData()[level]->getLeafIdxFromLevelIdx(levelIdxOut);

                    if (maxLevelCartIdx - minLevelCartIdx == 1 && levelCartDims[0] > 1 ) {
                        tranx.template data<double>()[minLevelCartIdx] = leafTrans.transmissibility(leafIdxIn, leafIdxOut);
                        continue; // skip other if clauses as they are false, last one needs some computation
                    }

                    if (maxLevelCartIdx - minLevelCartIdx == levelCartDims[0] && levelCartDims[1] > 1) {
                        trany.template data<double>()[minLevelCartIdx] = leafTrans.transmissibility(leafIdxIn, leafIdxOut);
                        continue; // skipt next if clause as it needs some computation
                    }

                    if ( maxLevelCartIdx - minLevelCartIdx == levelCartDims[0]*levelCartDims[1] ||
                         directVerticalNeighbors(levelCartDims,
                                                 levelCartToLevelCompressed[level],
                                                 minLevelCartIdx,
                                                 maxLevelCartIdx)) {
                        tranz.template data<double>()[minLevelCartIdx] = leafTrans.transmissibility(leafIdxIn, leafIdxOut);
                    }
                }
            } // end-for-loop-level
        } // end-if-Grid=CpGrid
    }
}

#endif // OPM_GRID_CPGRID_LGROUTPUTHELPERS_HEADER_INCLUDED
