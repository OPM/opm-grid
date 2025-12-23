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
#include <opm/grid/cpgrid/CartesianIndexMapper.hpp>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>
#include <opm/output/data/Cells.hpp>
#include <opm/output/data/Solution.hpp>
#include <opm/output/eclipse/RestartValue.hpp>

#include <algorithm>    // for std::min/max
#include <cstddef>      // for std::size_t
#include <unordered_map>
#include <utility>      // for std::move
#include <type_traits>  // for std::is_same_v
#include <vector>

template<class ScalarType>
struct MinChildrenData{
    void operator()(int level,
                    int levelIdx,
                    const std::vector<std::vector<ScalarType>>& levelVectors)
    {
        min = std::min(min, levelVectors[ level ][ levelIdx ]);
    }

    ScalarType getValue()
    {
        return min;
    }

    ScalarType min{};
};

template<class ScalarType>
struct MaxChildrenData{
    void operator()(int level,
                    int levelIdx,
                    const std::vector<std::vector<ScalarType>>& levelVectors)
    {
        max = std::max(max, levelVectors[ level ][ levelIdx ]);
    }

    ScalarType getValue()
    {
        return max;
    }

    ScalarType max{};
};

template<class ScalarType>
struct AverageChildrenData{
    void operator()(int level,
                    int levelIdx,
                    const std::vector<std::vector<ScalarType>>& levelVectors)
    {
        partialSum+= levelVectors[ level ][ levelIdx ];
        ++count;
    }

    ScalarType getValue()
    {
        return partialSum/count;
    }

    std::size_t count{};
    ScalarType partialSum{};
};

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

/// @breif Compute the parent data from the children data values
///
/// @param [in] levelVectors A collection of per-level data vectors.
///                          levelVectors[l][i] contains the data value for the element at level l
///                          with level index i. To extract the children data.
/// @param [in] element      A non-leaf element (parent cell). Its children data
///                          will be used to compute the average.
/// @param [in] grid         The grid from which the element is taken.
/// @return Parent data computed from children data values.
template <typename ScalarType>
ScalarType processChildrenData(const std::vector<std::vector<ScalarType>>& levelVectors,
                               const Dune::cpgrid::Entity<0>& element,
                               const Dune::CpGrid& grid);

/// @brief Populate level data vectors based on leaf vector, for a specific named data field.
///
/// @param [in]       grid
/// @param [in]       maxLevel
/// @param [in]       leafVector Containing one named data field (e.g., pressure, saturation).
/// @param [in]       toOutput_refinedLevels For level grids 1,2,..,maxLevel, a map to store
///                                          data in the order expected by outout files
///                                          (increasing level Cartesian indices).
/// @param [out]      levelVectors A collection of per-level data vectors.
///                                levelVectors[l][i] contains the data value for the element at level l
///                                with level index i.
template <typename ScalarType>
void populateDataVectorLevelGrids(const Dune::CpGrid& grid,
                                  int maxLevel,
                                  const std::vector<ScalarType>& leafVector,
                                  const std::vector<std::vector<int>>& toOutput_refinedLevels,
                                  std::vector<std::vector<ScalarType>>& levelVectors);

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
/// the average of its children values is assigned.
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
/// @param [out]   A vector of Opm::data::Solution objects, one for each refinement level
///                (from level 0 to grid.maxLevel()), where each entry contains data reordered
///                according to increasing level Cartesian indices for output.
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

/// @brief Constructs TRANS* values for level zero cells, that appear on the leaf.
///
/// @param [in] grid
/// @param [in] leafTrans
/// @param [out] tranx       Opm::data::CellData for TRANX
/// @param [out] trany       Opm::data::CellData for TRANY
/// @param [out] tranz       Opm::data::CellData for TRANZ
///                          tranx, trany, and tranz only contain TRANS* values,
///                          between level zero cells sharing an intersection, i.e.,
///                          not at the boundary of an LGR, or the grid domain.
/// @param [in] cartMapp     CartesianIndexMapper of the grid
/// @param [in] cartDims     Cartesian dimensions of the grid (coincide with level zero grid).
/// @param [in] directVerticalNeighbors
template <typename Grid, typename TransmissibilityType, typename CartesianMapper,
          typename IsNumAquCell, typename DirectVerticalNeighborsFunc>
void extractTransLevelZero(const Grid& grid,
                           const TransmissibilityType& leafTrans,
                           Opm::data::CellData& tranx,
                           Opm::data::CellData& trany,
                           Opm::data::CellData& tranz,
                           const CartesianMapper& cartMapp,
                           const std::array<int,3>& cartDims,
                           const IsNumAquCell& isNumAquCell,
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

template <typename ScalarType>
ScalarType Opm::Lgr::processChildrenData(const std::vector<std::vector<ScalarType>>& levelVectors,
                                         const Dune::cpgrid::Entity<0>& element,
                                         const Dune::CpGrid& grid)
{
    const auto& [level, level_indices] = grid.currentData()[element.level()]->getChildrenLevelAndIndexList(element.index());

    // For now, we compute field output properties as follows:
    // - int fields: use the maximum value among the children data
    // - double fields: use the average of the children data
    //
    // In the future, this could be extended to use the std::string field name
    // to apply field-specific computation logic.
    auto childrenDataFunc = [](){
        if constexpr (std::is_same_v<ScalarType, double>)
            return AverageChildrenData<ScalarType>{};
        else
            return MaxChildrenData<ScalarType>{};
    }();

    for (const auto& levelIdx : level_indices) {
        childrenDataFunc(level, levelIdx, levelVectors);
    }
    return childrenDataFunc.getValue();
}

template <typename ScalarType>
void Opm::Lgr::populateDataVectorLevelGrids(const Dune::CpGrid& grid,
                                            int maxLevel,
                                            const std::vector<ScalarType>& leafVector,
                                            const std::vector<std::vector<int>>& toOutput_refinedLevels,
                                            std::vector<std::vector<ScalarType>>& levelVectors)
{
    for (int level = 0; level <= maxLevel; ++level) {
        levelVectors[level].resize(grid.levelGridView(level).size(0));
    }
    // For level cells that appear in the leaf, extract the data value from leafVector
    // and assign it to the equivalent level cell.
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        levelVectors[element.level()][element.getLevelElem().index()] = leafVector[element.index()];
    }
    // Note that all cells from maxLevel have assigned values at this point.
    // Now, assign values for parent cells (for now, average of children values).
    if (maxLevel)  {
        for (int level = maxLevel-1; level >= 0; --level) {
            for (const auto& element : Dune::elements(grid.levelGridView(level))) {
                if (!element.isLeaf()) {
                    levelVectors[level][element.index()] = Opm::Lgr::processChildrenData(levelVectors,
                                                                                         element,
                                                                                         grid);
                }
            }
        }
    }
    // Use toOutput_levels to reorder in ascending level cartesian indices
    for (int level = 1; level<=maxLevel; ++level) { // exclude level zero (does not need reordering)
        levelVectors[level] = Opm::Lgr::reorderForOutput(levelVectors[level], toOutput_refinedLevels[level-1]);
    }
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
            toOutput_refinedLevels[level-1] = Opm::Lgr::mapLevelIndicesToCartesianOutputOrder(grid, levelCartMapp, level);
        }

        std::vector<Opm::data::Solution> dataSolutionLevels{};
        extractSolutionLevelGrids(grid,
                                  toOutput_refinedLevels,
                                  leafRestartValue.solution,
                                  dataSolutionLevels);



        for (int level = 0; level <= maxLevel; ++level) {
            restartValue_levels[level] = Opm::RestartValue(std::move(dataSolutionLevels[level]),
                                                           leafRestartValue.wells,
                                                           leafRestartValue.grp_nwrk,
                                                           leafRestartValue.aquifer,
                                                           level);
        }

        for (const auto& [rst_key, leafVector] : leafRestartValue.extra) {

            std::vector<std::vector<double>> levelVectors{};
            levelVectors.resize(maxLevel+1);

            Opm::Lgr::populateDataVectorLevelGrids<double>(grid,
                                                           maxLevel,
                                                           leafVector,
                                                           toOutput_refinedLevels,
                                                           levelVectors);

            for (int level = 0; level <= maxLevel; ++level) {
                restartValue_levels[level].addExtra(rst_key.key, rst_key.dim, std::move(levelVectors[level]));
            }
        }
    }
}

template <typename Grid, typename TransmissibilityType, typename CartesianMapper,
          typename IsNumAquCell, typename DirectVecticalNeighborsFunc>
void Opm::Lgr::extractTransLevelZero(const Grid& grid,
                                     const TransmissibilityType& leafTrans,
                                     Opm::data::CellData& tranx,
                                     Opm::data::CellData& trany,
                                     Opm::data::CellData& tranz,
                                     const CartesianMapper& cartMapp,
                                     const std::array<int,3>& cartDims,
                                     const IsNumAquCell& isNumAquCell,
                                     const DirectVecticalNeighborsFunc& directVerticalNeighbors)

{
    if constexpr (std::is_same_v<Grid, Dune::CpGrid>) {
        // To detect direct vertical neighboring cells:
        const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
        const auto levelCartToLevelCompressed = Opm::Lgr::levelCartesianToLevelCompressedMaps(grid, levelCartMapp);

        // Extract transmissibility values for intersections between leaf cells
        // belonging to level zero grid
        for (const auto& element : Dune::elements(grid.leafGridView())) {
            for (const auto& intersection : Dune::intersections(grid.leafGridView(), element)) {
                if (!intersection.neighbor())  // intersection is on the domain boundary
                    continue;
                if ((intersection.inside().level()>0) || (intersection.outside().level()>0))
                    continue; // for now, we only care about level zero cells

                const unsigned leafIdxIn = intersection.inside().index();
                const unsigned leafIdxOut = intersection.outside().index();

                // Leaf cells that belong to level zero grid, in particular,
                // they have diffent 'origin' cell in level zero. Therefore,
                // different Cartesian Index.
                const int cartIdxIn = cartMapp.cartesianIndex(leafIdxIn);
                const int cartIdxOut = cartMapp.cartesianIndex(leafIdxOut);

                if (isNumAquCell(cartIdxIn) || isNumAquCell(cartIdxOut)) {
                    // Connections involving numerical aquifers are always NNCs
                    // for the purpose of file output.  This holds even for
                    // connections between cells like (I,J,K) and (I+1,J,K)
                    // which are nominally neighbours in the Cartesian grid.
                    continue;
                }

                if (cartIdxIn > cartIdxOut)
                    continue; // we only need to handle each connection once.

                // We use min and max level zero Cartesian indices to distinguish
                // X, Y, and Z transmissibility values.

                if (cartIdxOut - cartIdxIn == 1 && cartDims[0] > 1 ) {
                    tranx.template data<double>()[cartIdxIn] = leafTrans.transmissibility(leafIdxIn, leafIdxOut);
                    continue;
                }

                if (cartIdxOut - cartIdxIn == cartDims[0] && cartDims[1] > 1) {
                    trany.template data<double>()[cartIdxIn] = leafTrans.transmissibility(leafIdxIn, leafIdxOut);
                    continue;
                }

                if (cartIdxOut - cartIdxIn == cartDims[0]*cartDims[1] ||
                    directVerticalNeighbors(cartDims,
                                            levelCartToLevelCompressed[/*level = */0],
                                            cartIdxIn,   // min
                                            cartIdxOut)) // max
                {
                    tranz.template data<double>()[cartIdxIn] = leafTrans.transmissibility(leafIdxIn, leafIdxOut);
                }
            } // end-intersection-loop
        } // end-element-loop
    } // end-if-Grid==CpGrid
}

#endif // OPM_GRID_CPGRID_LGROUTPUTHELPERS_HEADER_INCLUDED
