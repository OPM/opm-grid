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

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/common/version.hh>

#include <opm/grid/CpGrid.hpp>

#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>
#include <opm/grid/cpgrid/LgrHelpers.hpp>

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>
#include <opm/output/data/Solution.hpp>
#include <opm/output/data/Cells.hpp>

#include <array>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility> // for std::pair
#include <vector>

namespace Dune
{
namespace cpgrid
{
class CpGridData;
}
}

namespace Opm
{

template<class Grid> class LevelCartesianIndexMapper;

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


template <typename Scalar>
std::vector<Opm::data::Solution> extractDataForLevelGrids(const Dune::CpGrid& grid,
                                                          const Opm::data::Solution& leafSolution)
{
    int maxLevel = grid.maxLevel();

    // To restrict/create the level cell data, based on the leaf cells and the hierarchy
    std::vector<Opm::data::Solution> levelSolutions{};
    levelSolutions.resize(maxLevel+1);

    for (const auto& [name, leafCellData] : leafSolution)
    {
        const auto& leafVector = leafCellData.template data<Scalar>();
        const auto& messure =  leafCellData.dim;  // Opm::UnitSystem::measure;
        const auto& target = leafCellData.target; // Opm::data::TargetType>;

        if (leafVector.empty()) {
            continue;
        }

        std::vector<std::vector<Scalar>> levelVectors{};
        levelVectors.resize(maxLevel+1);

        const auto rubbish = -1;
        for (int level = 0; level <= maxLevel; ++level) {
            levelVectors[level].resize(grid.levelGridView(level).size(0), rubbish);
        }

        // For level cells that appear in the leaf, extract the data value from leafVector
        // and assign it the the equivalent level cell.
        // Notice that cells that vanished (parent cells) get the rubbish value.
        for (const auto& element : Dune::elements(grid.leafGridView())) {
            levelVectors[element.level()][element.getLevelElem().index()] = leafVector[element.index()];
        }

        // Reorder the containers in the order expected by outout files (increasing level Cartesian indices)
        const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
        for (int level = 0; level <= maxLevel; ++level) {

            const auto toOutput = Opm::Lgr::mapLevelIndicesToCartesianOutputOrder(grid, levelCartMapp, level);
            auto outputContainer = Opm::Lgr::reorderForOutput( levelVectors[level],
                                                               toOutput);

            levelSolutions[level].insert(name,
                                         messure,
                                         std::move(outputContainer),
                                         target);
        }
    }
    return levelSolutions;
}


} // namespace Lgr
} // namespace Opm


#endif // OPM_GRID_CPGRID_LGROUTPUTHELPERS_HEADER_INCLUDED
