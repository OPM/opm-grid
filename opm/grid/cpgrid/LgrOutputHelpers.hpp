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
