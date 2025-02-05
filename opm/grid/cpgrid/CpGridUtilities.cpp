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

#include <opm/grid/cpgrid/CpGridUtilities.hpp>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>

#include <array>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace Opm
{

std::tuple<std::vector<int>, std::unordered_map<int, int>, std::vector<std::array<int, 3>>>
lgrIJK(const Dune::CpGrid& grid, const std::string& lgr_name)
{
    // Check if lgr_name exists in lgr_names_
    const auto& lgr_names = grid.getLgrNameToLevel();
    auto it = lgr_names.find(lgr_name);
    if (it == lgr_names.end()) {
        OPM_THROW(std::runtime_error, "LGR name not found: " + lgr_name);
    }

    const auto level = it->second;
    const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapper(grid);
    const auto& levelView = grid.levelGridView(level);
    const auto numCells = levelView.size(0);

    std::vector<std::array<int, 3>> lgrIJK(numCells);
    std::vector<int> cellIdxToLgrCartesianIdx(numCells);
    std::unordered_map<int, int> lgrCartesianIdxToCellIdx;
    lgrCartesianIdxToCellIdx.reserve(numCells);

    // Iterate over (active) elements in the grid and populate the structures
    for (const auto& element : Dune::elements(grid.levelGridView(level))) {
        std::array<int, 3> ijk;
        levelCartMapper.cartesianCoordinate(element.index(), ijk, level);

        const int cellIndex = element.index();
        const int cartesianIdx = element.getLevelCartesianIdx();

        lgrIJK[cellIndex] = ijk;
        cellIdxToLgrCartesianIdx[cellIndex] = cartesianIdx;
        lgrCartesianIdxToCellIdx[cartesianIdx] = cellIndex;
    }

    return std::make_tuple(cellIdxToLgrCartesianIdx, lgrCartesianIdxToCellIdx, lgrIJK);
}

} // namespace Opm
