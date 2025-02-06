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

#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Opm
{

std::pair<std::unordered_map<int, int>, std::vector<std::array<int, 3>>>
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
    std::unordered_map<int, int> lgrCartesianIdxToCellIdx;
    lgrCartesianIdxToCellIdx.reserve(numCells);

    // Iterate over (active) elements in the grid and populate the structures
    for (const auto& element : Dune::elements(grid.levelGridView(level))) {
        std::array<int, 3> ijk;
        levelCartMapper.cartesianCoordinate(element.index(), ijk, level);

        const int cellIndex = element.index();
        const int cartesianIdx = element.getLevelCartesianIdx();

        lgrIJK[cellIndex] = ijk;
        lgrCartesianIdxToCellIdx[cartesianIdx] = cellIndex;
    }

    return std::make_pair(lgrCartesianIdxToCellIdx, lgrIJK);
}

std::vector<std::array<double, 6>> lgrCOORD(const Dune::CpGrid& grid,
                                            int level,
                                            const std::unordered_map<int, int>&  lgrCartesianIdxToCellIdx,
                                            const std::vector<std::array<int, 3>>& lgrIJK)
{
    const auto& levelGrid = *(grid.currentData()[level]);
    
    // Check not all cells are inactive
    const auto numCells = levelGrid.size(0);
    if (numCells == 0) {
        OPM_THROW(std::logic_error, "LGR in level " + std::to_string(level) + " has no active cells.\n");
    }

    // LGR dimensions
    const auto& lgr_dim = grid.currentData()[level]->logicalCartesianSize();
    const int nx = lgr_dim[0];
    const int ny = lgr_dim[1];
    const int nz = lgr_dim[2];

    // Initialize all pillars as inactive (setting COORD values to std::numeric_limits<double>::max()).
    std::vector<std::array<double,6>> lgrCOORD((nx+1)*(ny+1));
    for (auto& pillar : lgrCOORD) {
        pillar.fill(std::numeric_limits<double>::max());
    }

    // Map to determine min and max k per cell column (i, j) (min/max_k = 0, ..., nz-1).
    // Initialized as {nz, -1} to detect inactive cell columns.
    std::vector<std::array<int,2>> minMaxPerCellPillar(nx*ny, {nz, -1});

    // Compute the bottom and top k per cell pillar (i, j).
    for (const auto& ijk : lgrIJK) {
        int cell_pillar_idx = ijk[1] * nx + ijk[0];
        auto& minMax = minMaxPerCellPillar[cell_pillar_idx];

        minMax[0] = std::min(ijk[2], minMax[0]);
        minMax[1] = std::max(ijk[2], minMax[1]);
    }

    // Rewrite values for active pillars
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const int cell_pillar_idx = (j*nx) + i;

            // Get min/max k for pillar at (i,j)
            const auto& [bottom_k, top_k] = minMaxPerCellPillar[cell_pillar_idx];

            if ( (bottom_k == nz) || (top_k == -1)) {
                 continue; // no active pillar at (i,j)
            }

            const auto bottom_lgr_cartesian_idx = (bottom_k*nx*ny) + cell_pillar_idx;
            const auto top_lgr_cartesian_idx = (top_k*nx*ny) + cell_pillar_idx;

            const auto& bottomElemIdx = lgrCartesianIdxToCellIdx.at(bottom_lgr_cartesian_idx);
            const auto& topElemIdx = lgrCartesianIdxToCellIdx.at(top_lgr_cartesian_idx);

            const auto& bottomElem = Dune::cpgrid::Entity<0>(levelGrid, bottomElemIdx, true);
            const auto& topElem = Dune::cpgrid::Entity<0>(levelGrid, topElemIdx, true);

            Opm::processPillars(i,j, nx, topElem, bottomElem, lgrCOORD);
        }
    }
    return lgrCOORD;
}

void setPillarCoordinates(int i, int j, int nx,
                          int topCorner, int bottomCorner, int positionIdx,
                          const Dune::cpgrid::Entity<0>& topElem,
                          const Dune::cpgrid::Entity<0>& bottomElem,
                          std::vector<std::array<double, 6>>& lgrCOORD)
{
    // positionIdx (0,1,2, or 3) is used to distinguish the 4 pillars:
    //
    // pillar       positionIdx   positionIdx / 2   positionIdx % 2
    // (i,j)             0               0                  0
    // (i+1,j)           1               0                  1
    // (i,j+1)           2               1                  0
    // (i+1,j+1)         3               1                  1
    const int pillar = ((j + positionIdx / 2) * (nx + 1)) + (i + positionIdx % 2);

    // Top pillar's COORD values
    const auto& top_point = topElem.subEntity<3>(topCorner).geometry().center();
    std::copy(top_point.begin(), top_point.end(), lgrCOORD[pillar].begin());

    // Bottom pillar's COORD values
    const auto& bottom_point = bottomElem.subEntity<3>(bottomCorner).geometry().center();
    std::copy(bottom_point.begin(), bottom_point.end(), lgrCOORD[pillar].begin() + 3);
}

void processPillars(int i, int j, int nx,
                    const Dune::cpgrid::Entity<0>& topElem,
                    const Dune::cpgrid::Entity<0>& bottomElem,
                    std::vector<std::array<double, 6>>& lgrCOORD)
{
    // Recall that a cell has 8 corners:
    //        6 --- 7
    //       /     /   TOP FACE
    //      4 --- 5
    //        2 --- 3
    //       /     /   BOTTOM FACE
    //      0 --- 1

    // To take into account inactive cells, consider for each (i,j) column of cells, 4 pillars:
    // (i,j)     pillar associated with bottom element corner 0 and top element corner 4
    // (i+1,j)   pillar associated with bottom element corner 1 and top element corner 5
    // (i,j+1)   pillar associated with bottom element corner 2 and top element corner 6
    // (i+1,j+1) pillar associated with bottom element corner 3 and top element corner 7
    setPillarCoordinates(i, j, nx, 4 /*topCorner*/, 0 /*bottomCorner*/,  0 /*to select pillar (i,j)*/, topElem, bottomElem, lgrCOORD);
    setPillarCoordinates(i, j, nx, 5 /*topCorner*/, 1 /*bottomCorner*/,  1 /*to select pillar (i+1,j)*/, topElem, bottomElem, lgrCOORD);
    setPillarCoordinates(i, j, nx, 6 /*topCorner*/, 2 /*bottomCorner*/,  2 /*to select pillar (i,j+1)*/, topElem, bottomElem, lgrCOORD);
    setPillarCoordinates(i, j, nx, 7 /*topCorner*/, 3 /*bottomCorner*/,  3 /*to select pillar (i+1,j+1)*/, topElem, bottomElem, lgrCOORD);
}

} // namespace Opm
