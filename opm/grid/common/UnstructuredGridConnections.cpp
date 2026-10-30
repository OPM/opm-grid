/*
  Copyright 2026 Equinor ASA

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

#include <opm/grid/common/UnstructuredGridConnections.hpp>

#include <cstdlib>
#include <new>
#include <stdexcept>

namespace Dune::unstructured
{

namespace
{
// The tags of the sides of the cells of a corner-point grid: the sides in the three directions, the lower side first.
constexpr int topTag = 4;
constexpr int bottomTag = 5;
} // anonymous namespace

int connectCellsAcrossGaps(UnstructuredGrid& grid, const std::vector<std::pair<int, int>>& cells)
{
    if (grid.dimensions != 3 || grid.cell_facetag == nullptr) {
        throw std::invalid_argument("connectCellsAcrossGaps needs a three dimensional grid with the tags of the sides");
    }

    // The faces that are added to the faces of the cell below, with the face.
    std::vector<std::pair<int, int>> added;
    for (const auto& [above, below] : cells) {
        if (above < 0 || below < 0 || above >= grid.number_of_cells || below >= grid.number_of_cells || above == below) {
            continue;
        }
        for (auto pos = grid.cell_facepos[above]; pos < grid.cell_facepos[above + 1]; ++pos) {
            const int face = grid.cell_faces[pos];
            if (grid.cell_facetag[pos] == bottomTag && grid.face_cells[2 * face] == above && grid.face_cells[2 * face + 1] == -1) {
                grid.face_cells[2 * face + 1] = below;
                added.emplace_back(below, face);
                break;
            }
        }
    }
    if (added.empty()) {
        return 0;
    }

    // The faces of each cell are in cell_faces after each other, so the cells below get room for one more.
    const std::size_t oldSize = grid.cell_facepos[grid.number_of_cells];
    const std::size_t newSize = oldSize + added.size();
    auto* newPos = static_cast<grid_size_t*>(std::malloc((grid.number_of_cells + 1) * sizeof(grid_size_t)));
    auto* newFaces = static_cast<int*>(std::malloc(newSize * sizeof(int)));
    auto* newTags = static_cast<int*>(std::malloc(newSize * sizeof(int)));
    if (newPos == nullptr || newFaces == nullptr || newTags == nullptr) {
        std::free(newPos);
        std::free(newFaces);
        std::free(newTags);
        throw std::bad_alloc();
    }

    std::vector<std::vector<int>> facesAdded(grid.number_of_cells);
    for (const auto& [cell, face] : added) {
        facesAdded[cell].push_back(face);
    }
    grid_size_t out = 0;
    for (int c = 0; c < grid.number_of_cells; ++c) {
        newPos[c] = out;
        for (auto pos = grid.cell_facepos[c]; pos < grid.cell_facepos[c + 1]; ++pos) {
            newFaces[out] = grid.cell_faces[pos];
            newTags[out] = grid.cell_facetag[pos];
            ++out;
        }
        for (const int face : facesAdded[c]) {
            newFaces[out] = face;
            newTags[out] = topTag;
            ++out;
        }
    }
    newPos[grid.number_of_cells] = out;

    std::free(grid.cell_facepos);
    std::free(grid.cell_faces);
    std::free(grid.cell_facetag);
    grid.cell_facepos = newPos;
    grid.cell_faces = newFaces;
    grid.cell_facetag = newTags;
    return static_cast<int>(added.size());
}

} // namespace Dune::unstructured
