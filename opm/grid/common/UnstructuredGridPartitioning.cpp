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

#include <opm/grid/common/UnstructuredGridPartitioning.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace Dune::unstructured
{

namespace
{

// The cell on the other side of face f as seen from cell c, or -1.
int neighbor(const UnstructuredGrid& grid, const int face, const int cell)
{
    const int c0 = grid.face_cells[2 * face];
    const int c1 = grid.face_cells[2 * face + 1];
    return (c0 == cell) ? c1 : c0;
}

} // anonymous namespace

std::vector<std::vector<int>>
computeOverlapCells(const UnstructuredGrid& grid,
                    const std::vector<int>& cellPart,
                    const int numParts,
                    const int layers)
{
    const int numCells = grid.number_of_cells;
    if (static_cast<int>(cellPart.size()) != numCells) {
        OPM_THROW(std::invalid_argument,
                  "computeOverlapCells(): one partition number per cell is required.");
    }
    for (const int part : cellPart) {
        if (part < 0 || part >= numParts) {
            OPM_THROW(std::invalid_argument,
                      "computeOverlapCells(): partition number out of range.");
        }
    }

    std::vector<std::vector<int>> overlap(numParts);
    if (layers <= 0) {
        return overlap;
    }

    // Cells of each partition that were added in the previous layer.
    std::vector<std::vector<int>> frontier(numParts);
    std::vector<std::unordered_set<int>> included(numParts);

    // First layer: foreign cells across partition boundary faces.
    for (int face = 0; face < grid.number_of_faces; ++face) {
        const int c0 = grid.face_cells[2 * face];
        const int c1 = grid.face_cells[2 * face + 1];
        if (c0 < 0 || c1 < 0 || cellPart[c0] == cellPart[c1]) {
            continue;
        }
        if (included[cellPart[c1]].insert(c0).second) {
            frontier[cellPart[c1]].push_back(c0);
        }
        if (included[cellPart[c0]].insert(c1).second) {
            frontier[cellPart[c0]].push_back(c1);
        }
    }

    for (int layer = 1; layer < layers; ++layer) {
        for (int part = 0; part < numParts; ++part) {
            std::vector<int> next;
            for (const int cell : frontier[part]) {
                for (auto pos = grid.cell_facepos[cell]; pos < grid.cell_facepos[cell + 1]; ++pos) {
                    const int other = neighbor(grid, grid.cell_faces[pos], cell);
                    if (other < 0 || cellPart[other] == part) {
                        continue;
                    }
                    if (included[part].insert(other).second) {
                        next.push_back(other);
                    }
                }
            }
            frontier[part] = std::move(next);
        }
    }

    for (int part = 0; part < numParts; ++part) {
        overlap[part].assign(included[part].begin(), included[part].end());
        std::sort(overlap[part].begin(), overlap[part].end());
    }
    return overlap;
}

LocalGrid extractLocalGrid(const UnstructuredGrid& grid,
                           const std::vector<int>& cells)
{
    const int dim = grid.dimensions;

    std::unordered_map<int, int> globalToLocalCell;
    globalToLocalCell.reserve(cells.size());
    for (std::size_t i = 0; i < cells.size(); ++i) {
        const int cell = cells[i];
        if (cell < 0 || cell >= grid.number_of_cells) {
            OPM_THROW(std::invalid_argument,
                      "extractLocalGrid(): cell index out of range.");
        }
        if (!globalToLocalCell.emplace(cell, static_cast<int>(i)).second) {
            OPM_THROW(std::invalid_argument,
                      "extractLocalGrid(): duplicate cell index.");
        }
    }

    LocalGrid result;

    // Faces: all faces of the cells, in increasing global order.
    std::size_t numCellFaces = 0;
    for (const int cell : cells) {
        numCellFaces += grid.cell_facepos[cell + 1] - grid.cell_facepos[cell];
        for (auto pos = grid.cell_facepos[cell]; pos < grid.cell_facepos[cell + 1]; ++pos) {
            result.globalFace.push_back(grid.cell_faces[pos]);
        }
    }
    std::sort(result.globalFace.begin(), result.globalFace.end());
    result.globalFace.erase(std::unique(result.globalFace.begin(), result.globalFace.end()),
                            result.globalFace.end());

    // Nodes: all nodes of those faces, in increasing global order.
    std::size_t numFaceNodes = 0;
    for (const int face : result.globalFace) {
        numFaceNodes += grid.face_nodepos[face + 1] - grid.face_nodepos[face];
        for (auto pos = grid.face_nodepos[face]; pos < grid.face_nodepos[face + 1]; ++pos) {
            result.globalNode.push_back(grid.face_nodes[pos]);
        }
    }
    std::sort(result.globalNode.begin(), result.globalNode.end());
    result.globalNode.erase(std::unique(result.globalNode.begin(), result.globalNode.end()),
                            result.globalNode.end());

    std::unordered_map<int, int> globalToLocalFace;
    globalToLocalFace.reserve(result.globalFace.size());
    for (std::size_t i = 0; i < result.globalFace.size(); ++i) {
        globalToLocalFace.emplace(result.globalFace[i], static_cast<int>(i));
    }
    std::unordered_map<int, int> globalToLocalNode;
    globalToLocalNode.reserve(result.globalNode.size());
    for (std::size_t i = 0; i < result.globalNode.size(); ++i) {
        globalToLocalNode.emplace(result.globalNode[i], static_cast<int>(i));
    }

    const auto numFaces = result.globalFace.size();
    const auto numNodes = result.globalNode.size();

    result.grid.reset(allocate_grid(dim, cells.size(), numFaces, numFaceNodes,
                                    numCellFaces, numNodes));
    if (!result.grid) {
        OPM_THROW(std::runtime_error, "extractLocalGrid(): grid allocation failed.");
    }
    UnstructuredGrid& local = *result.grid;

    // Nodes.
    for (std::size_t n = 0; n < numNodes; ++n) {
        std::copy_n(grid.node_coordinates + dim * result.globalNode[n], dim,
                    local.node_coordinates + dim * n);
    }

    // Faces.
    local.face_nodepos[0] = 0;
    std::size_t nodePos = 0;
    for (std::size_t f = 0; f < numFaces; ++f) {
        const int gf = result.globalFace[f];
        for (auto pos = grid.face_nodepos[gf]; pos < grid.face_nodepos[gf + 1]; ++pos) {
            local.face_nodes[nodePos++] = globalToLocalNode.at(grid.face_nodes[pos]);
        }
        local.face_nodepos[f + 1] = nodePos;

        for (int side = 0; side < 2; ++side) {
            const int gc = grid.face_cells[2 * gf + side];
            const auto it = (gc < 0) ? globalToLocalCell.end() : globalToLocalCell.find(gc);
            local.face_cells[2 * f + side] = (it == globalToLocalCell.end()) ? -1 : it->second;
        }

        std::copy_n(grid.face_centroids + dim * gf, dim, local.face_centroids + dim * f);
        std::copy_n(grid.face_normals + dim * gf, dim, local.face_normals + dim * f);
        local.face_areas[f] = grid.face_areas[gf];
    }

    // Cells.
    const bool hasTags = (grid.cell_facetag != nullptr);
    if (!hasTags) {
        std::free(local.cell_facetag);
        local.cell_facetag = nullptr;
    }
    local.global_cell = static_cast<int*>(std::malloc(cells.size() * sizeof(int)));
    if (local.global_cell == nullptr) {
        OPM_THROW(std::runtime_error, "extractLocalGrid(): grid allocation failed.");
    }

    local.cell_facepos[0] = 0;
    std::size_t facePos = 0;
    for (std::size_t c = 0; c < cells.size(); ++c) {
        const int gc = cells[c];
        for (auto pos = grid.cell_facepos[gc]; pos < grid.cell_facepos[gc + 1]; ++pos) {
            local.cell_faces[facePos] = globalToLocalFace.at(grid.cell_faces[pos]);
            if (hasTags) {
                local.cell_facetag[facePos] = grid.cell_facetag[pos];
            }
            ++facePos;
        }
        local.cell_facepos[c + 1] = facePos;

        std::copy_n(grid.cell_centroids + dim * gc, dim, local.cell_centroids + dim * c);
        local.cell_volumes[c] = grid.cell_volumes[gc];
        local.global_cell[c] = grid.global_cell ? grid.global_cell[gc] : gc;
    }

    std::copy_n(grid.cartdims, 3, local.cartdims);

    return result;
}

} // namespace Dune::unstructured
