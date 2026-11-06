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

#if defined(HAVE_METIS) && HAVE_MPI
// METIS installed as part of Scotch needs some of its definitions.
#if IS_SCOTCH_METIS_HEADER
extern "C" {
#include <scotch.h>
}
#ifndef SCOTCH_METIS_VERSION
#define SCOTCH_METIS_VERSION 5
#endif
#endif
extern "C" {
#include <metis.h>
}
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <numeric>
#include <stdexcept>
#include <type_traits>
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


namespace
{

constexpr int packHeaderSize = 11;
// header: dimensions, number of cells, faces, nodes, face nodes and cell faces,
// whether there are face tags, the Cartesian dimensions, whether there is global_cell

} // anonymous namespace

PackedGrid packGrid(const UnstructuredGrid& grid,
                    const std::vector<char>& isOwner,
                    const std::vector<int>& globalCells)
{
    const int dim = grid.dimensions;
    const int nC = grid.number_of_cells;
    const int nF = grid.number_of_faces;
    const int nN = grid.number_of_nodes;
    const int nFN = static_cast<int>(grid.face_nodepos[nF]);
    const int nCF = grid.cell_facepos[nC];
    const bool hasTags = (grid.cell_facetag != nullptr);

    if (static_cast<int>(isOwner.size()) != nC || static_cast<int>(globalCells.size()) != nC) {
        OPM_THROW(std::invalid_argument, "packGrid(): one ownership flag and global index per cell is required.");
    }

    PackedGrid packed;
    auto& ints = packed.ints;
    ints.reserve(packHeaderSize + nFN + nF + 1 + 2 * nF + nCF + nC + 1 + (hasTags ? nCF : 0) + 3 * nC);
    ints.insert(ints.end(), {dim, nC, nF, nN, nFN, nCF, hasTags ? 1 : 0,
                             grid.cartdims[0], grid.cartdims[1], grid.cartdims[2],
                             grid.global_cell ? 1 : 0});
    ints.insert(ints.end(), grid.face_nodes, grid.face_nodes + nFN);
    ints.insert(ints.end(), grid.face_nodepos, grid.face_nodepos + nF + 1);
    ints.insert(ints.end(), grid.face_cells, grid.face_cells + 2 * nF);
    ints.insert(ints.end(), grid.cell_faces, grid.cell_faces + nCF);
    ints.insert(ints.end(), grid.cell_facepos, grid.cell_facepos + nC + 1);
    if (hasTags) {
        ints.insert(ints.end(), grid.cell_facetag, grid.cell_facetag + nCF);
    }
    if (grid.global_cell) {
        ints.insert(ints.end(), grid.global_cell, grid.global_cell + nC);
    }
    ints.insert(ints.end(), isOwner.begin(), isOwner.end());
    ints.insert(ints.end(), globalCells.begin(), globalCells.end());

    auto& doubles = packed.doubles;
    doubles.reserve(dim * (nN + 2 * nF + nC) + nF + nC);
    doubles.insert(doubles.end(), grid.node_coordinates, grid.node_coordinates + dim * nN);
    doubles.insert(doubles.end(), grid.face_centroids, grid.face_centroids + dim * nF);
    doubles.insert(doubles.end(), grid.face_normals, grid.face_normals + dim * nF);
    doubles.insert(doubles.end(), grid.face_areas, grid.face_areas + nF);
    doubles.insert(doubles.end(), grid.cell_centroids, grid.cell_centroids + dim * nC);
    doubles.insert(doubles.end(), grid.cell_volumes, grid.cell_volumes + nC);
    return packed;
}

OwnedGrid unpackGrid(const std::vector<int>& ints, const std::vector<double>& doubles)
{
    if (static_cast<int>(ints.size()) < packHeaderSize) {
        OPM_THROW(std::invalid_argument, "unpackGrid(): truncated data.");
    }
    const int dim = ints[0];
    const int nC = ints[1];
    const int nF = ints[2];
    const int nN = ints[3];
    const int nFN = ints[4];
    const int nCF = ints[5];
    const bool hasTags = ints[6] != 0;
    const bool hasGlobalCell = ints[10] != 0;

    const std::size_t expectedInts = packHeaderSize + nFN + (nF + 1) + 2 * nF + nCF + (nC + 1)
        + (hasTags ? nCF : 0) + (hasGlobalCell ? nC : 0) + 2 * nC;
    const std::size_t expectedDoubles = dim * (nN + 2 * nF + nC) + nF + nC;
    if (ints.size() != expectedInts || doubles.size() != expectedDoubles) {
        OPM_THROW(std::invalid_argument, "unpackGrid(): inconsistent data sizes.");
    }

    OwnedGrid result;
    result.grid.reset(allocate_grid(dim, nC, nF, nFN, nCF, nN));
    if (!result.grid) {
        OPM_THROW(std::runtime_error, "unpackGrid(): grid allocation failed.");
    }
    UnstructuredGrid& g = *result.grid;
    std::copy_n(ints.begin() + 7, 3, g.cartdims);

    auto in = ints.begin() + packHeaderSize;
    std::copy_n(in, nFN, g.face_nodes);          in += nFN;
    std::copy_n(in, nF + 1, g.face_nodepos);     in += nF + 1;
    std::copy_n(in, 2 * nF, g.face_cells);       in += 2 * nF;
    std::copy_n(in, nCF, g.cell_faces);          in += nCF;
    std::copy_n(in, nC + 1, g.cell_facepos);     in += nC + 1;
    if (hasTags) {
        std::copy_n(in, nCF, g.cell_facetag);    in += nCF;
    }
    else {
        std::free(g.cell_facetag);
        g.cell_facetag = nullptr;
    }
    if (hasGlobalCell) {
        g.global_cell = static_cast<int*>(std::malloc(nC * sizeof(int)));
        if (g.global_cell == nullptr) {
            OPM_THROW(std::runtime_error, "unpackGrid(): grid allocation failed.");
        }
        std::copy_n(in, nC, g.global_cell);      in += nC;
    }
    result.isOwner.assign(in, in + nC);                in += nC;
    result.globalCells.assign(in, in + nC);

    auto d = doubles.begin();
    std::copy_n(d, dim * nN, g.node_coordinates); d += dim * nN;
    std::copy_n(d, dim * nF, g.face_centroids);   d += dim * nF;
    std::copy_n(d, dim * nF, g.face_normals);     d += dim * nF;
    std::copy_n(d, nF, g.face_areas);             d += nF;
    std::copy_n(d, dim * nC, g.cell_centroids);   d += dim * nC;
    std::copy_n(d, nC, g.cell_volumes);
    return result;
}

OwnedGrid distributeGrid(const UnstructuredGrid* grid,
                         const std::vector<int>& cellPart,
                         const int overlapLayers,
                         const bool ownersFirst,
                         const Dune::Communication<Dune::MPIHelper::MPICommunicator>& comm)
{
    const int size = comm.size();
    const bool isRoot = (comm.rank() == 0);

    // Validate on the root and tell everybody, so that all ranks throw
    // instead of some of them waiting for data that never arrives.
    int invalid = 0;
    if (isRoot) {
        invalid = (grid == nullptr) || (static_cast<int>(cellPart.size()) != grid->number_of_cells);
        for (const int part : cellPart) {
            invalid = invalid || part < 0 || part >= size;
        }
    }
    comm.broadcast(&invalid, 1, 0);
    if (invalid) {
        OPM_THROW(std::invalid_argument,
                  "distributeGrid(): rank 0 needs the grid and a partition with one value "
                  "in [0, number of ranks) for each cell.");
    }

    // One packed grid per rank, only filled on the root.
    std::vector<PackedGrid> packed(isRoot ? size : 0);
    std::vector<std::vector<int>> cellsOfRank(isRoot ? size : 0);
    if (isRoot) {
        const auto overlap = computeOverlapCells(*grid, cellPart, size, overlapLayers);
        std::vector<std::vector<int>> owned(size);
        for (int c = 0; c < grid->number_of_cells; ++c) {
            owned[cellPart[c]].push_back(c);
        }
        for (int rank = 0; rank < size; ++rank) {
            std::vector<int> cells = owned[rank];
            cells.insert(cells.end(), overlap[rank].begin(), overlap[rank].end());
            if (!ownersFirst) {
                std::sort(cells.begin(), cells.end());
            }
            std::vector<char> isOwner(cells.size());
            for (std::size_t i = 0; i < cells.size(); ++i) {
                isOwner[i] = (cellPart[cells[i]] == rank);
            }
            packed[rank] = packGrid(*extractLocalGrid(*grid, cells).grid, isOwner, cells);
            cellsOfRank[rank] = std::move(cells);
        }
    }

    if (size == 1) {
        auto result = unpackGrid(packed[0].ints, packed[0].doubles);
        result.cellsOfRank = std::move(cellsOfRank);
        return result;
    }

    // Tell each rank how much data to expect and send it.
    std::vector<int> counts(isRoot ? 2 * size : 0);
    for (std::size_t rank = 0; rank < packed.size(); ++rank) {
        counts[2 * rank] = static_cast<int>(packed[rank].ints.size());
        counts[2 * rank + 1] = static_cast<int>(packed[rank].doubles.size());
    }
    std::array<int, 2> myCounts{};
    comm.scatter(counts.data(), myCounts.data(), 2, 0);

    auto scatterData = [&](auto member, const int myCount) {
        using Value = typename std::remove_reference_t<decltype((packed[0].*member)[0])>;
        std::vector<Value> send;
        std::vector<int> sendLen(isRoot ? size : 0), displ(isRoot ? size : 0);
        for (std::size_t rank = 0; rank < packed.size(); ++rank) {
            displ[rank] = static_cast<int>(send.size());
            sendLen[rank] = static_cast<int>((packed[rank].*member).size());
            send.insert(send.end(), (packed[rank].*member).begin(), (packed[rank].*member).end());
        }
        std::vector<Value> recv(myCount);
        comm.scatterv(send.data(), sendLen.data(), displ.data(), recv.data(), myCount, 0);
        return recv;
    };
    const auto ints = scatterData(&PackedGrid::ints, myCounts[0]);
    const auto doubles = scatterData(&PackedGrid::doubles, myCounts[1]);
    auto result = unpackGrid(ints, doubles);
    result.cellsOfRank = std::move(cellsOfRank);
    return result;
}

bool metisAvailable()
{
#if defined(HAVE_METIS) && HAVE_MPI
    return true;
#else
    return false;
#endif
}

std::vector<int> partitionWithMetis([[maybe_unused]] const UnstructuredGrid& grid,
                                    [[maybe_unused]] const int numParts,
                                    [[maybe_unused]] const double imbalanceTol)
{
#if defined(HAVE_METIS) && HAVE_MPI
    const int numCells = grid.number_of_cells;
    if (numParts < 1 || numParts > numCells) {
        OPM_THROW(std::invalid_argument,
                  "partitionWithMetis(): the number of parts must be between one and the number of cells.");
    }
    if (numParts == 1) {
        return std::vector<int>(numCells, 0);
    }

    // The graph of the cells. Cells that share more than one face get one edge, weighted by the number of faces.
#if IS_SCOTCH_METIS_HEADER
    using Index = SCOTCH_Num;
#elif defined(IDXTYPEWIDTH)
    using Index = ::idx_t;
#else
    using Index = int;
#endif
    std::vector<Index> xadj(numCells + 1, 0);
    std::vector<Index> adjncy;
    std::vector<Index> adjwgt;
    std::vector<int> neighbors;
    for (int cell = 0; cell < numCells; ++cell) {
        neighbors.clear();
        for (auto pos = grid.cell_facepos[cell]; pos < grid.cell_facepos[cell + 1]; ++pos) {
            const int other = neighbor(grid, grid.cell_faces[pos], cell);
            if (other >= 0 && other != cell) {
                neighbors.push_back(other);
            }
        }
        std::sort(neighbors.begin(), neighbors.end());
        for (std::size_t i = 0; i < neighbors.size();) {
            std::size_t j = i;
            while (j < neighbors.size() && neighbors[j] == neighbors[i]) {
                ++j;
            }
            adjncy.push_back(neighbors[i]);
            adjwgt.push_back(static_cast<Index>(j - i));
            i = j;
        }
        xadj[cell + 1] = static_cast<Index>(adjncy.size());
    }

    Index options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    // The imbalance is given in parts per thousand above one.
    options[METIS_OPTION_UFACTOR] = std::max<Index>(1, static_cast<Index>(std::lround((imbalanceTol - 1.0) * 1000.0)));

    Index numVertices = numCells;
    Index numConstraints = 1;
    Index numPartsMetis = numParts;
    Index objective = 0;
    std::vector<Index> part(numCells, 0);
    const int status = METIS_PartGraphKway(&numVertices, &numConstraints, xadj.data(), adjncy.data(),
                                           /* vwgt = */ nullptr, /* vsize = */ nullptr, adjwgt.data(),
                                           &numPartsMetis, /* tpwgts = */ nullptr, /* ubvec = */ nullptr,
                                           options, &objective, part.data());
    if (status != METIS_OK) {
        OPM_THROW(std::runtime_error, "partitionWithMetis(): METIS failed to partition the grid.");
    }
    return std::vector<int>(part.begin(), part.end());
#else
    OPM_THROW(std::runtime_error, "partitionWithMetis(): the grid library was built without METIS.");
#endif
}

} // namespace Dune::unstructured
