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
#ifndef OPM_UNSTRUCTUREDGRIDPARTITIONING_HEADER_INCLUDED
#define OPM_UNSTRUCTUREDGRIDPARTITIONING_HEADER_INCLUDED

#include <opm/grid/UnstructuredGrid.h>

#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <memory>
#include <vector>

/// \file
/// Grid-agnostic helpers for distributing an UnstructuredGrid.
///
/// These functions work directly on the plain-C UnstructuredGrid data and
/// involve no communication, so they can be tested in serial.  They are used
/// by PolyhedralGrid but deliberately do not depend on CpGrid.

namespace Dune::unstructured
{

/// Compute the overlap cells of every partition.
///
/// The first overlap layer of a partition consists of all cells that are not
/// owned by it but share a face with a cell it owns.  Layer k+1 consists of
/// the not yet included, not owned face neighbors of layer k.
///
/// \param grid       The global grid.
/// \param cellPart   For each cell the (zero-based) owning partition.
/// \param numParts   Number of partitions.
/// \param layers     Number of overlap layers (zero gives empty overlap).
/// \return For each partition the sorted global (compressed) indices of its
///         overlap cells.
std::vector<std::vector<int>>
computeOverlapCells(const UnstructuredGrid& grid,
                    const std::vector<int>& cellPart,
                    int numParts,
                    int layers);

/// Deleter for grids allocated by the UnstructuredGrid C API.
struct GridDeleter
{
    void operator()(UnstructuredGrid* grid) const
    {
        destroy_grid(grid);
    }
};

/// A grid extracted from a global grid, together with the mapping of its
/// entities to the entities of the global grid.
struct LocalGrid
{
    std::unique_ptr<UnstructuredGrid, GridDeleter> grid;
    /// Global face index of each local face.
    std::vector<int> globalFace;
    /// Global node index of each local node.
    std::vector<int> globalNode;
};

/// Extract the sub grid consisting of the given cells.
///
/// Local cell i is global cell cells[i], so the caller controls the local
/// cell order (for instance owner cells first, then overlap cells).  The
/// local faces are all faces of these cells in increasing global face order,
/// and the local nodes all their nodes in increasing global node order.
/// Faces whose neighbor is not part of the sub grid become boundary faces
/// (the missing neighbor is -1); the orientation of faces is retained.
///
/// global_cell of the local grid is inherited from the global grid, or is the
/// global compressed cell index if the global grid has no global_cell.
/// cartdims are copied.
///
/// \param grid   The global grid.
/// \param cells  Distinct global cell indices.
LocalGrid extractLocalGrid(const UnstructuredGrid& grid,
                           const std::vector<int>& cells);

/// A local grid in a flat form that can be sent between processes.
struct PackedGrid
{
    std::vector<int> ints;
    std::vector<double> doubles;
};

/// Flatten a grid.
///
/// \param grid        The grid to pack.
/// \param isOwner     For each cell whether the process owns it (1) or whether
///                    it is an overlap cell (0).
/// \param globalCells For each cell its index in the global grid.
PackedGrid packGrid(const UnstructuredGrid& grid,
                    const std::vector<char>& isOwner,
                    const std::vector<int>& globalCells);

/// A grid together with the ownership of its cells.
struct OwnedGrid
{
    std::unique_ptr<UnstructuredGrid, GridDeleter> grid;
    /// For each cell whether the process owns it (1) or whether it is an overlap cell (0).
    std::vector<char> isOwner;
    /// For each cell its index in the global grid.
    std::vector<int> globalCells;
    /// Only on rank 0 after distributeGrid(): for each rank the global index
    /// of each of its cells, in the order of the cells of that rank.
    std::vector<std::vector<int>> cellsOfRank;
};

/// Inverse of packGrid().
OwnedGrid unpackGrid(const std::vector<int>& ints, const std::vector<double>& doubles);

/// Distribute a grid from rank 0 to all ranks of a communicator.
///
/// Every rank receives the cells it owns according to the partitioning plus
/// the requested layers of overlap cells. This is a collective operation, and
/// all ranks throw if the input is invalid.
///
/// \param grid          The global grid, only used on rank 0 (may be null elsewhere).
/// \param cellPart      For each cell of the global grid the rank that owns
///                      it. Only used on rank 0. Values must be smaller than
///                      the size of the communicator.
/// \param overlapLayers Number of overlap layers.
/// \param ownersFirst   Order the owned cells before the overlap cells.
///                      Otherwise the cells are in increasing global order.
/// \param comm          The communicator.
OwnedGrid distributeGrid(const UnstructuredGrid* grid,
                         const std::vector<int>& cellPart,
                         int overlapLayers,
                         bool ownersFirst,
                         const Dune::Communication<Dune::MPIHelper::MPICommunicator>& comm);

} // namespace Dune::unstructured

#endif // OPM_UNSTRUCTUREDGRIDPARTITIONING_HEADER_INCLUDED
