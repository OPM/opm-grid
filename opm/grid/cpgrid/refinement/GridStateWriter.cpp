/*
  Copyright 2026 SINTEF Digital.

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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <opm/grid/cpgrid/refinement/GridStateWriter.hpp>

#include <opm/grid/cpgrid/Indexsets.hpp>
#include <opm/grid/cpgrid/PartitionTypeIndicator.hpp>

#include <memory>
#include <utility>

namespace Opm
{
namespace Refinement
{

void GridStateWriter::setLevel(Data& grid, int level)
{
    grid.levels_.level = level;
}

void GridStateWriter::setCellsPerDim(Data& grid, const std::array<int,3>& cellsPerDim)
{
    grid.levels_.cells_per_dim = cellsPerDim;
}

void GridStateWriter::setSubdivision(Data& grid,
                                     const std::array<Opm::Refinement::AxisSubdivision,3>& subdivision)
{
    grid.levels_.subdivision = subdivision;
}

void GridStateWriter::setParentRelations(Data& grid,
                                         std::vector<std::array<int,2>> childToParent,
                                         std::vector<int> idxInParent)
{
    grid.levels_.child_to_parent_cells = std::move(childToParent);
    grid.levels_.cell_to_idxInParentCell = std::move(idxInParent);
}

void GridStateWriter::setParentToChildren(Data& grid,
                                          std::vector<std::tuple<int, std::vector<int>>> parentToChildren)
{
    grid.levels_.parent_to_children_cells = std::move(parentToChildren);
}

void GridStateWriter::setLeafToLevel(Data& grid, std::vector<std::array<int,2>> leafToLevel)
{
    grid.levels_.leaf_to_level_cells = std::move(leafToLevel);
}

void GridStateWriter::setAquiferCells(Data& grid, std::vector<int> cells)
{
    grid.aquifer_cells_ = std::move(cells);
}

void GridStateWriter::setCornerHistory(Data& grid, std::vector<std::array<int,2>> cornerHistory)
{
    grid.levels_.corner_history = std::move(cornerHistory);
}

void GridStateWriter::setLogicalCartesianSize(Data& grid, const std::array<int,3>& size)
{
    grid.logical_cartesian_size_ = size;
}

void GridStateWriter::setGlobalCell(Data& grid, std::vector<int> globalCell)
{
    grid.global_cell_ = std::move(globalCell);
}

void GridStateWriter::setGlobalIdMapping(Data& grid,
                                         std::vector<int> cellIds,
                                         std::vector<int> faceIds,
                                         std::vector<int> pointIds)
{
    grid.global_id_set_->swap(cellIds, faceIds, pointIds);
}

void GridStateWriter::setIndexSet(Data& grid, std::size_t numCells, std::size_t numPoints)
{
    grid.index_set_ = std::make_unique<Dune::cpgrid::IndexSet>(numCells, numPoints);
}

void GridStateWriter::setRefinementMaxLevel(Data& grid, int maxLevel)
{
    grid.levels_.refinement_max_level = maxLevel;
}

void GridStateWriter::setCommunicator(Data& grid, const Data::Communication& comm)
{
    grid.ccobj_ = comm;
}

std::vector<std::array<int,8>>& GridStateWriter::cellToPoint(Data& grid)
{
    return grid.cell_to_point_;
}

Dune::cpgrid::OrientedEntityTable<0,1>& GridStateWriter::cellToFace(Data& grid)
{
    return grid.cell_to_face_;
}

Dune::cpgrid::OrientedEntityTable<1,0>& GridStateWriter::faceToCell(Data& grid)
{
    return grid.face_to_cell_;
}

Opm::SparseTable<int>& GridStateWriter::faceToPoint(Data& grid)
{
    return grid.face_to_point_;
}

Dune::cpgrid::EntityVariable<enum face_tag, 1>& GridStateWriter::faceTag(Data& grid)
{
    return grid.face_tag_;
}

Dune::cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>, 1>& GridStateWriter::faceNormals(Data& grid)
{
    return grid.face_normals_;
}

Dune::cpgrid::DefaultGeometryPolicy& GridStateWriter::geometry(Data& grid)
{
    return grid.geometry_;
}

const std::vector<std::array<int,2>>& GridStateWriter::childToParent(const Data& grid)
{
    return grid.levels_.child_to_parent_cells;
}

const std::vector<int>& GridStateWriter::idxInParent(const Data& grid)
{
    return grid.levels_.cell_to_idxInParentCell;
}

void GridStateWriter::setRetainedCornerPointInput(Data& grid,
                                                  std::shared_ptr<const RetainedCornerPointInput> input)
{
    grid.retained_cp_input_ = std::move(input);
}

std::shared_ptr<const RetainedCornerPointInput> GridStateWriter::retainedCornerPointInput(const Data& grid)
{
    return grid.retained_cp_input_;
}

const std::vector<char>& GridStateWriter::cellPartitionTypes(const Data& grid)
{
    return grid.partition_type_indicator_->cell_indicator_;
}

void GridStateWriter::setCellPartitionTypes(Data& grid, std::vector<char> types)
{
    grid.partition_type_indicator_->cell_indicator_ = std::move(types);
}

} // namespace Refinement
} // namespace Opm
