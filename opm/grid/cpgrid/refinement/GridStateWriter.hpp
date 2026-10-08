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
#ifndef OPM_GRID_REFINEMENT_GRIDSTATEWRITER_HEADER_INCLUDED
#define OPM_GRID_REFINEMENT_GRIDSTATEWRITER_HEADER_INCLUDED

#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/refinement/RetainedCornerPointInput.hpp>

#include <memory>

#include <opm/grid/cpgrid/refinement/RefinementRequest.hpp>

#include <array>
#include <tuple>
#include <vector>

namespace Opm
{
namespace Refinement
{

/// The single access point for multilevel/refinement state inside
/// CpGridData. Befriended by CpGridData so the refinement builder can stay
/// outside the grid classes without widening the friend surface: the
/// builder is the only writer of this state (docs/DESIGN-builder.md §5),
/// and every read/mutation of otherwise-private members goes through these
/// named operations.
struct GridStateWriter
{
    using Data = Dune::cpgrid::CpGridData;

    /// Set the level index a grid occupies in the hierarchy.
    static void setLevel(Data& grid, int level);

    /// Set the per-parent subdivision factors of a refined level grid.
    static void setCellsPerDim(Data& grid, const std::array<int,3>& cellsPerDim);

    /// Set a refined level's graded subdivision: where each of its refined
    /// columns sits inside its parent cell, per direction. cells_per_dim_ only
    /// describes a uniform box, so this is what geometryInFather() needs once
    /// the columns differ in size.
    static void setSubdivision(Data& grid,
                               const std::array<Opm::Refinement::AxisSubdivision,3>& subdivision);

    /// Set the parent relations of a refined level or leaf grid:
    /// childToParent[c] = {parent level, parent cell index} ({-1,-1} for
    /// cells without a parent), and idxInParent[c] = lattice index
    /// i + j*rx + k*rx*ry within the parent (-1 without a parent).
    static void setParentRelations(Data& grid,
                                   std::vector<std::array<int,2>> childToParent,
                                   std::vector<int> idxInParent);

    /// Set parent-to-children on the grid the parents live in:
    /// entry = {child level, child indices} or {-1, {}}.
    static void setParentToChildren(Data& grid,
                                    std::vector<std::tuple<int, std::vector<int>>> parentToChildren);

    /// Set leaf-to-level mapping on the leaf: entry = {level, level index}.
    static void setLeafToLevel(Data& grid, std::vector<std::array<int,2>> leafToLevel);

    /// Set the sorted cell indices of numerical aquifer cells.
    static void setAquiferCells(Data& grid, std::vector<int> cells);

    /// Set corner history: entry = {birth level, corner index there} or {-1,-1}.
    static void setCornerHistory(Data& grid, std::vector<std::array<int,2>> cornerHistory);

    static void setLogicalCartesianSize(Data& grid, const std::array<int,3>& size);
    static void setGlobalCell(Data& grid, std::vector<int> globalCell);
    static void setIndexSet(Data& grid, std::size_t numCells, std::size_t numPoints);
    static void setRefinementMaxLevel(Data& grid, int maxLevel);

    /// Set the collective communicator (ccobj_) a grid reports. Used by the
    /// refine-before-redistribute path: the new level grids and the leaf are
    /// assembled with a self-communicator (collective-free serial refinement),
    /// then restored to the world communicator here so the refined grid still
    /// reports the full parallel communicator and can be load balanced.
    static void setCommunicator(Data& grid, const Data::Communication& comm);

    /// Topology/geometry table access for assembling a grid (and reading
    /// the tables of source grids during assembly).
    static std::vector<std::array<int,8>>& cellToPoint(Data& grid);
    static Dune::cpgrid::OrientedEntityTable<0,1>& cellToFace(Data& grid);
    static Dune::cpgrid::OrientedEntityTable<1,0>& faceToCell(Data& grid);
    static Opm::SparseTable<int>& faceToPoint(Data& grid);
    static Dune::cpgrid::EntityVariable<enum face_tag, 1>& faceTag(Data& grid);
    static Dune::cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>, 1>& faceNormals(Data& grid);
    static Dune::cpgrid::DefaultGeometryPolicy& geometry(Data& grid);
    static const std::vector<std::array<int,2>>& childToParent(const Data& grid);
    static const std::vector<int>& idxInParent(const Data& grid);

    /// Corner-point input retained for refinement (nullptr unless the deck
    /// requested LGRs at construction time).
    static void setRetainedCornerPointInput(Data& grid,
                                            std::shared_ptr<const RetainedCornerPointInput> input);
    static std::shared_ptr<const RetainedCornerPointInput> retainedCornerPointInput(const Data& grid);

    /// Per-cell partition types (PartitionTypeIndicator::cell_indicator_).
    /// Empty means "not parallel; all interior". Used to give the leaf its
    /// partition types (coarse cells from level zero, refined cells from
    /// their parent) on a distributed grid.
    static const std::vector<char>& cellPartitionTypes(const Data& grid);
    static void setCellPartitionTypes(Data& grid, std::vector<char> types);

    /// Global ids of a distributed grid's cells, faces and points, by local
    /// index; replaces ids derived from local numbering.
    static void setGlobalIdMapping(Data& grid,
                                   std::vector<int> cellIds,
                                   std::vector<int> faceIds,
                                   std::vector<int> pointIds);
};

} // namespace Refinement
} // namespace Opm

#endif // OPM_GRID_REFINEMENT_GRIDSTATEWRITER_HEADER_INCLUDED
