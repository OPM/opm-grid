//===========================================================================
//
// File: CpGridLGR.cpp
//
// Created: 2025
//
// Author(s): Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
//===========================================================================

/*
  Copyright 2025 Equinor ASA.

  This file is part of The Open Porous Media project  (OPM).

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

#if HAVE_MPI
#include <opm/grid/utility/platform_dependent/disable_warnings.h>
#include "mpi.h"
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>
#endif

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#endif

#include "../CpGridLGR.hpp"
#include "LgrHelpers.hpp"
#include "ParentToChildrenCellGlobalIdHandle.hpp"
#include "NestedRefinementUtilities.hpp"
#include <opm/grid/common/CommunicationUtils.hpp>

#include <algorithm>
#include <iomanip>
#include <numeric>
#include <tuple>

namespace Dune
{

CpGridLGR::CpGridLGR()
    : CpGrid()
{
}

CpGridLGR::CpGridLGR(MPIHelper::MPICommunicator comm)
    : CpGrid(comm)
{
}

void CpGridLGR::computeGlobalCellLgr(const int& level, const std::array<int,3>& startIJK, std::vector<int>& global_cell_lgr)
{
    assert(level);
    for (const auto& element : elements(levelGridView(level))) {
        // Element belogns to an LGR, therefore has a father. Get IJK of the father in the level grid the father was born.
        // For CARFIN, parent cells belong to level 0.
        std::array<int,3> parentIJK = {0,0,0};
        currentData()[element.father().level()]->getIJK(element.father().index(), parentIJK);
        // Each parent cell has been refined in cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2] child cells.
        // element has certain 'position' inside its parent cell that can be described with 'IJK' indices, let's denote them by ijk,
        // where 0<= i < cells_per_dim[0], 0<= j < cells_per_dim[1], 0<= k < cells_per_dim[2].
        const auto& cells_per_dim = currentData()[level]->cells_per_dim_;
        //
        // Refined cell (here 'element') has "index in parent cell": k*cells_per_dim[0]*cells_per_dim[1] + j*cells_per_dim[0] + i
        // and it's stored in  cell_to_idxInParentCell_.
        auto idx_in_parent_cell =  currentData()[level]-> cell_to_idxInParentCell_[element.index()];
        // Find ijk.
        std::array<int,3> childIJK = Opm::Lgr::getIJK(idx_in_parent_cell, cells_per_dim);
        // The corresponding lgrIJK can be computed as follows:
        const std::array<int,3>& lgrIJK = { ( (parentIJK[0] - startIJK[0])*cells_per_dim[0] ) + childIJK[0],  // Shift parent index according to the startIJK of the LGR.
                                            ( (parentIJK[1] - startIJK[1])*cells_per_dim[1] ) + childIJK[1],
                                            ( (parentIJK[2] - startIJK[2])*cells_per_dim[2] ) + childIJK[2] };
        // Dimensions of the "patch of cells" formed when providing startIJK and endIJK for an LGR
        const auto& lgr_logical_cartesian_size = currentData()[level]->logical_cartesian_size_;
        global_cell_lgr[element.index()] = (lgrIJK[2]*lgr_logical_cartesian_size[0]*lgr_logical_cartesian_size[1]) + (lgrIJK[1]*lgr_logical_cartesian_size[0]) + lgrIJK[0];
    }
}

void CpGridLGR::computeGlobalCellLeafGridViewWithLgrs(std::vector<int>& global_cell_leaf)
{
    for (const auto& element: elements(leafGridView())) {
        // In the context of allowed nested refinement, we lookup for the oldest ancestor, belonging to level-zero-grid.
        auto ancestor = element.getOrigin();
        int origin_in_level_zero = ancestor.index();
        assert(origin_in_level_zero < currentData().front()->size(0));
        global_cell_leaf[element.index()] = currentData().front()-> global_cell_[origin_in_level_zero];
    }
}

std::vector<std::unordered_map<std::size_t, std::size_t>> CpGridLGR::mapLocalCartesianIndexSetsToLeafIndexSet() const
{
    std::vector<std::unordered_map<std::size_t, std::size_t>> localCartesianIdxSets_to_leafIdx(maxLevel()+1); // Plus level 0
    for (const auto& element : elements(leafGridView())) {
        const auto& global_cell_level = currentData()[element.level()]->globalCell()[element.getLevelElem().index()];
        localCartesianIdxSets_to_leafIdx[element.level()][global_cell_level] = element.index();
    }
    return localCartesianIdxSets_to_leafIdx;
}

std::vector<std::array<int,2>> CpGridLGR::mapLeafIndexSetToLocalCartesianIndexSets() const
{
    std::vector<std::array<int,2>> leafIdx_to_localCartesianIdxSets(currentLeafData().size(0));
    for (const auto& element : elements(leafGridView())) {
        const auto& global_cell_level = currentData()[element.level()]->globalCell()[element.getLevelElem().index()];
        leafIdx_to_localCartesianIdxSets[element.index()] = {element.level(), global_cell_level};
    }
    return leafIdx_to_localCartesianIdxSets;
}

void CpGridLGR::globalRefine (int refCount, bool throwOnFailure)
{
    if (refCount < 0) {
        OPM_THROW(std::logic_error, "Invalid argument. Provide a nonnegative integer for global refinement.");
    }

    // Throw if strict local refinement is detected.
    if(this->maxLevel()) {
        // For strict local refinement, these sizes are identical. Global refinement
        // results in a difference in at least one direction (x, y, or z).
        bool sameNX = currentLeafData().logicalCartesianSize()[0] == currentData().front()->logicalCartesianSize()[0];
        bool sameNY = currentLeafData().logicalCartesianSize()[1] == currentData().front()->logicalCartesianSize()[1];
        bool sameNZ = currentLeafData().logicalCartesianSize()[2] == currentData().front()->logicalCartesianSize()[2];
        if (sameNX && sameNY && sameNZ) {
            OPM_THROW(std::logic_error, "Global refinement of a mixed grid with coarse and refined cells is not supported yet.");
        }
    }
    if (refCount>0) {
        for (int refinedLevel = 0; refinedLevel < refCount; ++refinedLevel) {

            std::vector<int> assignRefinedLevel(current_data_->back()-> size(0));
            const auto& preAdaptMaxLevel = this ->maxLevel();
            std::vector<std::string> lgr_name_vec = { "GR" + std::to_string(preAdaptMaxLevel +1) };
            const std::array<int,3>& endIJK = currentLeafData().logicalCartesianSize();

            for(const auto& element: elements(this-> leafGridView())) {
                // Mark all the elements of the current leaf grid view for refinement
                if (mark(1, element, throwOnFailure))
                    assignRefinedLevel[element.index()] = preAdaptMaxLevel +1;
            }

            preAdapt();
            refineAndUpdateGrid(/* cells_per_dim_vec = */ {{2,2,2}}, assignRefinedLevel, lgr_name_vec, {{0,0,0}}, {endIJK});
            postAdapt();
        }
    }
}

template <int codim>

/// \brief Size of the ghost cell layer on the leaf level

/// \brief Size of the overlap on a given level

/// \brief Size of the ghost cell layer on a given level

//

/// \brief Get the number of faces.
/// \brief Get The number of vertices.

Dune::cpgrid::Intersection CpGridLGR::getParentIntersectionFromLgrBoundaryFace(const Dune::cpgrid::Intersection& intersection) const
{
    if ( intersection.neighbor()) {
        if ((intersection.inside().level() != intersection.outside().level())) {
            // one coarse and one refined neighboring cell
            /** Now, it could also be two refined cells. In that case, any of them will fit to search for the parent face */
            const auto& cellIn = intersection.inside();
            const auto& cellOut = intersection.outside();

            // Identify the coarse and the refined neighboring cell
            const auto coarseCell =  (cellIn.level() == 0) ? cellIn : cellOut;
            const auto refinedCell =  (coarseCell == cellIn) ? cellOut : cellIn;
            assert(coarseCell.level() != refinedCell.level());

            // Get parent cell of the refined cell
            const auto& parentCell = refinedCell.father();

            // Get the index inside and orientation from the leaf grid (refined) face
            const auto& intersectionIdxInInside = intersection.indexInInside();

            for(const auto& parentIntersection : intersections(this->levelGridView(0), parentCell)){
                // Get the inInsideIdx and orientation from the parent intersection
                const auto& parentIdxInInside = parentIntersection.indexInInside();
                if (parentIdxInInside == intersectionIdxInInside) {
                    return parentIntersection;
                }
            }
        }
        OPM_THROW(std::invalid_argument, "Parent intersection not found for face with index: " + std::to_string(intersection.id()) +
                  " and index in inside: " + std::to_string(intersection.indexInInside()));
    }
    OPM_THROW(std::invalid_argument, "Face is on the boundary of the grid");
}

void CpGridLGR::markElemAssignLevelDetectActiveLgrs(const std::vector<std::array<int,3>>& startIJK_vec,
                                                 const std::vector<std::array<int,3>>& endIJK_vec,
                                                 std::vector<int>& assignRefinedLevel,
                                                 std::vector<int>& lgr_with_at_least_one_active_cell)
{
    auto assignAndDetect = [this, &assignRefinedLevel, &lgr_with_at_least_one_active_cell](const cpgrid::Entity<0>& element, int level)
    {
        mark(1, element, /* throwOnFailure = */ true);
        assignRefinedLevel[element.index()] = level+1;
        // shifted since starting grid is level 0, and refined grids levels are >= 1.
        lgr_with_at_least_one_active_cell[level] = 1;
    };
    Opm::Lgr::computeOnLgrParents(*this, startIJK_vec, endIJK_vec, assignAndDetect);
}

void CpGridLGR::populateCellIndexSetRefinedGrid([[maybe_unused]] int level)
{
#if HAVE_MPI
    const auto& level_data = currentData()[level];
    const auto& level_global_id_set =  level_data->global_id_set_;
    auto& level_index_set =  level_data->cellIndexSet();
    // Clean up level cell index set - needed e.g. for synchronization of cell ids.
    level_index_set = ParallelIndexSet();

    level_index_set.beginResize();
    // ParallelIndexSet::LocalIndex( elemIdx, attribute /* owner or copy */, true/false)
    // The bool indicates whether the global index might exist on other processes with other attributes.
    // RemoteIndices::rebuild has a Boolean as a template parameter telling the method whether to ignore this
    // Boolean on the indices or not when building.
    //
    // For refined level grids, we check if the parent cell is fully interior. Then, its children won't be seen
    // by any other process. Therefore, the boolean is set to false.
    for(const auto& element : elements(levelGridView(level))) {
        if ( element.partitionType() == InteriorEntity) {

            // Check if it has an overlap neighbor
            bool parentFullyInterior = true;
            const auto& parent_cell = element.father();

            const auto& intersections = Dune::intersections(levelGridView(parent_cell.level()), parent_cell);
            for (const auto& intersection : intersections) {
                if ( intersection.neighbor() ) {
                    const auto& neighborPartitionType = intersection.outside().partitionType();
                    // To help detection of fully interior cells, i.e., without overlap neighbors
                    if (neighborPartitionType == OverlapEntity )  {
                        parentFullyInterior = false;
                        // Interior cells with overlap neighbor may appear in other processess -> true
                        level_index_set.add( level_global_id_set->id(element),
                                             ParallelIndexSet::LocalIndex(element.index(), AttributeSet(AttributeSet::owner), true));
                        // Store it only once
                        break;
                    }
                }
            }
            if(parentFullyInterior) { // Fully interior cells do not appear in any other process -> false
                level_index_set.add( level_global_id_set->id(element),
                                     ParallelIndexSet::LocalIndex(element.index(), AttributeSet(AttributeSet::owner), false));
            }
        }
        else { // overlap cell
            assert(element.partitionType() == OverlapEntity);
            level_index_set.add( level_global_id_set->id(element),
                                 ParallelIndexSet::LocalIndex(element.index(), AttributeSet(AttributeSet::copy), true));
        }
    }
    level_index_set.endResize();

    currentData()[level]->cellRemoteIndices().template rebuild<true>(); // true->ignore bool in ParallelIndexSet!
    // Check the size
    assert(static_cast<std::size_t>(level_data->cellIndexSet().size()) == static_cast<std::size_t>(level_data->size(0)) );
#endif
}

void CpGridLGR::populateCellIndexSetLeafGridView()
{
#if HAVE_MPI
    auto& leaf_index_set =  (*current_data_).back()->cellIndexSet();
    // Clean up leaf cell index set - needed e.g. for synchronization of cell ids.
    leaf_index_set = ParallelIndexSet();

    leaf_index_set.beginResize();

    // ParallelIndexSet::LocalIndex( elemIdx, attribute /* owner or copy */, true/false)
    // The bool indicates whether the global index might exist on other processes with other attributes.
    // RemoteIndices::rebuild has a Boolean as a template parameter telling the method whether to ignore this
    // Boolean on the indices or not when building.
    for(const auto& element : elements(leafGridView())) {
        const auto& elemPartitionType = element.getLevelElem().partitionType();
        if ( elemPartitionType == InteriorEntity) {
            // Check if it has an overlap neighbor
            bool isFullyInterior = true;
            for (const auto& intersection : intersections(leafGridView(), element)) {
                if ( intersection.neighbor() ) {
                    const auto& neighborPartitionType = intersection.outside().getLevelElem().partitionType();
                    // To help detection of fully interior cells, i.e., without overlap neighbors
                    if (neighborPartitionType == OverlapEntity )  {
                        isFullyInterior = false;
                        // Interior cells with overlap neighbor may appear in other processess -> false
                        leaf_index_set.add(globalIdSet().id(element),
                                           ParallelIndexSet::LocalIndex(element.index(), AttributeSet(AttributeSet::owner), true));
                        // Store it only once
                        break;
                    }
                }
            }
            if(isFullyInterior) { // Fully interior cells do not appear in any other process -> false
                leaf_index_set.add(globalIdSet().id(element),
                                   ParallelIndexSet::LocalIndex(element.index(), AttributeSet(AttributeSet::owner), false));
            }
        }
        else { // overlap cell
            assert(elemPartitionType == OverlapEntity);
            leaf_index_set.add(globalIdSet().id(element),
                               ParallelIndexSet::LocalIndex(element.index(), AttributeSet(AttributeSet::copy), true));
        }
    }
    leaf_index_set.endResize();

    (*current_data_).back()->cellRemoteIndices().template rebuild<true>(); // true->ignore bool in ParallelIndex!

#endif
}

void CpGridLGR::populateLeafGlobalIdSet()
{
    // Global id for the cells in leaf grid view
    std::vector<int> leafCellIds(current_data_->back()->size(0));
    for(const auto& element: elements(leafGridView())){
        // Notice that for level zero cells the global_id_set_ is given, for refined level grids was defined
        // under the assumption of each lgr being fully contained in the interior of a process.
        // Therefore, it is not needed here to distingish between owned and overlap cells.
        auto equivElem = element.getLevelElem();
        leafCellIds[element.index()] = (*current_data_)[element.level()]->global_id_set_->id(equivElem);
    }

    // Global id for the faces in leaf grid view. Empty vector (Entity<1> not supported for CpGrid).
    std::vector<int> leafFaceIds{};

    // Global id for the points in leaf grid view
    std::vector<int> leafPointIds(current_data_->back()->size(3));
    for(const auto& point : vertices(leafGridView())){
        const auto& level_pointLevelIdx = current_data_->back()->corner_history_[point.index()];
        assert(level_pointLevelIdx[0] != -1);
        assert(level_pointLevelIdx[1] != -1);
        const auto& pointLevelEntity =  cpgrid::Entity<3>(*( (*current_data_)[level_pointLevelIdx[0]]), level_pointLevelIdx[1], true);
        leafPointIds[point.index()] = (*current_data_)[level_pointLevelIdx[0]]->global_id_set_->id(pointLevelEntity);
    }

    current_data_->back()->global_id_set_->swap(leafCellIds, leafFaceIds, leafPointIds);
}

bool CpGridLGR::mark(int refCount, const cpgrid::Entity<0>& element, bool throwOnFailure)
{
    // Throw if element has a neighboring cell from a different level.
    // E.g., a coarse cell touching the boundary of an LGR, or
    // a refined cell with a coarser/finner neighboring cell.
    for (const auto& intersection : Dune::intersections(leafGridView(), element)){
        if (intersection.neighbor() && (intersection.outside().level() != element.level()) && (element.level()==0)) {
            // Refinement of cells at LGR boundaries is not supported, yet.
            if (throwOnFailure)
                OPM_THROW(std::invalid_argument, "Refinement of cells at LGR boundaries is not supported, yet.");
            else
                return false;
        }
    }
    // For serial run, mark elements also in the level they were born.
    std::optional<bool> mark0;
    if(currentData().size()>1) {
        // Mark element in its level
        mark0 = currentData()[element.level()] -> mark(refCount, element.getLevelElem(), throwOnFailure);
    }
    // Mark element (also in the serial run case) in leaf grid view. Note that if scatterGrid has been invoked, then
    // current_data_ == distributed_data_.
    bool mark1 = current_data_->back()->mark(refCount, element, throwOnFailure);
    // Consistency check between level view and current view (serial run case).
    if (mark0.value_or(mark1) != mark1)
        OPM_THROW(std::logic_error, "Inconsistent marking state between current view and level view.");
    return mark1;
}

int CpGridLGR::getMark(const cpgrid::Entity<0>& element) const
{
    return current_data_->back()->getMark(element);
}

bool CpGridLGR::preAdapt()
{
    // Check if elements in pre-adapt existing grids have been marked for refinment.
    // Serial run: currentData() = data_. Parallel run: currentData() = distributed_data_.
    bool isPreAdapted = false; // 0
    for (const auto& preAdaptGrid : currentData()) {
        isPreAdapted = std::max(isPreAdapted, preAdaptGrid -> preAdapt()); // could be 0 or 1
    }
    // If at least one process has marked elements, return true.
    return this->comm().max(isPreAdapted);
}

bool CpGridLGR::adapt()
{
    if(!preAdapt()) { // marked cells set can be empty
        return false; // the grid does not change at all.
    }

    const std::vector<std::array<int,3>>& cells_per_dim_vec = {{2,2,2}}; // Arbitrary chosen values.
    std::vector<int> assignRefinedLevel(current_data_->back()-> size(0));
    const auto& preAdaptMaxLevel = this ->maxLevel();

    int local_marked_elem_count = 0;
    for (int elemIdx = 0; elemIdx < current_data_->back()->size(0); ++elemIdx) {
        const auto& element = cpgrid::Entity<0>(*(current_data_->back()), elemIdx, true);
        assignRefinedLevel[elemIdx] = (this->getMark(element) == 1) ? (preAdaptMaxLevel +1) : 0;
        if (this->getMark(element) == 1) {
            ++local_marked_elem_count;
        }
    }

    std::vector<std::string> lgr_name_vec = { "LGR" + std::to_string(preAdaptMaxLevel +1) };

    auto global_marked_elem_count = comm().sum(local_marked_elem_count);
    auto global_cell_count_before_adapt = comm().sum(current_data_->back()-> size(0)); // Recall overlap cells are also marked
    // Check if its a global refinement
    bool is_global_refine = (global_marked_elem_count == global_cell_count_before_adapt);
    if (is_global_refine) { // parallel or sequential
        // Rewrite the lgr name (GR stands for GLOBAL REFINEMET)
        lgr_name_vec = { "GR" + std::to_string(preAdaptMaxLevel +1) };
        const std::array<int,3>& endIJK = currentLeafData().logicalCartesianSize();
        return this->refineAndUpdateGrid(cells_per_dim_vec, assignRefinedLevel, lgr_name_vec, {{0,0,0}}, {endIJK});
    }
    return this-> refineAndUpdateGrid(cells_per_dim_vec, assignRefinedLevel, lgr_name_vec);
}

bool CpGridLGR::refineAndUpdateGrid(const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                 const std::vector<int>& assignRefinedLevel,
                                 const std::vector<std::string>& lgr_name_vec,
                                 const std::vector<std::array<int,3>>& startIJK_vec,
                                 const std::vector<std::array<int,3>>& endIJK_vec)
{
    // To do: support coarsening.
    assert( static_cast<int>(assignRefinedLevel.size()) == currentLeafData().size(0));
    assert(cells_per_dim_vec.size() == lgr_name_vec.size());

    auto& data = currentData(); // data pointed by current_data_ (data_ or distributed_data_[if loadBalance() has been invoked before adapt()]).
    // Logical Cartesian Size before adapting the grid - to be used in case the entire grid will be refined.
    const auto& lcs =  data.back()->logical_cartesian_size_;

    bool isCARFIN = !startIJK_vec.empty();
    bool isGlobalRefine = isCARFIN; // One way of global-refine the grid is via startIJK and endIJK.
    for (int c = 0; c<3; ++c) {
        isGlobalRefine = isGlobalRefine && (startIJK_vec[0][c] == 0 ) && (endIJK_vec[0][c] == lcs[c]);
        if (!isGlobalRefine)
            break;
    }

    // Each marked element has its assigned level where its refined entities belong.
    const int& levels = cells_per_dim_vec.size();
    // Notice that "levels" represents also the total amount of new (after calling adapt) refined level grids.
    const int& preAdaptMaxLevel = this->maxLevel();
    // Copy corner history - needed to compute later ids, empty vector if the grid to be adapted is level 0 grid, or the grid has been distributed.
    const auto& preAdaptGrid_corner_history = (preAdaptMaxLevel>0) ? current_data_->back()->corner_history_ : std::vector<std::array<int,2>>();

    if (!global_id_set_ptr_) {
        global_id_set_ptr_ = std::make_shared<cpgrid::GlobalIdSet>(*data.back());

        for (int preAdaptLevel = 0; preAdaptLevel <= preAdaptMaxLevel; ++preAdaptLevel) {
            global_id_set_ptr_->insertIdSet(*data[preAdaptLevel]);
        }
    }

    // To determine if an LGR is not empty in a given process, we set
    // lgr_with_at_least_one_active_cell[in that level] to 1 if it contains
    // at least one active cell, and to 0 otherwise.
    std::vector<int> lgr_with_at_least_one_active_cell(levels);
    Opm::Lgr::detectActiveLgrs(*this, startIJK_vec, endIJK_vec, lgr_with_at_least_one_active_cell);

    // To store/build refined level grids.
    std::vector<std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>> refined_data_vec(levels, data);
    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>> refined_grid_ptr_vec(levels);

    std::vector<Dune::cpgrid::DefaultGeometryPolicy> refined_geometries_vec(levels);
    std::vector<std::vector<std::array<int,8>>> refined_cell_to_point_vec(levels);
    std::vector<cpgrid::OrientedEntityTable<0,1>> refined_cell_to_face_vec(levels);
    std::vector<Opm::SparseTable<int>> refined_face_to_point_vec(levels);
    std::vector<cpgrid::OrientedEntityTable<1,0>> refined_face_to_cell_vec(levels);

    // Mutable containers for refined corners, faces, cells, face tags, and face normals.
    std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>> refined_corners_vec(levels);
    std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>> refined_faces_vec(levels);
    std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>> refined_cells_vec(levels);
    std::vector<Dune::cpgrid::EntityVariableBase<enum face_tag>> mutable_refined_face_tags_vec(levels);
    typedef Dune::FieldVector<double,3> PointType;
    std::vector<Dune::cpgrid::EntityVariableBase<PointType>> mutable_refined_face_normals_vec(levels);

    std::vector<std::vector<int>> refined_global_cell_vec(levels);

    // To store adapted grid
    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& adaptedData = data;
#if HAVE_MPI
    auto adaptedGrid_ptr =
        std::make_shared<Dune::cpgrid::CpGridData>((*(data[0])).ccobj_, adaptedData);
#else
    // DUNE 2.7 is missing convertion to NO_COMM
    auto adaptedGrid_ptr = std::make_shared<Dune::cpgrid::CpGridData>(adaptedData);
#endif
    auto& adaptedGrid = *adaptedGrid_ptr;
    Dune::cpgrid::DefaultGeometryPolicy&                         adapted_geometries = adaptedGrid.geometry_;
    std::vector<std::array<int,8>>&                              adapted_cell_to_point = adaptedGrid.cell_to_point_;
    cpgrid::OrientedEntityTable<0,1>&                            adapted_cell_to_face = adaptedGrid.cell_to_face_;
    Opm::SparseTable<int>&                                       adapted_face_to_point = adaptedGrid.face_to_point_;
    cpgrid::OrientedEntityTable<1,0>&                            adapted_face_to_cell = adaptedGrid.face_to_cell_;
    cpgrid::EntityVariable<enum face_tag,1>&                     adapted_face_tags = adaptedGrid.face_tag_;
    cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>,1>& adapted_face_normals = adaptedGrid.face_normals_;
    // Mutable containers for adapted corners, faces, cells, face tags, and face normals.
    Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>& adapted_corners =
        *(adapted_geometries.geomVector(std::integral_constant<int,3>()));
    Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>& adapted_faces =
        *(adapted_geometries.geomVector(std::integral_constant<int,1>()));
    Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>& adapted_cells =
        *(adapted_geometries.geomVector(std::integral_constant<int,0>()));
    Dune::cpgrid::EntityVariableBase<enum face_tag>& mutable_face_tags = adapted_face_tags;
    typedef Dune::FieldVector<double,3> PointType;
    Dune::cpgrid::EntityVariableBase<PointType>& mutable_face_normals = adapted_face_normals;

    // Refine marked cells and provide marked-corner/face/cell - refined-corner/faces/cells relations.
    //
    // ------------------------ Marked elements parameters
    // -- markedElem_to_itsLgr :
    // Each marked element gets refined and we store this "auxiliary markedElementLGR", to later
    // build a unique level containing all the refined entities from all the marked elements.
    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData> > markedElem_to_itsLgr;
    markedElem_to_itsLgr.resize(current_data_->back()->size(0));
    // -- markedElem_count: Total amount of marked elements to be refined. It will be used to print grid info.
    int markedElem_count = 0;
    // -- cornerInMarkedElemWithEquivRefinedCorner :
    // For each corner from level zero, we store the marked elements where the corner appears and its equivalent
    // refined corner in  each auxiliary marked-element-lgr. Example: corner with index 5 appears in marked
    // elements 0 and 1, with refined equivalent corner indices 8 and 2 respectively. Then,
    // cornerInMarkedElemWithEquivRefinedCorner[5] = {{0, 8}, {1, 2}}.
    // For corners not appearing in any marked element, empty vector.
    std::vector<std::vector<std::array<int,2>>> cornerInMarkedElemWithEquivRefinedCorner;
    cornerInMarkedElemWithEquivRefinedCorner.resize(current_data_->back()->size(3) );
    // -- markedElemAndEquivRefinedCorner_to_corner :
    // To correctly build the level-refined and adapted-grid topology features, we need to keep track of the
    // corners that got replaced by equivalent refined corners, in each marked element where the corner appeared,
    // not only in its last appearance. The last appearance will be used to avoid repetition when storing.
    // Following the example above,
    // markedElemAndEquivRefinedCorner_to_corner[{0, 8}] = 5;
    // markedElemAndEquivRefinedCorner_to_corner[{1, 2}] = 5;
    std::map< std::array<int,2>, int > markedElemAndEquivRefinedCorn_to_corner;
    // -- faceInMarkedElemAndRefinedFaces :
    // For each face from level zero, we store the marked elements where the face appears (maximum 2 cells)
    // and its new-born refined faces from each auxiliary marked-element-lgr. Example: face with index 9
    // appears in marked elements 0 and 1. Then,
    // faceInMarkedElemAndRefinedFaces[9] = {{0, {refinedFace0_0, ..., refinedFaceN_0}},
    //                                       {1, {refinedFace0_1, ..., refinedFaceM_1}}}.
    // For faces not appearing in any marked element, empty vector.
    std::vector<std::vector<std::pair<int, std::vector<int>>>> faceInMarkedElemAndRefinedFaces;
    faceInMarkedElemAndRefinedFaces.resize(current_data_->back()->face_to_cell_.size());
    // ------------------------ Refined cells parameters
    // --- Refined cells and PreAdapt cells relations ---
    std::map<std::array<int,2>,std::array<int,2>> elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell;
    std::map<std::array<int,2>,std::array<int,2>> refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell;
    // Integer to count only REFINED cells (new-born refined cells from ANY marked element).
    std::vector<int> refined_cell_count_vec(levels, 0);
    // -- Parent-child relations --
    // Relation between the grids before adapt() and new refined cells (on the refined grid level - not in each individual lgr).
    std::vector<std::vector<std::tuple<int,std::vector<int>>>> preAdapt_parent_to_children_cells_vec(preAdaptMaxLevel +1);
    // ------------------------ Adapted cells parameters
    // --- Adapted cells and PreAdapt cells relations ---
    std::map<std::array<int,2>,int>           elemLgrAndElemLgrCell_to_adaptedCell;
    std::unordered_map<int,std::array<int,2>> adaptedCell_to_elemLgrAndElemLgrCell;
    // Integer to count adapted cells (mixed between cells from level0 (not involved in LGRs), and (new-born) refined cells).
    int cell_count = 0;
    // -- Some extra indices relations between preAdapt-grid and adapted-grid --
    // Relation between the grids before adapt() and leafview cell indices.
    std::vector<std::vector<int>> preAdapt_level_to_leaf_cells_vec(preAdaptMaxLevel +1);
    for (int preAdaptLevel = 0; preAdaptLevel < preAdaptMaxLevel +1; ++preAdaptLevel) {
        // Resize with the corresponding amount of cells of the preAdapt level. Deafualt {-1, empty vector} when the cell has no children.
        if ( (*data[preAdaptLevel]).parent_to_children_cells_.empty()){
            preAdapt_parent_to_children_cells_vec[preAdaptLevel].resize(data[preAdaptLevel]->size(0), std::make_pair(-1, std::vector<int>{}));
        }
        else {
            preAdapt_parent_to_children_cells_vec[preAdaptLevel] =  (*data[preAdaptLevel]).parent_to_children_cells_;
        }
        // Resize with the corresponding amount of cell of the preAdapt level. Dafualt -1 when the cell vanished and does not appear on the leaf grid view.
        // In entry 'level cell index', we store 'leafview cell index', or -1 when the cell vanished.
        preAdapt_level_to_leaf_cells_vec[preAdaptLevel].resize(data[preAdaptLevel]->size(0), -1);
    }
    //
    Opm::Lgr::refineAndProvideMarkedRefinedRelations( *this,
                                                      /* Marked elements parameters */
                                                      markedElem_to_itsLgr,
                                                      markedElem_count,
                                                      cornerInMarkedElemWithEquivRefinedCorner,
                                                      markedElemAndEquivRefinedCorn_to_corner,
                                                      faceInMarkedElemAndRefinedFaces,
                                                      /* Refined cells parameters */
                                                      elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell,
                                                      refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                                      refined_cell_count_vec,
                                                      assignRefinedLevel,
                                                      preAdapt_parent_to_children_cells_vec,
                                                      /* Adapted cells parameters */
                                                      elemLgrAndElemLgrCell_to_adaptedCell,
                                                      adaptedCell_to_elemLgrAndElemLgrCell,
                                                      cell_count,
                                                      preAdapt_level_to_leaf_cells_vec,
                                                      /* Additional parameters */
                                                      cells_per_dim_vec);

#if HAVE_MPI
    auto global_markedElem_count = comm().sum(markedElem_count);
    if ( global_markedElem_count == 0 ) {
        return false;
    }
#else
    if ( markedElem_count == 0 ) {
        return false;
    }
#endif

    // Update/define parent_to_children_cells_ and level_to_leaf_cells_ for all the existing level grids (level 0, 1, ..., preAdaptMaxLevel), before this call of adapt.
    for (int preAdaptLevel = 0; preAdaptLevel < preAdaptMaxLevel +1; ++preAdaptLevel) {
        (*data[preAdaptLevel]).parent_to_children_cells_ = preAdapt_parent_to_children_cells_vec[preAdaptLevel];
        (*data[preAdaptLevel]).level_to_leaf_cells_ =  preAdapt_level_to_leaf_cells_vec[preAdaptLevel];
    }

    // -- Child-parent relations --
    // refined_child_to_parent_cells_vec:   Refined child cells and their parents. Entry is {-1,-1} when cell has no father.
    //                                      Otherwise, {level parent cell, parent cell index}. Each entry represents a refined level.
    // refined_cell_to_idxInParentCell_vec: Each refined child cell has a unique index in its parent cell, to be used to build geometryInFather().
    //                                      Each entry represents a refined level.
    // adapted_child_to_parent_cells:       Adapted child cells and their parents. Entry is {-1,-1} when cell has no father. Otherwise, {level parent cell, parent cell index}
    // adapted_cell_to_idxInParentCell:     Each refined child cell has a unique index in its parent cell, to be used to build geometryInFather().
    //                                      When the cell has not been refined, -1.
    const auto& [refined_child_to_parent_cells_vec,
                 refined_cell_to_idxInParentCell_vec,
                 adapted_child_to_parent_cells,
                 adapted_cell_to_idxInParentCell] = Opm::Lgr::defineChildToParentAndIdxInParentCell(currentLeafData(),
                                                                                                    preAdaptMaxLevel,
                                                                                                    refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                                                                                    refined_cell_count_vec,adaptedCell_to_elemLgrAndElemLgrCell,
                                                                                                    cell_count);

    // -- Refined to Adapted cells and Adapted-cells to {level where the cell was born, cell index on that level} --
    // refined_level_to_leaf_cells_vec:  Relation between the refined grid and leafview cell indices.
    // leaf_to_level_cells:              Relation between an adapted cell and its equivalent cell coming either from pre-refined-leaf or the refined grid (level)
    const auto& [refined_level_to_leaf_cells_vec,
                 leaf_to_level_cells] = Opm::Lgr::defineLevelToLeafAndLeafToLevelCells(currentLeafData(),
                                                                                       preAdaptMaxLevel,
                                                                                       elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell,
                                                                                       refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                                                                       refined_cell_count_vec,
                                                                                       elemLgrAndElemLgrCell_to_adaptedCell,
                                                                                       adaptedCell_to_elemLgrAndElemLgrCell,
                                                                                       cell_count);

    // CORNERS
    // Stablish relationships between PreAdapt corners and refined or adapted ones ---
    //
    // --- Refined corners and PreAdapt corners relations ---
    std::map<std::array<int,2>,std::array<int,2>> elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner;
    std::map<std::array<int,2>,std::array<int,2>> refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner;
    std::map<std::array<int,2>,std::array<int,2>> vanishedRefinedCorner_to_itsLastAppearance;
    // Integer to count only refined corners.
    std::vector<int> refined_corner_count_vec(levels, 0);
    Opm::Lgr::identifyRefinedCornersPerLevel(currentLeafData(),
                                             preAdaptMaxLevel,
                                             /* Refined grid parameters */
                                             elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                             refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner,
                                             refined_corner_count_vec,
                                             vanishedRefinedCorner_to_itsLastAppearance,
                                             /* Additional parameters */
                                             markedElem_to_itsLgr,
                                             assignRefinedLevel,
                                             cornerInMarkedElemWithEquivRefinedCorner,
                                             faceInMarkedElemAndRefinedFaces,
                                             cells_per_dim_vec);

    // --- Adapted corners and PreAdapt corners relations ---
    std::map<std::array<int,2>,int>           elemLgrAndElemLgrCorner_to_adaptedCorner;
    std::unordered_map<int,std::array<int,2>> adaptedCorner_to_elemLgrAndElemLgrCorner;
    // Integer to count adapted corners (mixed between corners from pre-refined-leaf corners not involved in LGRs, and new-born refined corners).
    int corner_count = 0;
    Opm::Lgr::identifyLeafGridCorners(currentLeafData(),
                                      preAdaptMaxLevel,
                                      /* Adapted grid parameters */
                                      elemLgrAndElemLgrCorner_to_adaptedCorner,
                                      adaptedCorner_to_elemLgrAndElemLgrCorner,
                                      corner_count,
                                      /* Additional parameters */
                                      markedElem_to_itsLgr,
                                      assignRefinedLevel,
                                      cornerInMarkedElemWithEquivRefinedCorner,
                                      vanishedRefinedCorner_to_itsLastAppearance,
                                      faceInMarkedElemAndRefinedFaces,
                                      cells_per_dim_vec);

    // FACES
    // Stablish relationships between PreAdapt faces and refined or adapted ones ---
    // --- Refined faces and PreAdapt faces relations ---
    std::map<std::array<int,2>,std::array<int,2>> elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace;
    std::map<std::array<int,2>,std::array<int,2>> refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace;
    // Integer to count adapted faces (mixed between faces from level0 (not involved in LGRs), and (new-born) refined faces).
    std::vector<int> refined_face_count_vec(levels, 0);
    Opm::Lgr::identifyRefinedFacesPerLevel(currentLeafData(),
                                           preAdaptMaxLevel,
                                           elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                                           refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                                           refined_face_count_vec,
                                           markedElem_to_itsLgr,
                                           assignRefinedLevel,
                                           faceInMarkedElemAndRefinedFaces,
                                           cells_per_dim_vec);

    // --- Adapted faces and PreAdapt faces relations ---
    std::map< std::array<int,2>, int >           elemLgrAndElemLgrFace_to_adaptedFace;
    std::unordered_map< int, std::array<int,2> > adaptedFace_to_elemLgrAndElemLgrFace;
    // Integer to count adapted faces (mixed between faces from pre-refined-leaf faces not involved in LGRs and new-born refined faces).
    int face_count = 0;
    Opm::Lgr::identifyLeafGridFaces(currentLeafData(),
                                    preAdaptMaxLevel,
                                    elemLgrAndElemLgrFace_to_adaptedFace,
                                    adaptedFace_to_elemLgrAndElemLgrFace,
                                    face_count,
                                    markedElem_to_itsLgr,
                                    assignRefinedLevel,
                                    faceInMarkedElemAndRefinedFaces,
                                    cells_per_dim_vec);

    // Set refined level grids geometries
    // --- Refined corners  ---
    Opm::Lgr::populateRefinedCorners(refined_corners_vec,
                                     refined_corner_count_vec,
                                     markedElem_to_itsLgr,
                                     preAdaptMaxLevel,
                                     refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner);
    // --- Refined faces  ---
    Opm::Lgr::populateRefinedFaces(refined_faces_vec,
                                   mutable_refined_face_tags_vec,
                                   mutable_refined_face_normals_vec,
                                   refined_face_to_point_vec,
                                   refined_face_count_vec,
                                   refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                                   elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                   vanishedRefinedCorner_to_itsLastAppearance,
                                   markedElem_to_itsLgr,
                                   preAdaptMaxLevel,
                                   cornerInMarkedElemWithEquivRefinedCorner,
                                   markedElemAndEquivRefinedCorn_to_corner);
    // --- Refined cells  ---
    Opm::Lgr::populateRefinedCells(currentLeafData(),
                                   refined_cells_vec,
                                   refined_cell_to_point_vec,
                                   refined_global_cell_vec,
                                   refined_cell_count_vec,
                                   refined_cell_to_face_vec,
                                   refined_face_to_cell_vec,
                                   refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                   elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                                   faceInMarkedElemAndRefinedFaces,
                                   refined_geometries_vec,
                                   elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                   vanishedRefinedCorner_to_itsLastAppearance,
                                   markedElem_to_itsLgr,
                                   assignRefinedLevel,
                                   preAdaptMaxLevel,
                                   markedElemAndEquivRefinedCorn_to_corner,
                                   cornerInMarkedElemWithEquivRefinedCorner,
                                   cells_per_dim_vec);

    // Update leaf grid geometries
    // --- Adapted corners ---
    Opm::Lgr::populateLeafGridCorners(currentLeafData(),
                                      adapted_corners,
                                      corner_count,
                                      markedElem_to_itsLgr,
                                      adaptedCorner_to_elemLgrAndElemLgrCorner);
    // --- Adapted faces ---
    Opm::Lgr::populateLeafGridFaces(currentLeafData(),
                                    adapted_faces,
                                    mutable_face_tags,
                                    mutable_face_normals,
                                    adapted_face_to_point,
                                    face_count,
                                    adaptedFace_to_elemLgrAndElemLgrFace,
                                    elemLgrAndElemLgrCorner_to_adaptedCorner,
                                    vanishedRefinedCorner_to_itsLastAppearance,
                                    markedElem_to_itsLgr,
                                    assignRefinedLevel,
                                    markedElemAndEquivRefinedCorn_to_corner,
                                    cornerInMarkedElemWithEquivRefinedCorner,
                                    cells_per_dim_vec,
                                    preAdaptMaxLevel);
    // --- Adapted cells ---
    Opm::Lgr::populateLeafGridCells(currentLeafData(),
                                    adapted_cells,
                                    adapted_cell_to_point,
                                    cell_count,
                                    adapted_cell_to_face,
                                    adapted_face_to_cell,
                                    adaptedCell_to_elemLgrAndElemLgrCell,
                                    elemLgrAndElemLgrFace_to_adaptedFace,
                                    faceInMarkedElemAndRefinedFaces,
                                    adapted_geometries,
                                    elemLgrAndElemLgrCorner_to_adaptedCorner,
                                    vanishedRefinedCorner_to_itsLastAppearance,
                                    markedElem_to_itsLgr,
                                    assignRefinedLevel,
                                    markedElemAndEquivRefinedCorn_to_corner,
                                    cornerInMarkedElemWithEquivRefinedCorner,
                                    cells_per_dim_vec,
                                    preAdaptMaxLevel);

    for (int level = 0; level < levels; ++level) {
        const int refinedLevelGridIdx = level + preAdaptMaxLevel +1;
#if HAVE_MPI
        refined_grid_ptr_vec[level] = std::make_shared<Dune::cpgrid::CpGridData>((*(data[0])).ccobj_, refined_data_vec[level]);
#else
        // DUNE 2.7 is missing convertion to NO_COMM
        refined_grid_ptr_vec[level] = std::make_shared<Dune::cpgrid::CpGridData>(refined_data_vec[level]);
#endif
        // Store refined grid
        if ((level == 0) && (preAdaptMaxLevel>0)) { // Overwrite the leaf-grid-view with the first new-refined-level-grid
            data[preAdaptMaxLevel+1] = refined_grid_ptr_vec[level];
        }
        else {
            data.push_back(refined_grid_ptr_vec[level]);
        }

        Dune::cpgrid::DefaultGeometryPolicy&  refinedLevel_geometries = (*data[refinedLevelGridIdx]).geometry_;
        // Mutable containers for adapted corners, faces, cells, face tags, and face normals.
        Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>& level_corners =
            *(refinedLevel_geometries.geomVector(std::integral_constant<int,3>()));
        Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>& level_faces =
            *(refinedLevel_geometries.geomVector(std::integral_constant<int,1>()));
        Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>& level_cells =
            *(refinedLevel_geometries.geomVector(std::integral_constant<int,0>()));

        level_corners.swap(refined_corners_vec[level]);
        level_faces.swap(refined_faces_vec[level]);
        level_cells.swap(refined_cells_vec[level]);

        (*data[refinedLevelGridIdx]).cell_to_point_.swap(refined_cell_to_point_vec[level]);
        (*data[refinedLevelGridIdx]).cell_to_face_.swap(refined_cell_to_face_vec[level]);

        (*data[refinedLevelGridIdx]).face_to_point_.swap(refined_face_to_point_vec[level]);
        (*data[refinedLevelGridIdx]).face_to_cell_.swap(refined_face_to_cell_vec[level]);

        cpgrid::EntityVariable<enum face_tag,1>& level_face_tags =   (*data[refinedLevelGridIdx]).face_tag_;
        Dune::cpgrid::EntityVariableBase<enum face_tag>& level_mutable_face_tags = level_face_tags;
        level_mutable_face_tags.swap(mutable_refined_face_tags_vec[level]);

        cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>,1>&  level_face_normals =   (*data[refinedLevelGridIdx]).face_normals_;
        Dune::cpgrid::EntityVariableBase<PointType>& level_mutable_face_normals = level_face_normals;
        level_mutable_face_normals.swap(mutable_refined_face_normals_vec[level]);

        // Further Refined grid Attributes
        //
        // Populate some attributes of the level LGR
        (*data[refinedLevelGridIdx]).level_data_ptr_ = &(this -> currentData());
        (*data[refinedLevelGridIdx]).level_ = refinedLevelGridIdx;
        this -> lgr_names_[lgr_name_vec[level]] = refinedLevelGridIdx; // {"name_lgr", level}
        (*data[refinedLevelGridIdx]).child_to_parent_cells_ = refined_child_to_parent_cells_vec[level];
        (*data[refinedLevelGridIdx]).cell_to_idxInParentCell_ = refined_cell_to_idxInParentCell_vec[level];
        (*data[refinedLevelGridIdx]).level_to_leaf_cells_ =  refined_level_to_leaf_cells_vec[level];
        (*data[refinedLevelGridIdx]).index_set_ = std::make_unique<cpgrid::IndexSet>(data[refinedLevelGridIdx]->size(0),
                                                                                     data[refinedLevelGridIdx]->size(3));
        (*data[refinedLevelGridIdx]).refinement_max_level_ = levels + preAdaptMaxLevel;
        // Determine the amount of cells per direction, per parent cell, of the corresponding LGR.
        (*data[refinedLevelGridIdx]).cells_per_dim_ = cells_per_dim_vec[level];
        // TO DO: This new code for refinement do not assume Cartesian Shape. How does logical_cartesian_size_ should be defined then?
        // When the refined level grid has been originated from a block of cells, then its logical Cartesian size
        // corresponds to the inner product between cells_per_dim_vec[level] and the dimension of the block (amount of cells in each direction).
        // In the case of a block of cells, e.g., when CARFIN keyword is used, we need the following:
        if (isCARFIN) {
            const auto& blockDim = Opm::Lgr::getPatchDim(startIJK_vec[level], endIJK_vec[level]);
            (*data[refinedLevelGridIdx]).logical_cartesian_size_ = { cells_per_dim_vec[level][0]*blockDim[0],
                                                                     cells_per_dim_vec[level][1]*blockDim[1],
                                                                     cells_per_dim_vec[level][2]*blockDim[2] };
        }
        else {
            (*data[refinedLevelGridIdx]).logical_cartesian_size_ = (*data[0]).logical_cartesian_size_;
            (*data[refinedLevelGridIdx]).global_cell_.swap(refined_global_cell_vec[level]);
        }
        // One alternative definition for logical_cartesian_size_ in the case where the marked elements for refinement do not form a block of cells,
        // therefore, are not associated with the keyword CARFIN, is to imagine that we put all the marked elements one next to the other, along
        // the x-axis. Then, the "imaginary" logical Cartesian size of the refined level grid would be
        // { (# marked elemnts)x cells_per_dim_vec[level][0], cells_per_dim_vec[level][1], cells_per_dim_vec[level][2]}.
        /** To do: how the definition of refined level grids logical_cartesian_size_ affects LookUpData class (and LookUpCartesianData)*/
    }

    // Store adapted grid
    data.push_back(adaptedGrid_ptr);

    // Further Adapted  grid Attributes
    (*data[levels + preAdaptMaxLevel +1]).child_to_parent_cells_ = adapted_child_to_parent_cells;
    (*data[levels + preAdaptMaxLevel +1]).cell_to_idxInParentCell_ = adapted_cell_to_idxInParentCell;
    (*data[levels + preAdaptMaxLevel +1]).leaf_to_level_cells_ =  leaf_to_level_cells;
    (*data[levels + preAdaptMaxLevel +1]).index_set_ = std::make_unique<cpgrid::IndexSet>(data[levels + preAdaptMaxLevel +1]->size(0),
                                                                                          data[levels + preAdaptMaxLevel +1]->size(3));
    (*data[levels + preAdaptMaxLevel +1]).refinement_max_level_ = levels + preAdaptMaxLevel;

    if (isGlobalRefine) {
        assert(cells_per_dim_vec.size() == 1);
        (*data[levels + preAdaptMaxLevel +1]).logical_cartesian_size_ =  { lcs[0]*cells_per_dim_vec[0][0],
                                                                           lcs[1]*cells_per_dim_vec[0][1],
                                                                           lcs[2]*cells_per_dim_vec[0][2] };
    }
    else {
        (*data[levels + preAdaptMaxLevel +1]).logical_cartesian_size_ =  (*data[0]).logical_cartesian_size_;
    }

    // When the refinement is determined by startIJK and endIJK values, the LGR has a (local) Cartesian size.
    // Therefore, each refined cell belonging to the LGR can be associated with a (local) IJK and its (local) Cartesian index.
    // If the LGR has NXxNYxNZ dimension, then the Cartesian indices take values
    // k*NN*NY + j*NX + i, where i<NX, j<Ny, k<NZ.
    // This index is stored in <refined-level-grid>.global_cell_[ refined cell index (~element.index()) ] =  k*NN*NY + j*NX + i.
    if (isCARFIN) {
        for (int level = 0; level < levels; ++level) {
            const int refinedLevelGridIdx = level + preAdaptMaxLevel +1;
            std::vector<int> global_cell_lgr(data[refinedLevelGridIdx]->size(0));
            computeGlobalCellLgr(refinedLevelGridIdx, startIJK_vec[level], global_cell_lgr);
            (*data[refinedLevelGridIdx]).global_cell_.swap(global_cell_lgr);
        }
    }

    std::vector<int> global_cell_leaf( data[levels + preAdaptMaxLevel +1]->size(0));
    computeGlobalCellLeafGridViewWithLgrs(global_cell_leaf);
    (*data[levels + preAdaptMaxLevel +1]).global_cell_.swap(global_cell_leaf);

    updateCornerHistoryLevels(cornerInMarkedElemWithEquivRefinedCorner,
                              elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                              adaptedCorner_to_elemLgrAndElemLgrCorner,
                              corner_count,
                              preAdaptGrid_corner_history,
                              preAdaptMaxLevel,
                              levels);

    // Insert the new id sets into the grid global_id_set_ptr_
    for (int level = 0; level < levels; ++level) {
        const int refinedLevelGridIdx = level + preAdaptMaxLevel +1;
        this->global_id_set_ptr_->insertIdSet(*data[refinedLevelGridIdx]);
    }
    this->global_id_set_ptr_->insertIdSet(*data.back());

    // Only for parallel runs
    // - Define global ids for refined level grids (level 1, 2, ..., maxLevel)
    // - Define GlobalIdMapping (cellMapping, faceMapping, pointMapping required per level)
    // - Define ParallelIndex for overlap cells and their neighbors
    if(comm().size()>1) {
        globalIdsPartitionTypesLgrAndLeafGrids(assignRefinedLevel,
                                               cells_per_dim_vec,
                                               lgr_with_at_least_one_active_cell);
    }

    // Print total amount of cells on the adapted grid
    Opm::OpmLog::info(std::to_string(markedElem_count) + " elements have been marked (in " + std::to_string(comm().rank()) + " rank).\n");
    Opm::OpmLog::info(std::to_string(levels)  + " (new) refined level grid(s) (in " + std::to_string(comm().rank()) + " rank).\n");
    Opm::OpmLog::info(std::to_string(cell_count)  + " total cells on the leaf grid view (in " + std::to_string(comm().rank()) + " rank).\n");

    return (preAdaptMaxLevel < this->maxLevel()); // true if at least one entity was refined
}

void CpGridLGR::postAdapt()
{
    // - Resize with the new amount of cells on the leaf grid view
    // - Set marks equal to zero (representing 'doing nothing')
    current_data_ ->back()-> postAdapt();
}

void CpGridLGR::globalIdsPartitionTypesLgrAndLeafGrids([[maybe_unused]] const std::vector<int>& assignRefinedLevel,
                                                    [[maybe_unused]] const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                                    [[maybe_unused]] const std::vector<int>& lgr_with_at_least_one_active_cell)
{
#if HAVE_MPI
    // Prediction min cell and point global ids per process
    //
    // Predict how many new cells/points (born in refined level grids) need new globalIds, so we can assign unique
    // new ids ( and anticipate the maximum).
    // The grid is already refined according to the LGR specification
    // At this point, neither cell_index_set_ nor partition_type_indicator_ are populated.
    // Refined level grid cells:
    //    1. Inherit their partition type from their parent cell (i.e., element.father().partitionType()).
    //    2. Assign global ids only for interior cells.
    //    3. Communicate the already assigned cell global ids from interior to overlap refined cells.
    // Refined level grid points/vertices:
    // There are 4 partition types: interior, border, front, overlap. This classification requires that both
    // cell and face partition types are already defined, not available yet for refined level grids.
    //    1. Assign for all partition type points a 'candidate of global id' (unique in each process).
    //       Except the points that coincide with a point from level zero.
    //    2. Re-write the values for points that are corners of overlap refined cells, via communication.
    // Under the assumption of LGRs fully-interior, no communication is needed. In the general case, communication will be used
    // to populate overlap cell/point global ids on the refined level grids.
    /** Warning: due to the overlap layer size (equal to 1) cells that share corners or edges (not faces) with interior cells
        are not included/seen by the process. This, in some cases, ends up in multiple ids for the same point. */

    int min_globalId_cell_in_proc = 0;
    int min_globalId_point_in_proc = 0;
    Opm::Lgr::predictMinCellAndPointGlobalIdPerProcess(*this,
                                                       assignRefinedLevel,
                                                       cells_per_dim_vec,
                                                       lgr_with_at_least_one_active_cell,
                                                       min_globalId_cell_in_proc,
                                                       min_globalId_point_in_proc);

    // Only for level 1,2,.., maxLevel grids.
    // For each level, define the local-to-global maps for cells and points (for faces: empty).
    // 1) Assignment of new global ids is done only for owned cells and non-overlap points.
    // 2) For overlap cells and points: communicate. Not needed under the assumption of fully interior LGRs.
    std::vector<std::vector<int>> localToGlobal_cells_per_level(cells_per_dim_vec.size());
    std::vector<std::vector<int>> localToGlobal_points_per_level(cells_per_dim_vec.size());
    // Ignore faces - empty vectors.
    std::vector<std::vector<int>> localToGlobal_faces_per_level(cells_per_dim_vec.size());

    Opm::Lgr::assignCellIdsAndCandidatePointIds(*this,
                                                localToGlobal_cells_per_level,
                                                localToGlobal_points_per_level,
                                                min_globalId_cell_in_proc,
                                                min_globalId_point_in_proc,
                                                cells_per_dim_vec);

    const auto& parent_to_children = current_data_->front()->parent_to_children_cells_;
    ParentToChildrenCellGlobalIdHandle parentToChildrenGlobalId_handle(parent_to_children, localToGlobal_cells_per_level);
    currentData().front()->communicate(parentToChildrenGlobalId_handle,
                                       Dune::InteriorBorder_All_Interface,
                                       Dune::ForwardCommunication );

    // After assigning global IDs to points in refined-level grids, a single point may have
    // a "unique" global ID in each local leaf grid view for every process to which it belongs.
    // To ensure true uniqueness, since global IDs must be distinct across the global leaf view
    // and consistent across each refined-level grid, we will rewrite the entries in
    // localToGlobal_points_per_level.
    //
    // This correction is done using cell_to_point_ across all refined cells through
    // communication: gathering the 8 corner points of each interior cell and scattering the
    // 8 corner points of overlapping cells, for all child cells of a parent cell in level zero grid.
    //
    // Design decision: Why we communicate via level zero grid instead of in each refined level grid.
    // The reason is that how children are stored (the ordering) in parent_to_children_cells_
    // is always the same, accross all processes.
    // Even though the ordering of the corners in cell_to_point_ is the same accross all processes,
    // this may not be enough to correctly overwrite the "winner" point global ids for refined cells.
    //
    /** Current approach avoids duplicated point ids when
     // 1. the LGR is distributed in P_{i_0}, ..., P_{i_n}, with n+1 < comm().size(),
     // AND
     // 2. there is no coarse cell seen by a process P with P != P_{i_j}, j = 0, ..., n.
     // Otherwise, there will be duplicated point ids.
     //
     // Reason: neighboring cells that only share corners (not faces) are NOT considered in the
     // overlap layer of the process.*/
    Opm::Lgr::selectWinnerPointIds(*this,
                                   localToGlobal_points_per_level,
                                   parent_to_children,
                                   cells_per_dim_vec);

    for (std::size_t level = 1; level < cells_per_dim_vec.size()+1; ++level) {
        // Global id set for each (refined) level grid.
        if(lgr_with_at_least_one_active_cell[level-1]>0) { // Check if LGR is active in currect process.
            (*current_data_)[level]->global_id_set_->swap(localToGlobal_cells_per_level[level-1],
                                                          localToGlobal_faces_per_level[level-1],
                                                          localToGlobal_points_per_level[level-1]);
        }
    }

    for (std::size_t level = 1; level < cells_per_dim_vec.size()+1; ++level) {

        populateCellIndexSetRefinedGrid(level);
        // Compute the partition type for cell
        currentData()[level]->computeCellPartitionType();
        // Compute the partition type for point
        currentData()[level]->computePointPartitionType();
        // Now we can compute the communication interface.
        currentData()[level]->computeCommunicationInterfaces(currentData()[level]->size(3));
    }

    populateLeafGlobalIdSet();

    // Insert the new id sets into the grid global_id_set_ptr_
    for (std::size_t level = 0; level < cells_per_dim_vec.size()+1; ++level) {
        this->global_id_set_ptr_->insertIdSet(*(*current_data_)[level]);
    }
    this->global_id_set_ptr_->insertIdSet(*currentData().back());

    populateCellIndexSetLeafGridView();

    // Compute the partition type for cell
    (*current_data_).back()->computeCellPartitionType();

    // Compute the partition type for point
    (*current_data_).back()->computePointPartitionType();

    // Now we can compute the communication interface.
    current_data_->back()->computeCommunicationInterfaces(current_data_->back()->size(3));
    assert(static_cast<std::size_t>(current_data_->back()->cellIndexSet().size()) == static_cast<std::size_t>(current_data_->back()->size(0)) );
#endif
}

void CpGridLGR::getFirstChildGlobalIds([[maybe_unused]] std::vector<int>& parentToFirstChildGlobalIds)
{
#if HAVE_MPI
    switchToGlobalView();
    Opm::Lgr::getFirstChildGlobalIds(*this, parentToFirstChildGlobalIds);
#endif
}

void CpGridLGR::syncDistributedGlobalCellIds()
{
#if HAVE_MPI
    std::vector<int> parentToFirstChildGlobalIds;
    getFirstChildGlobalIds(parentToFirstChildGlobalIds);

    auto size = comm().max(parentToFirstChildGlobalIds.size());

    switchToDistributedView();

    parentToFirstChildGlobalIds.resize(size);
    comm().broadcast(parentToFirstChildGlobalIds.data(), parentToFirstChildGlobalIds.size(), 0);

    const int maxLevel = this->maxLevel();

    // Preallocate syncCellIds (and vertexIds, which will NOT be synchronized)
    std::vector<std::vector<int>> syncCellIds(maxLevel);
    std::vector<std::vector<int>> vertexIds(maxLevel);
    for (int level = 1; level <= maxLevel; ++level) {
        syncCellIds[level-1].resize(currentData()[level]->size(0));
        vertexIds[level-1].resize(currentData()[level]->size(3));
    }

    const auto& globalIdSet = this->globalIdSet();

    // Populate syncCellIds and vertexIds
    for (int level = 1; level <= maxLevel; ++level) {
        const auto& elements = Dune::elements(levelGridView(level));
        for (const auto& element : elements) {
            const int parent_globalId = globalIdSet.id(element.father());
            const int idx_in_parent = element.getIdxInParentCell();
            const int first_child_id = parentToFirstChildGlobalIds[parent_globalId];
            const int new_elem_globalId = first_child_id + idx_in_parent;

            syncCellIds[element.level()-1][element.index()] = new_elem_globalId;
        }

        for (const auto& vertex : Dune::vertices(levelGridView(level))){
            vertexIds[level-1][vertex.index()] = globalIdSet.id(vertex);
        }
    }

    // Re-assign new cell global ids for all refined level grids
    std::vector<int> faceIds; // empty for all
    for (int level = 1; level <= maxLevel; ++level) {
        if(currentData()[level]->size(0)) { // Check if LGR is active in currect process.
            currentData()[level]->global_id_set_->swap(syncCellIds[level-1],
                                                       faceIds,
                                                       vertexIds[level-1]);

            populateCellIndexSetRefinedGrid(level);
            // Insert the new id sets into the grid global_id_set_ptr_
            this->global_id_set_ptr_->insertIdSet(*(*current_data_)[level]);
        }
    }

    populateLeafGlobalIdSet();
    this->global_id_set_ptr_->insertIdSet(*currentData().back());

    assert(static_cast<std::size_t>(current_data_->back()->cellIndexSet().size()) == static_cast<std::size_t>(current_data_->back()->size(0)) );
#endif
}

void CpGridLGR::addLgrsUpdateLeafView(const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                   const std::vector<std::array<int,3>>& startIJK_vec,
                                   const std::vector<std::array<int,3>>& endIJK_vec,
                                   const std::vector<std::string>& lgr_name_vec,
                                   const std::vector<std::string>& lgr_parent_grid_name_vec)
{
    // For parallel run, level zero grid is stored in distributed_data_[0]. If CpGridLGR::scatterGrid has been invoked,
    // then current_data_ == distributed_data_.
    // For serial run, level zero grid is stored in data_[0]. In this case, current_data_ == data_.
    // Note: currentData() returns data_ (if grid is not distributed) or distributed_data_ otherwise.

    Opm::Lgr::validStartEndIJKs(startIJK_vec, endIJK_vec);

    // If no parent grid name vector has been provided, then default "GLOBAL" for all (new) level grids.
    std::vector<std::string> parent_grid_names = lgr_parent_grid_name_vec;
    if (parent_grid_names.size() == 0){ // No parent grid name given->default "GLOBAL" parent grid
        parent_grid_names.resize(cells_per_dim_vec.size(), "GLOBAL");
    }

    // Sizes of provided vectors (number of subivisions per cells and lgrs name) should coincide.
    bool matchingSizeHasFailed = false;
    if ( (cells_per_dim_vec.size() != startIJK_vec.size()) ||
         (lgr_name_vec.size() != startIJK_vec.size()) ||
         (parent_grid_names.size() != startIJK_vec.size())) {
        matchingSizeHasFailed = true;
    }
    matchingSizeHasFailed = comm().max(matchingSizeHasFailed);
    if (matchingSizeHasFailed) {
        OPM_THROW(std::invalid_argument, "Sizes of provided vectors with subdivisions per cell and LGR names need to match.");
    }

    // Discard LGRs whose subdivisions do not trigger actual refinement, i.e., cells_per_dim_ = {1,1,1}
    const auto [filtered_cells_per_dim_vec,
                filtered_startIJK_vec,
                filtered_endIJK_vec,
                filtered_lgr_name_vec,
                filtered_lgr_parent_grid_name_vec] = Opm::Lgr::excludeFakeSubdivisions(cells_per_dim_vec,
                                                                                       startIJK_vec,
                                                                                       endIJK_vec,
                                                                                       lgr_name_vec,
                                                                                       parent_grid_names);
    if (filtered_cells_per_dim_vec.size() == 0) { // if all LGRs expect 1 child per direction, then no refinement will be done.
        return;
    }

    if (!Opm::areParentGridsAvailableBeforeTheirLgrs(getLgrNameToLevel(),
                                                     filtered_lgr_name_vec,
                                                     filtered_lgr_parent_grid_name_vec)) {
        OPM_THROW(std::invalid_argument, "Parent grid (name) must exist before its LGRs.");
    }

    // Refinement proceeds in steps. In each step, for every parent grid:
    //   1. Gather the data of its child LGRs.
    //   2. Create the corresponding LGRs.
    //   3. Update the leaf grid view.
    //
    // To achieve this, we first collect the set of unique parent grid names
    // (avoiding duplicates).
    std::set<std::string> non_repeated_parent_grid_names(filtered_lgr_parent_grid_name_vec.begin(),
                                                         filtered_lgr_parent_grid_name_vec.end());
    int tmp_maxLevel = this->maxLevel();

    for (const auto& parent_grid_name : non_repeated_parent_grid_names) {
        //   1. Gather the data of its child LGRs.
        auto [cells_per_dim_vec_parent_grid,
              startIJK_vec_parent_grid,
              endIJK_vec_parent_grid,
              lgr_name_vec_parent_grid] =  Opm::filterLgrDataPerParentGridName(filtered_cells_per_dim_vec,
                                                                               filtered_startIJK_vec,
                                                                               filtered_endIJK_vec,
                                                                               filtered_lgr_name_vec,
                                                                               filtered_lgr_parent_grid_name_vec,
                                                                               parent_grid_name);

        int parent_grid_index = getLgrNameToLevel().at(parent_grid_name);
        // Determine the assigned level for the refinement of each marked cell
        std::vector<int> assignRefinedLevel(currentLeafData().size(0));

        // Compatibility of numbers of subdivisions of neighboring LGRs".
        // The method compatibleSubdivision returns a bool. We convert it into an int since MPI within DUNE does not support bool directly.
        int compatibleSubdivisions = Opm::Lgr::compatibleSubdivisions(filtered_cells_per_dim_vec,
                                                                      filtered_startIJK_vec,
                                                                      filtered_endIJK_vec,
                                                                      currentData()[parent_grid_index]->logicalCartesianSize());
        compatibleSubdivisions = comm().min(compatibleSubdivisions); // 0 when at least one process returns false (un-compatible subdivisions).
        if(!compatibleSubdivisions) {
            if (comm().rank()==0){
                OPM_THROW(std::logic_error, "Subdivisions of neighboring LGRs sharing at least one face do not coincide. Not suppported yet.");
            }
            else{
                OPM_THROW_NOLOG(std::logic_error, "Subdivisions of neighboring LGRs sharing at least one face do not coincide. Not suppported yet.");
            }
        }

        // To determine if an LGR is not empty in a given process, for each
        // parent grid, we set at_least_one_active_parent[in that level] to 1
        // if it contains at least one active cell, and to 0 otherwise.
        std::vector<int> at_least_one_active_parent(startIJK_vec_parent_grid.size());

        // Find out which (ACTIVE) elements belong to the block cells defined by startIJK and endIJK values.
        for(const auto& element: elements(this->leafGridView())) {
            std::array<int,3> ijk;
            // Skip the element if its level is not equal to parent_grid_index
            // Note: Inconsistency with DUNE Grid interface. element.level()
            //       returns the index to access the level grid where the entity was born.
            //       This means that, in general,
            //       elment.level() != element.father().level() + 1/
            //       It can happen that | element.level() - element.father().level()| >1.
            if (parent_grid_index != element.level())
                continue;
            currentData()[element.level()]->getIJK(element.getLevelElem().index(), ijk);

            for (std::size_t level = 0; level < startIJK_vec_parent_grid.size(); ++level) {
                bool belongsToLevel = true;
                for (int c = 0; c < 3; ++c) {
                    belongsToLevel = belongsToLevel && ( (ijk[c] >= startIJK_vec_parent_grid[level][c]) && (ijk[c] < endIJK_vec_parent_grid[level][c]) );
                    if (!belongsToLevel)
                        break;
                }
                if(belongsToLevel) {
                    this->mark(1, element, /* throwOnFailure = */ true);
                    at_least_one_active_parent[level] = 1;
                    assignRefinedLevel[element.index()] = tmp_maxLevel + level +1;
                }
            }
        }
        tmp_maxLevel = tmp_maxLevel + startIJK_vec_parent_grid.size(); // Update the maxLevel

        //   2. Create the corresponding LGRs. and  3. Update the leaf grid view.
        refineAndUpdateGrid(cells_per_dim_vec_parent_grid,
                            assignRefinedLevel,
                            lgr_name_vec_parent_grid,
                            startIJK_vec_parent_grid,
                            endIJK_vec_parent_grid);

        int non_empty_lgrs = 0;
        for (std::size_t level = 0; level < startIJK_vec_parent_grid.size(); ++level){
            // Do not throw if all cells of an LGR are inactive in a parallel run (The process might not 'see' those cells.)
            if (at_least_one_active_parent[level]) {
                ++non_empty_lgrs;
            }
            if ((comm().max(at_least_one_active_parent[level]) == 0) && (comm().rank() == 0)) {
                Opm::OpmLog::warning(lgr_name_vec_parent_grid[level]+ " contains only inactive cells (in rank " + std::to_string(comm().rank()) +").\n");
            }
        }

        // Notice that in a parallel run, non_empty_lgrs represents the local active lgrs, i.e. the lgrs containing active cells which also belong
        // to the current process. Considered per parent grid.
        auto globalActiveLgrs_currentParentGrid = comm().sum(non_empty_lgrs);
        if((globalActiveLgrs_currentParentGrid == 0) && (comm().rank() == 0)) {
            Opm::OpmLog::warning("All the LGRs with parent grid " + parent_grid_name + " contain only inactive cells.\n");
        }
    }
}

void CpGridLGR::autoRefine(const std::array<int,3>& nxnynz)
{
    // Refinement factors must be odd and positive.
    for (const auto& nd : nxnynz) {
        if (nd<=0 || (nd%2==0)) {
            OPM_THROW(std::invalid_argument, "Refinement factor must be odd and positive.\n");
        }
    }
    const auto endIJK = this->logicalCartesianSize();
    addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {nxnynz},
                          /* startIJK_vec = */ {{0,0,0}},
                          /* endIJK_vec = */ {endIJK},
                          /* lgr_name_vec = */ {"GLOBAL_REFINED"});
}

const std::map<std::string,int>& CpGridLGR::getLgrNameToLevel() const{
    return lgr_names_;
}

void CpGridLGR::updateCornerHistoryLevels(const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                       const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                       const std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner,
                                       const int& corner_count,
                                       const std::vector<std::array<int,2>>& preAdaptGrid_corner_history,
                                       const int& preAdaptMaxLevel,
                                       const int& newLevels)
{
    for (int level = preAdaptMaxLevel+1; level < preAdaptMaxLevel + newLevels+1; ++level) {
        currentData()[level]->corner_history_.resize( currentData()[level] ->size(3), std::array<int,2>({-1,-1}));
    }
    // corner_history_ for levels 0, level 1, ..., preAdapt-maxLevel (maximum level before calling (again) adapt) should be already populated
    // corner_history_[ new corner ] = {-1,-1}
    // corner_history_[ corner equivalent to a corner in a previous level ] = { level where the corner was born, its index in that level grid}.
    for (std::size_t corner = 0; corner < cornerInMarkedElemWithEquivRefinedCorner.size(); ++corner) {
        if (!cornerInMarkedElemWithEquivRefinedCorner[corner].empty()) {
            const auto& [refinedLevel, refinedCorner] = elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.at(cornerInMarkedElemWithEquivRefinedCorner[corner].back());
            currentData()[refinedLevel]->corner_history_[refinedCorner] = preAdaptGrid_corner_history.empty() ? std::array<int,2>{{0, static_cast<int>(corner)}} :  preAdaptGrid_corner_history[corner];
        }
    }

    // corner_history_ leaf grid view
    for ( int leafCorner = 0; leafCorner < corner_count; ++leafCorner){
        currentData().back()->corner_history_.resize(corner_count);
        const auto& [elemLgr, elemLgrCorner] = adaptedCorner_to_elemLgrAndElemLgrCorner.at(leafCorner);
        if (elemLgr != -1) {
            const auto& [refinedLevel, refinedCorner] = elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.at({elemLgr, elemLgrCorner});
            currentData().back()->corner_history_[leafCorner] = { refinedLevel, refinedCorner};
        }
        else {
            currentData().back()->corner_history_[leafCorner] =  preAdaptGrid_corner_history.empty() ? std::array<int,2>{{0, elemLgrCorner}} : preAdaptGrid_corner_history[elemLgrCorner];
        }
    }
}

} // namespace Dune
