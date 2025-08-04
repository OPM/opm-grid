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

#ifndef OPM_GRID_CPGRID_LGRHELPERS_HEADER_INCLUDED
#define OPM_GRID_CPGRID_LGRHELPERS_HEADER_INCLUDED


#include <dune/grid/common/mcmgmapper.hh>
#include <dune/common/version.hh>


#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>


#include <algorithm>
#include <array>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>


namespace Opm
{


/// --------------- Auxiliary methods to support Adaptivity (begin) ---------------

/// @brief Refine each marked element and establish relationships between corners, faces, and cells marked for refinement,
///        with the refined corners, refined faces, and refined cells.
///
/// --- Marked elements parameters ---
/// @param [out] markedElem_to_itsLgr:                      Each marked element gets refined and we store this "auxiliary markedElementLGR", to later
///                                                         build each refined level grid containing all the refined entities assigned for that grid.
/// @param [out] markedElem_count:                          Total amount of marked elements to be refined. It will be used to print grid info.
/// @param [out] cornerInMarkedElemWithEquivRefinedCorner:  For each corner from level zero, we store the marked elements where the corner appears and its equivalent
///                                                         refined corner in  each auxiliary marked-element-lgr. Example: corner with index 5 appears in marked
///                                                         elements 0 and 1, with refined equivalent corner indices 8 and 2 respectively. Then,
///                                                         cornerInMarkedElemWithEquivRefinedCorner[5] = {{0, 8}, {1, 2}}.
///                                                         For corners not appearing in any marked element, empty vector.
/// @param [out] markedElemAndEquivRefinedCorner_to_corner: To correctly build the level-refined and adapted-grid topology features, we need to keep track of the
///                                                         corners that got replaced by equivalent refined corners, in each marked element where the corner appeared,
///                                                         not only in its last appearance. The last appearance will be used to avoid repetition when storing.
///                                                         Following the example above,
///                                                         markedElemAndEquivRefinedCorner_to_corner[{0, 8}] = 5;
///                                                         markedElemAndEquivRefinedCorner_to_corner[{1, 2}] = 5;
/// @param [out] faceInMarkedElemAndRefinedFaces:           For each face from level zero, we store the marked elements where the face appears (maximum 2 cells)
///                                                         and its new-born refined faces from each auxiliary marked-element-lgr. Example: face with index 9
///                                                         appears in marked elements 0 and 1. Then,
///                                                         faceInMarkedElemAndRefinedFaces[9] = {{0, {refinedFace0_0, ..., refinedFaceN_0}},
///                                                                                               {1, {refinedFace0_1, ..., refinedFaceM_1}}}.
///                                                         For faces not appearing in any marked element, empty vector.
/// --- Refined cells parameters ---
/// @param [out] elemLgrAndElemLgrCell_to_refinedLevelAdRefinedCell:  Each marked element has been refined in its "own elemLgr". Refined entities should be stored in
///                                                                   the corresponding assigned refined level grid. To keep track of the cell index relation,
///                                                                   associate each
///                                                                   { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement } with
///                                                                   { refined level grid assigned for the marked element, refined cell index in refined level grid }.
/// @param [out] refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell: Each marked element has been assigned to certain refined level grid. To keep track of the "inverse"
///                                                                   cell index relation, associate each
///                                                                   { refined level grid assigned for the marked element, refined cell index in refined level grid }
///                                                                   with { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement }.
/// @param [out] refined_cell_count_vec:                              Total amount of refined cells, per level (i.e. in each refined level grid).
/// @param [in] assignRefinedLevel:                                   Each marked element can be assigned to certain refined level grid. This vector has entries 0 for
///                                                                   non marked elements, and the corresponding integer representing a refined level grid for marked
///                                                                   elements.
/// @param [out] preAdapt_parent_to_children_cells_vec:               Parent cells and their refined children. Entry is {-1, {}} when cell has no children. Othewise,
///                                                                   {refined grid level where children were born, {child0, child1, ...}}
///                                                                   Each vector entry represents an existing level grid before calling adapt.
/// --- Adapted cells parameters ---
/// @param [out] elemLgrAndElemLgrCell_to_adaptedCell:                Each marked element has been refined in its "own elemLgr". Refined entities should be also stored in
///                                                                   the corresponding leaf grid view (or adapted grid). To keep track of the cell index relation,
///                                                                   associate each
///                                                                   { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement } with
///                                                                   refined cell index inthe leaf grid view (or adapted grid).
/// @param [out] adaptedCell_to_elemLgrAndElemLgrCell:                Each marked element has been refined in its "own elemLgr". Refined entities should be also stored in
///                                                                   the corresponding leaf grid view (or adapted grid). To keep track of the "inverse" cell index
///                                                                   relation, associate the refined cell index inthe leaf grid view (or adapted grid) with
///                                                                   { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement }.
/// @param [out] cell_count:                                          Total amount of cells on the leaf grid view (or adapted grid).
/// @param [out] preAdapt_level_to_leaf_cells_vec:                    For each existing grid before calling adapt, we stablish the index relation between preAdapt cells
///                                                                   and cells on the leaf grid view (or adapted cells).-1 means that the cell vanished.
/// --- Additional parameters ---
/// @param [in] cells_per_dim_vec:                                    For each set of marked elements for refinement, that will belong to a same
///                                                                   refined level grid, number of (refined) cells in each direction that each
///                                                                   parent cell should be refined to.
void refineAndProvideMarkedRefinedRelations(const Dune::CpGrid& grid,/* Marked elements parameters */
                                            std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                            int& markedElem_count,
                                            std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                            std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                                            std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                            /* Refined cells parameters */
                                            std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCell_to_refinedLevelAdRefinedCell,
                                            std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                            std::vector<int>& refined_cell_count_vec,
                                            const std::vector<int>& assignRefinedLevel,
                                            std::vector<std::vector<std::tuple<int,std::vector<int>>>>& preAdapt_parent_to_children_cells_vec,
                                            /* Adapted cells parameters */
                                            std::map<std::array<int,2>,int>& elemLgrAndElemLgrCell_to_adaptedCell,
                                            std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                                            int& cell_count,
                                            std::vector<std::vector<int>>& preAdapt_level_to_leaf_cells_vec,
                                            /* Additional parameters */
                                            const std::vector<std::array<int,3>>& cells_per_dim_vec);

/// @brief  Define child-parent relations from the new refined cells of the new refined level grids to its parent cells (belonging to pre-existing grid,
///         before adapting the grid/before updating the leaf grid view). Define the index in parent cell (-1 when cell has no parent).
///
/// @param [in] refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell: Each marked element has been assigned to certain refined level grid. To keep track of the "inverse"
///                                                                  cell index relation, associate each
///                                                                  { refined level grid assigned for the marked element, refined cell index in refined level grid }
///                                                                  with { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement }.
/// @param [in] refined_cell_count_vec:                              Total amount of refined cells, per level (i.e. in each refined level grid).
/// @param [in] adaptedCell_to_elemLgrAndElemLgrCell:                Each marked element has been refined in its "own elemLgr". Refined entities should be also stored in
///                                                                  the corresponding leaf grid view (or adapted grid). To keep track of the "inverse" cell index
///                                                                  relation, associate the refined cell index inthe leaf grid view (or adapted grid) with
///                                                                  { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement }.
/// @param [in] cell_count:                                          Total amount of cells on the leaf grid view (or adapted grid).
///
/// @return refined_child_to_parent_cells_vec:   Refined child cells and their parents. Entry is {-1,-1} when cell has no father. Otherwise,
///                                              {level parent cell, parent cell index}. Each vector entry represents a refined level grid.
///         refined_cell_to_idxInParentCell_vec: Each refined child cell has a unique index in its parent cell, to be used to build geometryInFather().
///                                              Each vector entry represents a refined level grid.
///         adapted_child_to_parent_cell:        Refined child cells and their parents. Entry is {-1,-1} when cell has no father. Otherwise,
///                                              {level parent cell, parent cell index}
///         adapted_cell_to_idxInParentCell:     Each refined child cell has a unique index in its parent cell, to be used to build geometryInFather(). -1 when has no father.
std::tuple< std::vector<std::vector<std::array<int,2>>>,
            std::vector<std::vector<int>>,
            std::vector<std::array<int,2>>,
            std::vector<int>>
defineChildToParentAndIdxInParentCell( const Dune::CpGrid& grid,
                                       const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                       const std::vector<int>& refined_cell_count_vec,
                                       const std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                                       const int& cell_count);

/// @brief Define refined level grid cells indices and leaf grid view (or adapted grid) cells indices relations. Namely, level_to_leaf_cells_ for each new
///        refined level grid, and leaf_to_level_cells_ for the updated leaf grid view.
///
/// @param [in] elemLgrAndElemLgrCell_to_refinedLevelAdRefinedCell:  Each marked element has been refined in its "own elemLgr". Refined entities should be stored in
///                                                                  the corresponding assigned refined level grid. To keep track of the cell index relation, we
///                                                                  associate each
///                                                                  { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement } with
///                                                                  { refined level grid assigned for the marked element, refined cell index in refined level grid }.
/// @param [in] refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell: Each marked element has been assigned to certain refined level grid. To keep track of the "inverse"
///                                                                  cell index relation, associate each
///                                                                  { refined level grid assigned for the marked element, refined cell index in refined level grid }
///                                                                  with { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement }.
/// @param [in] refined_cell_count_vec:                              Total amount of refined cells, per level (i.e. in each refined level grid).
/// @param [in] elemLgrAndElemLgrCell_to_adaptedCell:                Each marked element has been refined in its "own elemLgr". Refined entities should be also stored in
///                                                                  the corresponding leaf grid view (or adapted grid). To keep track of the cell index relation,
///                                                                  associate each
///                                                                  { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement } with
///                                                                  refined cell index inthe leaf grid view (or adapted grid).
/// @param [in] adaptedCell_to_elemLgrAndElemLgrCell:                Each marked element has been refined in its "own elemLgr". Refined entities should be also stored in
///                                                                  the corresponding leaf grid view (or adapted grid). To keep track of the "inverse" cell index
///                                                                  relation, associate the refined cell index inthe leaf grid view (or adapted grid) with
///                                                                  { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement }.
/// @param [in] cell_count:                                          Total amount of cells on the leaf grid view (or adapted grid).
///
/// @return refined_level_to_leaf_cells_vec:                         refined_level_to_leaf_cells_vec[ levelGridIdx ] [ cell idx in that level grid ] = equivalent leaf cell idx
///         leaf_to_level_cells:                                     leaf_to_level_cells[ leaf cell idx ] = {level where cell was born, cell idx on that level}
std::pair<std::vector<std::vector<int>>, std::vector<std::array<int,2>>>
defineLevelToLeafAndLeafToLevelCells(const Dune::CpGrid& grid,
                                     const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell,
                                     const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                     const std::vector<int>& refined_cell_count_vec,
                                     const std::map<std::array<int,2>,int>& elemLgrAndElemLgrCell_to_adaptedCell,
                                     const std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                                     const int& cell_count);


}

#endif // OPM_GRID_CPGRID_LGRHELPERS_HEADER_INCLUDED
