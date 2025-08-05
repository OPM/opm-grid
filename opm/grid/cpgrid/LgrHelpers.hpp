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
/// --------------- Auxiliary methods to support refinement ---------------

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

/// @brief Define various corner relations. 1. refined corners from auxiliary single marked element refinement to its corresponding refined level grid, and vice versa.
///                                         2. refined corners from single-cell-refinements that vanish in the "storing only once each entity process". To avoid repetition,
///                                            we store such corners in their "last apperance". We keep track of all the appearances since that is needed for correctly
///                                            define CpGridData attributes such as cell_to_point_ and face_to_point_.
///
/// @param [out] elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner: Each marked element has been refined in its "own elemLgr". Refined corners should be stored in
///                                                                       the corresponding assigned refined level grid. To keep track of the corner index relation, we
///                                                                       associate each
///                                                                       { marked element index ("elemLgr"), refined corner index in the auxiliary single-cell-refinement } with
///                                                                       { refined level grid assigned for the marked element, refined corner index in refined level grid }.
/// @param [out] refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner: Each marked element has been refined in its "own elemLgr". Refined corners should be stored in
///                                                                       the corresponding assigned refined level grid. To keep track of the corner index relation, we
///                                                                       associate each
///                                                                       { refined level grid assigned for the marked element, refined corner index in refined level grid } with
///                                                                       { marked element index ("elemLgr"), refined corner index in the auxiliary single-cell-refinement }.
/// @param [out] refined_corner_count_vec:                                Total amount of refined corners, per level (each vector entry corresponds to a refined level grid).
/// @param [out] vanishedRefinedCorner_to_itsLastAppearance:              A refined corner might appear in several single-cell-refinements, we store it only in its last
///                                                                       appearance, but keep track of the vanishing. Example, a corner appears in total 3
///                                                                       single-cell-refinements, with indices { elemLgr1, elemLgr1Corner }, { elemLgr2, elemLgr2Corner },
///                                                                       and { elemLgr3, elemLgr3Corner }. Then, for X = 1, and X=2, we store
///                                                                       vanishedRefinedCorner_to_itsLastAppearance[{elemLgrX, elemLgrXCorner}] = {elemLgr3, elemLgr3Corner}.
/// @param [in] markedElem_to_itsLgr
/// @param [in] assignRefinedLevel
/// @param [in] cornerInMarkedElemWithEquivRefinedCorner
/// @param [in] faceInMarkedElemAndRefinedFaces
/// @param [in] cells_per_dim_vec
void identifyRefinedCornersPerLevel(const Dune::CpGrid& grid,
                                    std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                    std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner,
                                    std::vector<int>& refined_corner_count_vec,
                                    std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                    const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                    const std::vector<int>& assignRefinedLevel,
                                    const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                    const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                    const std::vector<std::array<int,3>>& cells_per_dim_vec);

/// @brief Determine if a refined corner is located in the interior of the single-cell-refinement.
///
/// @param [in] cells_per_dim:    Total children cells in each direction (x-,y-, and z-direction) of the single-cell-refinement.
/// @param [in] cornerIdxInLgr:   Corner index in the single-cell-refinement.
bool isRefinedCornerInInteriorLgr(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr);

/// @brief Get the ijk index of a refined corner, given its corner index of a single-cell-refinement.
///
/// Given a single-cell, we refine it in {nx, ny, nz} refined children cells (per direction). Then, this single-cell-refinement
/// has in total (nx +1)(ny +1)(nz +1) refined corners. Each of this corners has an ijk value associated since they are stored
/// (following order defined in Geometry::refine) with the index (j*(nx+1)(nz+1)) + (i(nz+1)) + k, where i=0,...,nx, j=0,...,ny,
/// and k=0,...,nz. This method returns the ijk, given the cornerIdxInLgr = 0,...,(nx +1)(ny +1)(nz +1).
///
/// @param [in] cells_per_dim:    Total children cells in each direction (x-,y-, and z-direction) of the single-cell-refinement.
/// @param [in] cornerIdxInLgr:   Corner index in the single-cell-refinement.
std::array<int,3> getRefinedCornerIJK(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr);

/// @brief Determine if a new refined corner is located on an edge of the parent cell. In particular, it's on the boundary of
///        the single-cell-refinement, and does not coincide with  a preAdapt-existing corner.
///
/// @param [in] cells_per_dim:    Total children cells in each direction (x-,y-, and z-direction) of the single-cell-refinement.
/// @param [in] cornerIdxInLgr:   Corner index in the single-cell-refinement.
bool newRefinedCornerLiesOnEdge(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr);

/// @brief Get the parent faces that containes the edge where the new refined corner lies on.
///
/// @param [in] cells_per_dim:    Total children cells in each direction (x-,y-, and z-direction) of the single-cell-refinement.
/// @param [in] cornerIdxInLgr:   Corner index in the single-cell-refinement.
/// @param [in] elemLgr:          Cell index from starting grid, that has been refined into a single-cell-refinement.
std::array<int,2> getParentFacesAssocWithNewRefinedCornLyingOnEdge(const Dune::CpGrid& grid,
                                                                   const std::array<int,3>& cells_per_dim,
                                                                   int cornerIdxInLgr,
                                                                   int elemLgr);


/// @brief Determine if a refined corner is located on the boundary of the single-cell-refinement, and does not coincide with
///        a preAdapt-existing corner.
///
/// @param [in] cells_per_dim:    Total children cells in each direction (x-,-, and z-direction) of the single-cell-refinement.
/// @param [in] cornerIdxInLgr:   Corner index in the single-cell-refinement.
bool isRefinedNewBornCornerOnLgrBoundary(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr);

/// @brief Get the parent face where the new refined corner lays on.
///
/// @param [in] cells_per_dim:    Total children cells in each direction (x-,y-, and z-direction) of the single-cell-refinement.
/// @param [in] cornerIdxInLgr:   Corner index in the single-cell-refinement.
/// @param [in] elemLgr:          Cell index from starting grid, that has been refined into a single-cell-refinement.
int getParentFaceWhereNewRefinedCornerLiesOn(const Dune::CpGrid& grid,
                                             const std::array<int,3>& cells_per_dim,
                                             int cornerIdxInLgr,
                                             int elemLgr);

/// @brief A refined corner appears in two single-cell-refinements. Given the corner index in the first single-cell-refinement, compute the
///         corner index in the neighboring single-cell-refinement.
///
/// @param [in] cells_per_dim_lgr1:    Total children cells in each direction (x-,y-, and z-direction) of the elemLgr1 single-cell-refinement.
/// @param [in] cornerIdxInLgr1:       Corner index in the elemLgr1 single-cell-refinement.
/// @param [in] cells_per_dim_lgr2:    Total children cells in each direction (x-,y-, and z-direction) of the elemLgr2 single-cell-refinement.
int replaceLgr1CornerIdxByLgr2CornerIdx(const std::array<int,3>& cells_per_dim_lgr1,
                                        int cornerIdxLgr1,
                                        const std::array<int,3>& cells_per_dim_lgr2);

/// @brief A new refined corner lays on an edge and appears in at least two single-cell-refinements. Given the corner index in one single-cell-refinement, compute the
///        corner index in a neighboring single-cell-refinement.
///
/// @param [in] cells_per_dim_lgr1:    Total children cells in each direction (x-,y-, and z-direction) of the elemLgr1 single-cell-refinement.
/// @param [in] cornerIdxInLgr1:       Corner index in the elemLgr1 single-cell-refinement.
/// @param [in] parentFaceLastAppearanceIdx: Parent face index where the refined corner appears for last time.
/// @param [in] cells_per_dim_lgr2:    Total children cells in each direction (x-,y-, and z-direction) of the elemLgr2 single-cell-refinement.
int replaceLgr1CornerIdxByLgr2CornerIdx(const Dune::CpGrid& grid,
                                        const std::array<int,3>& cells_per_dim_lgr1, int cornerIdxLgr1, int elemLgr1, int parentFaceLastAppearanceIdx,
                                        const std::array<int,3>& cells_per_dim_lgr2);

/// @brief Identify corners that appear on the leaf grid view.
///        Define various corner relations. preAdapt or refined corners from auxiliary single marked element refinement to the leaf grid view (or adapted grid), and vice versa.
///
/// @param [out] elemLgrAndElemLgrCorner_to_adaptedCorner: Each marked element has been refined in its "own elemLgr". When the element has not been refined, elemLgr == -1.
///                                                        To keep track of the corner index relation, we associate each
///                                                        { marked element index ("elemLgr"),  corner index in the auxiliary single-cell-refinement }, or
///                                                        { -1 ("elemLgr"),  corner index in the starting grid}, with
///                                                        corner index in the leaf grid view .
/// @param [out] adaptedCorner_to_elemLgrAndElemLgrCorner: Each marked element has been refined in its "own elemLgr". Refined corners should be stored in
///                                                        the corresponding assigned refined level grid. To keep track of the cell index relation, we
///                                                        associate each corner index in the leaf grid view (or adapted grid) with
///                                                        { marked element index ("elemLgr"), refined corner index in the auxiliary single-cell-refinement }, or
///                                                        { -1 ("elemLgr"),  corner index in the starting grid},
/// @param [out] corner_count:                             Total amount of corners on the leaf grid view (or adapted grid).
/// @param [in] markedElem_to_itsLgr
/// @param [in] assignRefinedLevel
/// @param [in] cornerInMarkedElemWithEquivRefinedCorner
/// @param [in] vanishedRefinedCorner_to_itsLastAppearance
/// @param [in] faceInMarkedElemAndRefinedFaces
/// @param [in] cells_per_dim_vec
void identifyLeafGridCorners(const Dune::CpGrid& grid,
                             std::map<std::array<int,2>,int>& elemLgrAndElemLgrCorner_to_adaptedCorner,
                             std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner,
                             int& corner_count,
                             const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                             const std::vector<int>& assignRefinedLevel,
                             const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                             std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                             const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                             const std::vector<std::array<int,3>>& cells_per_dim_vec);

/// @brief Define relations between single-cell-refinement faces and refined level faces.
///
/// @param [out] elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace: Each marked element has been refined in its "own elemLgr". Refined faces should be stored in
///                                                                   the corresponding assigned refined level grid. To keep track of the face index relation, we
///                                                                   associate each
///                                                                   { marked element index ("elemLgr"), refined face index in the auxiliary single-cell-refinement } with
///                                                                   { refined level grid assigned for the marked element, refined face index in refined level grid }.
/// @param [out] refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace: Each marked element has been refined in its "own elemLgr". Refined faces should be stored in
///                                                                   the corresponding assigned refined level grid. To keep track of the face index relation, we
///                                                                   associate each
///                                                                   { refined level grid assigned for the marked element, refined face index in refined level grid } with
///                                                                   { marked element index ("elemLgr"), refined face index in the auxiliary single-cell-refinement }.
/// @param [out] refined_face_count_vec:                              Total amount of refined corners, per level (each vector entry corresponds to a refined level grid).
/// @param [in] markedElem_to_itsLgr
/// @param [in] assignRefinedLevel
/// @param [in] faceInMarkedElemAndRefinedFaces
/// @param [in] cells_per_dim_vec
void identifyRefinedFacesPerLevel(const Dune::CpGrid& grid,
                                  std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                                  std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                                  std::vector<int>& refined_face_count_vec,
                                  const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                  const std::vector<int>& assignRefinedLevel,
                                  const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                  const std::vector<std::array<int,3>>& cells_per_dim_vec);

/// @brief Identify faces that appear on the leaf grid view.
///        Define various face relations. preAdapt or refined faces from auxiliary single marked element refinement to the leaf grid view (or adapted grid), and vice versa.
///
/// @param [out] elemLgrAndElemLgrFace_to_adaptedFace: Each marked element has been refined in its "own elemLgr". When the element has not been refined, elemLgr == -1.
///                                                    To keep track of the face index relation, we associate each
///                                                    { marked element index ("elemLgr"),  face index in the auxiliary single-cell-refinement }, or
///                                                    { -1 ("elemLgr"),  face index in the starting grid }, with
///                                                    face index in the leaf grid view .
/// @param [out] adaptedFace_to_elemLgrAndElemLgrFace: Each marked element has been refined in its "own elemLgr". Refined corners should be stored in
///                                                    the corresponding assigned refined level grid. To keep track of the face index relation, we associate each
///                                                    face index in the leaf grid view (or adapted grid) with
///                                                    { marked element index ("elemLgr"), refined corner index in the auxiliary single-cell-refinement }.
/// @param [out] face_count:                           Total amount of faces on the leaf grid view (or adapted grid).
/// @param [in] markedElem_to_itsLgr
/// @param [in] assignRefinedLevel
/// @param [in] faceInMarkedElemAndRefinedFaces
/// @param [in] cells_per_dim_vec
void identifyLeafGridFaces(const Dune::CpGrid& grid,
                           std::map<std::array<int,2>,int>& elemLgrAndElemLgrFace_to_adaptedFace,
                           std::unordered_map<int,std::array<int,2>>& adaptedFace_to_elemLgrAndElemLgrFace,
                           int& face_count,
                           const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                           const std::vector<int>& assignRefinedLevel,
                           const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                           const std::vector<std::array<int,3>>& cells_per_dim_vec);

/// @brief Get the ijk index of a reined face, given its corner index of a single-cell-refinement.
///
/// Given a single-cell, we refine it in {nx, ny, nz} refined children cells (per direction). Then, this single-cell-refinement
/// has in total ((nx+1)*ny*nz) + (nx*(ny+1)*nz) + (nx*ny*(nz+1)) refined faces. Each of this faces has an ijk value associated since
/// they are stored  with the index (following order defined in Geometry::refine):
/// K_FACES  (k*nx*ny) + (j*nx) + i
/// I_FACES  (nx*ny*(nz+1)) + (i*ny*nz) + (k*ny) + j
/// J_FACES   (nx*ny*(nz+1)) + ((nx+1)*ny*nz) + (j*nx*nz) + (i*nz) + k
///  where i=0,...,nx-1, j=0,...,ny-1, and k=0,...,nz-1. This method returns the corresponding ijk
/// given a faceIdxInLgr = 0,...,((nx+1)nynz) + (nx*(ny+1)nz) + (nxny*(nz+1))
///
/// @param [in] cells_per_dim:    Total children cells in each direction (x-,y-, and z-direction) of the single-cell-refinement.
/// @param [in] faceIdxInLgr:     Face index in the single-cell-refinement.
/// @param [in] elemLgr_ptr:      Pointer to the single-cell-refinement grid.
std::array<int,3> getRefinedFaceIJK(const std::array<int,3>& cells_per_dim, int faceIdxInLgr,
                                    const std::shared_ptr<Dune::cpgrid::CpGridData>& elemLgr_ptr);

/// @brief Determine if a refined face is located in the interior of the single-cell-refinement.
///
/// @param [in] cells_per_dim:    Total children cells in each direction (x-,y-, and z-direction) of the single-cell-refinement.
/// @param [in] faceIdxInLgr:     Face index in the single-cell-refinement.
/// @param [in] elemLgr_ptr:      Pointer to the single-cell-refinement grid.
bool isRefinedFaceInInteriorLgr(const std::array<int,3>& cells_per_dim, int faceIdxInLgr,
                                const std::shared_ptr<Dune::cpgrid::CpGridData>& elemLgr_ptr);

/// @brief Determine if a refined face is located on the boundary of the single-cell-refinement.
///
/// @param [in] cells_per_dim:    Total children cells in each direction (x-,y-, and z-direction) of the single-cell-refinement.
/// @param [in] faceIdxInLgr:     Face index in the single-cell-refinement.
/// @param [in] elemLgr_ptr:      Pointer to the single-cell-refinement grid.
bool isRefinedFaceOnLgrBoundary(const std::array<int,3>& cells_per_dim, int faceIdxInLgr,
                                const std::shared_ptr<Dune::cpgrid::CpGridData>& elemLgr_ptr);


/// @brief Get the parent face index where the new refined face lays on.
///
/// @param [in] cells_per_dim:    Total children cells in each direction (x-,y-, and z-direction) of the single-cell-refinement.
/// @param [in] faceIdxInLgr:     Face index in the single-cell-refinement.
/// @param [in] elemLgr_ptr:      Pointer to the elemLgr single-cell-refinement grid.
/// @param [in] elemLgr:          Cell index from starting grid, that has been refined into a single-cell-refinement.
int getParentFaceWhereNewRefinedFaceLiesOn(const Dune::CpGrid& grid,
                                           const std::array<int,3>& cells_per_dim, int faceIdxInLgr,
                                           const std::shared_ptr<Dune::cpgrid::CpGridData>& elemLgr_ptr,
                                           int elemLgr);



/// @brief Define the corners (geometry) for each refined level grid.
void populateRefinedCorners(std::vector<Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>>& refined_corners_vec,
                            const std::vector<int>& refined_corner_count_vec,
                            const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                            const int& preAdaptMaxLevel,
                            const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner);

/// @brief Define the faces, face tags, face normarls, and face_to_point_, for each refined level grid.
void populateRefinedFaces(std::vector<Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<2,3>>>& refined_faces_vec,
                          std::vector<Dune::cpgrid::EntityVariableBase<enum face_tag>>& mutable_refined_face_tags_vec,
                          std::vector<Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>>& mutable_refine_face_normals_vec,
                          std::vector<Opm::SparseTable<int>>& refined_face_to_point_vec,
                          const std::vector<int>& refined_face_count_vec,
                          const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                          const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                          const std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                          const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                          const int& preAdaptMaxLevel,
                          const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                          const std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner);


/// @brief Define the cells, cell_to_point_, global_cell_, cell_to_face_, face_to_cell_, for each refined level grid.
void populateRefinedCells(const Dune::CpGrid& grid,
                          std::vector<Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<3,3>>>& refined_cells_vec,
                          std::vector<std::vector<std::array<int,8>>>& refined_cell_to_point_vec,
                          std::vector<std::vector<int>>& refined_global_cell_vec,
                          const std::vector<int>& refined_cell_count_vec,
                          std::vector<Dune::cpgrid::OrientedEntityTable<0,1>>& refined_cell_to_face_vec,
                          std::vector<Dune::cpgrid::OrientedEntityTable<1,0>>& refined_face_to_cell_vec,
                          const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                          const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                          const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                          const std::vector<Dune::cpgrid::DefaultGeometryPolicy>& refined_geometries_vec,
                          const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                          const std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                          const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                          const std::vector<int>& assignRefinedLevel,
                          const int& preAdaptMaxLevel,
                          const std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                          const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                          const std::vector<std::array<int,3>>&  cells_per_dim_vec);


/// @brief A new refined face lays on the boudndary of a single-cell-refinement appears in at most two single-cell-refinements. Given the face index in one
///        single-cell-refinement, compute the face index in a neighboring single-cell-refinement.
///
/// @param [in] cells_per_dim_lgr1:    Total children cells in each direction (x-,y-, and z-direction) of the elemLgr1 single-cell-refinement.
/// @param [in] faceIdxInLgr1:         Face index in the elemLgr1 single-cell-refinement.
/// @param [in] elemLgr1_ptr:          Pointer to the elemLgr1 single-cell-refinement grid.
/// @param [in] cells_per_dim_lgr2:    Total children cells in each direction (x-,y-, and z-direction) of the elemLgr2 single-cell-refinement.
int replaceLgr1FaceIdxByLgr2FaceIdx(const std::array<int,3>& cells_per_dim_lgr1, int faceIdxInLgr1,
                                    const std::shared_ptr<Dune::cpgrid::CpGridData>& elemLgr1_ptr,
                                    const std::array<int,3>& cells_per_dim_lgr2);


/// @brief Set geometrical and topological attributes for each refined level grid.
void setRefinedLevelGridsGeometries( const Dune::CpGrid& grid,
                                     /* Refined corner arguments */
                                     std::vector<Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>>& refined_corners_vec,
                                     const std::vector<int>& refined_corner_count_vec,
                                     /* Refined face arguments */
                                     std::vector<Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<2,3>>>& refined_faces_vec,
                                     std::vector<Dune::cpgrid::EntityVariableBase<enum face_tag>>& mutable_refined_face_tags_vec,
                                     std::vector<Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>>& mutable_refine_face_normals_vec,
                                     std::vector<Opm::SparseTable<int>>& refined_face_to_point_vec,
                                     const std::vector<int>& refined_face_count_vec,
                                     /* Refined cell argumets */
                                     std::vector<Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<3,3>>>& refined_cells_vec,
                                     std::vector<std::vector<std::array<int,8>>>& refined_cell_to_point_vec,
                                     std::vector<std::vector<int>>& refined_global_cell_vec,
                                     const std::vector<int>& refined_cell_count_vec,
                                     std::vector<Dune::cpgrid::OrientedEntityTable<0,1>>& refined_cell_to_face_vec,
                                     std::vector<Dune::cpgrid::OrientedEntityTable<1,0>>& refined_face_to_cell_vec,
                                     /* Auxiliary arguments */
                                     const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                     const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                                     const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner,
                                     const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                     const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                                     const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                     const std::vector<Dune::cpgrid::DefaultGeometryPolicy>& refined_geometries_vec,
                                     const std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                     const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                     const std::vector<int>& assignRefinedLevel,
                                     const int& preAdaptMaxLevel,
                                     const std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                                     const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                     const std::vector<std::array<int,3>>&  cells_per_dim_vec);


}


#endif // OPM_GRID_CPGRID_LGRHELPERS_HEADER_INCLUDED
