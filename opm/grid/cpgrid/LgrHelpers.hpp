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

void markVanishedCorner(const std::array<int,2>& vanished,
                               const std::array<int,2>& lastAppearance,
                        std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance);

void processInteriorCorners(int elemIdx, int shiftedLevel,
                            const std::shared_ptr<Dune::cpgrid::CpGridData>& lgr,
                            int& corner_count,
                            std::map<std::array<int,2>,int>& elemLgrAndElemLgrCorner_to_adaptedCorner,
                            std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner,
                            const std::vector<std::array<int,3>>& cells_per_dim_vec);


void processEdgeCorners(int elemIdx, int shiftedLevel,
                        const std::shared_ptr<Dune::cpgrid::CpGridData>& lgr,
                        int& corner_count,
                        std::map<std::array<int,2>,int>& elemLgrAndElemLgrCorner_to_adaptedCorner,
                        std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner,
                        std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                        const Dune::CpGrid& grid,
                        const std::vector<int>& assignRefinedLevel,
                        const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                        const std::vector<std::array<int,3>>& cells_per_dim_vec);

void processBoundaryCorners(int elemIdx, int shiftedLevel,
                            const std::shared_ptr<Dune::cpgrid::CpGridData>& lgr,
                            int& corner_count,
                            std::map<std::array<int,2>,int>& elemLgrAndElemLgrCorner_to_adaptedCorner,
                            std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner,
                            std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                            const Dune::CpGrid& grid,
                            const std::vector<int>& assignRefinedLevel,
                             const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                            const std::vector<std::array<int,3>>& cells_per_dim_vec);

// Small helper: insert bidirectional mapping and increment counter
void insertBidirectional(std::map<std::array<int,2>,std::array<int,2>>& a_to_b,
                         std::map<std::array<int,2>,std::array<int,2>>& b_to_a,
                         const std::array<int,2>& keyA,
                         const std::array<int,2>& keyB,
                         int& counter);

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


// Generic helper: insert bidirectional mapping and increment counter
void insertBidirectionalArrayToInt( std::map<std::array<int,2>,int>& a_to_b,
                                    std::unordered_map<int,std::array<int,2>>& b_to_a,
                                    const std::array<int,2>& keyA, const int& keyB,
                                    int& counter);

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


/// @brief Define the corners (gemotry) for the leaf grid view (or adapted grid).
void populateLeafGridCorners(const Dune::CpGrid& grid,
                             Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& adapted_corners,
                             const int& corners_count,
                             const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                             const std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner);

/// @brief Define the faces, face tags, face normarls, and face_to_point_, for the leaf grid view.
void populateLeafGridFaces(const Dune::CpGrid& grid,
                           Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<2,3>>& adapted_faces,
                           Dune::cpgrid::EntityVariableBase<enum face_tag>& mutable_face_tags,
                           Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& mutable_face_normals,
                           Opm::SparseTable<int>& adapted_face_to_point,
                           const int& face_count,
                           const std::unordered_map<int,std::array<int,2>>& adaptedFace_to_elemLgrAndElemLgrFace,
                           const std::map<std::array<int,2>,int>& elemLgrAndElemLgrCorner_to_adaptedCorner,
                           const std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                           const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                           const std::vector<int>& assignRefinedLevel,
                           const std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                           const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                           const std::vector<std::array<int,3>>& cells_per_dim_vec,
                           const int& preAdaptMaxLevel);

/// @brief Define the cells, cell_to_point_, cell_to_face_, face_to_cell_, for the leaf grid view (or adapted grid).
void populateLeafGridCells(const Dune::CpGrid& grid,
                           Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<3,3>>& adapted_cells,
                           std::vector<std::array<int,8>>& adapted_cell_to_point,
                           const int& cell_count,
                           Dune::cpgrid::OrientedEntityTable<0,1>& adapted_cell_to_face,
                           Dune::cpgrid::OrientedEntityTable<1,0>& adapted_face_to_cell,
                           const std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                           const std::map<std::array<int,2>,int>& elemLgrAndElemLgrFace_to_adaptedFace,
                           const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                           const Dune::cpgrid::DefaultGeometryPolicy& adapted_geometries,
                           const std::map<std::array<int,2>,int>& elemLgrAndElemLgrCorner_to_adaptedCorner,
                           const std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                           const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                           const std::vector<int>& assignRefinedLevel,
                           const std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                           const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                           const std::vector<std::array<int,3>>& cells_per_dim_vec,
                           const int& preAdaptMaxLevel);

/// @brief Define geometrical and topological attributes for the leaf grid view (or adapted grid).
void updateLeafGridViewGeometries( const Dune::CpGrid& grid,
                                   /* Leaf grid View Corners arguments */
                                   Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& adapted_corners,
                                   const int& corner_count,
                                   /* Leaf grid View Faces arguments */
                                   Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<2,3>>& adapted_faces,
                                   Dune::cpgrid::EntityVariableBase<enum face_tag>& mutable_face_tags,
                                   Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& mutable_face_normals,
                                   Opm::SparseTable<int>& adapted_face_to_point,
                                   const int& face_count,
                                   /* Leaf grid View Cells argumemts  */
                                   Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<3,3>>& adapted_cells,
                                   std::vector<std::array<int,8>>& adapted_cell_to_point,
                                   const int& cell_count,
                                   Dune::cpgrid::OrientedEntityTable<0,1>& adapted_cell_to_face,
                                   Dune::cpgrid::OrientedEntityTable<1,0>& adapted_face_to_cell,
                                   /* Auxiliary arguments */
                                   const std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner,
                                   const std::unordered_map<int,std::array<int,2>>& adaptedFace_to_elemLgrAndElemLgrFace,
                                   const std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                                   const std::map<std::array<int,2>,int>& elemLgrAndElemLgrFace_to_adaptedFace,
                                   const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                   const Dune::cpgrid::DefaultGeometryPolicy& adapted_geometries,
                                   const std::map<std::array<int,2>,int>& elemLgrAndElemLgrCorner_to_adaptedCorner,
                                   const std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                   const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                   const std::vector<int>& assignRefinedLevel,
                                   const std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                                   const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                   const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                   const int& preAdaptMaxLevel);

/// @brief Auxilliary function to compute one or more properties on selected block of parent cells.
///
/// @param [in] startIJK_vec    Vector of ijk values denoting the start of each block of cells selected for refinement.
/// @param [in] endIJK_vec      Vector of ijk values denoting the end of each block of cells selected for refinement.
/// @param [in] function        Lambda expression/function that computes the desired properties for each parent cell.
/// The full definition needs to be in the header so that the compiler can instantiate it when needed.
template<class T>
void computeOnLgrParents(const Dune::CpGrid& grid,
                         const std::vector<std::array<int,3>>& startIJK_vec,
                         const std::vector<std::array<int,3>>& endIJK_vec,
                         T func)
{
    // Find out which (ACTIVE) elements belong to the block cells defined by startIJK and endIJK values.
    for(const auto& element: Dune::elements(grid.leafGridView())) {
        std::array<int,3> ijk;
        grid.getIJK(element.index(), ijk);
        for (std::size_t level = 0; level < startIJK_vec.size(); ++level) {
            bool belongsToLevel = true;
            for (int c = 0; c < 3; ++c) {
                belongsToLevel = belongsToLevel && ( (ijk[c] >= startIJK_vec[level][c]) && (ijk[c] < endIJK_vec[level][c]) );
                if (!belongsToLevel)
                    break;
            }
            if(belongsToLevel) {
                func(element, level);
            }
        }
    }
}


/// @brief Detect active LGRs in each process.
///
/// Given blocks of cells selected for refinement on a level zero distributed grid, detect which LGRs are active
/// in each process.
///
/// @param [in] startIJK_vec    Vector of ijk values denoting the start of each block of cells selected for refinement.
/// @param [in] endIJK_vec      Vector of ijk values denoting the end of each block of cells selected for refinement.
/// @param [out] lgr_with_at_least_one_active_cell Determine if an LGR is not empty in a given process, we set
///                                                lgr_with_at_least_one_active_cell[in that level] to 1 if it contains
///                                                at least one active cell, and to 0 otherwise.
void detectActiveLgrs(const Dune::CpGrid& grid,
                      const std::vector<std::array<int,3>>& startIJK_vec,
                      const std::vector<std::array<int,3>>& endIJK_vec,
                      std::vector<int>& lgr_with_at_least_one_active_cell);

/// @brief Predict minimum cell and point global ids per process.
///
/// Predict how many new cells/points (born in refined level grids) need new globalIds, so we can assign unique
/// new ids ( and anticipate the maximum). At this point, the grid is already refined according to the LGR specification.
///
/// @param [in] assignRefinedLevel   Assign level for the refinement of each marked cell. Example: refined element from
///                                  LGR1 have level 1, refined element rfom LGR2 have level 2, etc.
/// @param [in] cells_per_dim_vec    Total child cells in each direction (x-,y-, and z-direction) per block of cells.
/// @param [in] lgr_with_at_least_one_active_cell  Determine if an LGR is not empty in a given process:
///                                                lgr_with_at_least_one_active_cell[level] = 1 if it contains
///                                                at least one active cell in the current process, and 0 otherwise.
/// @param [out] min_globalId_cell_in_proc
/// @param [out] min_globalId_point_in_proc
void predictMinCellAndPointGlobalIdPerProcess(const Dune::CpGrid& grid,
                                              const std::vector<int>& assignRefinedLevel,
                                              const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                              const std::vector<int>& lgr_with_at_least_one_active_cell,
                                              int& min_globalId_cell_in_proc,
                                              int& min_globalId_point_in_proc);

/// @brief Assign cell global ids of new born cell from refined level grids. Assign 'candidate' point global ids
///        for points in refined level grids.
///
/// @param [out] localToGlobal_cells_per_level    Relation local element.index() to assigned cell global id.
/// @param [out] localToGlobal_points_per_level   Relation local point.index() to assigned 'candidate' global id.
/// @param [in] min_globalId_cell_in_proc         Minimum cell global id per process.
/// @param [in] min_globalId_point_in_proc        Minimum point global id per process.
/// @param [in] cells_per_dim_vec                 Total child cells in each direction (x-,y-, and z-direction) per block of cells.
void assignCellIdsAndCandidatePointIds( const Dune::CpGrid& grid,
                                        std::vector<std::vector<int>>& localToGlobal_cells_per_level,
                                        std::vector<std::vector<int>>& localToGlobal_points_per_level,
                                        int min_globalId_cell_in_proc,
                                        int min_globalId_point_in_proc,
                                        const std::vector<std::array<int,3>>& cells_per_dim_vec);

/// @brief Select and re-write point global ids.
///
/// After assigning global IDs to points in refined-level grids, a single point may have
/// "multiple unique" global IDs, one in each process to which it belongs.
/// To reduce the unnucesary id assigments, since global IDs must be distinct across the global leaf view
/// and consistent across each refined-level grid, we will rewrite the entries in
/// localToGlobal_points_per_level. Using cell_to_point_ across all refined cells through
/// communication: gathering the 8 corner points of each interior cell and scattering the
/// 8 corner points of overlapping cells, for all child cells of a parent cell in level zero grid.
///
/// @param [out] localToGlobal_points_per_level   Relation local point.index() to assigned 'candidate' global id.
/// @param [in] parent_to_children                The communication step is based on level zero grid, via the relation parent-children-cells.
/// @param [in] cells_per_dim_vec                 Total child cells in each direction (x-,y-, and z-direction) per block of cells.
void selectWinnerPointIds(const Dune::CpGrid& grid,
                          std::vector<std::vector<int>>&  localToGlobal_points_per_level,
                          const std::vector<std::tuple<int,std::vector<int>>>& parent_to_children,
                          const std::vector<std::array<int,3>>& cells_per_dim_vec);

/// @brief Retrieves the global ids of the first child for each parent cell in the grid.
///
/// If a cell has no children, its entry is set to -1, indicating an invalid id.
///
/// @param[out] parentToFirstChildGlobalIds A vector that will be filled with the first child global IDs.
///                                         The vector is resized to match the number of parent cells.
void getFirstChildGlobalIds(const Dune::CpGrid& grid,
                            std::vector<int>& parentToFirstChildGlobalIds);

/// @brief Extract Cartesian index triplet (i,j,k) given an index between 0 and NXxNYxNZ -1
///    where NX, NY, and NZ is the total amoung of cells in each direction x-,y-,and z- respectively.
///
/// @param [in] idx      Integer between 0 and cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]-1
/// @param [in] cells_per_dim
/// @return Cartesian index triplet.
std::array<int,3> getIJK(int idx_in_parent_cell, const std::array<int,3>& cells_per_dim);

/// @brief Check startIJK and endIJK of each patch of cells to be refined are valid, i.e.
///        startIJK and endIJK vectors have the same size and, startIJK < endIJK coordenate by coordenate.
///
/// @param [in]  startIJK_vec       Vector of Cartesian triplet indices where each patch starts.
/// @param [in]  endIJK_vec         Vector of Cartesian triplet indices where each patch ends.
///                                 Last cell part of the lgr will be {endIJK_vec[patch][0]-1, ..., endIJK_vec[patch][2]-1}.
void validStartEndIJKs(const std::vector<std::array<int,3>>& startIJK_vec,
                       const std::vector<std::array<int,3>>& endIJK_vec);

/// @brief Compute patch boundary face indices (Cartesian grid required).
///
/// @param [in]  startIJK  Cartesian triplet index where the patch starts.
/// @param [in]  endIJK    Cartesian triplet index where the patch ends.
///                        Last cell part of the lgr will be {endijk[0]-1, ... endIJK[2]-1}.
///
/// @return patch_boundary_faces
std::array<std::vector<int>,6> getBoundaryPatchFaces(const std::array<int,3>& startIJK,
                                                     const std::array<int,3>& endIJK,
                                                     const std::array<int,3>& grid_dim);

/// @brief Compute amount of cells in each direction of a patch of cells. (Cartesian grid required).
///
/// @param [in]  startIJK  Cartesian triplet index where the patch starts.
/// @param [in]  endIJK    Cartesian triplet index where the patch ends.
///                        Last patch cell Cartesian triplet is {endijk[0]-1, ... endIJK[2]-1}.
///
/// @return patch_dim Patch dimension {#cells in x-direction, #cells in y-direction, #cells in z-direction}.
std::array<int,3> getPatchDim(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK);

/// @brief Determine if a finite amount of patches (of cells) share a face.
///
/// @param [in]  startIJK_vec  Vector of Cartesian triplet indices where each patch starts.
/// @param [in]  endIJK_vec    Vector of Cartesian triplet indices where each patch ends.
///                            Last patch Cartesian triplet is {endIJK_vec[<patch>][0]-1, ... ,endIJK_vec[<patch>][2]-1}.
bool patchesShareFace(const std::vector<std::array<int,3>>& startIJK_vec,
                      const std::vector<std::array<int,3>>& endIJK_vec,
                      const std::array<int,3>& grid_dim);

int sharedFaceTag(const std::vector<std::array<int,3>>& startIJK_2Patches,
                  const std::vector<std::array<int,3>>& endIJK_2Patches,
                  const std::array<int,3>& grid_dim);

/// @brief Filter out LGR entries that do not result in any actual refinement.
///
/// This function removes entries where the number of subdivisions in each direction is 0
/// (i.e., cells_per_dim is equal to {1, 1, 1}) which would result in no grid refinement.
///
/// A warning is logged for each excluded LGR name.
///
/// @param [in] startIJK_vec          Vector of Cartesian triplet indices where each patch starts.
/// @param [in] endIJK_vec            Vector of Cartesian triplet indices where each patch ends.
/// @param [in] lgr_name_vec          Names (std::string) for the LGRs/levels.
///
/// @return A tuple containing a bool and the filtered vectors:
///         - allUndesired         True if all LGRs have cells_per_dim_ = {1,1,1}
///         - cells_per_dim_vec
///         - startIJK_vec
///         - endIJK_vec
///         - lgr_name_vec
std::tuple<bool,
           std::vector<std::array<int,3>>,
           std::vector<std::array<int,3>>,
           std::vector<std::array<int,3>>,
           std::vector<std::string>>
filterUndesiredNumberOfSubdivisions(const std::vector<std::array<int, 3>>& cells_per_dim_vec,
                                    const std::vector<std::array<int, 3>>& startIJK_vec,
                                    const std::vector<std::array<int, 3>>& endIJK_vec,
                                    const std::vector<std::string>& lgr_name_vec);

}


#endif // OPM_GRID_CPGRID_LGRHELPERS_HEADER_INCLUDED
