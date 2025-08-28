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

#include <array>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility> // for std::pair
#include <vector>

namespace Dune
{
namespace cpgrid
{
class CpGridData;
}
}

namespace Opm
{
namespace Lgr
{
/// --------------- Auxiliary methods to support refinement ---------------

/// @brief Refines all marked elements and establishes mappings between
///        original (corners, faces, cells) and their refined/adapted counterparts.
///
/// --- Marked element data ---
/// @param [out] markedElem_to_itsLgr   Auxiliary LGR for each marked element, later used to build refined level grids.
/// @param [out] markedElem_count       Total number of marked elements (used in grid info).
/// @param [out] cornerInMarkedElemWithEquivRefinedCorner
///     Maps each corner (from starting grid) to the marked elements it belongs to and its equivalent refined corner indices.
/// @param [out] markedElemAndEquivRefinedCorn_to_corner
///     Reverse mapping: (marked element idx, marked-element-refined corner idx) -> original corner idx.
/// @param [out] faceInMarkedElemAndRefinedFaces
///     For each original face, stores marked element indices where it appears (<=2) and their refined face indices.
///
/// --- Refined cell data ---
/// @param [out] elemLgrAndElemLgrCell_to_refinedLevelAdRefinedCell
///     Maps (marked element idx, marked-element-refined cell idx) -> (refined level grid, refined level cell idx).
/// @param [out] refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell
///     Inverse mapping: (refined level grid, refined cell idx) -> (marked element idx, merked-element-refined cell idx).
/// @param [out] refined_cell_count_vec  Number of refined cells per refined level grid.
/// @param [in]  assignRefinedLevel      For each element: 0 = not marked, >0 = assigned refined level.
/// @param [out] preAdapt_parent_to_children_cells_vec
///     For each pre-adapt/refine grid cell: children in refined grid, or {-1,{}} if none.
///
/// --- Adapted cell data ---
/// @param [out] elemLgrAndElemLgrCell_to_adaptedCell
///     Maps (marked element idx, marked-element-refined cell idx) -> leaf/adapted grid cell.
/// @param [out] adaptedCell_to_elemLgrAndElemLgrCell
///     Inverse mapping: adapted cell idx -> (marked element idx, marked-element-refined cell idx).
/// @param [out] cell_count              Total number of cells in adapted (leaf) grid.
/// @param [out] preAdapt_level_to_leaf_cells_vec
///     Maps pre-adapt grid cells to adapted cells (−1 if vanished).
///
/// --- Additional ---
/// @param [in]  cells_per_dim_vec       Refinement factors per dimension for each refined level grid.
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

/// @brief Establish child–parent relations for refined cells:
///        - Maps each refined cell (in new refined grids) to its parent cell in the pre-adapt grid(s).
///        - Assigns each child cell an index within its parent (-1 if no parent).
///
/// @param [in] refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell
///     Maps (refined level grid, refined cell idx)->(marked element idx, marked-element-refined cell idx).
/// @param [in] refined_cell_count_vec
///     Number of refined cells per refined level grid.
/// @param [in] adaptedCell_to_elemLgrAndElemLgrCell
///     Maps adapted cell idx -> (marked element idx, marked-element-refined cell idx).
/// @param [in] cell_count
///     Total number of cells in the adapted (leaf) grid.
///
/// @return
/// - refined_child_to_parent_cells_vec : For each refined level grid, maps child cell->{parent level, parent index}, or {−1,−1} if none.
/// - refined_cell_to_idxInParentCell_vec: For each refined level grid, child cell’s index within parent (for geometryInFather), -1 if none.
/// - adapted_child_to_parent_cell       : Maps adapted child cell idx->{parent level, parent index}, or {−1,−1} if none.
/// - adapted_cell_to_idxInParentCell    : Index of adapted child cell within parent, -1 if none.
std::tuple< std::vector<std::vector<std::array<int,2>>>,
            std::vector<std::vector<int>>,
            std::vector<std::array<int,2>>,
            std::vector<int>>
defineChildToParentAndIdxInParentCell(const Dune::cpgrid::CpGridData& current_data,
                                      int preAdaptMaxLevel,
                                      const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                      const std::vector<int>& refined_cell_count_vec,
                                      const std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                                      const int& cell_count);

/// @brief Define index mappings between refined level grids and the leaf (adapted) grid:
///        - level_to_leaf_cells: maps each cell idx in a refined level grid to its leaf cell idx
///        - leaf_to_level_cells: maps each leaf cell idx back to its originating level cell idx
///
/// @param [in] elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell
///     Maps (marked element idx, marked-element-refined cell idx) -> (refined level grid, refined cell idx).
/// @param [in] refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell
///     Inverse mapping: (refined level grid, refined cell) -> (marked element idx, marked-element-refined cell idx).
/// @param [in] refined_cell_count_vec
///     Number of refined cells per level.
/// @param [in] elemLgrAndElemLgrCell_to_adaptedCell
///     Maps (marked element idx, marked-element-refined cell idx) -> leaf (adapted) cell idx.
/// @param [in] adaptedCell_to_elemLgrAndElemLgrCell
///     Inverse mapping: leaf (adapted) cell -> (marked element idx, marked-element-refined cell idx).
/// @param [in] cell_count
///     Total number of cells in the leaf (adapted) grid.
///
/// @return
/// - refined_level_to_leaf_cells_vec : For each refined level grid, maps cell idx -> leaf cell idx
/// - leaf_to_level_cells             : For each leaf cell, {level grid idx, cell idx in that level}
std::pair<std::vector<std::vector<int>>, std::vector<std::array<int,2>>>
defineLevelToLeafAndLeafToLevelCells(const Dune::cpgrid::CpGridData& current_data,
                                     int preAdaptMaxLevel,
                                     const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell,
                                     const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                     const std::vector<int>& refined_cell_count_vec,
                                     const std::map<std::array<int,2>,int>& elemLgrAndElemLgrCell_to_adaptedCell,
                                     const std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                                     const int& cell_count);

/// @brief Define refined corner relations:
///        1. Map corners from single marked-element refinements to their refined level grids (and inverse).
///        2. Handle corners that appear in multiple single-cell refinements but are stored only once
///           (last appearance kept, earlier ones mapped to it).
///
/// @param [out] elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner
///     Maps (marked element idx, marked-element-refined corner idx) -> (refined level grid, refined corner idx).
/// @param [out] refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner
///     Inverse mapping: (refined level grid, refined corner idx) -> (marked element idx, marked-element-refined corner idx).
/// @param [out] refined_corner_count_vec
///     Number of refined corners per refined level grid.
/// @param [out] vanishedRefinedCorner_to_itsLastAppearance
///     For corners present in multiple refinements: map each earlier occurrence -> its last appearance.
///
/// @param [in] markedElem_to_itsLgr
/// @param [in] assignRefinedLevel
/// @param [in] cornerInMarkedElemWithEquivRefinedCorner
/// @param [in] faceInMarkedElemAndRefinedFaces
/// @param [in] cells_per_dim_vec
void identifyRefinedCornersPerLevel(const Dune::cpgrid::CpGridData& current_data,
                                    int preAdaptMaxLevel,
                                    std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                    std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner,
                                    std::vector<int>& refined_corner_count_vec,
                                    std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                    const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                    const std::vector<int>& assignRefinedLevel,
                                    const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                    const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                    const std::vector<std::array<int,3>>& cells_per_dim_vec);

/// @brief Check if a refined corner lies in the interior of a single-cell refinement.
///
/// @param [in] cells_per_dim  Number of child cells in {x, y, z} directions.
/// @param [in] cornerIdxInLgr Corner index within the single-cell refinement.
/// @return true if the corner is interior, false otherwise.
bool isRefinedCornerInInteriorLgr(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr);

/// @brief Compute the {i,j,k} index of a refined corner from its linear index
///        in a single-cell refinement.
///
/// A cell refined into {nx, ny, nz} children has (nx+1)(ny+1)(nz+1) corners.
/// Corners are ordered as:
///   idx = j*(nx+1)(nz+1) + i*(nz+1) + k,
/// with i in [0,nx], j in [0,ny], k in [0,nz].
/// This function converts cornerIdxInLgr (0..(nx+1)(ny+1)(nz+1)-1) to {i,j,k}.
///
/// @param [in] cells_per_dim  Number of child cells in {x,y,z} directions.
/// @param [in] cornerIdxInLgr Corner index in the single-cell refinement.
/// @return {i,j,k} index of the corner.
std::array<int,3> getRefinedCornerIJK(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr);

/// @brief Check if a refined corner lies on a parent-cell edge.
///        Specifically, on the boundary of the single-cell refinement,
///        but not coinciding with a pre-adapt corner.
///
/// @param [in] cells_per_dim  Number of child cells in {x,y,z} directions.
/// @param [in] cornerIdxInLgr Corner index in the single-cell refinement.
/// @return true if the corner is on an edge, false otherwise.
bool newRefinedCornerLiesOnEdge(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr);

/// @brief Get the parent faces that contain the edge on which a new refined corner lies.
///
/// @param [in] cells_per_dim  Number of child cells in {x,y,z} directions.
/// @param [in] cornerIdxInLgr Corner index in the single-cell refinement.
/// @param [in] elemLgr        Parent cell index (from the original grid) refined into a single-cell refinement.
/// @return Indices of the two parent faces sharing the edge.
std::array<int,2> getParentFacesAssocWithNewRefinedCornLyingOnEdge(const Dune::cpgrid::CpGridData& current_data,
                                                                   const std::array<int,3>& cells_per_dim,
                                                                   int cornerIdxInLgr,
                                                                   int elemLgr);

/// @brief Check if a refined corner lies on the boundary of the single-cell refinement
///        and is not a pre-adapt (original) corner.
///
/// @param [in] cells_per_dim  Number of child cells in {x,y,z} directions.
/// @param [in] cornerIdxInLgr Corner index in the single-cell refinement.
/// @return true if the corner is on the boundary and new, false otherwise.
bool isRefinedNewBornCornerOnLgrBoundary(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr);

// @brief Get the parent face containing the new refined corner.
///
/// @param [in] cells_per_dim  Number of child cells in {x, y, z} directions.
/// @param [in] cornerIdxInLgr Corner index in the single-cell refinement.
/// @param [in] elemLgr        Parent cell index (from the original grid) refined into a single-cell refinement.
/// @return Index of the parent face containing the corner.
int getParentFaceWhereNewRefinedCornerLiesOn(const Dune::cpgrid::CpGridData& current_data,
                                             const std::array<int,3>& cells_per_dim,
                                             int cornerIdxInLgr,
                                             int elemLgr);

/// @brief Map a refined corner from one single-cell refinement to a neighboring refinement.
///
/// Given a corner index in the first single-cell refinement, compute the corresponding corner index
/// in the second single-cell refinement.
///
/// @param [in] cells_per_dim_lgr1 Number of child cells in {x,y,z} directions for the first refinement.
/// @param [in] cornerIdxLgr1      Corner index in the first single-cell refinement.
/// @param [in] cells_per_dim_lgr2 Number of child cells in {x,y,z} directions for the second refinement.
/// @return Corresponding corner index in the second single-cell refinement.
int replaceLgr1CornerIdxByLgr2CornerIdx(const std::array<int,3>& cells_per_dim_lgr1,
                                        int cornerIdxLgr1,
                                        const std::array<int,3>& cells_per_dim_lgr2);

/// @brief Map a new refined corner on an edge from one single-cell refinement to a neighboring refinement.
///
/// Given a corner index in one refinement, compute its index in a neighboring refinement. The corner
/// lies on an edge and may appear in multiple refinements.
///
/// @param [in] cells_per_dim_lgr1           Number of child cells in {x,y,z} directions for the first refinement.
/// @param [in] cornerIdxLgr1                Corner index in the first single-cell refinement.
/// @param [in] elemLgr1                     Parent cell index of the first refinement.
/// @param [in] parentFaceLastAppearanceIdx  Parent face index where the corner appears last.
/// @param [in] cells_per_dim_lgr2           Number of child cells in {x,y,z} directions for the second refinement.
/// @return Corresponding corner index in the second single-cell refinement.
int replaceLgr1CornerIdxByLgr2CornerIdx(const Dune::cpgrid::CpGridData& current_data,
                                        const std::array<int,3>& cells_per_dim_lgr1,
                                        int cornerIdxLgr1, int elemLgr1, int parentFaceLastAppearanceIdx,
                                        const std::array<int,3>& cells_per_dim_lgr2);

/// @brief Identify corners on the leaf (adapted) grid and establish corner mappings.
///
/// Maps pre-adapt and refined corners from single-element refinements to the leaf grid,
/// and vice versa.
///
/// @param [out] elemLgrAndElemLgrCorner_to_adaptedCorner  
///     Maps (marked element idx, marked-element-refined corner idx) or {-1, corner in original grid} -> leaf grid corner idx.
/// @param [out] adaptedCorner_to_elemLgrAndElemLgrCorner  
///     Inverse mapping: leaf grid corner idx -> (marked element idx, marked-element-refined corner idx) or {-1, original corner idx}.
/// @param [out] corner_count  Total number of corners in the leaf (adapted) grid.
/// @param [in] markedElem_to_itsLgr  
/// @param [in] assignRefinedLevel  
/// @param [in] cornerInMarkedElemWithEquivRefinedCorner  
/// @param [in] vanishedRefinedCorner_to_itsLastAppearance  
/// @param [in] faceInMarkedElemAndRefinedFaces  
/// @param [in] cells_per_dim_vec
void identifyLeafGridCorners(const Dune::cpgrid::CpGridData& current_data,
                             int preAdaptMaxLevel,
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
                        const Dune::cpgrid::CpGridData& current_data,
                        int preAdaptMaxLevel,
                        const std::vector<int>& assignRefinedLevel,
                        const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                        const std::vector<std::array<int,3>>& cells_per_dim_vec);

void processBoundaryCorners(int elemIdx, int shiftedLevel,
                            const std::shared_ptr<Dune::cpgrid::CpGridData>& lgr,
                            int& corner_count,
                            std::map<std::array<int,2>,int>& elemLgrAndElemLgrCorner_to_adaptedCorner,
                            std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner,
                            std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                            const Dune::cpgrid::CpGridData& current_data,
                            int preAdaptMaxLevel,
                            const std::vector<int>& assignRefinedLevel,
                            const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                            const std::vector<std::array<int,3>>& cells_per_dim_vec);

// To insert bidirectional mapping and increment counter
void insertBidirectional(std::map<std::array<int,2>,std::array<int,2>>& a_to_b,
                         std::map<std::array<int,2>,std::array<int,2>>& b_to_a,
                         const std::array<int,2>& keyA,
                         const std::array<int,2>& keyB,
                         int& counter,
                         bool useFullKeyB = false);

/// @brief Define mappings between single-cell-refinement faces and refined level faces.
///
/// Maps faces from each marked element’s single-cell refinement to its assigned refined level grid,
/// and provides the inverse mapping.
///
/// @param [out] elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace  
///     Maps (marked element idx, marked-element-refined face idx)->(refined level grid, refined face idx).
/// @param [out] refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace  
///     Inverse mapping: (refined level grid, refined face idx)->(marked element idx, marked-element-refined face idx).
/// @param [out] refined_face_count_vec  Number of refined faces per refined level grid.
/// @param [in] markedElem_to_itsLgr
/// @param [in] assignRefinedLevel
/// @param [in] faceInMarkedElemAndRefinedFaces
/// @param [in] cells_per_dim_vec
void identifyRefinedFacesPerLevel(const Dune::cpgrid::CpGridData& current_data,
                                  int preAdaptMaxLevel,
                                  std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                                  std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                                  std::vector<int>& refined_face_count_vec,
                                  const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                  const std::vector<int>& assignRefinedLevel,
                                  const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                  const std::vector<std::array<int,3>>& cells_per_dim_vec);

/// @brief Identify faces on the leaf (adapted) grid and establish face mappings.
///
/// Maps pre-adapt and refined faces from single-element refinements to the leaf grid,
/// and provides the inverse mapping.
///
/// @param [out] elemLgrAndElemLgrFace_to_adaptedFace
///     Maps (marked element idx, marked-element-refined face idx) or {-1, face idx in original grid} -> leaf grid face index.
/// @param [out] adaptedFace_to_elemLgrAndElemLgrFace
///     Inverse mapping: leaf grid face -> (marked element idx, marked-element-refined face idx) or {-1, original face idx}.
/// @param [out] face_count  Total number of faces in the leaf (adapted) grid.
/// @param [in] markedElem_to_itsLgr
/// @param [in] assignRefinedLevel
/// @param [in] faceInMarkedElemAndRefinedFaces
/// @param [in] cells_per_dim_vec
void identifyLeafGridFaces(const Dune::cpgrid::CpGridData& current_data,
                           int preAdaptMaxLevel,
                           std::map<std::array<int,2>,int>& elemLgrAndElemLgrFace_to_adaptedFace,
                           std::unordered_map<int,std::array<int,2>>& adaptedFace_to_elemLgrAndElemLgrFace,
                           int& face_count,
                           const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                           const std::vector<int>& assignRefinedLevel,
                           const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                           const std::vector<std::array<int,3>>& cells_per_dim_vec);

/// @brief Compute the {i,j,k} index of a refined face from its linear index in a single-cell refinement.
///
/// A single cell refined into {nx, ny, nz} children has
/// ((nx+1)*ny*nz) + (nx*(ny+1)*nz) + (nx*ny*(nz+1)) faces, stored in the order defined by Geometry::refine:
/// - K_FACES: (k*nx*ny) + (j*nx) + i  
/// - I_FACES: (nx*ny*(nz+1)) + (i*ny*nz) + (k*ny) + j  
/// - J_FACES: (nx*ny*(nz+1)) + ((nx+1)*ny*nz) + (j*nx*nz) + (i*nz) + k  
/// where i=0..nx-1, j=0..ny-1, k=0..nz-1.  
/// This function converts faceIdxInLgr to its {i,j,k} index.
///
/// @param [in] cells_per_dim  Number of child cells in {x,y,z} directions.
/// @param [in] faceIdxInLgr   Face index in the single-cell refinement.
/// @param [in] elemLgr_ptr    Pointer to the single-cell refinement grid.
/// @return {i,j,k} index of the face
std::array<int,3> getRefinedFaceIJK(const std::array<int,3>& cells_per_dim, int faceIdxInLgr,
                                    const std::shared_ptr<Dune::cpgrid::CpGridData>& elemLgr_ptr);

/// @brief Check if a refined face lies in the interior of a single-cell refinement.
///
/// @param [in] cells_per_dim  Number of child cells in {x,y,z} directions.
/// @param [in] faceIdxInLgr   Face index in the single-cell refinement.
/// @param [in] elemLgr_ptr    Pointer to the single-cell refinement grid.
/// @return true if the face is interior, false otherwise.
bool isRefinedFaceInInteriorLgr(const std::array<int,3>& cells_per_dim, int faceIdxInLgr,
                                const std::shared_ptr<Dune::cpgrid::CpGridData>& elemLgr_ptr);

/// @brief Check if a refined face lies on the boundary of a single-cell refinement.
///
/// @param [in] cells_per_dim  Number of child cells in {x,y,z} directions.
/// @param [in] faceIdxInLgr   Face index in the single-cell refinement.
/// @param [in] elemLgr_ptr    Pointer to the single-cell refinement grid.
/// @return true if the face is on the boundary, false otherwise.
bool isRefinedFaceOnLgrBoundary(const std::array<int,3>& cells_per_dim, int faceIdxInLgr,
                                const std::shared_ptr<Dune::cpgrid::CpGridData>& elemLgr_ptr);

/// @brief Get the parent face containing a new refined face.
///
/// @param [in] cells_per_dim  Number of child cells in {x,y,z} directions.
/// @param [in] faceIdxInLgr   Face index in the single-cell refinement.
/// @param [in] elemLgr_ptr    Pointer to the single-cell refinement grid.
/// @param [in] elemLgr        Parent cell index from the original grid.
/// @return Index of the parent face containing the refined face.
int getParentFaceWhereNewRefinedFaceLiesOn(const Dune::cpgrid::CpGridData& current_data,
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
void populateRefinedCells(const Dune::cpgrid::CpGridData& current_data,
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

/// @brief Map a refined boundary face from one single-cell refinement to a neighboring refinement.
///
/// A refined face on the boundary may appear in up to two single-cell refinements. Given its index
/// in the first refinement, this function computes the corresponding face index in the neighboring refinement.
///
/// @param [in] cells_per_dim_lgr1  Number of child cells in {x,y,z} directions for the first refinement.
/// @param [in] faceIdxInLgr1       Face index in the first single-cell refinement.
/// @param [in] elemLgr1_ptr        Pointer to the first single-cell refinement grid.
/// @param [in] cells_per_dim_lgr2  Number of child cells in {x,y,z} directions for the second refinement.
/// @return Corresponding face index in the second single-cell refinement.
int replaceLgr1FaceIdxByLgr2FaceIdx(const std::array<int,3>& cells_per_dim_lgr1, int faceIdxInLgr1,
                                    const std::shared_ptr<Dune::cpgrid::CpGridData>& elemLgr1_ptr,
                                    const std::array<int,3>& cells_per_dim_lgr2);

/// @brief Set geometrical and topological attributes for each refined level grid.
void setRefinedLevelGridsGeometries(const Dune::cpgrid::CpGridData& current_data,
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

/// @brief Detect active local refinement grids (LGRs) on each process.
///
/// For blocks of cells selected for refinement on a level-zero distributed grid, this function
/// marks which LGRs contain at least one active cell on the current process.
///
/// @param [in] startIJK_vec  Start indices {i,j,k} of each refinement block.
/// @param [in] endIJK_vec    End indices {i,j,k} of each refinement block.
/// @param [out] lgr_with_at_least_one_active_cell
///     For each level, set to 1 if the LGR contains at least one active cell, 0 otherwise.
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

} // namespace Lgr
} // namespace Opm


#endif // OPM_GRID_CPGRID_LGRHELPERS_HEADER_INCLUDED
