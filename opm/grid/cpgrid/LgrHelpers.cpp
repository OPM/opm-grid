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
#include "config.h"


#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/cpgrid/LgrHelpers.hpp>
#include <opm/grid/cpgrid/ParentToChildCellToPointGlobalIdHandle.hpp>

#include <algorithm>
#include <array>
#include <memory>
#include <numeric>
#include <vector>

namespace Opm {

void refineAndProvideMarkedRefinedRelations(const Dune::CpGrid& grid, /* Marked elements parameters */
                                            std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                            int& markedElem_count,
                                            std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                            std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                                            std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                            /* Refined cells parameters */
                                            std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell,
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
                                            const std::vector<std::array<int,3>>& cells_per_dim_vec)
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.

    // Each marked element for refinement (mark equal to 1), will be refined individuality, creating its own Lgr. The element index will
    // be also used to identify its lgr. Even though, in the end, all the refined entities will belong to a unique level grid.
    // For this reason, we associate "-1" with those elements that are not involved in any refinement and will appear
    // as "coarse" cells in the leaf-grid-view (adapted-grid).

    // Max level before calling adapt.
    const int& preAdaptMaxLevel = grid.maxLevel();

    for (const auto& element : Dune::elements(grid.leafGridView())) {
        // When the element is marked with 0 ("doing nothing"), it will appear in the adapted grid with same geometrical features (center, volume).
        if (grid.getMark(element) ==  0) {
            elemLgrAndElemLgrCell_to_adaptedCell[{-1, element.index()}] = cell_count;
            adaptedCell_to_elemLgrAndElemLgrCell[cell_count] = {-1, element.index()};
            cell_count +=1;
            preAdapt_level_to_leaf_cells_vec[element.level()][element.getLevelElem().index()] = cell_count;
        }

        // When the element is marked for refinement, we also mark its corners and faces
        // since they will get replaced by refined ones.
        if (grid.getMark(element) ==  1) {
            markedElem_count +=1;
            const auto& markedElemLevel = assignRefinedLevel[element.index()];
            assert(markedElemLevel > preAdaptMaxLevel);
            // Shift the markedElemRefinedLevel to access data containers
            const auto& shiftedLevel = markedElemLevel - preAdaptMaxLevel-1;
            // Build auxiliary LGR for the refinement of this element
            const auto& [elemLgr_ptr,
                         parentCorners_to_equivalentRefinedCorners,
                         parentFace_to_itsRefinedFaces,
                         parentCell_to_itsRefinedCells,
                         refinedFace_to_itsParentFace,
                         refinedCell_to_itsParentCell]
                = grid.currentData().back()->refineSingleCell(cells_per_dim_vec[shiftedLevel], element.index());
            markedElem_to_itsLgr[ element.index() ] = elemLgr_ptr;

            const auto& childrenCount = cells_per_dim_vec[shiftedLevel][0]*cells_per_dim_vec[shiftedLevel][1]*cells_per_dim_vec[shiftedLevel][2];
            std::vector<int> refinedChildrenList(childrenCount);

            for (int refinedCell = 0; refinedCell < childrenCount; ++refinedCell) {

                elemLgrAndElemLgrCell_to_adaptedCell[{element.index(), refinedCell}] = cell_count;
                adaptedCell_to_elemLgrAndElemLgrCell[cell_count] = {element.index(), refinedCell};
                cell_count +=1;

                elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell[{element.index(), refinedCell}] = { markedElemLevel, refined_cell_count_vec[shiftedLevel]};
                refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell.
                    insert_or_assign(std::array{markedElemLevel, refined_cell_count_vec[shiftedLevel]},
                                     std::array{element.index(), refinedCell});
                refinedChildrenList[refinedCell] = refined_cell_count_vec[shiftedLevel];
                refined_cell_count_vec[shiftedLevel] +=1;

            }

            preAdapt_parent_to_children_cells_vec[element.level()][element.getLevelElem().index()] = std::make_pair( markedElemLevel, refinedChildrenList);
            for (const auto& [markedCorner, lgrEquivCorner] : parentCorners_to_equivalentRefinedCorners) {
                cornerInMarkedElemWithEquivRefinedCorner[markedCorner].push_back({element.index(), lgrEquivCorner});
                markedElemAndEquivRefinedCorn_to_corner[ {element.index(), lgrEquivCorner}] = markedCorner;
            }
            for (const auto& [markedFace, itsRefinedFaces] : parentFace_to_itsRefinedFaces) {
                faceInMarkedElemAndRefinedFaces[markedFace].push_back({element.index(), itsRefinedFaces});
            }
        } // end-if-elemMark==1
    } // end-elem-for-loop
}

std::tuple<std::vector<std::vector<std::array<int,2>>>, std::vector<std::vector<int>>, std::vector<std::array<int,2>>, std::vector<int>>
defineChildToParentAndIdxInParentCell(const Dune::CpGrid& grid,
                                      const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                      const std::vector<int>& refined_cell_count_vec,
                                      const std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                                      const int& cell_count)
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.

    // ------------------------ Refined grid parameters
    // Refined child cells and their parents. Entry is {-1,-1} when cell has no father. Otherwise, {level parent cell, parent cell index}
    // Each entry represents a refined level.
    std::vector<std::vector<std::array<int,2>>> refined_child_to_parent_cells_vec(refined_cell_count_vec.size());
    // Each entry represents a refined level. Each refined child cell has a unique index in its parent cell, to be used to build geometryInFather().
    std::vector<std::vector<int>> refined_cell_to_idxInParentCell_vec(refined_cell_count_vec.size());
    // ------------------------ Adapted grid parameters
    // Adapted child cells and their parents. Entry is {-1,-1} when cell has no father. Otherwise, {level parent cell, parent cell index}
    std::vector<std::array<int,2>> adapted_child_to_parent_cells;
    std::vector<int> adapted_cell_to_idxInParentCell;

    adapted_child_to_parent_cells.resize(cell_count, std::array<int,2>{-1,-1}); //  Entry is {-1,-1} when cell has no father, otherwise, {level parent cell, parent cell index}.
    adapted_cell_to_idxInParentCell.resize(cell_count, -1);

    // Rewrite only the entries of adapted cells that have a parent cell
    for (int cell = 0; cell < cell_count; ++cell) {
        // For the element with index "cell" in the adapted grid (or updated-leaf-grid-view),
        // - if "cell" is a new born refined cell (i.e. born in the current adapt-call), we get the parent cell from the preAdapt-leaf-grid-view ("current_view_data_").
        //   In this case, "preAdapt_parent_or_elem" represents "preAdapt_parent".
        // - if "cell" is either a coarse cell or a refined cell that was born in a preAdapt-refined-level-grid, we get the equivalent cell in the
        //   preAdapt-leaf-grid-view ("current_view_data_"). In this case, "preAdapt_parent_or_elem" represents "preAdapt_elem".
        const auto& [elemLgr, elemLgrCell] = adaptedCell_to_elemLgrAndElemLgrCell.at(cell);
        // Get the element of either a parent cell of a new born refined cell with index "cell" or an equivalent cell, in the preAdapt-leaf-grid-view ("current_view_data_").
        const auto& preAdapt_parent_or_elem = Dune::cpgrid::Entity<0>(*grid.currentData().back(), ((elemLgr != -1) ? elemLgr : elemLgrCell), true);
        if (elemLgr != -1) { // "cell" is a new born refined cell
            adapted_child_to_parent_cells[cell] = {preAdapt_parent_or_elem.level(), preAdapt_parent_or_elem.getLevelElem().index()};
            adapted_cell_to_idxInParentCell[cell] = elemLgrCell;
        }
        else {// "cell" is either a coarse cell or a refined cell that was born in a preAdapt-refined-level-grid
            // Only populate the entries of refined cells that were born in preAdapt-refined-level-grids.
            if (preAdapt_parent_or_elem.hasFather()) {
                adapted_child_to_parent_cells[cell] =  {preAdapt_parent_or_elem.father().level(), preAdapt_parent_or_elem.father().index() };
                adapted_cell_to_idxInParentCell[cell] = preAdapt_parent_or_elem.getLevelElem().getIdxInParentCell();
            }
        }
    }

    // Max level before calling adapt.
    const int& preAdaptMaxLevel = grid.maxLevel();
    for (std::size_t shiftedLevel = 0; shiftedLevel < refined_cell_count_vec.size(); ++shiftedLevel) {
        refined_child_to_parent_cells_vec[shiftedLevel].resize(refined_cell_count_vec[shiftedLevel]);
        refined_cell_to_idxInParentCell_vec[shiftedLevel].resize(refined_cell_count_vec[shiftedLevel]);
        const int& level =  shiftedLevel + preAdaptMaxLevel +1;
        // Every refined cell has a parent cell
        for (int cell = 0; cell < refined_cell_count_vec[shiftedLevel]; ++cell) {
            const auto& [elemLgr, elemLgrCell] = refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell.at({level, cell});
            // (elemLgr == -1) means that the adapted cell is equivalent to a cell from the starting grid (current_view_data_)
            assert(elemLgr != -1);
            // Search for the level where the parent cell was born, and its index in that level grid.
            // Notice that elemLgr is a cell index in current_view_data_ (the starting grid where elements got marked for (further) refinement).
            const auto& element = Dune::cpgrid::Entity<0>(*grid.currentData().back(), elemLgr, true);  // elemLgr == parent cell index in starting grid.
            refined_child_to_parent_cells_vec[shiftedLevel][cell] = {element.level(), element.getLevelElem().index()};
            refined_cell_to_idxInParentCell_vec[shiftedLevel][cell] = elemLgrCell;
        }
    }

    return std::make_tuple<std::vector<std::vector<std::array<int,2>>>, std::vector<std::vector<int>>,
                           std::vector<std::array<int,2>>, std::vector<int>>(std::move(refined_child_to_parent_cells_vec),
                                                                             std::move(refined_cell_to_idxInParentCell_vec),
                                                                             std::move(adapted_child_to_parent_cells),
                                                                             std::move(adapted_cell_to_idxInParentCell));

}

std::pair<std::vector<std::vector<int>>, std::vector<std::array<int,2>>>
defineLevelToLeafAndLeafToLevelCells(const Dune::CpGrid& grid,
                                     const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell,
                                     const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                     const std::vector<int>& refined_cell_count_vec,
                                     const std::map<std::array<int,2>,int>& elemLgrAndElemLgrCell_to_adaptedCell,
                                     const std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                                     const int& cell_count)
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.

    // -- Refined to Adapted cells and Adapted-cells to {level where the cell was born, cell index on that level} --
    // Relation between the refined grid and leafview cell indices.
    std::vector<std::vector<int>> refined_level_to_leaf_cells_vec(refined_cell_count_vec.size());
    // Relation between an adapted cell and its equivalent cell coming either from current_view_data_ or from the refined grid (level)
    std::vector<std::array<int,2>> leaf_to_level_cells;
    leaf_to_level_cells.resize(cell_count);

    // Max level before calling adapt.
    const int& preAdaptMaxLevel = grid.maxLevel();

    // -- Adapted to {level, cell index in that level}  --
    for (int cell = 0; cell < cell_count; ++cell) {
        const auto& [elemLgr, elemLgrCell] = adaptedCell_to_elemLgrAndElemLgrCell.at(cell);
        // elemLgr == -1 means that this adapted cell is equivalent to a cell from the starting grid. So we need to find out the level where that equivalent
        // cell was born, as well as its cell index in that level.
        if (elemLgr == -1) {
            const auto& element = Dune::cpgrid::Entity<0>(*grid.currentData().back(), elemLgrCell, true);
            leaf_to_level_cells[cell] = { element.level(), element.getLevelElem().index()};
        }
        else {
            leaf_to_level_cells[cell] = elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell.at({elemLgr, elemLgrCell});
        }
    }
    // -- Refined to adapted cells --
    for (std::size_t shiftedLevel = 0; shiftedLevel < refined_cell_count_vec.size(); ++shiftedLevel){
        refined_level_to_leaf_cells_vec[shiftedLevel].resize(refined_cell_count_vec[shiftedLevel]);
        for (int cell = 0; cell < refined_cell_count_vec[shiftedLevel]; ++cell) {
            refined_level_to_leaf_cells_vec[shiftedLevel][cell] =
                elemLgrAndElemLgrCell_to_adaptedCell.at(refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell.at({static_cast<int>(shiftedLevel) + preAdaptMaxLevel +1, cell}));
        }
    }
    return std::make_pair<std::vector<std::vector<int>>, std::vector<std::array<int,2>>>(std::move(refined_level_to_leaf_cells_vec), std::move(leaf_to_level_cells));
}

void identifyRefinedCornersPerLevel(const Dune::CpGrid& grid,
                                    std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                    std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner,
                                    std::vector<int>& refined_corner_count_vec,
                                    std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                    const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                    const std::vector<int>& assignRefinedLevel,
                                    const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                    const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                    const std::vector<std::array<int,3>>& cells_per_dim_vec)
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.

    // Max level before calling adapt.
    const int& preAdaptMaxLevel = grid.maxLevel();

    // Step 1. Replace the corners from the preAdapt grid involved in LGR by the equivalent ones, born in LGRs.
    //         In this case, we avoid repetition considering the last appearance of the preAdapt corner
    //         in the LGRs.
    for (int corner = 0; corner < grid.currentData().back()->size(3); ++corner) {
        if (!cornerInMarkedElemWithEquivRefinedCorner[corner].empty()) {
            // corner involved in refinement, so we search for it in one LGR (the last one where it appears)
            // Get the lgr corner that replaces the marked corner from level zero.
            // Note: Recall that lgr coincides with the marked element index from the preAdapt grid that got refined.
            //       Since the container is a map, the lgr and the lgr corner index correspond to the last
            //       appearance of the marked corner (from the starting grid - where elements got marked).
            const auto& [lastAppearanceLgr, lastAppearanceLgrCorner] = cornerInMarkedElemWithEquivRefinedCorner[corner].back();

            const auto& lastAppearanceLgrLevel = assignRefinedLevel[lastAppearanceLgr];
            assert(lastAppearanceLgrLevel>0);
            // To access containers with refined level grid information
            const auto& shiftedLevel = lastAppearanceLgrLevel - preAdaptMaxLevel -1;

            elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.
                insert_or_assign(std::array{lastAppearanceLgr, lastAppearanceLgrCorner},
                                 std::array{lastAppearanceLgrLevel, refined_corner_count_vec[shiftedLevel]});
            refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner.
                insert_or_assign(std::array{lastAppearanceLgrLevel, refined_corner_count_vec[shiftedLevel]},
                                 std::array{lastAppearanceLgr, lastAppearanceLgrCorner});
            refined_corner_count_vec[shiftedLevel] +=1;

            if (cornerInMarkedElemWithEquivRefinedCorner[corner].size()>1) {
                for (const auto& [elemLgr, elemLgrCorner] : cornerInMarkedElemWithEquivRefinedCorner[corner]) {
                    const auto& elemLgrLevel = assignRefinedLevel[elemLgr];
                    if (elemLgrLevel != lastAppearanceLgrLevel) {
                        const auto& shiftedElemLgrLevel = elemLgrLevel - preAdaptMaxLevel -1;
                        elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner[{elemLgr, elemLgrCorner}] = {elemLgrLevel, refined_corner_count_vec[shiftedElemLgrLevel]};
                        refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner.
                            insert_or_assign(std::array{elemLgrLevel, refined_corner_count_vec[shiftedElemLgrLevel]},
                                             std::array{elemLgr, elemLgrCorner});
                        refined_corner_count_vec[shiftedElemLgrLevel] +=1;
                    }
                }
            }
        }
    } // end corner-forloop

    for (int elemIdx = 0; elemIdx < grid.currentData().back()->size(0); ++elemIdx) {
        if (markedElem_to_itsLgr.at(elemIdx)!= nullptr) {
            const auto& level = assignRefinedLevel[elemIdx];
            assert(level>0);
            // To access containers with refined level grid information
            const auto& shiftedLevel = level - preAdaptMaxLevel -1;
            for (int corner = 0; corner < markedElem_to_itsLgr.at(elemIdx) ->size(3); ++corner) {
                // Discard marked corners. Store (new born) refined corners

                // INTERIOR
                if (isRefinedCornerInInteriorLgr(cells_per_dim_vec[shiftedLevel], corner)) { // It's a refined interior corner, so we store it.
                    // In this case, the corner is a new born refined corner that does not
                    // coincide with any corner from the GLOBAL grid (level 0). Therefore,
                    // it has to be stored.
                    elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner[{elemIdx, corner}] = {level, refined_corner_count_vec[shiftedLevel]};
                    refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner.
                        insert_or_assign(std::array{level, refined_corner_count_vec[shiftedLevel]},
                                         std::array{elemIdx, corner});
                    refined_corner_count_vec[shiftedLevel] +=1;
                }

                // LYING ON EDGES
                //
                // Refined corners lying on edges - Refined edge has a 'coarse' parent edge (line between 2 corners of the parent cell)
                // To avoid repetition, we distinguish the case where the refined corner lies on an edge of its parent cell.
                // We detect the two coarse faces involved (Notice that the extremes of the parent cell have been stored previously).
                // When the marked faces appears only once, we store the corner now. Otherwise, we store the refined corner on its
                // last appearence associated with one of these parent faces, taking also into account the elemLgr. For example, when
                // the refined corners lie on an edge connecting I_FACE false and K_FACE true of the parent cell, let's say iFaceIdx,
                // kFaceIdx, with each of those faces appearing twice (maximum) :
                // iFaceIdx appearing in current "elem" and elemLgr1
                // kFaceIdx appearing in current "elem" and elemLgr2
                // Then, we take the max(elemLgr1, elemLgr2) and store the refined corner only if this maximum equals elem.
                if (newRefinedCornerLiesOnEdge(cells_per_dim_vec[shiftedLevel], corner)) {
                    const auto& markedFacesTouchingEdge = getParentFacesAssocWithNewRefinedCornLyingOnEdge(grid,
                                                                                                           cells_per_dim_vec[shiftedLevel],
                                                                                                           corner,
                                                                                                           elemIdx);
                    const auto& [markedFace1, markedFace2] = markedFacesTouchingEdge;

                    int lastAppearanceMarkedFace1 = faceInMarkedElemAndRefinedFaces[markedFace1].back().first; // elemLgr1
                    int lastAppearanceMarkedFace2 = faceInMarkedElemAndRefinedFaces[markedFace2].back().first; // elemLgr2

                    int maxLastAppearance = std::max(lastAppearanceMarkedFace1, lastAppearanceMarkedFace2);
                    int faceAtMaxLastAppearance = (maxLastAppearance == lastAppearanceMarkedFace1) ? markedFace1 : markedFace2;

                    // Save the relationship between the vanished refined corner and its last appearance
                    const auto& maxLastAppearanceLevel = assignRefinedLevel[maxLastAppearance];
                    const auto& maxLastAppearanceLevelShifted = assignRefinedLevel[maxLastAppearance] - preAdaptMaxLevel -1;

                    bool atLeastOneFaceAppearsTwice = (faceInMarkedElemAndRefinedFaces[markedFace1].size()>1) ||
                        (faceInMarkedElemAndRefinedFaces[markedFace2].size()>1);
                    if (atLeastOneFaceAppearsTwice && (maxLastAppearance != elemIdx)) {
                        const auto& neighboringLgrCornerIdx = replaceLgr1CornerIdxByLgr2CornerIdx(grid,
                                                                                                  cells_per_dim_vec[shiftedLevel],
                                                                                                  corner, elemIdx, faceAtMaxLastAppearance,
                                                                                                  cells_per_dim_vec[maxLastAppearanceLevelShifted]);
                        vanishedRefinedCorner_to_itsLastAppearance[{elemIdx, corner}] = {maxLastAppearance, neighboringLgrCornerIdx};
                        // Notice that, when we use these container to locate vanished corners, we might need a while-loop,
                        // since {elem, corner} leads to {lastMaxAppearance, neighboringLgrCornerIdx}, which can also vanish.
                        // So we need something like:
                        // if (elemLgrAndElemLgrCorner_to_adapted/refinedCorner.count({elem, corner}) == 0)
                        //    int updateElemLgr =  vanishedRefinedCorner_to_itsLastAppearance[{elem, corner}][0];
                        //    int updateElemLgrCorner =  vanishedRefinedCorner_to_itsLastAppearance[{elem, corner}][1];
                        //     while (elemLgrAndElemLgrCorner_to_adapted/refinedCorner.count({updateElemLgr, updateElemLgCorner}) == 0)
                        //        int tempElemLgr =  updateElemLgr;
                        //        int tempElemLgrCorner =  updateElemLgrCorner;
                        //        updateElemLgr =  vanishedRefinedCorner_to_itsLastAppearance[{ tempElemLgr ,  tempElemLgrCorner}][0];
                        //        updateElemLgrCorner =  vanishedRefinedCorner_to_itsLastAppearance[{ tempElemLgr ,  tempElemLgrCorner}][1];
                        // Then, use the lastest update to search for the corner in teh refined/adapted grid (which would be the one that
                        // gives elemLgrAndElemLgrCorner_to_adapted/refinedCorner.count({updateElemLgr, updateElemLgCorner}) == 1).
                    }
                    if ((maxLastAppearance == elemIdx) || (level!= maxLastAppearanceLevel)) {
                        // Store the refined corner in its last appearence - to avoid repetition.
                        elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.
                            insert_or_assign(std::array{elemIdx, corner},
                                             std::array{level, refined_corner_count_vec[shiftedLevel]});
                        refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner.
                            insert_or_assign(std::array{level,refined_corner_count_vec[shiftedLevel]},
                                             std::array{elemIdx, corner});
                        refined_corner_count_vec[shiftedLevel] +=1;
                    }
                }

                // LYING ON BOUNDARY LGR - NOT ON AN EDGE - NOT COINCIDING WITH A MARKED CORNER
                //
                // If the refined corner lies on the boundary of the LGR, e.i., it was born on one of the faces
                // of the marked element that got refined, then, we have two cases:
                // - the marked face appears only in one marked element -> then, we store this corner now.
                // - the marked face appears twice (maximum times) in two marked elements -> we store it later.
                if ( isRefinedNewBornCornerOnLgrBoundary(cells_per_dim_vec[shiftedLevel], corner) &&
                     !newRefinedCornerLiesOnEdge(cells_per_dim_vec[shiftedLevel], corner)) {
                    // Get the index of the marked face where the refined corner was born.
                    const auto& markedFace = getParentFaceWhereNewRefinedCornerLiesOn(grid, cells_per_dim_vec[shiftedLevel], corner, elemIdx);
                    // check how many times marked face appearn
                    // Get the last LGR (marked element) where the marked face appeared.
                    int lastLgrWhereMarkedFaceAppeared = faceInMarkedElemAndRefinedFaces[markedFace].back().first;

                    const auto& lastLgrLevel = assignRefinedLevel[lastLgrWhereMarkedFaceAppeared ];
                    const auto& lastLgrLevelShifted = lastLgrLevel - preAdaptMaxLevel -1;
                    // Save the relationship between the vanished refined corner and its last appearance
                    if ((faceInMarkedElemAndRefinedFaces[markedFace].size()>1) && (lastLgrWhereMarkedFaceAppeared != elemIdx)) {
                        const auto& neighboringLgrCornerIdx = replaceLgr1CornerIdxByLgr2CornerIdx(cells_per_dim_vec[shiftedLevel], corner,
                                                                                                  cells_per_dim_vec[lastLgrLevelShifted]);
                        vanishedRefinedCorner_to_itsLastAppearance[{elemIdx, corner}] = {lastLgrWhereMarkedFaceAppeared, neighboringLgrCornerIdx};
                    }

                    if ((lastLgrWhereMarkedFaceAppeared == elemIdx) || (lastLgrLevel != level)) {
                        // Store the refined corner in its last appearence - to avoid repetition.
                        elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.
                            insert_or_assign(std::array{elemIdx, corner},
                                             std::array{level, refined_corner_count_vec[shiftedLevel]});
                        refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner.
                            insert_or_assign(std::array{level, refined_corner_count_vec[shiftedLevel]},
                                             std::array{elemIdx, corner});
                        refined_corner_count_vec[shiftedLevel] +=1;
                    }
                }
            } // end-corner-for-loop
        } // end-if-nullptr
    } // end-elem-for-loop
}


bool isRefinedCornerInInteriorLgr(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr)
{
    assert(cells_per_dim[0]>0);
    assert(cells_per_dim[1]>0);
    assert(cells_per_dim[2]>0);

    const auto& ijk = getRefinedCornerIJK(cells_per_dim, cornerIdxInLgr);
    return ((ijk[0]%cells_per_dim[0] > 0) &&  (ijk[1]%cells_per_dim[1]>0) && (ijk[2]%cells_per_dim[2]>0));
}



std::array<int,3> getRefinedCornerIJK(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr)
{
    const auto& total_corners = (cells_per_dim[0] +1)*(cells_per_dim[1]+1)*(cells_per_dim[2]+1);
    if (cornerIdxInLgr >= total_corners) {
        OPM_THROW(std::logic_error, "Invalid corner index from single-cell-refinement.\n");
    }
    // Order defined in Geometry::refine
    //  (j*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (i*(cells_per_dim[2]+1)) + k
    std::array<int,3> ijk;
    ijk[2] = cornerIdxInLgr % (cells_per_dim[2] +1);
    cornerIdxInLgr -= ijk[2];
    cornerIdxInLgr /= (cells_per_dim[2] +1);
    ijk[0] = cornerIdxInLgr % (cells_per_dim[0]+1);
    cornerIdxInLgr -=ijk[0];
    ijk[1] = cornerIdxInLgr / (cells_per_dim[0]+1);
    return ijk;
}

bool newRefinedCornerLiesOnEdge(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr)
{
    const auto& ijk = getRefinedCornerIJK(cells_per_dim, cornerIdxInLgr);
    // Edges laying on bottom face
    bool isNewBornOnEdge01 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == 0) && (ijk[2] == 0);
    bool isNewBornOnEdge23 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == cells_per_dim[1]) && ( ijk[2] == 0);
    bool isNewBornOnEdge02 = (ijk[0] == 0) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == 0);
    bool isNewBornOnEdge13 = (ijk[0] == cells_per_dim[0]) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == 0);

    // Edges connecting bottom and top faces
    bool isNewBornOnEdge04 = (ijk[0] == 0) && (ijk[1] == 0) && (ijk[2] % cells_per_dim[2] != 0);
    bool isNewBornOnEdge26 = (ijk[0] == 0) && (ijk[1] == cells_per_dim[1]) && (ijk[2] % cells_per_dim[2] != 0);
    bool isNewBornOnEdge15 = (ijk[0] == cells_per_dim[0]) && (ijk[1] == 0) && (ijk[2] % cells_per_dim[2] != 0);
    bool isNewBornOnEdge37 = (ijk[0] == cells_per_dim[0]) && (ijk[1] == cells_per_dim[1]) && (ijk[2] % cells_per_dim[2] != 0);

    // Edges laying on top face
    bool isNewBornOnEdge45 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == 0) && (ijk[2] == cells_per_dim[2]);
    bool isNewBornOnEdge67 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == cells_per_dim[1]) && ( ijk[2] == cells_per_dim[2]);
    bool isNewBornOnEdge46 = (ijk[0] == 0) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == cells_per_dim[2]);
    bool isNewBornOnEdge57 = (ijk[0] == cells_per_dim[0]) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == cells_per_dim[2]);

    bool isOnEdge = isNewBornOnEdge01 || isNewBornOnEdge23 || isNewBornOnEdge02 || isNewBornOnEdge13 ||
        isNewBornOnEdge04 || isNewBornOnEdge26 || isNewBornOnEdge15 || isNewBornOnEdge37 ||
        isNewBornOnEdge45 || isNewBornOnEdge67 || isNewBornOnEdge46 || isNewBornOnEdge57;

    return isOnEdge;
}

std::array<int,2> getParentFacesAssocWithNewRefinedCornLyingOnEdge(const Dune::CpGrid& grid,
                                                                   const std::array<int,3>& cells_per_dim,
                                                                   int cornerIdxInLgr,
                                                                   int elemLgr)
{
    assert(newRefinedCornerLiesOnEdge(cells_per_dim, cornerIdxInLgr));

    const auto& parentCell_to_face = grid.cellFaceRow(elemLgr);
    if(parentCell_to_face.size()>6){
        OPM_THROW(std::logic_error, "The associted parent cell has more than six faces. Refinment/Adaptivity not supported yet.");
    }
    // Corners Order defined in Geometry::refine  (j*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (i*(cells_per_dim[2]+1)) + k
    const auto& ijk = getRefinedCornerIJK(cells_per_dim, cornerIdxInLgr);
    // Edges laying on bottom face
    bool isNewBornOnEdge01 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == 0) && (ijk[2] == 0);
    bool isNewBornOnEdge23 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == cells_per_dim[1]) && ( ijk[2] == 0);
    bool isNewBornOnEdge02 = (ijk[0] == 0) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == 0);
    bool isNewBornOnEdge13 = (ijk[0] == cells_per_dim[0]) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == 0);

    // Edges connecting bottom and top faces
    bool isNewBornOnEdge04 = (ijk[0] == 0) && (ijk[1] == 0) && (ijk[2] % cells_per_dim[2] != 0);
    bool isNewBornOnEdge26 = (ijk[0] == 0) && (ijk[1] == cells_per_dim[1]) && (ijk[2] % cells_per_dim[2] != 0);
    bool isNewBornOnEdge15 = (ijk[0] == cells_per_dim[0]) && (ijk[1] == 0) && (ijk[2] % cells_per_dim[2] != 0);
    bool isNewBornOnEdge37 = (ijk[0] == cells_per_dim[0]) && (ijk[1] == cells_per_dim[1]) && (ijk[2] % cells_per_dim[2] != 0);

    // Edges laying on top face
    bool isNewBornOnEdge45 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == 0) && (ijk[2] == cells_per_dim[2]);
    bool isNewBornOnEdge67 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == cells_per_dim[1]) && ( ijk[2] == cells_per_dim[2]);
    bool isNewBornOnEdge46 = (ijk[0] == 0) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == cells_per_dim[2]);
    bool isNewBornOnEdge57 = (ijk[0] == cells_per_dim[0]) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == cells_per_dim[2]);

    std::vector<int> auxFaces;
    auxFaces.reserve(2);

    for (const auto& face : parentCell_to_face) {
        const auto& faceTag = grid.currentData().back()->faceTag(face.index());
        // Add I_FACE false
        bool addIfalse = isNewBornOnEdge02 || isNewBornOnEdge04 || isNewBornOnEdge26 || isNewBornOnEdge46;
        if( addIfalse && (faceTag == 0)  && (!face.orientation()))  {
            auxFaces.push_back(face.index());
        }
        // Add J_FACE false
        bool addJfalse = isNewBornOnEdge01 || isNewBornOnEdge04 || isNewBornOnEdge15 || isNewBornOnEdge45;
        if( addJfalse && (faceTag == 1)  && (!face.orientation()))  {
            auxFaces.push_back(face.index());
        }
        // Add K_FACE false
        bool addKfalse = isNewBornOnEdge01 ||  isNewBornOnEdge13 || isNewBornOnEdge23 || isNewBornOnEdge02;
        if( addKfalse && (faceTag == 2) && (!face.orientation())) {
            auxFaces.push_back(face.index());
        }
        // Add I_FACE true
        bool addItrue = isNewBornOnEdge13 || isNewBornOnEdge15 || isNewBornOnEdge37 || isNewBornOnEdge57;
        if( addItrue && (faceTag == 0)  && (face.orientation()))  {
            auxFaces.push_back(face.index());
        }
        // Add J_FACE true
        bool addJtrue = isNewBornOnEdge23|| isNewBornOnEdge26 || isNewBornOnEdge37 || isNewBornOnEdge67;
        if( addJtrue && (faceTag == 1)  && (face.orientation()))  {
            auxFaces.push_back(face.index());
        }
        // Add K_FACE true
        bool addKtrue = isNewBornOnEdge45 || isNewBornOnEdge67 || isNewBornOnEdge46 || isNewBornOnEdge57;
        if(addKtrue && (faceTag == 2) && (face.orientation())) {
            auxFaces.push_back(face.index());
        }
    }
    return {auxFaces[0], auxFaces[1]};
}

bool isRefinedNewBornCornerOnLgrBoundary(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr)
{
    const auto& ijk = getRefinedCornerIJK(cells_per_dim, cornerIdxInLgr);
    bool isOnParentCell_I_FACEfalse_and_newBornCorn = ( (ijk[0] == 0) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_I_FACEtrue_and_newBornCorn = ( (ijk[0] == cells_per_dim[0]) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_J_FACEfalse_and_newBornCorn = ( (ijk[1] == 0) && ((ijk[0] % cells_per_dim[0] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_J_FACEtrue_and_newBornCorn = ( (ijk[1] == cells_per_dim[1]) && ((ijk[0] % cells_per_dim[0] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_K_FACEfalse_and_newBornCorn = ( (ijk[2] == 0) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[0] % cells_per_dim[0] !=0) ));
    bool isOnParentCell_K_FACEtrue_and_newBornCorn = ( (ijk[2] == cells_per_dim[2]) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[0] % cells_per_dim[0] !=0) ));
    bool isOnParentCell_I_FACE = isOnParentCell_I_FACEfalse_and_newBornCorn || isOnParentCell_I_FACEtrue_and_newBornCorn;
    bool isOnParentCell_J_FACE = isOnParentCell_J_FACEfalse_and_newBornCorn || isOnParentCell_J_FACEtrue_and_newBornCorn;
    bool isOnParentCell_K_FACE = isOnParentCell_K_FACEfalse_and_newBornCorn || isOnParentCell_K_FACEtrue_and_newBornCorn;
    return (isOnParentCell_I_FACE || isOnParentCell_J_FACE || isOnParentCell_K_FACE);
}

int getParentFaceWhereNewRefinedCornerLiesOn(const Dune::CpGrid& grid,
                                             const std::array<int,3>& cells_per_dim,
                                             int cornerIdxInLgr, int elemLgr)
{
    assert(isRefinedNewBornCornerOnLgrBoundary(cells_per_dim, cornerIdxInLgr));

    const auto& parentCell_to_face = grid.cellFaceRow(elemLgr);
    if(parentCell_to_face.size()>6){
        OPM_THROW(std::logic_error, "The associted parent cell has more than six faces. Refinment/Adaptivity not supported yet.");
    }
    const auto& ijk = getRefinedCornerIJK(cells_per_dim, cornerIdxInLgr);

    bool isOnParentCell_I_FACEfalse_and_newBornCorn = ( (ijk[0] == 0) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_I_FACEtrue_and_newBornCorn = ( (ijk[0] == cells_per_dim[0]) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_J_FACEfalse_and_newBornCorn = ( (ijk[1] == 0) && ((ijk[0] % cells_per_dim[0] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_J_FACEtrue_and_newBornCorn = ( (ijk[1] == cells_per_dim[1]) && ((ijk[0] % cells_per_dim[0] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_K_FACEfalse_and_newBornCorn = ( (ijk[2] == 0) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[0] % cells_per_dim[0] !=0) ));
    bool isOnParentCell_K_FACEtrue_and_newBornCorn = ( (ijk[2] == cells_per_dim[2]) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[0] % cells_per_dim[0] !=0) ));

    for (const auto& face : parentCell_to_face) {
        const auto& faceTag = grid.currentData().back()->faceTag(face.index());
        if (isOnParentCell_I_FACEfalse_and_newBornCorn && (faceTag == 0) && !face.orientation()) { // I_FACE false
            return face.index();
        }
        if (isOnParentCell_I_FACEtrue_and_newBornCorn && (faceTag == 0) && face.orientation()) { // I_FACE true
            return face.index();
        }
        if (isOnParentCell_J_FACEfalse_and_newBornCorn && (faceTag == 1) && !face.orientation()) { // J_FACE false
            return face.index();
        }
        if (isOnParentCell_J_FACEtrue_and_newBornCorn && (faceTag == 1) && face.orientation()) { // J_FACE true
            return face.index();
        }
        if (isOnParentCell_K_FACEfalse_and_newBornCorn && (faceTag == 2) && !face.orientation()) { // K_FACE false
            return face.index();
        }
        if (isOnParentCell_K_FACEtrue_and_newBornCorn && (faceTag == 2) && face.orientation()) { // K_FACE true
            return face.index();
        }
    }
    OPM_THROW(std::logic_error, "Cannot find parent face index where new refined corner lays on.");
}


int replaceLgr1CornerIdxByLgr2CornerIdx(const std::array<int,3>& cells_per_dim_lgr1,
                                        int cornerIdxLgr1,
                                        const std::array<int,3>& cells_per_dim_lgr2)
{
    const auto& ijkLgr1 = getRefinedCornerIJK(cells_per_dim_lgr1, cornerIdxLgr1);
    // Order defined in Geometry::refine
    // (j*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (i*(cells_per_dim[2]+1)) + k

    // On a parallel run, no symmetry between neighboring elements should be assumed. Therefore, all the six cases
    // (i = 0, cells_per_dim[0], j = 0, cells_per_dim[1], and k = 0, cells_per_dim[2]) have to be taken into account.
    // On a serial run, it would be enough to consider i = cells_per_dim[0], j = cells_per_dim[1], and k = cells_per_dim[2].
    // To cover all possible scenarios, serial and parallel, we consider the six cases.

    if (ijkLgr1[0] == cells_per_dim_lgr1[0]) { // same j, k, but i = 0
        return   (ijkLgr1[1]*(cells_per_dim_lgr2[0]+1)*(cells_per_dim_lgr2[2]+1)) + ijkLgr1[2];
    }
    if (ijkLgr1[1] == cells_per_dim_lgr1[1]) { // same i,k, but j = 0
        return  (ijkLgr1[0]*(cells_per_dim_lgr2[2]+1)) + ijkLgr1[2];
    }
    if (ijkLgr1[2] == cells_per_dim_lgr1[2]) { // same i,j, but k = 0
        return  (ijkLgr1[1]*(cells_per_dim_lgr2[0]+1)*(cells_per_dim_lgr2[2]+1)) + (ijkLgr1[0]*(cells_per_dim_lgr2[2]+1));
    }
    if (ijkLgr1[0] == 0) { // same j,k, but i = cells_per_dim[0]
        return   (ijkLgr1[1]*(cells_per_dim_lgr2[0]+1)*(cells_per_dim_lgr2[2]+1)) + (cells_per_dim_lgr2[0]*(cells_per_dim_lgr2[2]+1))+ ijkLgr1[2];
    }
    if (ijkLgr1[1] == 0) { // same i,k, but j = cells_per_Dim[1]
        return  (cells_per_dim_lgr2[1]*(cells_per_dim_lgr2[0]+1)*(cells_per_dim_lgr2[2]+1)) + (ijkLgr1[0]*(cells_per_dim_lgr2[2]+1)) + ijkLgr1[2];
    }
    if (ijkLgr1[2] == 0) { // same i,j, but k = cells_per_dim[2]
        return  (ijkLgr1[1]*(cells_per_dim_lgr2[0]+1)*(cells_per_dim_lgr2[2]+1)) + (ijkLgr1[0]*(cells_per_dim_lgr2[2]+1)) + cells_per_dim_lgr2[2];
    }
    else {
        const auto& message = "Cannot convert corner index from one LGR to its neighboring LGR.";
        OPM_THROW(std::logic_error, message);
    }
}


int replaceLgr1CornerIdxByLgr2CornerIdx(const Dune::CpGrid& grid,
                                        const std::array<int,3>& cells_per_dim_lgr1,
                                        int cornerIdxLgr1,
                                        int elemLgr1,
                                        int parentFaceLastAppearanceIdx,
                                        const std::array<int,3>& cells_per_dim_lgr2)
{
    assert(newRefinedCornerLiesOnEdge(cells_per_dim_lgr1, cornerIdxLgr1));
    const auto& faces = getParentFacesAssocWithNewRefinedCornLyingOnEdge(grid, cells_per_dim_lgr1, cornerIdxLgr1, elemLgr1);
    assert( (faces[0] == parentFaceLastAppearanceIdx) || (faces[1] == parentFaceLastAppearanceIdx));

    const auto& ijkLgr1 = getRefinedCornerIJK(cells_per_dim_lgr1, cornerIdxLgr1);
    const auto& parentCell_to_face = grid.cellFaceRow(elemLgr1);

    if(parentCell_to_face.size()>6){
        const auto& message = "The associated parent cell has more than six faces. Refinement/Adaptivity not supported yet.";
        if (grid.comm().rank() == 0){
            OPM_THROW(std::logic_error, message);
        }
        else{
            OPM_THROW_NOLOG(std::logic_error, message);
        }
    }

    // Order defined in Geometry::refine
    //  (j*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (i*(cells_per_dim[2]+1)) + k

    for (const auto& face : parentCell_to_face) {
        const auto& faceTag = grid.currentData().back()->faceTag(face.index());
        if (parentFaceLastAppearanceIdx == face.index()) {
            if ( face.orientation() ){
                if (faceTag == 0) { // I_FACE true. The same new born refined corner will have equal j and k, but i == 0.
                    return  (ijkLgr1[1]*(cells_per_dim_lgr2[0]+1)*(cells_per_dim_lgr2[2]+1))  + ijkLgr1[2];
                }
                if (faceTag == 1) {// J_FACE true. The same new born refined corner will have equal i and k, but j == 0.
                    return  (ijkLgr1[0]*(cells_per_dim_lgr2[2]+1)) + ijkLgr1[2];
                }
                if (faceTag == 2) {// K_FACE true. The same new born refined corner will have equal  i and j, but k == 0.
                    return  (ijkLgr1[1]*(cells_per_dim_lgr2[0]+1)*(cells_per_dim_lgr2[2]+1)) + (ijkLgr1[0]*(cells_per_dim_lgr2[2]+1));
                }
            }
            if(!face.orientation()) {
                if (faceTag == 0) {// I_FACE false. The same new born refined corner will have equal values of j and k, but i == cells_per_dim[0].
                    return  (ijkLgr1[1]*(cells_per_dim_lgr2[0]+1)*(cells_per_dim_lgr2[2]+1)) +
                        (cells_per_dim_lgr2[0]*(cells_per_dim_lgr2[2]+1)) +ijkLgr1[2];
                }
                if (faceTag == 1) {// J_FACE false. The same new born refined corner will have equal  i and k, but j == cells_per_dim[1].
                    return   (cells_per_dim_lgr2[1]*(cells_per_dim_lgr2[0]+1)*(cells_per_dim_lgr2[2]+1))
                        + (ijkLgr1[0]*(cells_per_dim_lgr2[2]+1)) + ijkLgr1[2];
                }
                if (faceTag == 2) {// K_FACE false.  The same new born refined corner will have equal  i and j, but k == cells_per_dim[2].
                    return  (ijkLgr1[1]*(cells_per_dim_lgr2[0]+1)*(cells_per_dim_lgr2[2]+1)) + (ijkLgr1[0]*(cells_per_dim_lgr2[2]+1))
                        + cells_per_dim_lgr2[2];
                }
            }
        }
    }
    const auto& message = "Cannot convert corner index from one LGR to its neighboring LGR.";
    if (grid.comm().rank() == 0){
        OPM_THROW(std::logic_error, message);
    }
    else{
        OPM_THROW_NOLOG(std::logic_error, message);
    }
}


void identifyLeafGridCorners(const Dune::CpGrid& grid,
                             std::map<std::array<int,2>,int>& elemLgrAndElemLgrCorner_to_adaptedCorner,
                             std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner,
                             int& corner_count,
                             const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                             const std::vector<int>& assignRefinedLevel,
                             const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                             std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                             const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                             const std::vector<std::array<int,3>>& cells_per_dim_vec)
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.
    
    // Step 1. Select/store the corners from the starting grid not involved in any (new) LGR.
    //         Replace the corners from level zero involved in LGR by the equivalent ones, born in LGRs.
    //         In this case, we avoid repetition considering the last appearance of the level zero corner
    //         in the LGRs.
    for (int corner = 0; corner < grid.currentData().back()->size(3); ++corner) {
        if (cornerInMarkedElemWithEquivRefinedCorner[corner].empty()) { // save it
            // Note: Since we are associating each LGR with its parent cell index, and this index can take
            //       the value 0, we will represent the grid before adapt is called with the value -1.
            elemLgrAndElemLgrCorner_to_adaptedCorner[{-1, corner}] = corner_count;
            adaptedCorner_to_elemLgrAndElemLgrCorner[corner_count] = {-1, corner};
            corner_count +=1;
        }
        else { // corner involved in refinement, so we search it in one LGR (the last one where it appears)
            assert(!cornerInMarkedElemWithEquivRefinedCorner[corner].empty());
            // Get the lgr corner that replaces the marked corner from level zero.
            // Note: Recall that lgr coincides with the marked element index from level 0 that got refined.
            //       Since the container is a map, the lgr and the lgr corner index correspond to the last
            //       appearance of the marked corner (from level 0).
            const auto& [lastAppearanceLgr, lastAppearanceLgrCorner] = cornerInMarkedElemWithEquivRefinedCorner[corner].back();
            // Build the relationships between adapted corner and level corner, for future search due topology aspects.
            elemLgrAndElemLgrCorner_to_adaptedCorner[{lastAppearanceLgr, lastAppearanceLgrCorner}] = corner_count;
            adaptedCorner_to_elemLgrAndElemLgrCorner[corner_count] = {lastAppearanceLgr, lastAppearanceLgrCorner};
            corner_count +=1;
        }
    } // end corner-forloop

    // Max level before calling adapt.
    const int& preAdaptMaxLevel = grid.maxLevel();

    for (int elemIdx = 0; elemIdx < grid.currentData().back()->size(0); ++elemIdx) {
        if (markedElem_to_itsLgr[elemIdx]!= nullptr) {
            const auto& level = assignRefinedLevel[elemIdx];
            assert(level>0);
            // To access containers with refined level grid information
            const auto& shiftedLevel = level - preAdaptMaxLevel -1;
            for (int corner = 0; corner < markedElem_to_itsLgr[elemIdx] ->size(3); ++corner) {
                // Discard marked corners. Store (new born) refined corners

                // INTERIOR
                if (isRefinedCornerInInteriorLgr(cells_per_dim_vec[shiftedLevel], corner)) { // It's a refined interior corner, so we store it.
                    // In this case, the corner is a new born refined corner that does not
                    // coincide with any corner from the GLOBAL grid (level 0). Therefore,
                    // it has to be stored.
                    elemLgrAndElemLgrCorner_to_adaptedCorner[{elemIdx, corner}] = corner_count;
                    adaptedCorner_to_elemLgrAndElemLgrCorner[corner_count] = {elemIdx, corner};
                    corner_count += 1;
                }

                // LYING ON EDGES
                //
                // Refined corners lying on edges - Refined edge has a 'coarse' parent edge (line between 2 corners of the parent cell)
                // To avoid repetition, we distinguish the case where the refined corner lies on an edge of its parent cell.
                // We detect the two coarse faces involved (Notice that the extremes of the parent cell have been stored previously).
                // When the marked faces appeares only once, we store the corner now. Otherwise, we store the refined corner on its
                // last apparance assocaited to one of these parent faces, taking also into account the elemLgr. For example, when
                // the refined corners lies on an edge connecting I_FACE false and K_FACE true of the parent cell, let's say iFaceIdx,
                // kFaceIdx, with each of those faces appearing twice (maximum) :
                // iFaceIdx appearing in current "elem" and elemLgr1
                // kFaceIdx appearing in current "elem" and elemLgr2
                // Then, we take the max(elemLgr1, elemLgr2) and store the refined corner only if this maximum equals elem.
                if (newRefinedCornerLiesOnEdge(cells_per_dim_vec[shiftedLevel], corner)) {
                    const auto& markedFacesTouchingEdge = getParentFacesAssocWithNewRefinedCornLyingOnEdge(grid,
                                                                                                           cells_per_dim_vec[shiftedLevel],
                                                                                                           corner,
                                                                                                           elemIdx);
                    const auto& [markedFace1, markedFace2] = markedFacesTouchingEdge;

                    int lastAppearanceMarkedFace1 = faceInMarkedElemAndRefinedFaces[markedFace1].back().first; // elemLgr1
                    int lastAppearanceMarkedFace2 = faceInMarkedElemAndRefinedFaces[markedFace2].back().first; // elemLgr2

                    int maxLastAppearance = std::max(lastAppearanceMarkedFace1, lastAppearanceMarkedFace2); // max(elemLgr1, elemLgr2)
                    int faceAtMaxLastAppearance = (maxLastAppearance == lastAppearanceMarkedFace1) ? markedFace1 : markedFace2;

                    int maxLastAppearanceLevel = assignRefinedLevel[maxLastAppearance];


                    // Save the relationship between the vanished refined corner and its last appearance
                    bool atLeastOneFaceAppearsTwice = (faceInMarkedElemAndRefinedFaces[markedFace1].size()>1) ||
                        (faceInMarkedElemAndRefinedFaces[markedFace2].size()>1);
                    if ((atLeastOneFaceAppearsTwice && (maxLastAppearance != elemIdx)) || (maxLastAppearanceLevel != level)) {
                        const int maxLastAppearanceLevelShifted = assignRefinedLevel[maxLastAppearance] -
                                                                  preAdaptMaxLevel - 1;
                        const auto& neighboringLgrCornerIdx = replaceLgr1CornerIdxByLgr2CornerIdx(grid,
                                                                                                  cells_per_dim_vec[shiftedLevel],
                                                                                                  corner, elemIdx, faceAtMaxLastAppearance,
                                                                                                  cells_per_dim_vec[maxLastAppearanceLevelShifted]);
                        vanishedRefinedCorner_to_itsLastAppearance[{elemIdx, corner}] = {maxLastAppearance, neighboringLgrCornerIdx};

                        // Notice that, when we use this container to locate vanished corners, we might need a while-loop,
                        // since {elem, corner} leads to {lastMaxAppearance, neighboringLgrCornerIdx}, which can also vanish.
                        // So we need something like:
                        // if (elemLgrAndElemLgrCorner_to_adapted/refinedCorner.count({elem, corner}) == 0)
                        //    int updateElemLgr =  vanishedRefinedCorner_to_itsLastAppearance[{elem, corner}][0];
                        //    int updateElemLgrCorner =  vanishedRefinedCorner_to_itsLastAppearance[{elem, corner}][1];
                        //     while (elemLgrAndElemLgrCorner_to_adapted/refinedCorner.count({updateElemLgr, updateElemLgCorner}) == 0)
                        //        int tempElemLgr =  updateElemLgr;
                        //        int tempElemLgrCorner =  updateElemLgrCorner;
                        //        updateElemLgr =  vanishedRefinedCorner_to_itsLastAppearance[{ tempElemLgr ,  tempElemLgrCorner}][0];
                        //        updateElemLgrCorner =  vanishedRefinedCorner_to_itsLastAppearance[{ tempElemLgr ,  tempElemLgrCorner}][1];
                        // Then, use the lastest update to search for the corner in teh refined/adapted grid (which would be the one that
                        // gives elemLgrAndElemLgrCorner_to_adapted/refinedCorner.count({updateElemLgr, updateElemLgCorner}) == 1).
                    }
                    if (maxLastAppearance == elemIdx) {
                        // Store the refined corner in its last appearence - to avoid repetition.
                        elemLgrAndElemLgrCorner_to_adaptedCorner[{elemIdx, corner}] = corner_count;
                        adaptedCorner_to_elemLgrAndElemLgrCorner[corner_count] = {elemIdx, corner};
                        corner_count += 1;
                    }
                }

                // LYING ON BOUNDARY LGR - NOT ON AN EDGE - NOT COINCIDING WITH A MARKED CORNER
                //
                // If the refined corner lies on the boundary of the LGR, e.i., it was born on one of the faces
                // of the marked element that got refined, then, we have two cases:
                // - the marked face appears only in one marked element -> then, we store this corner now.
                // - the marked face appears twice (maximum times) in two marked elements -> we store it later.
                if ( isRefinedNewBornCornerOnLgrBoundary(cells_per_dim_vec[shiftedLevel], corner) &&
                     !newRefinedCornerLiesOnEdge(cells_per_dim_vec[shiftedLevel], corner)) {
                    // Get the index of the marked face where the refined corner was born.
                    const auto& markedFace = getParentFaceWhereNewRefinedCornerLiesOn(grid,
                                                                                      cells_per_dim_vec[shiftedLevel], corner, elemIdx);
                    // check how many times marked face appearn
                    // Get the last LGR (marked element) where the marked face appeared.
                    int lastLgrWhereMarkedFaceAppeared = faceInMarkedElemAndRefinedFaces[markedFace].back().first;

                    const auto& lastLgrLevel = assignRefinedLevel[lastLgrWhereMarkedFaceAppeared ];
                    const auto& lastLgrLevelShifted = lastLgrLevel - preAdaptMaxLevel -1;

                    // Save the relationship between the vanished refined corner and its last appearance
                    if ((faceInMarkedElemAndRefinedFaces[markedFace].size()>1) && (lastLgrWhereMarkedFaceAppeared != elemIdx)) {
                        const auto& neighboringLgrCornerIdx = replaceLgr1CornerIdxByLgr2CornerIdx(cells_per_dim_vec[shiftedLevel], corner,
                                                                                                  cells_per_dim_vec[lastLgrLevelShifted]);
                        vanishedRefinedCorner_to_itsLastAppearance[{elemIdx, corner}] = {lastLgrWhereMarkedFaceAppeared, neighboringLgrCornerIdx};
                    }

                    if (lastLgrWhereMarkedFaceAppeared == elemIdx) {
                        // Store the refined corner in its last appearence - to avoid repetition.
                        elemLgrAndElemLgrCorner_to_adaptedCorner[{elemIdx, corner}] = corner_count;
                        adaptedCorner_to_elemLgrAndElemLgrCorner[corner_count] = {elemIdx, corner};
                        corner_count += 1;
                    }
                }
            } // end-corner-for-loop
        } // end-if-nullptr
    } // end-elem-for-loop
}

void identifyRefinedFacesPerLevel(const Dune::CpGrid& grid,
                                  std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                                  std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                                  std::vector<int>& refined_face_count_vec,
                                  const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                  const std::vector<int>& assignRefinedLevel,
                                  const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                  const std::vector<std::array<int,3>>& cells_per_dim_vec)
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.

    // Max level before calling adapt.
    const int& preAdaptMaxLevel = grid.maxLevel();

    // Step 1. Add the LGR faces, for each LGR
    for (int elem = 0; elem < grid.currentData().back()->size(0); ++elem) {
        if (markedElem_to_itsLgr[elem]!=nullptr)  {
            const auto& level = assignRefinedLevel[elem];
            assert(level>0);
            // To access containers with refined level grid information
            const auto& shiftedLevel = level - preAdaptMaxLevel -1;
            for (int face = 0; face < markedElem_to_itsLgr[elem] ->numFaces(); ++face) {
                // Discard marked faces. Store (new born) refined faces
                bool isNewRefinedFaceOnLgrBoundary = isRefinedFaceOnLgrBoundary(cells_per_dim_vec[shiftedLevel], face, markedElem_to_itsLgr[elem]);
                if (!isNewRefinedFaceOnLgrBoundary) { // It's a refined interior face, so we store it
                    // In this case, the face is a new born refined face that does not
                    // have any "parent face" from the GLOBAL grid (level 0).
                    elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace[{elem, face}] = {level, refined_face_count_vec[shiftedLevel]};
                    refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace.
                        insert_or_assign(std::array{level, refined_face_count_vec[shiftedLevel]},
                                         std::array{elem, face});
                    refined_face_count_vec[shiftedLevel] +=1;
                }
                // If the refined face lays on the boundary of the LGR, e.i., it was born on one of the faces
                // of the marked element that got refined, then, we have two cases:
                // - the marked face appears only in one marked element -> then, we store this face now.
                // - the marked face appears twice (maximum times) in two marked elements -> we store it later.
                else {
                    // Get the index of the marked face where the refined corner was born.
                    int markedFace = getParentFaceWhereNewRefinedFaceLiesOn(grid,
                                                                            cells_per_dim_vec[shiftedLevel],
                                                                            face,
                                                                            markedElem_to_itsLgr[elem],
                                                                            elem);
                    assert(!faceInMarkedElemAndRefinedFaces[markedFace].empty());
                    // Get the last LGR (marked element) where the marked face appeared.
                    int lastLgrWhereMarkedFaceAppeared = faceInMarkedElemAndRefinedFaces[markedFace].back().first;
                    if (lastLgrWhereMarkedFaceAppeared == elem) {
                        // Store the refined face in its last appearence - to avoid repetition.
                        elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace.
                            insert_or_assign(std::array{elem, face},
                                             std::array{level, refined_face_count_vec[shiftedLevel]});
                        refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace.
                            insert_or_assign(std::array{level, refined_face_count_vec[shiftedLevel]},
                                             std::array{elem, face});
                        refined_face_count_vec[shiftedLevel] +=1;
                    }
                    if(faceInMarkedElemAndRefinedFaces[markedFace].size()>1) { // maximum size is 2
                        const auto& firstMarkedElem = faceInMarkedElemAndRefinedFaces[markedFace][0].first;
                        const auto& firstMarkedElemLevel = assignRefinedLevel[firstMarkedElem];
                        if (firstMarkedElemLevel != level) {
                            const auto& shiftedFirstMarkedElemLevel = firstMarkedElemLevel - preAdaptMaxLevel -1;
                            elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace.
                                insert_or_assign(std::array{firstMarkedElem, face},
                                                 std::array{firstMarkedElemLevel,
                                                            refined_face_count_vec[shiftedFirstMarkedElemLevel]});
                            refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace.
                                insert_or_assign(std::array{firstMarkedElemLevel,
                                                            refined_face_count_vec[shiftedFirstMarkedElemLevel]},
                                    std::array{firstMarkedElem, face});
                            refined_face_count_vec[shiftedFirstMarkedElemLevel] +=1;
                        }
                    }
                }
            } // end-face-for-loop
        } // end-if-nullptr
    } // end-elem-for-loop
}

void identifyLeafGridFaces(const Dune::CpGrid& grid,
                           std::map<std::array<int,2>,int>& elemLgrAndElemLgrFace_to_adaptedFace,
                           std::unordered_map<int,std::array<int,2>>& adaptedFace_to_elemLgrAndElemLgrFace,
                           int& face_count,
                           const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                           const std::vector<int>& assignRefinedLevel,
                           const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                           const std::vector<std::array<int,3>>& cells_per_dim_vec)
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.

    // Max level before calling adapt.
    const int& preAdaptMaxLevel = grid.maxLevel(); 
    
    // Step 1. Add the LGR faces, for each LGR
    for (int elem = 0; elem < grid.currentData().back()->size(0); ++elem) {
        if (markedElem_to_itsLgr[elem]!=nullptr)  {
            const auto& level = assignRefinedLevel[elem];
            assert(level>0);
            // To access containers with refined level grid information
            const auto& shiftedLevel = level - preAdaptMaxLevel -1;
            for (int face = 0; face < markedElem_to_itsLgr[elem] -> numFaces(); ++face) {
                // Discard marked faces. Store (new born) refined faces
                bool isNewRefinedFaceOnLgrBoundary = isRefinedFaceOnLgrBoundary(cells_per_dim_vec[shiftedLevel], face, markedElem_to_itsLgr[elem]);
                if (!isNewRefinedFaceOnLgrBoundary) { // It's a refined interior face, so we store it
                    //  if (isRefinedFaceInInteriorLgr(cells_per_dim, face, markedElem_to_itsLgr[elem])) {
                    // In this case, the face is a new born refined face that does not
                    // have any "parent face" from the GLOBAL grid (level 0).
                    elemLgrAndElemLgrFace_to_adaptedFace[{elem, face}] = face_count;
                    adaptedFace_to_elemLgrAndElemLgrFace[face_count] = {elem, face};
                    face_count += 1;
                }
                // If the refined face lays on the boundary of the LGR, e.i., it was born on one of the faces
                // of the marked element that got refined, then, we have two cases:
                // - the marked face appears only in one marked element -> then, we store this face now.
                // - the marked face appears twice (maximum times) in two marked elements -> we store it later.
                else {
                    // Get the index of the marked face where the refined corner was born.
                    int markedFace = getParentFaceWhereNewRefinedFaceLiesOn(grid,
                                                                            cells_per_dim_vec[shiftedLevel],
                                                                            face,
                                                                            markedElem_to_itsLgr[elem],
                                                                            elem);
                    assert(!faceInMarkedElemAndRefinedFaces[markedFace].empty());
                    // Get the last LGR (marked element) where the marked face appeared.
                    int lastLgrWhereMarkedFaceAppeared = faceInMarkedElemAndRefinedFaces[markedFace].back().first;
                    if (lastLgrWhereMarkedFaceAppeared == elem) {
                        // Store the refined face in its last appearence - to avoid repetition.
                        elemLgrAndElemLgrFace_to_adaptedFace[{elem, face}] = face_count;
                        adaptedFace_to_elemLgrAndElemLgrFace[face_count] = {elem, face};
                        face_count += 1;
                    }
                }
            } // end-face-for-loop
        } // end-if-nullptr
    } // end-elem-for-loop
    // Step 2. Select/store the faces from the starting grid (where cells got marked) not involved in any LGR.
    //         Replace the faces from level zero involved in LGR by the equivalent ones, born in LGRs.
    //         In this case, we avoid repetition considering the last appearance of the level zero corner
    //         in the LGRs.
    for (int face = 0; face < grid.currentData().back()->numFaces(); ++face) {
        if (faceInMarkedElemAndRefinedFaces[face].empty()) { // save it
            // Note: Since we are associating each LGR with its parent cell index, and this index can take
            //       the value 0, we will represent the current_view_data_ with the value -1
            elemLgrAndElemLgrFace_to_adaptedFace[{-1, face}] = face_count;
            adaptedFace_to_elemLgrAndElemLgrFace[face_count] = {-1, face};
            face_count +=1;
        }
    } // end face-forloop
}

std::array<int,3> getRefinedFaceIJK(const std::array<int,3>& cells_per_dim,
                                    int faceIdxInLgr,
                                    const std::shared_ptr<Dune::cpgrid::CpGridData>& elemLgr_ptr)
{
    // Order defined in Geometry::refine
    // K_FACES  (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i
    // I_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1))
    //           + (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j
    // J_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2] +1))
    //                    + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2])
    //                    + (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2]) + k
    const auto& i_faces =  (cells_per_dim[0] +1)*cells_per_dim[1]*cells_per_dim[2];
    const auto& j_faces =  cells_per_dim[0]*(cells_per_dim[1]+1)*cells_per_dim[2];
    const auto& k_faces =  cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1);

    if (faceIdxInLgr >= i_faces + j_faces + k_faces) {
        OPM_THROW(std::logic_error, "Invalid face index from single-cell-refinement.\n");
    }
    
    const auto& faceTag =  elemLgr_ptr ->faceTag(faceIdxInLgr);
    std::array<int,3> ijk;
    switch (faceTag) {
    case I_FACE:
        faceIdxInLgr -= (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1));
        // faceIdxInLgr =  (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j
        ijk[1] = faceIdxInLgr % cells_per_dim[1];
        faceIdxInLgr -= ijk[1]; // (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1])
        faceIdxInLgr /= cells_per_dim[1]; // (i*cells_per_dim[2]) + k
        ijk[2] = faceIdxInLgr % cells_per_dim[2];
        faceIdxInLgr -=ijk[2]; // i*cells_per_dim[2]
        ijk[0] = faceIdxInLgr / cells_per_dim[2];
        break;
    case J_FACE:
        faceIdxInLgr -=  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2] +1))
            + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2]);
        // faceIdxInLgr =  (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2]) + k
        ijk[2] = faceIdxInLgr % cells_per_dim[2];
        faceIdxInLgr -= ijk[2]; // (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2])
        faceIdxInLgr /= cells_per_dim[2]; // (j*cells_per_dim[0]) + i
        ijk[0] = faceIdxInLgr % cells_per_dim[0];
        faceIdxInLgr -=ijk[0]; // j*cells_per_dim[0]
        ijk[1] = faceIdxInLgr / cells_per_dim[0];
        break;
    case K_FACE:
        //  (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i
        ijk[0] = faceIdxInLgr % cells_per_dim[0];
        faceIdxInLgr -= ijk[0]; // (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0])
        faceIdxInLgr /= cells_per_dim[0]; // (k*cells_per_dim[1]) + j
        ijk[1] = faceIdxInLgr % cells_per_dim[1];
        faceIdxInLgr -=ijk[1]; // k*cells_per_dim[1]
        ijk[2] = faceIdxInLgr / cells_per_dim[1];
        break;
    default:
        OPM_THROW(std::logic_error, "FaceTag is not I, J, or K!");
    }
    return ijk;
}

bool isRefinedFaceInInteriorLgr(const std::array<int,3>& cells_per_dim, int faceIdxInLgr, const std::shared_ptr<Dune::cpgrid::CpGridData>& elemLgr_ptr)
{

    int refined_k_faces = cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1);
    int refined_i_faces = (cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2];

    bool isKface = (faceIdxInLgr < refined_k_faces);
    bool isIface = (faceIdxInLgr >= refined_k_faces) && (faceIdxInLgr < refined_k_faces + refined_i_faces);
    bool isJface = (faceIdxInLgr >= refined_k_faces + refined_i_faces);

    const auto& ijk = getRefinedFaceIJK(cells_per_dim, faceIdxInLgr, elemLgr_ptr);
    return ((ijk[0]%cells_per_dim[0] > 0 && isIface) ||  (ijk[1]%cells_per_dim[1]>0 && isJface) || (ijk[2]%cells_per_dim[2]>0 && isKface));
}


bool isRefinedFaceOnLgrBoundary(const std::array<int,3>& cells_per_dim, int faceIdxInLgr,
                                const std::shared_ptr<Dune::cpgrid::CpGridData>& elemLgr_ptr)
{
    const auto& ijk = getRefinedFaceIJK(cells_per_dim, faceIdxInLgr, elemLgr_ptr);

    int refined_k_faces = cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1);
    int refined_i_faces = (cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2];

    bool isKface = (faceIdxInLgr < refined_k_faces);
    bool isIface = (faceIdxInLgr >= refined_k_faces) && (faceIdxInLgr < refined_k_faces + refined_i_faces);
    bool isJface = (faceIdxInLgr >= refined_k_faces + refined_i_faces);

    bool isOnParentCell_I_FACE = isIface && (ijk[0] % cells_per_dim[0] == 0) && (ijk[1]<cells_per_dim[1]) && (ijk[2]<cells_per_dim[2]);
    bool isOnParentCell_J_FACE = isJface && (ijk[1] % cells_per_dim[1] == 0) && (ijk[0]<cells_per_dim[0]) && (ijk[2]<cells_per_dim[2]);
    bool isOnParentCell_K_FACE = isKface && (ijk[2] % cells_per_dim[2] == 0) && (ijk[0]<cells_per_dim[0]) && (ijk[1]<cells_per_dim[1]);

    return (isOnParentCell_I_FACE || isOnParentCell_J_FACE || isOnParentCell_K_FACE);
}

int getParentFaceWhereNewRefinedFaceLiesOn(const Dune::CpGrid& grid,
                                           const std::array<int,3>& cells_per_dim,
                                           int faceIdxInLgr,
                                           const std::shared_ptr<Dune::cpgrid::CpGridData>& elemLgr_ptr,
                                           int elemLgr)
{
    assert(isRefinedFaceOnLgrBoundary(cells_per_dim, faceIdxInLgr, elemLgr_ptr));
    const auto& ijk = getRefinedFaceIJK(cells_per_dim, faceIdxInLgr, elemLgr_ptr);
    const auto& parentCell_to_face = grid.cellFaceRow(elemLgr);
    // cell_to_face_ [ element ] = { I false, I true, J false, J true, K false, K true } if current_view_data_ is level zero

    if(parentCell_to_face.size()>6){
        const auto& message = "The associated parent cell has more than six faces. Refinement/Adaptivity not supported yet.";
        if (grid.comm().rank() == 0){
            OPM_THROW(std::logic_error, message);
        }
        else{
            OPM_THROW_NOLOG(std::logic_error, message);
        }
    }

    // Order defined in Geometry::refine (to be used for distinguishing if faceIdxInLgr is K, I, or J face)
    //
    // K_FACES  (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i
    // I_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1))
    //           + (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j
    // J_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2] +1))
    //                    + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2])
    //                    + (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2]) + k
    int refined_k_faces = cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1);
    int refined_i_faces = (cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2];
    int refined_j_faces = cells_per_dim[0]*(cells_per_dim[1]+1)*cells_per_dim[2];

    assert( faceIdxInLgr < refined_k_faces + refined_i_faces + refined_j_faces);

    for (const auto& face : parentCell_to_face) {
        const auto& faceTag =  grid.currentData().back()->faceTag(face.index());
        if (faceIdxInLgr <  refined_k_faces ) { // It's a K_FACE
            if ((ijk[2] == 0) && (faceTag == 2) && !face.orientation()) { // {K_FACE, false}
                return face.index();
            }
            if ((ijk[2] == cells_per_dim[2]) && (faceTag == 2) && face.orientation()) { // {K_FACE, true}
                return face.index();
            }
        }
        if ((faceIdxInLgr >= refined_k_faces) && (faceIdxInLgr < refined_k_faces + refined_i_faces)) { // It's I_FACE
            if ((ijk[0] == 0) && (faceTag == 0) && !face.orientation()) { // {I_FACE, false}
                return face.index();
            }
            if ((ijk[0] == cells_per_dim[0]) && (faceTag == 0) && face.orientation()) { // {I_FACE, true}
                return face.index();
            }
        }
        if (faceIdxInLgr >= refined_k_faces + refined_i_faces) {// It's J_FACE
            if ((ijk[1] == 0) && (faceTag == 1) && !face.orientation()) { // {J_FACE, false}
                return face.index();
            }
            if ((ijk[1] == cells_per_dim[1]) && (faceTag == 1) && face.orientation()) { // {J_FACE, true}
                return face.index();
            }
        }
    }
    const auto& message = "Cannot find index of parent face where the new refined face lies on.";
    if (grid.comm().rank() == 0){
        OPM_THROW(std::logic_error, message);
    }
    else{
        OPM_THROW_NOLOG(std::logic_error, message);
    }
}

void populateRefinedCorners(std::vector<Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>>& refined_corners_vec,
                            const std::vector<int>& refined_corner_count_vec,
                            const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                            const int& preAdaptMaxLevel,
                            const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner)
{
    for (std::size_t shiftedLevel = 0; shiftedLevel < refined_corner_count_vec.size(); ++shiftedLevel) {
        refined_corners_vec[shiftedLevel].resize(refined_corner_count_vec[shiftedLevel]);
        for (int corner = 0; corner < refined_corner_count_vec[shiftedLevel]; ++corner) {
            const auto& [elemLgr, elemLgrCorner] = refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner.at({static_cast<int>(shiftedLevel) + preAdaptMaxLevel +1,corner});
            refined_corners_vec[shiftedLevel][corner] =  markedElem_to_itsLgr[elemLgr] -> getGeometry().geomVector(std::integral_constant<int,3>()) -> get(elemLgrCorner);
        }
    }
}


void populateRefinedFaces(std::vector<Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<2,3>>>& refined_faces_vec,
                          std::vector<Dune::cpgrid::EntityVariableBase<enum face_tag>>& mutable_refined_face_tags_vec,
                          std::vector<Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>>& mutable_refined_face_normals_vec,
                          std::vector<Opm::SparseTable<int>>& refined_face_to_point_vec,
                          const std::vector<int>& refined_face_count_vec,
                          const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                          const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                          const std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                          const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                          const int& preAdaptMaxLevel,
                          const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                          const std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner)
{
    for (std::size_t shiftedLevel = 0; shiftedLevel < refined_face_count_vec.size(); ++shiftedLevel) {

        // Store the refined faces
        refined_faces_vec[shiftedLevel].resize(refined_face_count_vec[shiftedLevel]);
        mutable_refined_face_tags_vec[shiftedLevel].resize(refined_face_count_vec[shiftedLevel]);
        mutable_refined_face_normals_vec[shiftedLevel].resize(refined_face_count_vec[shiftedLevel]);

        // Auxiliary integer to count all the points in refined_face_to_point.
        int refined_num_points = 0;
        // Auxiliary vector to store refined_face_to_point with non consecutive indices.
        std::vector<std::vector<int>> aux_refined_face_to_point;
        aux_refined_face_to_point.resize(refined_face_count_vec[shiftedLevel]);
        for (int face = 0; face < refined_face_count_vec[shiftedLevel]; ++face) {

            const auto& [elemLgr, elemLgrFace] = refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace.at({static_cast<int>(shiftedLevel) + preAdaptMaxLevel +1,face});
            const auto& elemLgrFaceEntity =  Dune::cpgrid::EntityRep<1>(elemLgrFace, true);

            // Get the face geometry.
            refined_faces_vec[shiftedLevel][face] = (*(markedElem_to_itsLgr.at(elemLgr)->getGeometry().geomVector(std::integral_constant<int,1>())))[elemLgrFaceEntity];
            // Get the face tag.
            mutable_refined_face_tags_vec[shiftedLevel][face] = markedElem_to_itsLgr.at(elemLgr)->faceTag(elemLgrFace);
            // Get the face normal.
            mutable_refined_face_normals_vec[shiftedLevel][face] = markedElem_to_itsLgr.at(elemLgr)->faceNormals(elemLgrFace);
            // Get face_to_point_ before adapting - we need to replace the level corners by the adapted ones.
            //Opm::SparseTable<int>::mutable_row_type
            const auto& preAdapt_face_to_point = markedElem_to_itsLgr.at(elemLgr)->faceToPoint(elemLgrFace);
            // Add the amount of points to the count num_points.
            refined_num_points += preAdapt_face_to_point.size();

            // Face_to_point
            for (std::size_t corn = 0; corn < preAdapt_face_to_point.size(); ++corn) {
                const auto& elemLgrCorn = preAdapt_face_to_point[corn];
                std::size_t refinedCorn = 0; // It'll be rewritten.
                // Corner is stored in adapted_corners
                if (auto candidate = elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.find({elemLgr, elemLgrCorn});
                    candidate != elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.end())  {
                    refinedCorn = candidate->second[1];
                }
                else{
                    // Corner might have vanished - Search its equivalent lgr-corner in that case -
                    // last lgr where the corner appears -
                    std::array<int,2> lastAppearanceLgr_lgrEquivCorner = {0, 0}; // It'll get rewritten.
                    if(auto corner_candidate = markedElemAndEquivRefinedCorn_to_corner.find({elemLgr, elemLgrCorn});
                       corner_candidate != markedElemAndEquivRefinedCorn_to_corner.end()) {
                        lastAppearanceLgr_lgrEquivCorner = cornerInMarkedElemWithEquivRefinedCorner[corner_candidate->second].back();
                    }
                    else {
                        // To locate vanished corners, we need a while-loop, since {elemLgr, elemLgrcorner} leads to
                        // {neighboringElemLgr, neighboringElemLgrCornerIdx}, which might have also vanished.
                        // Then, use the lastest appearance of the current corner, meaning, the first (and unique one - by construction) that
                        // gives elemLgrAndElemLgrCorner_to_refinedCorner.count(lastAppearanceLgr_lgrCorner) == 1).
                        // This corner lies on the area occupied by a coarse face that got refined and belonged to two marked elements.
                        // Get the index of this corner with respect to the greatest marked element index, using find instead of count.
                        lastAppearanceLgr_lgrEquivCorner = vanishedRefinedCorner_to_itsLastAppearance.at({elemLgr, elemLgrCorn});
                        while (elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.find(lastAppearanceLgr_lgrEquivCorner) ==
                               elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.end()) {
                            const auto& tempLgr_lgrCorner  = lastAppearanceLgr_lgrEquivCorner;
                            lastAppearanceLgr_lgrEquivCorner =  vanishedRefinedCorner_to_itsLastAppearance.at(tempLgr_lgrCorner);
                        }
                    }
                    refinedCorn =  elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.at(lastAppearanceLgr_lgrEquivCorner)[1];
                }
                aux_refined_face_to_point[face].push_back(refinedCorn);
            }
        }
        // Refined face_to_point.
        refined_face_to_point_vec[shiftedLevel].reserve(refined_face_count_vec[shiftedLevel], refined_num_points);
        for (int face = 0; face < refined_face_count_vec[shiftedLevel]; ++face) {
            refined_face_to_point_vec[shiftedLevel].appendRow(aux_refined_face_to_point[face].begin(), aux_refined_face_to_point[face].end());
        }
    }
}

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
                          const std::vector<std::array<int,3>>&  cells_per_dim_vec)
{   
    // --- Refined cells ---
    for (std::size_t shiftedLevel = 0; shiftedLevel < refined_cell_count_vec.size(); ++shiftedLevel) {

        refined_cells_vec[shiftedLevel].resize(refined_cell_count_vec[shiftedLevel]);
        refined_cell_to_point_vec[shiftedLevel].resize(refined_cell_count_vec[shiftedLevel]);
        refined_global_cell_vec[shiftedLevel].resize(refined_cell_count_vec[shiftedLevel]);

        const auto& allLevelCorners = refined_geometries_vec[shiftedLevel].geomVector(std::integral_constant<int,3>());

        for (int cell = 0; cell < refined_cell_count_vec[shiftedLevel]; ++cell) {

            const auto& [elemLgr, elemLgrCell] = refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell.at({static_cast<int>(shiftedLevel) + preAdaptMaxLevel +1, cell});
            assert(elemLgr >-1);
            const auto& elemLgrCellEntity = Dune::cpgrid::EntityRep<0>(elemLgrCell, true);
            // Auxiliary cell_to_face
            std::vector<Dune::cpgrid::EntityRep<1>> aux_refined_cell_to_face;

            // global_cell_ values of refined cells from level grid created via adapt() (i.e., without parameters startIJK, endIJK
            // delimiting a block of parent cells) inherit the global_cell_ value of its parent cell.
            refined_global_cell_vec[shiftedLevel][cell] = grid.currentData().back()->globalCell()[elemLgr];
            
            // Get pre-adapt corners of the cell that will be replaced with leaf view ones.
            const auto& preAdapt_cell_to_point = markedElem_to_itsLgr.at(elemLgr)->cellToPoint(elemLgrCell);
            // Get pre-adapt faces of the cell that will be replaced with leaf view ones.
            const auto& preAdapt_cell_to_face = markedElem_to_itsLgr.at(elemLgr)->cellToFace(elemLgrCell);

            // Cell to point.
            for (int corn = 0; corn < 8; ++corn) {
                int refinedCorn;
                const auto& preAdaptCorn = preAdapt_cell_to_point[corn];
                auto refined_candidate = elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.find({elemLgr, preAdaptCorn});
                if (refined_candidate == elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.end()) {
                    // Corner might have vanished - Search its equivalent lgr-corner in that case -
                    // last lgr where the corner appears -
                    std::array<int,2> lastAppearanceLgr_lgrCorner = {0, 0}; // It'll be rewritten
                    if(auto candidate = markedElemAndEquivRefinedCorn_to_corner.find({elemLgr, preAdaptCorn}); candidate != markedElemAndEquivRefinedCorn_to_corner.end()) {
                        lastAppearanceLgr_lgrCorner = cornerInMarkedElemWithEquivRefinedCorner[candidate->second].back();
                    }
                    else {
                        // To locate vanished corners, we need a while-loop, since {elemLgr, elemLgrcorner} leads to
                        // {neighboringElemLgr, neighboringElemLgrCornerIdx}, which might have also vanished.
                        // Then, use the lastest appearance of the current corner, meaning, the first (and unique one - by construction) that
                        // gives elemLgrAndElemLgrCorner_to_adaptedCorner.count( lastAppearanceLgr_lgrCorner ) == 1).
                        // This corner lies on the area occupied by a coarse face that got refined and belonged to two marked elements.
                        // Get the index of this corner with respect to the greatest marked element index.
                        lastAppearanceLgr_lgrCorner = vanishedRefinedCorner_to_itsLastAppearance.at({elemLgr, preAdaptCorn});
                        while (elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.find(lastAppearanceLgr_lgrCorner) == elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.end()) {
                            const auto& tempLgr_lgrCorner =   lastAppearanceLgr_lgrCorner;
                            lastAppearanceLgr_lgrCorner =  vanishedRefinedCorner_to_itsLastAppearance.at(tempLgr_lgrCorner);
                        }
                    }
                    refinedCorn = elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.at(lastAppearanceLgr_lgrCorner)[1];
                }
                // Corner is stored in adapted_corners
                else {
                    refinedCorn = refined_candidate->second[1];
                }
                refined_cell_to_point_vec[shiftedLevel][cell][corn] = refinedCorn;
            } // end-cell_to_point

            // Cell to face.
            for (const auto& face : preAdapt_cell_to_face) {
                const auto& preAdaptFace = face.index();
                int refinedFace = 0; // It'll be rewritten.
                // Face might have vanished - Search its refined lgr-children faces in that case -
                // last lgr where the face appears
                auto face_candidate = elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace.find({elemLgr, preAdaptFace});
                if (face_candidate == elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace.end()) {
                    // Get the index of the marked face where the refined face was born.
                    const auto& markedFace = getParentFaceWhereNewRefinedFaceLiesOn(grid,
                                                                                    cells_per_dim_vec[shiftedLevel],
                                                                                    preAdaptFace,
                                                                                    markedElem_to_itsLgr[elemLgr], elemLgr);
                    // Get the last LGR (marked element) where the marked face appeared.
                    const int& lastLgrWhereMarkedFaceAppeared = faceInMarkedElemAndRefinedFaces[markedFace].back().first;
                    const auto& lastAppearanceLgrEquivFace = replaceLgr1FaceIdxByLgr2FaceIdx(cells_per_dim_vec[shiftedLevel],
                                                                                             preAdaptFace,
                                                                                             markedElem_to_itsLgr[elemLgr],
                                                                                             cells_per_dim_vec[assignRefinedLevel[lastLgrWhereMarkedFaceAppeared] - preAdaptMaxLevel -1]);
                    refinedFace = elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace.at({lastLgrWhereMarkedFaceAppeared, lastAppearanceLgrEquivFace})[1];
                    aux_refined_cell_to_face.push_back({refinedFace, face.orientation()});
                }
                // Face is stored in adapted_faces
                else {
                    refinedFace = face_candidate->second[1];
                    aux_refined_cell_to_face.push_back({refinedFace, face.orientation()});
                }
            } // end-cell_to_face
            // Refined cell to face.
            refined_cell_to_face_vec[shiftedLevel].appendRow(aux_refined_cell_to_face.begin(), aux_refined_cell_to_face.end());
            // Get the cell geometry.
            const auto& elemLgrGeom =  (*( markedElem_to_itsLgr.at(elemLgr)->getGeometry().geomVector(std::integral_constant<int,0>())))[elemLgrCellEntity];

            // Create a pointer to the first element of "refined_cell_to_point" (required as the fourth argement to construct a Geometry<3,3> type object).
            int* indices_storage_ptr = refined_cell_to_point_vec[shiftedLevel][cell].data();
            refined_cells_vec[shiftedLevel][cell] = Dune::cpgrid::Geometry<3,3>(elemLgrGeom.center(), elemLgrGeom.volume(), allLevelCorners, indices_storage_ptr);
        } // refined_cells
        // Refined face to cell.
        refined_cell_to_face_vec[shiftedLevel].makeInverseRelation(refined_face_to_cell_vec[shiftedLevel]);
    } // end-shiftedLevel-for-loop
}

int replaceLgr1FaceIdxByLgr2FaceIdx(const std::array<int,3>& cells_per_dim_lgr1,
                                    int faceIdxInLgr1,
                                    const std::shared_ptr<Dune::cpgrid::CpGridData>& elemLgr1_ptr,
                                    const std::array<int,3>& cells_per_dim_lgr2)
{
    const auto& ijkLgr1 = getRefinedFaceIJK(cells_per_dim_lgr1, faceIdxInLgr1, elemLgr1_ptr);
    // lgr1 represents an element index < lgr2 (neighboring cells sharing a face with lgr1-element)
    // Order defined in Geometry::refine
    // K_FACES (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i
    // I_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1))
    //           + (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j
    // J_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2] +1))
    //                    + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2])
    //                    + (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2]) + k
    const int& kFacesLgr2 = cells_per_dim_lgr2[0]*cells_per_dim_lgr2[1]*(cells_per_dim_lgr2[2]+1);
    const int& iFacesLgr2 = ((cells_per_dim_lgr2[0]+1)*cells_per_dim_lgr2[1]*cells_per_dim_lgr2[2]);

    const auto& face_tag = elemLgr1_ptr-> faceTag(faceIdxInLgr1);

    if (face_tag == I_FACE) {
        assert( cells_per_dim_lgr1[1] == cells_per_dim_lgr2[1]);
        assert( cells_per_dim_lgr1[2] == cells_per_dim_lgr2[2]);
        if (ijkLgr1[0] == cells_per_dim_lgr1[0]) { // same j,k, but i = 0
            return  kFacesLgr2 + (ijkLgr1[2]*cells_per_dim_lgr2[1]) + ijkLgr1[1];
        }
        else { // same j,k, but i = cells_per_dim[0]
            return  kFacesLgr2 + (cells_per_dim_lgr2[0]*cells_per_dim_lgr2[1]*cells_per_dim_lgr2[2])
                + (ijkLgr1[2]*cells_per_dim_lgr2[1]) + ijkLgr1[1];
        }
    }
    if (face_tag == J_FACE) {
        assert( cells_per_dim_lgr1[0] == cells_per_dim_lgr2[0]);
        assert( cells_per_dim_lgr1[2] == cells_per_dim_lgr2[2]);
        if (ijkLgr1[1] == cells_per_dim_lgr1[1]) { // same i,k, but j = 0
            return kFacesLgr2 + iFacesLgr2 + (ijkLgr1[0]*cells_per_dim_lgr2[2]) + ijkLgr1[2];
        }
        else { // same i,k, but j = cells_per_dim[1]
            return kFacesLgr2 + iFacesLgr2 + (cells_per_dim_lgr2[1]*cells_per_dim_lgr2[0]*cells_per_dim_lgr2[2])
                + (ijkLgr1[0]*cells_per_dim_lgr2[2]) + ijkLgr1[2];
        }
    }
    if (face_tag == K_FACE) {
        assert( cells_per_dim_lgr1[0] == cells_per_dim_lgr2[0]);
        assert( cells_per_dim_lgr1[1] == cells_per_dim_lgr2[1]);
        if (ijkLgr1[2] == cells_per_dim_lgr1[2]) { // same i,j, but k = 0
            return  (ijkLgr1[1]*cells_per_dim_lgr2[0]) + ijkLgr1[0];
        }
        else{ // same i, j, but k = cells_per_dim[2]
            return  (cells_per_dim_lgr2[2]*cells_per_dim_lgr2[0]*cells_per_dim_lgr2[1])
                + (ijkLgr1[1]*cells_per_dim_lgr2[0]) + ijkLgr1[0];
        }
    }
    OPM_THROW(std::logic_error,  "Cannot convert face index from one LGR to its neighboring LGR.");
}

void setRefinedLevelGridsGeometries(const Dune::CpGrid& grid,
                                    /* Refined corner arguments */
                                    std::vector<Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>>& refined_corners_vec,
                                    const std::vector<int>& refined_corner_count_vec,
                                    /* Refined face arguments */
                                    std::vector<Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<2,3>>>& refined_faces_vec,
                                    std::vector<Dune::cpgrid::EntityVariableBase<enum face_tag>>& mutable_refined_face_tags_vec,
                                    std::vector<Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>>& mutable_refined_face_normals_vec,
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
                                    const std::vector<std::array<int,3>>&  cells_per_dim_vec)
{
    // --- Refined corners  ---
    populateRefinedCorners(refined_corners_vec,
                           refined_corner_count_vec,
                           markedElem_to_itsLgr,
                           preAdaptMaxLevel,
                           refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner);
    // --- Refined faces  ---
    populateRefinedFaces(refined_faces_vec,
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
    populateRefinedCells(grid,
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

}

void populateLeafGridCorners(const Dune::CpGrid& grid,
                             Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& adapted_corners,
                             const int& corner_count,
                             const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                             const std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner)
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.
    adapted_corners.resize(corner_count);
    for (int corner = 0; corner < corner_count; ++corner) {
        const auto& [elemLgr, elemLgrCorner] = adaptedCorner_to_elemLgrAndElemLgrCorner.at(corner);
        // Note: Since we are associating each LGR with its parent cell index, and this index can take
        //       the value 0, we represent the current_view_data_ with the value -1
        adapted_corners[corner] = ((elemLgr == -1) ?
                                   *grid.currentData().back() :
                                   *markedElem_to_itsLgr[elemLgr]).getGeometry().geomVector(std::integral_constant<int,3>())-> get(elemLgrCorner);
    }
}

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
                           const int& preAdaptMaxLevel)
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.
    
    adapted_faces.resize(face_count);
    mutable_face_tags.resize(face_count);
    mutable_face_normals.resize(face_count);

    // Auxiliary integer to count all the points in leaf_face_to_point.
    int num_points = 0;
    // Auxiliary vector to store face_to_point with non consecutive indices.
    std::vector<std::vector<int>> aux_face_to_point;
    aux_face_to_point.resize(face_count);
    for (int face = 0; face < face_count; ++face) {

        const auto& [elemLgr, elemLgrFace] = adaptedFace_to_elemLgrAndElemLgrFace.at(face);
        // Note: Since we are associating each LGR with its parent cell index, and this index can take
        //       the value 0, with the value -1 we refer to current_view_data_ (grid where the elements have been marked).
        const auto& elemLgrFaceEntity =  Dune::cpgrid::EntityRep<1>(elemLgrFace, true);
        const auto& grid_or_elemLgr_data = (elemLgr == -1) ? *grid.currentData().back() : *markedElem_to_itsLgr[elemLgr];

        // Get the face geometry.
        adapted_faces[face] = (*(grid_or_elemLgr_data.getGeometry().geomVector(std::integral_constant<int,1>())))[elemLgrFaceEntity];
        // Get the face tag.
        mutable_face_tags[face] = grid_or_elemLgr_data.faceTag(elemLgrFace);
        // Get the face normal.
        mutable_face_normals[face] = grid_or_elemLgr_data.faceNormals(elemLgrFace);
        // Get face_to_point_ before adapting - we need to replace the level corners by the adapted ones.
        const auto& preAdapt_face_to_point = grid_or_elemLgr_data.faceToPoint(elemLgrFace);
        // Add the amount of points to the count num_points.
        num_points += preAdapt_face_to_point.size();

        for (std::size_t corn = 0; corn < preAdapt_face_to_point.size(); ++corn) {
            std::size_t adaptedCorn = 0; // It'll get rewritten.
            const auto& elemLgrCorn = preAdapt_face_to_point[corn];
            // Corner is stored in adapted_corners
            if (auto candidate = elemLgrAndElemLgrCorner_to_adaptedCorner.find({elemLgr, elemLgrCorn}); candidate !=  elemLgrAndElemLgrCorner_to_adaptedCorner.end()) {
                adaptedCorn = candidate->second;
            }
            else{
                // Corner might have vanished - Search its equivalent lgr-corner in that case -
                // last lgr where the corner appears -
                std::array<int,2> lastAppearanceLgr_lgrEquivCorner = {0, 0}; // It'll get rewritten.
                if ((elemLgr ==-1) && (!cornerInMarkedElemWithEquivRefinedCorner[elemLgrCorn].empty())) {
                    lastAppearanceLgr_lgrEquivCorner = cornerInMarkedElemWithEquivRefinedCorner[elemLgrCorn].back();
                }
                if (elemLgr>-1) { // It represents the refinement of the element with index "elemIdx == elemLgr"
                    if (auto corner_candidate = markedElemAndEquivRefinedCorn_to_corner.find({elemLgr, elemLgrCorn}); corner_candidate != markedElemAndEquivRefinedCorn_to_corner.end()) {
                        lastAppearanceLgr_lgrEquivCorner =  cornerInMarkedElemWithEquivRefinedCorner[corner_candidate->second].back();
                    }
                    else {
                        const auto& shiftedLevel = assignRefinedLevel[elemLgr] - preAdaptMaxLevel -1; // Assigned level > preAdapt maxLevel
                        bool isNewRefinedCornInInteriorLgr = isRefinedCornerInInteriorLgr(cells_per_dim_vec[shiftedLevel], elemLgrCorn);
                        assert(!isNewRefinedCornInInteriorLgr);
                        // To locate vanished corners, we need a while-loop, since {elemLgr, elemLgrcorner} leads to
                        // {neighboringElemLgr, neighboringElemLgrCornerIdx}, which might have also vanished.
                        // Then, use the lastest appearance of the current corner, meaning, the first (and unique one - by construction) that
                        // gives elemLgrAndElemLgrCorner_to_adaptedCorner.count(lastAppearanceLgr_lgrEquivCorner) == 1).
                        // This corner lies on the area occupied by a coarse face that got refined and belonged to two marked elements.
                        // Get the index of this corner with respect to the greatest marked element index, using find instead of count.
                        lastAppearanceLgr_lgrEquivCorner = vanishedRefinedCorner_to_itsLastAppearance.at({elemLgr, elemLgrCorn});
                        while (elemLgrAndElemLgrCorner_to_adaptedCorner.find( lastAppearanceLgr_lgrEquivCorner ) == elemLgrAndElemLgrCorner_to_adaptedCorner.end()) {
                            const auto& tempLgr_lgrCorner = lastAppearanceLgr_lgrEquivCorner;
                            lastAppearanceLgr_lgrEquivCorner =  vanishedRefinedCorner_to_itsLastAppearance.at(tempLgr_lgrCorner);
                        }
                    }
                } // end-if(elemLgr>-1)
                adaptedCorn =   elemLgrAndElemLgrCorner_to_adaptedCorner.at(lastAppearanceLgr_lgrEquivCorner);
            }
            aux_face_to_point[face].push_back(adaptedCorn);
        }
    } // end-adapted-faces
    // Adapted/Leaf-Grid-View face_to_point.
    adapted_face_to_point.reserve(face_count, num_points);
    for (int face = 0; face < face_count; ++face) {
        adapted_face_to_point.appendRow(aux_face_to_point[face].begin(), aux_face_to_point[face].end());
    }
}

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
                           const int& preAdaptMaxLevel)
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.
    
    // --- Adapted cells ---
    // Store the adapted cells. Main difficulty: to lookup correctly the indices of the corners and faces of each cell.
    adapted_cells.resize(cell_count);
    adapted_cell_to_point.resize(cell_count);
    for (int cell = 0; cell < cell_count; ++cell) {

        const auto& [elemLgr, elemLgrCell] = adaptedCell_to_elemLgrAndElemLgrCell.at(cell);
        const auto& elemLgrCellEntity =  Dune::cpgrid::EntityRep<0>(elemLgrCell, true);

        // Auxiliary cell_to_face
        std::vector<Dune::cpgrid::EntityRep<1>> aux_cell_to_face;

        const auto& allCorners = adapted_geometries.geomVector(std::integral_constant<int,3>());
        const auto& grid_or_elemLgr_data = (elemLgr == -1) ? *grid.currentData().back() : *markedElem_to_itsLgr.at(elemLgr);

        // Get the cell geometry.
        const auto& cellGeom = (*(grid_or_elemLgr_data.getGeometry().geomVector(std::integral_constant<int,0>()) ) )[elemLgrCellEntity];
        // Get pre-adapt corners of the cell that will be replaced with leaf view ones.
        const auto& preAdapt_cell_to_point  =  grid_or_elemLgr_data.cellToPoint(elemLgrCell);
        // Get pre-adapt faces of the cell that will be replaced with leaf view ones.
        const auto& preAdapt_cell_to_face =  grid_or_elemLgr_data.cellToFace(elemLgrCell);

        // Cell to point.
        for (int corn = 0; corn < 8; ++corn) {
            int adaptedCorn = 0; // It'll get rewritten.
            const auto& preAdaptCorn = preAdapt_cell_to_point[corn];
            auto adapted_candidate =  elemLgrAndElemLgrCorner_to_adaptedCorner.find({elemLgr, preAdaptCorn});
            if ( adapted_candidate == elemLgrAndElemLgrCorner_to_adaptedCorner.end() ) {
                // Corner might have vanished - Search its equivalent lgr-corner in that case -
                // last lgr where the corner appears -
                std::array<int,2> lastAppearanceLgr_lgrEquivCorner = {0, 0}; // It'll get rewritten.
                if ((elemLgr ==-1) && (!cornerInMarkedElemWithEquivRefinedCorner[preAdaptCorn].empty())) {
                    lastAppearanceLgr_lgrEquivCorner = cornerInMarkedElemWithEquivRefinedCorner[preAdaptCorn].back();
                }
                if (elemLgr > -1) {
                    // Corner might have vanished - Search its equivalent lgr-corner in that case -
                    // last lgr where the corner appears -
                    if(auto candidate = markedElemAndEquivRefinedCorn_to_corner.find({elemLgr, preAdaptCorn}); candidate != markedElemAndEquivRefinedCorn_to_corner.end()) {
                        lastAppearanceLgr_lgrEquivCorner = cornerInMarkedElemWithEquivRefinedCorner[candidate->second].back();
                    }
                    else {
                        // To locate vanished corners, we need a while-loop, since {elemLgr, elemLgrcorner} leads to
                        // {neighboringElemLgr, neighboringElemLgrCornerIdx}, which might have also vanished.
                        // Then, use the lastest appearance of the current corner, meaning, the first (and unique one - by construction) that
                        // gives elemLgrAndElemLgrCorner_to_adaptedCorner.count( lastAppearanceLgr_lgrEquivCorner ) == 1).
                        // This corner lies on the area occupied by a coarse face that got refined and belonged to two marked elements.
                        // Get the index of this corner with respect to the greatest marked element index.
                        lastAppearanceLgr_lgrEquivCorner = vanishedRefinedCorner_to_itsLastAppearance.at({elemLgr, preAdaptCorn});
                        while (elemLgrAndElemLgrCorner_to_adaptedCorner.find(lastAppearanceLgr_lgrEquivCorner) == elemLgrAndElemLgrCorner_to_adaptedCorner.end()) {
                            const auto& tempLgr_lgrCorner =  lastAppearanceLgr_lgrEquivCorner;
                            lastAppearanceLgr_lgrEquivCorner =  vanishedRefinedCorner_to_itsLastAppearance.at(tempLgr_lgrCorner);
                        }
                    }
                }
                adaptedCorn =  elemLgrAndElemLgrCorner_to_adaptedCorner.at(lastAppearanceLgr_lgrEquivCorner);
            }
            // Corner is stored in adapted_corners
            else {
                assert(adapted_candidate != elemLgrAndElemLgrCorner_to_adaptedCorner.end());
                adaptedCorn = adapted_candidate->second;
            }
            adapted_cell_to_point[cell][corn] = adaptedCorn;
        } // end-cell_to_point

        // Cell to face.
        for (const auto& face : preAdapt_cell_to_face) {
            const auto& preAdaptFace = face.index();
            // Face is stored in adapted_faces
            if (auto candidate = elemLgrAndElemLgrFace_to_adaptedFace.find({elemLgr, preAdaptFace}); candidate != elemLgrAndElemLgrFace_to_adaptedFace.end()) {
                const int adaptedFace = candidate->second;
                aux_cell_to_face.push_back({adaptedFace, face.orientation()});
            }
            else{
                // Face might have vanished - Search its refined lgr-children faces in that case -
                // last lgr where the face appears
                if (elemLgr ==-1) { // Coarse face got replaced by its children - from the last appearance of the marked face.
                    assert(!faceInMarkedElemAndRefinedFaces[preAdaptFace].empty());
                    const auto& [lastAppearanceLgr, lastAppearanceLgrFaces] = faceInMarkedElemAndRefinedFaces[preAdaptFace].back();
                    for (const auto& refinedFace : lastAppearanceLgrFaces) {
                        const int adaptedFace = elemLgrAndElemLgrFace_to_adaptedFace.at({lastAppearanceLgr, refinedFace});
                        aux_cell_to_face.push_back({adaptedFace, face.orientation()});
                    }
                }
                if (elemLgr>-1) { // Refined face vanished and its equivalent refined face from a neighboring lgr got stored.
                    // Get shifted level
                    const auto& shiftedLevel = assignRefinedLevel[elemLgr] - preAdaptMaxLevel -1; // Assigned level > preAdapt maxLevel
                    // Get the index of the marked face where the refined face was born.
                    const auto& markedFace = getParentFaceWhereNewRefinedFaceLiesOn(grid,
                                                                                    cells_per_dim_vec[shiftedLevel],
                                                                                    preAdaptFace,
                                                                                    markedElem_to_itsLgr[elemLgr], elemLgr);
                    // Get the last LGR (marked element) where the marked face appeared.
                    const auto& lastLgrWhereMarkedFaceAppeared = faceInMarkedElemAndRefinedFaces[markedFace].back().first;
                    const auto& lastAppearanceLgrEquivFace = replaceLgr1FaceIdxByLgr2FaceIdx(cells_per_dim_vec[shiftedLevel], preAdaptFace, markedElem_to_itsLgr[elemLgr],
                                                                                             cells_per_dim_vec[assignRefinedLevel[lastLgrWhereMarkedFaceAppeared] - preAdaptMaxLevel -1]);
                    const int adaptedFace = elemLgrAndElemLgrFace_to_adaptedFace.at({lastLgrWhereMarkedFaceAppeared, lastAppearanceLgrEquivFace});
                    aux_cell_to_face.push_back({adaptedFace, face.orientation()});
                }
            }
        } // end-cell_to_face
        // Adapted/Leaf-grid-view cell to face.
        adapted_cell_to_face.appendRow(aux_cell_to_face.begin(), aux_cell_to_face.end());

        // Create a pointer to the first element of "adapted_cell_to_point" (required as the fourth argement to construct a Geometry<3,3> type object).
        int* indices_storage_ptr = adapted_cell_to_point[cell].data();
        adapted_cells[cell] = Dune::cpgrid::Geometry<3,3>(cellGeom.center(), cellGeom.volume(), allCorners, indices_storage_ptr);
    } // adapted_cells

    // Adapted/Leaf-grid-view face to cell.
    adapted_cell_to_face.makeInverseRelation(adapted_face_to_cell);
}

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
                                   const int& preAdaptMaxLevel)
{
    // --- Adapted corners ---
    populateLeafGridCorners(grid,
                            adapted_corners,
                            corner_count,
                            markedElem_to_itsLgr,
                            adaptedCorner_to_elemLgrAndElemLgrCorner);
    // --- Adapted faces ---
    populateLeafGridFaces(grid,
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
    populateLeafGridCells(grid,
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
}

void detectActiveLgrs(const Dune::CpGrid& grid,
                      const std::vector<std::array<int,3>>& startIJK_vec,
                              const std::vector<std::array<int,3>>& endIJK_vec,
                              std::vector<int>& lgr_with_at_least_one_active_cell)
{
    auto markLgr = [&lgr_with_at_least_one_active_cell]([[maybe_unused]] const Dune::cpgrid::Entity<0>& element, int level)
    {
        // shifted since starting grid is level 0, and refined grids levels are >= 1.
        lgr_with_at_least_one_active_cell[level] = 1;
    };
    computeOnLgrParents(grid, startIJK_vec, endIJK_vec, markLgr);
}



void predictMinCellAndPointGlobalIdPerProcess([[maybe_unused]] const Dune::CpGrid& grid,
                                                      [[maybe_unused]] const std::vector<int>& assignRefinedLevel,
                                                      [[maybe_unused]] const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                                      [[maybe_unused]] const std::vector<int>& lgr_with_at_least_one_active_cell,
                                                      [[maybe_unused]] int& min_globalId_cell_in_proc,
                                                      [[maybe_unused]] int& min_globalId_point_in_proc)
{
#if HAVE_MPI
    // Maximum global id from level zero. (Then, new entities get global id values greater than max_globalId_levelZero).
    // Recall that only cells and points are taken into account; faces are ignored (do not have any global id).
    auto max_globalId_levelZero = grid.comm().max(grid.currentData().front()->globalIdSet().getMaxGlobalId());

    // Predict how many new cell ids per process are needed.
    std::vector<std::size_t> cell_ids_needed_by_proc(grid.comm().size());
    std::size_t local_cell_ids_needed = 0;
    for ( const auto& element : Dune::elements( grid.levelGridView(0), Dune::Partitions::interior) ) {
        // Get old mark (from level zero). After calling adapt, all marks are set to zero.
        bool hasBeenMarked = grid.currentData().front()->getMark(element) == 1;
        if ( hasBeenMarked ) {
            const auto& level = assignRefinedLevel[element.index()];
            // Shift level (to level -1) since cells_per_dim_vec stores number of subdivisions in each direction (xyz)
            // per parent cell, per level, starting from level 1, ..., maxLevel.
            local_cell_ids_needed += cells_per_dim_vec[level-1][0]*cells_per_dim_vec[level-1][1]*cells_per_dim_vec[level-1][2];
        }
    }
    grid.comm().allgather(&local_cell_ids_needed, 1, cell_ids_needed_by_proc.data());

    // Overestimate ('predict') how many new point ids per process are needed.
    // Assign for all partition type points a 'candidate of global id' (unique in each process).
    std::vector<std::size_t> point_ids_needed_by_proc(grid.comm().size());
    std::size_t local_point_ids_needed = 0;
    for (std::size_t level = 1; level < cells_per_dim_vec.size()+1; ++level){
        if(lgr_with_at_least_one_active_cell[level-1]>0) {
            // Amount of local_point_ids_needed might be overestimated.
            for (const auto& point : Dune::vertices(grid.levelGridView(level))){
                // If point coincides with an existing corner from level zero, then it does not need a new global id.
                if ( grid.currentData()[level]->cornerHistorySize() ) { //  !(*current_data_)[level]->corner_history_.empty()
                    const auto& bornLevel_bornIdx =  grid.currentData()[level]->getCornerHistory(point.index());
                    if (bornLevel_bornIdx[0] == -1)  { // Corner is new-> it needs a new(candidate) global id
                        local_point_ids_needed += 1;
                    }
                }
            }
        }
    }
    grid.comm().allgather(&local_point_ids_needed, 1, point_ids_needed_by_proc.data());

    auto expected_max_globalId_cell = std::accumulate(cell_ids_needed_by_proc.begin(),
                                                      cell_ids_needed_by_proc.end(),
                                                      max_globalId_levelZero + 1);
    min_globalId_cell_in_proc = std::accumulate(cell_ids_needed_by_proc.begin(),
                                                cell_ids_needed_by_proc.begin()+grid.comm().rank(),
                                                max_globalId_levelZero + 1);
    min_globalId_point_in_proc = std::accumulate(point_ids_needed_by_proc.begin(),
                                                 point_ids_needed_by_proc.begin()+ grid.comm().rank(),
                                                 expected_max_globalId_cell+1);
#endif
}

void assignCellIdsAndCandidatePointIds( const Dune::CpGrid& grid,
                                        std::vector<std::vector<int>>& localToGlobal_cells_per_level,
                                        std::vector<std::vector<int>>& localToGlobal_points_per_level,
                                        int min_globalId_cell_in_proc,
                                        int min_globalId_point_in_proc,
                                        const std::vector<std::array<int,3>>& cells_per_dim_vec )
{
    for (std::size_t level = 1; level < cells_per_dim_vec.size()+1; ++level) {
        localToGlobal_cells_per_level[level-1].resize(grid.currentData()[level]-> size(0));
        localToGlobal_points_per_level[level-1].resize(grid.currentData()[level]-> size(3));
        // Notice that in general, (*current_data_)[level]-> size(0) != local owned cells/points.

        // Global ids for cells (for owned cells)
        for (const auto& element : Dune::elements(grid.levelGridView(level))) {
            // At his point, partition_type_indicator_ of refined level grids is not set. However, for refined cells,
            // element.partitionType() returns the partition type of the parent cell. Therefore, all child cells of
            // an interior/overlap parent cell are also interior/overlap.
            if (element.partitionType() == Dune::InteriorEntity) {
                localToGlobal_cells_per_level[level - 1][element.index()] = min_globalId_cell_in_proc;
                ++min_globalId_cell_in_proc;
            }
        }
        for (const auto& point : Dune::vertices(grid.levelGridView(level))) {
            // Check if it coincides with a corner from level zero. In that case, no global id is needed.
            const auto& bornLevel_bornIdx =  grid.currentData()[level]->getCornerHistory(point.index());
            if (bornLevel_bornIdx[0] != -1)  { // Corner in the refined grid coincides with a corner from level 0.
                // Therefore, search and assign the global id of the previous existing equivalent corner.
                const auto& equivPoint = Dune::cpgrid::Entity<3>(*(grid.currentData()[bornLevel_bornIdx[0]]), bornLevel_bornIdx[1], true);
                localToGlobal_points_per_level[level-1][point.index()] = grid.currentData().front()->globalIdSet().id( equivPoint );
            }
            else {
                // Assign CANDIDATE global id to (all partition type) points that do not coincide with
                // any corners from level zero.
                // TO DO after invoking this method: make a final decision on a unique id for points of refined level grids,
                // via a communication step.
                localToGlobal_points_per_level[level-1][point.index()] = min_globalId_point_in_proc;
                ++min_globalId_point_in_proc;
            }
        }
    }
}

void selectWinnerPointIds([[maybe_unused]] const Dune::CpGrid& grid,
                          [[maybe_unused]] std::vector<std::vector<int>>&  localToGlobal_points_per_level,
                          [[maybe_unused]] const std::vector<std::tuple<int,std::vector<int>>>& parent_to_children,
                          [[maybe_unused]] const std::vector<std::array<int,3>>& cells_per_dim_vec)
{
#if HAVE_MPI
    // To store cell_to_point_ information of all refined level grids.
    std::vector<std::vector<std::array<int,8>>> level_cell_to_point(cells_per_dim_vec.size());
    // To decide which "candidate" point global id wins, the rank is stored. The smallest ranks wins,
    // i.e., the other non-selected candidates get rewritten with the values from the smallest (winner) rank.
    std::vector<std::vector<int>> level_winning_ranks(cells_per_dim_vec.size());

    for (std::size_t level = 1; level < cells_per_dim_vec.size()+1; ++level) {

        level_cell_to_point[level -1] = grid.currentData()[level]->cellToPoint();
        // Set std::numeric_limits<int>::max() to make sure that, during communication, the rank of the interior cell
        // wins (int between 0 and comm().size()).
        level_winning_ranks[level-1].resize(grid.currentData()[level]->size(3), std::numeric_limits<int>::max());

        for (const auto& element : Dune::elements(grid.levelGridView(level))) {
            // For interior cells, rewrite the rank value - later used in "point global id competition".
            if (element.partitionType() == Dune::InteriorEntity) {
                for (const auto& corner : grid.currentData()[level]->cellToPoint(element.index())) {
                    int rank = grid.comm().rank();
                    level_winning_ranks[level -1][corner] = rank;
                }
            }
        }
    }
    ParentToChildCellToPointGlobalIdHandle parentToChildCellToPointGlobalId_handle(grid.comm(),
                                                                                   parent_to_children,
                                                                                   level_cell_to_point,
                                                                                   level_winning_ranks,
                                                                                   localToGlobal_points_per_level);
    grid.currentData().front()->communicate(parentToChildCellToPointGlobalId_handle,
                                            Dune::InteriorBorder_All_Interface,
                                            Dune::ForwardCommunication );
#endif
}

void getFirstChildGlobalIds(const Dune::CpGrid& grid,
                            std::vector<int>& parentToFirstChildGlobalIds)
{

    const auto& data = grid.currentData();
    const auto& parentToChildrenBeforeLoadBalance = data[0]->getParentToChildren();
    const auto& globalIdSet = grid.globalIdSet();

    parentToFirstChildGlobalIds.resize(data[0]->size(0), -1); // Initialize with -1 (invalid value, for non parent cells).

    const auto& elements = Dune::elements(grid.levelGridView(0));
    for (const auto& element : elements) {
        const auto& [level, children] = parentToChildrenBeforeLoadBalance[element.index()];

        if (!children.empty()) {
            const auto& levelData = *data[level];
            const auto& first_child = Dune::cpgrid::Entity<0>(levelData, children[0], true);

            // Rewrite parent global id entry with first child global id
            parentToFirstChildGlobalIds[globalIdSet.id(element)] = globalIdSet.id(first_child);
        }
    }
}

/// @brief Extract Cartesian index triplet (i,j,k) given an index between 0 and NXxNYxNZ -1
///    where NX, NY, and NZ is the total amoung of cells in each direction x-,y-,and z- respectively.
///
/// @param [in] idx      Integer between 0 and cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]-1
/// @param [in] cells_per_dim
/// @return Cartesian index triplet.
std::array<int,3> getIJK(int idx_in_parent_cell, const std::array<int,3>& cells_per_dim) 
{
    // idx = k*cells_per_dim_[0]*cells_per_dim_[1] + j*cells_per_dim_[0] + i
    // with 0<= i < cells_per_dim_[0], 0<= j < cells_per_dim_[1], 0<= k <cells_per_dim_[2].
    assert(cells_per_dim[0]);
    assert(cells_per_dim[1]);
    assert(cells_per_dim[2]);

    std::array<int,3> ijk = {0,0,0};
    ijk[0] = idx_in_parent_cell % cells_per_dim[0]; idx_in_parent_cell /= cells_per_dim[0];
    ijk[1] = idx_in_parent_cell % cells_per_dim[1];
    ijk[2] = idx_in_parent_cell /cells_per_dim[1];
    return ijk;
}

void validStartEndIJKs(const std::vector<std::array<int,3>>& startIJK_vec,
                       const std::vector<std::array<int,3>>& endIJK_vec)
{
    if (startIJK_vec.size() == endIJK_vec.size()) {
        for (unsigned int patch = 0; patch < startIJK_vec.size(); ++patch) {
            bool validPatch = true;
            for (int c = 0; c < 3; ++c) {
                // valid startIJK and endIJK for each patch
                validPatch = validPatch && (startIJK_vec[patch][c] < endIJK_vec[patch][c]);
                if (!validPatch) {
                    OPM_THROW(std::invalid_argument, "Invalid IJK-indices in LGR"+std::to_string(patch +1)+", end I/J/K need to be larger than start I/J/K.\n");
                }
            }
        }
    }
    else {
        OPM_THROW(std::invalid_argument, "Sizes of provided vectors with start and end IJK need to match.\n");
    }
}

}
