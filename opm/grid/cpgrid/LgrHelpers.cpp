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


}
