//===========================================================================
//
// File: CpGridLGR.hpp
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

#ifndef OPM_CPGRIDLGR_HEADER
#define OPM_CPGRIDLGR_HEADER

#include <opm/grid/CpGrid.hpp>

namespace Dune
{

    ////////////////////////////////////////////////////////////////////////
    //
    //   CpGridLGR
    //
    ////////////////////////////////////////////////////////////////////////

    /// \brief [<em> extends CpGrid with LGR (Local Grid Refinement) functionality </em>]
    ///
    /// CpGridLGR inherits from CpGrid and adds all methods related to
    /// Local Grid Refinement (LGR), including adaptivity (mark, adapt, etc.),
    /// LGR creation and management, and related utility methods.
    class CpGridLGR : public CpGrid
    {
    public:
        /// Default constructor
        CpGridLGR();

        explicit CpGridLGR(MPIHelper::MPICommunicator comm);

        /// \brief  Refine the grid refCount times using the default refinement rule.
        ///         This behaves like marking all elements for refinement and then calling preAdapt, adapt and postAdapt.
        ///         The state after globalRefine is comparable to the state after postAdapt.
        /// @param [in] refCount Refinement level
        /// @param [in] throwOnFailure If true, the function will throw an exception if the marking of any entity is invalid.
        void globalRefine (int refCount, bool throwOnFailure = false);

        /// @brief Create a grid out of a coarse one and (at most) 2 refinements(LGRs) of selected block-shaped disjoint patches
        ///        of cells from that coarse grid.
        ///
        /// Level0 refers to the coarse grid, assumed to be this-> data_[0]. Level1 and level2 refer to the LGRs (stored in this->data_[1]
        /// data_[2]). LeafView (stored in this-> data_[3]) is built with the level0-entities which weren't involded in the
        /// refinenment, together with the new born entities created in level1 and level2.
        /// Old-corners and old-faces (from coarse grid) lying on the boundary of the patches, get replaced by new-born-equivalent corners
        /// and new-born-faces.
        ///
        /// @param [in] cells_per_dim_vec         Vector of Number of (refined) cells in each direction that each
        ///                                       parent cell should be refined to.
        /// @param [in] startIJK_vec              Vector of Cartesian triplet indices where each patch starts.
        /// @param [in] endIJK_vec                Vector of Cartesian triplet indices where each patch ends.
        ///                                       Last cell part of each patch(lgr) will be
        ///                                       {endIJK_vec[<patch-number>][0]-1, ..., endIJK_vec[<patch-number>][2]-1}.
        /// @param [in] lgr_name_vec              Names (std::string) for the LGRs/levels.
        /// @param [in] lgr_parent_grid_name_vec  Names (std::string) for the LGRs/levels parent grids.
        ///                                       A parent grid may have more than one child LGR.
        void addLgrsUpdateLeafView(const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                   const std::vector<std::array<int,3>>& startIJK_vec,
                                   const std::vector<std::array<int,3>>& endIJK_vec,
                                   const std::vector<std::string>& lgr_name_vec,
                                   const std::vector<std::string>& lgr_parent_grid_name_vec = std::vector<std::string>{});

        /// @brief Global refine the grid with different refinement factors in each direction.
        ///
        ///        Related to AUTOREF keyword. Each refinement factor must be odd and positive.
        ///        The restriction on odd factors might be related to placing the wells in the central
        ///        refined cell, or row of cells, depending on the well direction.
        /// @param [in] nxnynz      Refinement factors in x-, y-, and z-direction.
        void autoRefine(const std::array<int,3>& nxnynz);

        // @brief TO BE DONE
        const std::map<std::string,int>& getLgrNameToLevel() const;

        // @brief Return parent (coarse) intersection (face) of a refined face on the leaf grid view, whose neighboring cells
        //        are two: one coarse cell (equivalent to its origin cell from level 0), and one refined cell
        //        from certain LGR.
        //        Used in Transmissibility_impl.hpp
        Dune::cpgrid::Intersection getParentIntersectionFromLgrBoundaryFace(const Dune::cpgrid::Intersection& intersection) const;

        /// --------------- Adaptivity (begin) ---------------
        /// @brief Mark entity for refinement (or coarsening).
        ///
        /// Refinement on CpGrid is partially supported for Cartesian grids, with the keyword CARFIN.
        /// Nested refinement is not supported yet, so the so-called "host grid" is defined by default
        /// equal to the GLOBAL Grid (level zero). Therefore, we mark elements (from the GLOBAL grid)
        /// for refinement.
        /// This only works for entities of codim 0.
        /// In distributed grids, element markings are synchronized across processes.
        /// If an element is marked by one process, all other processes that share that element
        /// will also receive and apply the same mark.
        ///
        /// @param [in] refCount   To mark the element for
        ///                        - refinement, refCount == 1
        ///                        - doing nothing, refCount == 0
        ///                        - coarsening, refCount == -1 (not applicable yet)
        /// @param [in] element    Entity<0>. Currently, an element from the GLOBAL grid (level zero).
        /// @param [in] throwOnFailure If true, the function will throw an exception if the marking is invalid.
        /// @return true, if marking was succesfull.
        ///         false, if marking was not possible.
        bool mark(int refCount, const cpgrid::Entity<0>& element, bool throwOnFailure = false);

        /// @brief Return refinement mark for entity.
        ///
        /// @return refinement mark (1,0,-1)  Currently, only 1 (refinement), or 0 (doing nothing).
        int getMark(const cpgrid::Entity<0>& element) const;

        /// @brief Set mightVanish flags for elements that will be refined in the next adapt() call
        ///        Need to be called after elements have been marked for refinement.
        bool preAdapt();

        /// @brief Triggers the grid refinement process.
        ///        Returns true if the grid has changed, false otherwise.
        bool adapt();

        /// @brief Triggers the grid refinement process, allowing to select diffrent refined level grids.
        ///
        /// @param [in] cells_per_dim_vec    For each set of marked elements for refinement, that will belong to a same
        ///                                  refined level grid, number of (refined) cells in each direction that each
        ///                                  parent cell should be refined to.
        /// @param [in] assignRefinedLevel   Vector with size equal to total amount of cells of the starting grid where
        ///                                  the marked elements belong. In each entry, the refined level grid where the
        ///                                  refined entities of the (parent) marked element should belong is stored.
        /// @param [in] lgr_name_vec         Each refined level grid name, e.g. {"LGR1", "LGR2"}.
        /// @param [in] startIJK_vec         Default empty vector. When isCARFIN, the starting ijk Cartesian index of each
        ///                                  block of cells to be refined.
        /// @param [in] endIJK_vec           Default empty vector. When isCARFIN, the final ijk Cartesian index of each
        ///                                  block of cells to be refined.
        bool refineAndUpdateGrid(const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                 const std::vector<int>& assignRefinedLevel,
                                 const std::vector<std::string>& lgr_name_vec,
                                 const std::vector<std::array<int,3>>& startIJK_vec = std::vector<std::array<int,3>>{},
                                 const std::vector<std::array<int,3>>& endIJK_vec = std::vector<std::array<int,3>>{});

        /// @brief Clean up refinement markers - set every element to the mark 0 which represents 'doing nothing'
        void postAdapt();
        /// --------------- Adaptivity (end) ---------------

        /// @brief Synchronizes cell global ids across processes after load balancing.
        ///
        /// LGRs (Local Grid Refinements) can be added either in the undistributed view first and then in the distributed view,
        /// or vice versa. This method ensures consistency by rewriting the global cell ids in the distributed view
        /// using the corresponding ids from the undistributed view.
        void syncDistributedGlobalCellIds();

        /// @brief Compute for each level grid, a map from the global_cell_[ cell index in level grid ] to the leaf index of the equivalent cell
        ///        on the leaf grid view.
        ///        Notice that cells that vanished and do not appear on the leaf grid view will not be considered.
        ///        global_cell_[ cell index in level grid ] coincide with (local) Cartesian Index.
        std::vector<std::unordered_map<std::size_t, std::size_t>> mapLocalCartesianIndexSetsToLeafIndexSet() const;

        /// @brief Reverse map: from leaf index cell to { level, local/level Cartesian index of the cell }
        std::vector<std::array<int,2>> mapLeafIndexSetToLocalCartesianIndexSets() const;

    private:
        void updateCornerHistoryLevels(const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                       const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                       const std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner,
                                       const int& corner_count,
                                       const std::vector<std::array<int,2>>& preAdaptGrid_corner_history,
                                       const int& preAdaptMaxLevel,
                                       const int& newLevels);

        void globalIdsPartitionTypesLgrAndLeafGrids(const std::vector<int>& assignRefinedLevel,
                                                    const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                                    const std::vector<int>& lgr_with_at_least_one_active_cell);

        /// @brief Retrieves the global ids of the first child for each parent cell in the grid.
        ///
        /// If a cell has no children, its entry is set to -1, indicating an invalid id.
        ///
        /// @param[out] parentToFirstChildGlobalIds A vector that will be filled with the first child global IDs.
        ///                                         The vector is resized to match the number of parent cells.
        void getFirstChildGlobalIds([[maybe_unused]] std::vector<int>& parentToFirstChildGlobalIds);

        /// @brief For refined level grids created based on startIJK and endIJK values, compute the "local ijk/Cartesian index" within the LGR.
        ///
        /// It's confusing that this "localIJK" is stored in CpGridData member global_cell_. Potential explanation: level zero grid is also called
        /// "GLOBAL" grid.
        /// Example: a level zero grid with dimension 4x3x3, an LGR with startIJK = {1,2,2}, endIJK = {3,3,3}, and cells_per_dim = {2,2,2}.
        /// Then the dimension of the LGR is (3-1)*2 x (3-2)*2 x (3-2)* 2 = 4x2x2 = 16. Therefore the global_cell_lgr minimim value should be 0,
        /// the maximum should be 15.
        /// To invoke this method, each refined level grid must have 1. logical_cartesian_size_, 2. cell_to_idxInParentCell_, and 3. cells_per_dim_
        /// already populated.
        ///
        /// @param [in] level    Grid index where LGR is stored
        /// @param [out] global_cell_lgr
        void computeGlobalCellLgr(const int& level, const std::array<int,3>& startIJK, std::vector<int>& global_cell_lgr);

        /// @brief For a leaf grid with with LGRs, we assign the global_cell_ values of either the parent cell or the equivalent cell from
        ///        level zero.
        ///        For nested refinement, we lookup the oldest ancestor, from level zero.
        void computeGlobalCellLeafGridViewWithLgrs(std::vector<int>& global_cell_leaf);

        /// @brief Mark selected elements, assign them their corresponding level, and detect active LGRs.
        ///
        /// Given blocks of cells selected for refinement, mark selected elements and assign them their corresponding
        /// (refined) level (grid). When level zero grid is distributed before refinement, detect which LGRs are active
        /// in each process.
        ///
        /// @param [in] startIJK_vec    Vector of ijk values denoting the start of each block of cells selected for refinement.
        /// @param [in] endIJK_vec      Vector of ijk values denoting the end of each block of cells selected for refinement.
        /// @param [out] assignRefinedLevel   Assign level for the refinement of each marked cell. Example: refined element from
        ///                                   LGR1 have level 1, refined element rfom LGR2 have level 2, etc.
        /// @param [out] lgr_with_at_least_one_active_cell Determine if an LGR is not empty in a given process, we set
        ///                                                lgr_with_at_least_one_active_cell[in that level] to 1 if it contains
        ///                                                at least one active cell, and to 0 otherwise.
        void markElemAssignLevelDetectActiveLgrs(const std::vector<std::array<int,3>>& startIJK_vec,
                                                  const std::vector<std::array<int,3>>& endIJK_vec,
                                                  std::vector<int>& assignRefinedLevel,
                                                  std::vector<int>& lgr_with_at_least_one_active_cell);

        /// @brief For a grid whose level zero has been distributed and then locally refined, populate the cell_index_set_ of each refined level grid.
        void populateCellIndexSetRefinedGrid(int level);

        /// @brief For a grid whose level zero has been distributed and then locally refined, populate the cell_index_set_ of the leaf grid view.
        void populateCellIndexSetLeafGridView();

        /// @brief For a grid whose level zero has been distributed and then locally refined, populate the global_id_set_ of the leaf grid view.
        void populateLeafGlobalIdSet();

        /** @brief To get the level given the lgr-name. Default, {"GLOBAL", 0}. */
        std::map<std::string,int> lgr_names_ = {{"GLOBAL", 0}};

    }; // end Class CpGridLGR

} // end namespace Dune

#endif // OPM_CPGRIDLGR_HEADER
