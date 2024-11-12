//===========================================================================
//
// File: RefinedCellToPointGlobalIdHandle.hpp
//
// Created: October 29 2024
//
// Author(s): Antonella Ritorto <antonella.ritorto@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2024 TBD
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

#ifndef OPM_REFINEDCELLTOPOINTGLOBALIDHANDLE_HEADER
#define OPM_REFINEDCELLTOPOINTGLOBALIDHANDLE_HEADER


#include <opm/grid/cpgrid/Entity.hpp>

#include <array>
#include <vector>


namespace
{
#if HAVE_MPI

/// \brief Handle for correcting assignment of global ids of refined points.
///
struct RefinedCellToPointGlobalIdHandle {
    //   - The container used for gather and scatter contains "may-be-rewritten-global ids" for points of the refined level grids (LGRs).
    //     Access is done with the level and local index of a refined cell, and the local indices of its 8 corners.
    //     level_point_global_ids[ element.level() -1 ] [ corner ] = may-be-rewritten-point-global-id,
    //     when corner belongs to
    //     level_cell_to_point_[ element.level() -1 ][ element.index() ] = { corner0, corner1, ..., corner7 }.
    //     The shift 'element.level() -1' is due to not taken into account level zero grid.
    //   - We use the number of entries in the scatter method. That number is the number of corners when sending.
    //   - We skip non-refined cells, i.e., coarse cells equivalent to level zero cells.

    using DataType = std::pair<int,int>;

    /// \param level_cell_to_point     Map from all refined level grids, cell_to_point_
    ///                                level_cell_to_point_[ element.level()-1 ][ element.index() ] = { corner0, corner1, ..., corner7 }.
    /// \param level_point_global_ids  A container that for a refined level grid contains all global point ids.
    ///                                level_point_global_ids[ element.level()-1 ][ corner ] = its global id
    ///                                when corner belongs to
    ///                                level_cell_to_point_[ element.level() -1 ][ element.index() ] = { corner0, corner1, ..., corner7 }.
    RefinedCellToPointGlobalIdHandle(const std::vector<std::vector<std::array<int,8>>>& level_cell_to_point,
                                     std::vector<std::vector<int>>& level_point_global_ids,
                                     std::vector<std::vector<DataType>>& level_point_global_ids_with_ranks)
        : level_cell_to_point_(level_cell_to_point)
        , level_point_global_ids_(level_point_global_ids)
        , level_point_global_ids_with_ranks_(level_point_global_ids_with_ranks)
    {
    }

    // We do not take into account cells from level zero.
    bool fixedSize(std::size_t, std::size_t)
    {
        return false;
    }
    // Only communicate values attached to refined cells.
    bool contains(std::size_t, std::size_t codim)
    {
        return codim == 0;
    }
    // Communicate variable size: 8 corners per refined cell; skip coarse cells from level zero.
    template <class T> // T = Entity<0>
    std::size_t size(const T& element)
    {
        // Skip values that are not interior
        if ( (element.getLevelElem().partitionType() != Dune::InteriorEntity)  || (element.level() == 0) )
            return 1;
        return level_cell_to_point_[element.level()-1][element.getLevelElem().index()].size(); // each refined cell has 8 corners.
    }

    // Gather global ids of points of an interior refined cell
    template <class B, class T> // T = Entity<0>
    void gather(B& buffer, const T& element)
    {
        // Skip values that are non interior or coarse cells.
        if ( (element.getLevelElem().partitionType() != Dune::InteriorEntity)  || (element.level() == 0) ) {
            DataType tmp;
            buffer.write(tmp);
            return;
        }
        // Store/rewritte global ids in the buffer when the element is an interior refined cell.
        for (const auto& corner : level_cell_to_point_[element.level()-1][element.getLevelElem().index()]) {
            // const auto& [corner_id, rank] = level_point_global_ids_with_ranks_[element.level()-1][corner];
            // buffer.write(level_point_global_ids_[element.level()-1][corner]);
            buffer.write(level_point_global_ids_with_ranks_[element.level()-1][corner]); // [corner_id, rank]
        }
    }

    // Scatter global ids of overlap refined cells.
    template <class B, class T> // T = Entity<0>
    void scatter(B& buffer, const T& element, std::size_t num_corners)
    {
        // Read all values to advance the pointer used by the buffer to the correct index.
        // Skip interior refined cells and overlap coarse cells.
        if  ( (element.getLevelElem().partitionType() == Dune::InteriorEntity) || (element.level() == 0) ) {
            // Read all values to advance the pointer used by the buffer
            // to the correct index
            for (std::size_t corner = 0; corner < num_corners;  ++corner) {
                DataType tmp;
                buffer.read(tmp);
            }
        }
        else { // Overlap refined cell.
            // Read and store the values in the correct location directly.
            for (const auto& corner: level_cell_to_point_[element.level()-1][element.getLevelElem().index()]){
                auto cornerId_rank  = level_point_global_ids_with_ranks_[element.level()-1][corner];
                DataType tmp;
                buffer.read(tmp);
                if (cornerId_rank.second < tmp.second) { // re-write
                    level_point_global_ids_with_ranks_[element.level()-1][corner] = tmp;
                    level_point_global_ids_[element.level()-1][corner] = tmp.first;
                }
                else {
                    // auto target_entry = level_point_global_ids_[element.level()-1][corner];
                    //  buffer.read(cornerId_rank);
                }
              
            }
        }
    }

private:
    const std::vector<std::vector<std::array<int,8>>>& level_cell_to_point_;
    std::vector<std::vector<int>>& level_point_global_ids_;
    std::vector<std::vector<DataType>>& level_point_global_ids_with_ranks_;
    
};
#endif // HAVE_MPI
} // namespace
#endif
