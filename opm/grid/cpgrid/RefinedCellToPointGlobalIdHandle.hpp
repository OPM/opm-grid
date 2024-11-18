//===========================================================================
//
// File: RefinedCellToPointGlobalIdHandle.hpp
//
// Created: October 29 2024
//
// Author(s): Antonella Ritorto <antonella.ritorto@opm-op.com>
//            Markus Blatt <markus.blatt@opm-op.com>
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

//#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/CommunicationUtils.hpp>
#include <opm/grid/cpgrid/Entity.hpp>

#include <array>
#include <limits>
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

    using DataType = int;

    /// \param level_cell_to_point     Map from all refined level grids, cell_to_point_
    ///                                level_cell_to_point_[ element.level()-1 ][ element.index() ] = { corner0, corner1, ..., corner7 }.
    /// \param level_point_global_ids  A container that for a refined level grid contains all global point ids.
    ///                                level_point_global_ids[ element.level()-1 ][ corner ] = its global id
    ///                                when corner belongs to
    ///                                level_cell_to_point_[ element.level() -1 ][ element.index() ] = { corner0, corner1, ..., corner7 }.
    /** documentation will be fixed soon */
    RefinedCellToPointGlobalIdHandle(const Dune::CpGrid::Communication& comm,
                                     const std::vector<std::array<int,8>>& cell_to_point,
                                     std::vector<DataType>& point_global_ids,
                                     std::vector<DataType>& winning_ranks)
        : comm_(comm)
        , cell_to_point_(cell_to_point)
        , point_global_ids_(point_global_ids)
        , winning_ranks_(winning_ranks)
    {
    }

    // 
    bool fixedSize(std::size_t, std::size_t)
    {
        return false;
    }
    // Only communicate values attached to refined cells.
    bool contains(std::size_t, std::size_t codim)
    {
        return codim == 0;
    }
    // Communicate variable size: for interior refined cells 9( 8 corners per refined cell plus rank)
    template <class T> // T = Entity<0>
    std::size_t size(const T& element)
    {
        if (element.partitionType() == Dune::InteriorEntity ) {
            return 9; //1+cell_to_point_[element.index()].size(); // rank, each refined cell has 8 corners.
        }
        else {
            return 1;
        }
    }

    // Gather global ids of the 8 corners of an interior refined cell
    template <class B, class T> // T = Entity<0>
    void gather(B& buffer, const T& element)
    {
        // Skip values that are non interior.
        if ( element.partitionType() != Dune::InteriorEntity ) {
            buffer.write(42); // For overlap cells, size() is equal to 1. 
            return;
        }
        // Store/rewritte global ids in the buffer when the element is an interior refined cell.
        // Write the rank, for example via the "corner 0" of cell_to_point_ of current element.
        buffer.write(winning_ranks_[ cell_to_point_[ element.index() ][0] ]);
        for (const auto& corner : cell_to_point_[element.index()]) {
            buffer.write(point_global_ids_[corner]);
        }
    }

    // Scatter global ids of overlap refined cells.
    template <class B, class T> // T = Entity<0>
    void scatter(B& buffer, const T& element, [[maybe_unused]] std::size_t num_corners_plus_one)
    {
        // Read all values (9 in total) to advance the pointer used by the buffer to the correct index.
        // Skip interior refined cells.
        if  ( element.partitionType() == Dune::InteriorEntity ) {
            // Read all values to advance the pointer used by the buffer
            // to the correct index
            DataType tmp_rank;
            buffer.read(tmp_rank);
            for (std::size_t corner = 0; corner < 8;  ++corner) {
                DataType tmp_id;
                buffer.read(tmp_id);
                }
        }
        else { // Overlap refined cell.
            // Read and store the values in the correct location directly.
            DataType tmp_rank;
            buffer.read(tmp_rank);
                for (const auto& corner: cell_to_point_[element.index()]){
                    auto& min_rank = winning_ranks_[corner];
                    if (tmp_rank < min_rank) {
                         min_rank = tmp_rank; 
                         auto& target_entry = point_global_ids_[corner];
                         // DataType min_rank_target_entity;
                        buffer.read(target_entry);
                        // Rewrite
                        //point_global_ids_[corner] = min_rank_target_entity;
                    } else {
                        DataType rubbish;
                        buffer.read(rubbish);
                    }
                }
            }
    }

private:
    const Dune::CpGrid::Communication& comm_;
    const std::vector<std::array<int,8>>& cell_to_point_;
    std::vector<DataType>& point_global_ids_;
    std::vector<DataType>& winning_ranks_;
};
#endif // HAVE_MPI
} // namespace
#endif

