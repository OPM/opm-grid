//===========================================================================
//
// File: ParentToChildCellToPointGlobalIdHandle.hpp
//
// Created: November 19 2024
//
// Author(s): Antonella Ritorto <antonella.ritorto@opm-op.com>
//            Markus Blatt      <markus.blatt@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2024 Equinor ASA
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

#ifndef OPM_PARENTTOCHILDCELLTOPOINTGLOBALIDHANDLE_HEADER
#define OPM_PARENTTOCHILDCELLTOPOINTGLOBALIDHANDLE_HEADER

#include <opm/grid/common/CommunicationUtils.hpp>
#include <opm/grid/cpgrid/Entity.hpp>

#include <array>
#include <tuple>
#include <vector>


namespace
{
#if HAVE_MPI

/// \brief Handle for assignment of point global ids of refined cells.
struct ParentToChildCellToPointGlobalIdHandle {
    //   - The container used for gather and scatter contains "candidates" point global ids for interior elements of the refined level grids (LGRs).
    //     Access is done with the local index of its parent cell index and its parent cell list of children.
    //     level_point_global_ids[ level-1 ][ child_cell_local_index [ corner ] ] = "candidate" point global id,
    //     when child_cell_local_index belongs to the children_local_index_list:
    //     parent_to_children_[ element.index() ] =  parent_to_children_[ element.index() ] = { level, children_local_index_list }
    //     and corner = 0, ...,7.
    //     To decide which "candidate" point global id wins, we use the rank. The smallest ranks wins,
    //     i.e., the other non-selected candidates get rewritten with the values from the smallest (winner) rank.
    //   - In the scatter method, the "winner" rank and the 8 point global ids of each number of children) get rewritten.

    using DataType = int;

    /// \param comm                    Communication
    /// \param parent_to_children      Map from parent index to all children, and the level they are stored.
    ///                                parent_to_children_[ element.index() ] = { level, children_list local indices }
    /// \param level_cell_to_point
    /// \param level_point_global_ids  A container that for the elements of a level contains all candidate point global ids.
    /// \param level_winning_ranks
    /// \param level_point_global_ids
    ParentToChildCellToPointGlobalIdHandle(const Dune::CpGrid::Communication& comm,
                                           const std::vector<std::tuple<int, std::vector<int>>>& parent_to_children,
                                           const std::vector<std::vector<std::array<int,8>>>& level_cell_to_point,
                                           std::vector<std::vector<DataType>>& level_winning_ranks,
                                           std::vector<std::vector<DataType>>& level_point_global_ids)
    : comm_(comm)
    , parent_to_children_(parent_to_children)
    , level_cell_to_point_(level_cell_to_point)
    , level_winning_ranks_(level_winning_ranks)
    , level_point_global_ids_(level_point_global_ids)
    {}

    bool fixedSize(std::size_t, std::size_t)
    {
        // Not every cell has children. When they have children, the amount might vary.
        return false;
    }

    bool contains(std::size_t, std::size_t codim)
    {
        // Only communicate values attached to cells.
        return codim == 0;
    }

    template <class T> // T = Entity<0>
    std::size_t size(const T& element)
    {
        // Communicate variable size: 1 (rank) + (8* amount of child cells) from an interior parent cell from level zero grid.
        // Skip values that are not interior, or have no children (in that case, 'invalid' level = -1)
        const auto& [level, children] = parent_to_children_[element.index()];
        // [Bug in dune-common] VariableSizeCommunicator will deadlock if a process attempts to send a message of size zero.
        // This can happen if the size method returns zero for all entities that are shared with another process.
        // Therefore, when skipping cells without children or for overlap cells, we set the size to 1.
        if ( (element.partitionType() != Dune::InteriorEntity) || (level == -1))
            return 1;
        return 1 + ( 8*children.size()); // rank + 8 "winner" point global ids per child cell
    }

    // Gather global ids of child cells of a coarse interior parent cell
    template <class B, class T> // T = Entity<0>
    void gather(B& buffer, const T& element)
    {
        // Skip values that are not interior, or have no children (in that case, 'invalid level' = -1)
        const auto& [level, children] = parent_to_children_[element.index()];
        // [Bug in dune-common] VariableSizeCommunicator will deadlock if a process tries to send a message with size zero.
        // To avoid this, for cells without children or for overlap cells, we set the size to 1 and write a single DataType
        // value (e.g., '42').
        if ( (element.partitionType() != Dune::InteriorEntity) || (level==-1)) {
            buffer.write(42);
            return;
        }
        // Store the children's corner global ids in the buffer when the element is interior and has children.
        // Write the rank first, for example via the "corner 0" of cell_to_point_ of the first child:
        // First child: children[0]
        // First corner of first child:  level_cell_to_point_[ level -1 ][children[0]] [0]
        buffer.write( comm_.rank() );
        for (const auto& child : children)
            for (const auto& corner : level_cell_to_point_[level -1][child])
                buffer.write(level_point_global_ids_[level-1][corner]);
    }

    // Scatter global ids of child cells of a coarse overlap parent cell
    template <class B, class T> // T = Entity<0>
    void scatter(B& buffer, const T& element, std::size_t size)
    {
        const auto& [level, children] = parent_to_children_[element.index()];
        // Read all values to advance the pointer used by the buffer to the correct index.
        // (Skip overlap-cells-without-children and interior-cells).
        if ( ( (element.partitionType() == Dune::OverlapEntity) && (level==-1) ) || (element.partitionType() == Dune::InteriorEntity ) ) {
            // Read all values to advance the pointer used by the buffer
            // to the correct index
            for (std::size_t received_int = 0; received_int < size;  ++received_int) {
                DataType tmp;
                buffer.read(tmp);
            }
        }
        else { // Overlap cell with children.
            // Read and store the values in the correct location directly.
            // The order of the children is the same on each process.
            assert(children.size()>0);
            assert(level>0);
            // Read and store the values in the correct location directly.
            DataType tmp_rank;
            buffer.read(tmp_rank);
            for (const auto& child : children) {
                for (const auto& corner : level_cell_to_point_[level -1][child]) {
                    auto& min_rank = level_winning_ranks_[level-1][corner];
                    // Rewrite the rank (smaller rank wins)
                    if (tmp_rank < min_rank) {
                        min_rank = tmp_rank;
                        auto& target_entry = level_point_global_ids_[level-1][corner];
                        buffer.read(target_entry);
                    } else {
                        DataType rubbish;
                        buffer.read(rubbish);
                    }
                }
            }
        }
    }

private:
    const Dune::CpGrid::Communication& comm_;
    const std::vector<std::tuple<int, std::vector<int>>>& parent_to_children_;
    const std::vector<std::vector<std::array<int,8>>>& level_cell_to_point_;
    std::vector<std::vector<DataType>>& level_winning_ranks_;
    std::vector<std::vector<DataType>>& level_point_global_ids_;
};
#endif // HAVE_MPI
} // namespace
#endif
