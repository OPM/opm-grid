//===========================================================================
//
// File: ParentToChildrenCellGlobalIdHandle.hpp
//
// Created: October 24 2024
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

#ifndef OPM_PARENTTOCHILDRENCELLGLOBALIDHANDLE_HEADER
#define OPM_PARENTTOCHILDRENCELLGLOBALIDHANDLE_HEADER


#include <opm/grid/cpgrid/Entity.hpp>

#include <tuple>
#include <vector>


namespace
{
#if HAVE_MPI

/// \brief Handle for assignment of global ids of refined cells (cells from refined level grids). 
struct ParentToChildrenCellGlobalIdHandle {
    //   - The container used for gather and scatter contains global ids for interior elements of the refined level grids (LGRs).
    //     Access is done with the local index of a refined cell (child_cell_local_index),
    //     via its parent cell index and its parent cell list of children.
    //     level_cell_global_ids[ level-1 ][ child_cell_local_index ] = global id of the child.
    //     when child_cell_local_index belongs to the children_local_index_list:
    //     parent_to_children_[ element.index() ] =  parent_to_children_[ element.index() ] = { level, children_local_index_list }
    //   - We use the number of entries in the scatter method. That number is the number of children when sending.
    //     (we still assume that the number entries is the same as the number of children of the element).

    using DataType = int;

    /// \param parent_to_children Map from parent index to all children, and the level they are stored.
    ///                           parent_to_children_[ element.index() ] = { level, children_list }
    /// \param level_cell_global_ids   A container that for the elements of a level contains all ids.
    ParentToChildrenCellGlobalIdHandle(const std::vector<std::tuple<int, std::vector<int>>>& parent_to_children,
                                       std::vector<std::vector<DataType>>& level_cell_global_ids)
        : parent_to_children_(parent_to_children)
        , level_cell_global_ids_(level_cell_global_ids)
    {
    }

    // Not every cell has children. When they have children, the amount might vary.
    bool fixedSize(std::size_t, std::size_t)
    {
        return false;
    }
    // Only communicate values attached to cells.
    bool contains(std::size_t, std::size_t codim)
    {
        return codim == 0;
    }
    // Communicate variable size: amount of child cells from a parent cell from level zero grid.
    template <class T> // T = Entity<0>
    std::size_t size(const T& element)
    {
        // Skip values that are not interior, or have no children (in that case, 'invalid' level = -1)
        const auto& [level, children] = parent_to_children_[element.index()];
        if ( (element.partitionType() != Dune::InteriorEntity) || (level == -1))
            return 0;
        return children.size();
    }

    // Gather global ids of child cells of a coarse interior parent cell
    template <class B, class T> // T = Entity<0>
    void gather(B& buffer, const T& element)
    {
        // Skip values that are not interior, or have no children (in that case, 'invalid level' = -1)
        const auto& [level, children] = parent_to_children_[element.index()];
        if ( (element.partitionType() != Dune::InteriorEntity) || (level==-1)) {
            return;
        }
        // Store the children's global ids in the buffer when the element is interior and has children.
        for (const auto& child : children) 
            // Shift level-> level-1 since level_global_ids_ stores only refined level grids global ids.
            // level_cell_global_ids_[0] corresponds to level 1, ..., level_cell_global_ids_[ maxLevel -1 ] to maxLevel grid.
            buffer.write(level_cell_global_ids_[level-1][child]);
    }

    // Scatter global ids of child cells of a coarse overlap parent cell
    template <class B, class T> // T = Entity<0>
    void scatter(B& buffer, const T& element, std::size_t num_children)
    {
        const auto& [level, children] = parent_to_children_[element.index()];
        if ( (element.partitionType() == Dune::OverlapEntity) && (level==-1) ) {
            // Read all values to advance the pointer used by the buffer
            // to the correct index
            for (std::size_t child = 0; child < num_children;  ++child) {
                DataType tmp;
                buffer.read(tmp);
            }
        }
        else {
            // Read and store the values in the correct location directly.
            // Careful, we assume that the order of the children is the same on
            // each process.
           for (const auto& child : children) {
                // Shift level-> level-1 since level_cell_global_ids_ stores only refined level grids cell global ids.
                // level_cell_global_ids_[0] corresponds to level 1, ..., level_cell_global_ids_[ maxLevel -1 ] to maxLevel grid.
                auto& target_entry = level_cell_global_ids_[level-1][child];
                buffer.read(target_entry);
            }
        }
    }

private:
const std::vector<std::tuple<int, std::vector<int>>>& parent_to_children_;
std::vector<std::vector<DataType>>& level_cell_global_ids_;
};
#endif // HAVE_MPI
} // namespace
#endif
