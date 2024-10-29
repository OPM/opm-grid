//===========================================================================
//
// File: ParentToChildrenCellGlobalIdHandle.hpp
//
// Created: October 24 2024
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

#ifndef OPM_PARENTTOCHILDRENCELLGLOBALIDHANDLE_HEADER
#define OPM_PARENTTOCHILDRENCELLGLOBALIDHANDLE_HEADER


#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/cpgrid/Indexsets.hpp>

#include <vector>


namespace
{
#if HAVE_MPI

/// \brief Handle for assignment of global ids of refined cells.
///       Idea: each refined cell 'element' has a parent in level zero and certain 'storage index'
///       in the children list of its parent cell, where
///       (level zero grid).parent_to_children_[ element.father().index() ] = { level, children_list }.
///       Assuming parent_to_children_cells_globalsIds maps [global id of parent cell, global ids of children]
///       where only interior cells are taken into account, then we can gather/scatter to get the global id of
///       a refined overlap cell.
struct ParentToChildrenCellGlobalIdHandle {
    using DataType = std::vector<int>;

    ParentToChildrenCellGlobalIdHandle(std::vector<std::vector<int>> parent_to_children_cells_globalIds,
                                       const Dune::cpgrid::LevelGlobalIdSet& globalIdSet_levelZero)
        : gather_parent_to_children_cells_globalIds_(parent_to_children_cells_globalIds)
        , scatter_parent_to_children_cells_globalIds_(parent_to_children_cells_globalIds)
        , globalIdSet_levelZero_(globalIdSet_levelZero)
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
    // Communicate variable size: amount of children cells from a parent cell from level zero grid.
    template <class T> // T = Entity<0>
    std::size_t size(const T& element)
    {
        const auto& globalId = globalIdSet_levelZero_.id(element);
        if (!gather_parent_to_children_cells_globalIds_[globalId].empty()) {
            const auto& children_list = gather_parent_to_children_cells_globalIds_[globalId];
            std::cout << children_list.size() << std::endl;
            return children_list.size(); // zero for non-parent cells
        } else {
            return 0;
        }
    }

    // Gather global ids of children cells of a coarse interior cell
    template <class B, class T> // T = Entity<0>
    void gather(B& buffer, const T& element)
    {
        // Gather only for interior cells.
        assert(element.partitionType() == Dune::InteriorEntity);
        // Get global id of the cell.
        const auto& globalId = globalIdSet_levelZero_.id(element);
        buffer.write(gather_parent_to_children_cells_globalIds_[globalId]);
    }

    // Scatter global ids of children cells of a coarse overlap cell
    template <class B, class T> // T = Entity<0>
    void scatter(B& buffer, const T& element, std::size_t)
    {
        // Scatter only for overlap cells.
        const auto& globalId = globalIdSet_levelZero_.id(element);
        DataType tmp = gather_parent_to_children_cells_globalIds_[globalId];
        buffer.read(tmp);


        if (element.partitionType() == Dune::OverlapEntity) {
            scatter_parent_to_children_cells_globalIds_[globalId] = tmp;
        }
    }

private:
    const std::vector<std::vector<int>>& gather_parent_to_children_cells_globalIds_;
    std::vector<std::vector<int>>& scatter_parent_to_children_cells_globalIds_;
    const Dune::cpgrid::LevelGlobalIdSet& globalIdSet_levelZero_;
};
#endif // HAVE_MPI
} // namespace
#endif
