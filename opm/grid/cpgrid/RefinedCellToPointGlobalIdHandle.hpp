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

#include <unordered_map>
#include <vector>


namespace
{
#if HAVE_MPI

/// \brief Handle for assignment of global ids of refined points....
///
struct RefinedCellToPointGlobalIdHandle {
    using DataType = std::vector<std::size_t>;

    RefinedCellToPointGlobalIdHandle(
        std::unordered_map<int, std::vector<std::size_t>> refined_cell_to_point_global_ids,
        const std::vector<int>& local_to_global_refined_cells)
        : gather_refined_cell_to_point_global_ids_(refined_cell_to_point_global_ids)
        , scatter_refined_cell_to_point_global_ids_(refined_cell_to_point_global_ids)
        , local_to_global_refined_cells_(local_to_global_refined_cells)
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
        const auto& global_id = local_to_global_refined_cells_[element.index()];
        // check at
        if (gather_refined_cell_to_point_global_ids_.find(global_id)
            != gather_refined_cell_to_point_global_ids_.end()) {

            const auto& cell_to_point = gather_refined_cell_to_point_global_ids_.at(global_id);
            return cell_to_point.size(); // zero for overlap refined cells
        } else {
            return 0;
        }
    }

    // Gather global ids of points of an interior refined cell
    template <class B, class T> // T = Entity<0>
    void gather(B& buffer, const T& element)
    {
        // Gather only for interior cells.
        assert(element.partitionType() == Dune::InteriorEntity);
        // Get global id of the cell.
        const auto& global_id = local_to_global_refined_cells_[element.index()];
        if (gather_refined_cell_to_point_global_ids_.find(global_id)
            != gather_refined_cell_to_point_global_ids_.end()) {
            buffer.write(gather_refined_cell_to_point_global_ids_.at(global_id));
        }
    }

    // Scatter global ids of children cells of a coarse overlap cell
    template <class B, class T> // T = Entity<0>
    void scatter(B& buffer, const T& element, std::size_t)
    {
        // Scatter only for overlap cells.
        const auto& global_id = local_to_global_refined_cells_[element.index()];
        if (gather_refined_cell_to_point_global_ids_.find(global_id)
            != gather_refined_cell_to_point_global_ids_.end()) {
            DataType tmp = gather_refined_cell_to_point_global_ids_.at(global_id);
            buffer.read(tmp);
        
            if (element.partitionType() == Dune::OverlapEntity) {
                scatter_refined_cell_to_point_global_ids_[global_id] = tmp;
            }
        }
    }

private:
    const std::unordered_map<int, std::vector<std::size_t>>& gather_refined_cell_to_point_global_ids_;
    std::unordered_map<int, std::vector<std::size_t>>& scatter_refined_cell_to_point_global_ids_;
    const std::vector<int>& local_to_global_refined_cells_;
};
#endif // HAVE_MPI
} // namespace
#endif
