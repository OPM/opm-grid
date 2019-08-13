//===========================================================================
//
// File: GridPartitioning.hpp
//
// Created: Mon Sep  7 10:09:13 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010, 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.
  Copyright 2013 Dr. Markus Blatt - HPC-Simulation-Software & Services

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

#ifndef OPM_GRIDPARTITIONING_HEADER
#define OPM_GRIDPARTITIONING_HEADER

#include <vector>
#include <array>
#include <set>

namespace Dune
{

    class CpGrid;

    struct OrderByFirst
    {
        bool operator()(const std::pair<int,int>& o, const std::pair<int,int>& v)
        {
            return o.first < v.first;
        }
    };

    /// Partition a CpGrid based on (ijk) coordinates, with splitting to ensure that each partition is connected.
    /// @param[in] grid the grid to partition
    /// @param[in] initial_split the number of parts in which to partition the grid, in each cardinal direction.
    ///                          Their product is the expected number of partitions produced.
    /// @param[out] num_part the resulting number of partitions. This may be lower than expected,
    ///                      because of inactive cells, or higher than expected,
    ///                      because of splits to ensure connectedness.
    /// @param[out] cell_part a vector containing, for each cell, its partition number
    void partition(const CpGrid& grid,
                   const std::array<int, 3>& initial_split,
                   int& num_part,
                   std::vector<int>& cell_part,
                   bool recursive = false,
                   bool ensureConnectivity = true);

/// \brief Adds a layer of overlap cells to a partitioning.
/// \param[in] grid The grid that is partitioned.
/// \param[in] cell_part a vector containing each cells partition number.
/// \param[out] cell_overlap a vector of sets that contains for each cell all
///             the partition numbers that it is an overlap cell of.
/// \param[in] mypart The partition number of the processor.
/// \param[in] all Whether to compute the overlap for all partions or just the
///            one associated by mypart.
     void addOverlapLayer(const CpGrid& grid,
                          const std::vector<int>& cell_part,
                          std::vector<std::set<int> >& cell_overlap,
                          int mypart, int overlapLayers, bool all=false);

/// \brief Marks cells of my_rank as either overlap or ghost/overlap cells and creates a natural ordering of local cells.
/// \param[in]  overlap a vector of sets that contains for each cell all the partition numbers that it is an overlap cell of.
/// \param[in]  cell_part a vector containing each cells partition number.
/// \param[in]  my_rank integer rank id of process.
/// \param[out] naturalOrder vector that on return has size equal to number of cells on process my_rank. Contains
///             global id of the local cells.
/// \param[out] partitionType vector of ints with size equal to total number of cells. On return it contains the partition
///             type of cell with global id gid in partitionType[gid]. There are three types: exterior=0, interior=2 and overlap=1.
    void findInteriorAndOverlapCells(std::vector<std::set<int>>& overlap, const std::vector<int>& cell_part, int my_rank,
                                     std::vector<int>& naturalOrder, std::vector<int>& partitionType);

/// \brief Reorder the local cells, so that ghost/overlap cells come after interior cells
/// \param[in] naturalOrder vector that contains global id of the local cells.
/// \param[in] pType vector that contains partition type of cell with global id gid in pType[gid]. 
/// \param[in] reorderMethod if 0 do not reorder. If 1 reorder.
/// \return vector containing new ordering of local cells.
    std::vector<int> reorderLocalCells(const std::vector<int>& naturalOrder, const std::vector<int>& pType,
                                       int reorderMethod);

} // namespace Dune


#endif // OPM_GRIDPARTITIONING_HEADER
