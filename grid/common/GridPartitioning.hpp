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
  Copyright 2009 SINTEF ICT, Applied Mathematics.
  Copyright 2009 Statoil ASA.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPENRS_GRIDPARTITIONING_HEADER
#define OPENRS_GRIDPARTITIONING_HEADER

#include <vector>
#include <boost/array.hpp>

namespace Dune
{

    class CpGrid;

    /// Partition a CpGrid based on (ijk) coordinates, with splitting to ensure that each partition is connected.
    /// @param[in] grid the grid to partition
    /// @param[in] initial_split the number of parts in which to partition the grid, in each cardinal direction.
    ///                          Their product is the expected number of partitions produced.
    /// @param[out] num_part the resulting number of partitions. This may be lower than expected,
    ///                      because of inactive cells, or higher than expected,
    ///                      because of splits to ensure connectedness.
    /// @param[out] cell_part a vector containing, for each cell, its partition number
    void partition(const CpGrid& grid,
		   const boost::array<int, 3>& initial_split,
		   int& num_part,
		   std::vector<int>& cell_part,
		   bool recursive = false);

} // namespace Dune


#endif // OPENRS_GRIDPARTITIONING_HEADER
