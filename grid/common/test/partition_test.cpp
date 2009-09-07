//===========================================================================
//
// File: partition_test.cpp
//
// Created: Mon Sep  7 12:14:13 2009
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

#include "../GridPartitioning.hpp"
#include <dune/grid/CpGrid.hpp>
#include <cstdlib>

using namespace Dune;

int main()
{
    // Make grid.
    CpGrid g;
    array<int, 3> dims = {{ 2, 3, 4 }};
    array<double, 3> sizes = {{ 1.0, 1.0, 1.0 }};
    g.createCartesian(dims, sizes);

    // Partition.
    boost::array<int, 3> split = {{ 2, 2, 2 }};
    int num_part = -1;
    std::vector<int> cell_part;
    partition(g, split, num_part, cell_part);

    // Check.
    if (num_part != 8) {
	return EXIT_FAILURE;
    }
    int cell_part_correct[] = { 0, 0, 4, 4, 1, 1, 5, 5,
				0, 0, 4, 4, 1, 1, 5, 5,
				2, 2, 6, 6, 3, 3, 7, 7 };
    if (!std::equal(cell_part.begin(), cell_part.end(), cell_part_correct)) {
	return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
