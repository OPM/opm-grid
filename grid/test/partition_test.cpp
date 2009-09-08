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

#include "config.h"
#include <dune/grid/common/GridPartitioning.hpp>
#include <dune/common/param/ParameterGroup.hpp>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/CpGrid.hpp>
#include <dune/solvers/common/ReservoirPropertyCapillary.hpp>
#include <dune/solvers/common/setupGridAndProps.hpp>
#include <cstdlib>

using namespace Dune;

int testPartition()
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


int main(int argc, char** argv)
{
    if (argc == 1) {
	// Running in test mode.
	return testPartition();
    }

    // Running in normal mode. Make grid.
    parameter::ParameterGroup param(argc, argv);
    CpGrid grid;
    ReservoirPropertyCapillary<3> res_prop;
    setupGridAndProps(param, grid, res_prop);

    // Partition.
    boost::array<int, 3> split = {{ param.getDefault("sx", 1), 
				    param.getDefault("sy", 1),
				    param.getDefault("sz", 1) }};
    int num_part = -1;
    std::vector<int> cell_partition;
    partition(grid, split, num_part, cell_partition);
    std::cout << "Grid with " << cell_partition.size()
	      << " cells was split into " << num_part
	      << " partitions." << std::endl;

    // Output.
    VTKWriter<CpGrid::LeafGridView> vtkwriter(grid.leafView());
    vtkwriter.addCellData(cell_partition, "partition");
    std::string fname = param.get<std::string>("filename");
    vtkwriter.write(fname.substr(0, fname.size() - 7) + "-partition", VTKOptions::ascii);

}
