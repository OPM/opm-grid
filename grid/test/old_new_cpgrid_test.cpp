//===========================================================================
//
// File: old_new_cpgrid_test.cpp
//
// Created: Fri Jun  25 10:33:18 2009
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

#include "../CpGrid.hpp"
#include <dune/grid/io/file/vtk/vtkwriter.hh>

using namespace Dune;

int main(int argc, char** argv)
{
    if (argc < 2) {
	std::cerr << "Usage: " << argv[0]
		  << "[eclipse grdecl filename ] "
		  << "[z-tolerance] "
		  << "[...]" << std::endl;
	exit(EXIT_FAILURE);
    }

    // readEclipseFormat
    CpGrid grid;
    const std::string filename = argv[1];
    double z_tolerance = 0.0;
    if (argc > 2) {
	z_tolerance = atof(argv[2]);
    }
    std::cout << "z-tolerance = " << z_tolerance << std::endl;
    grid.readEclipseFormat(filename, z_tolerance);

    Dune::VTKWriter<CpGrid::LeafGridView> vtkwriter(grid.leafView());
    vtkwriter.write("dune.vtk", Dune::VTKOptions::ascii);

    return 0;
}
