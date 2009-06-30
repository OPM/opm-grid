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

using namespace Dune;

int main(int argc, char** argv)
{
    if (argc < 2) {
	std::cerr << "Usage: " << argv[0]
			  << "[eclipse grdecl filename ] "
			  << "[...]" << std::endl;
	exit(EXIT_FAILURE);
    }

    CpGrid grid;
    const std::string filename = argv[1];
    grid.readEclipseFormat(filename, 0.0);

//     CpGrid slf_grid;
//     const std::string slf_grid_prefix = argv[2];
//     slf_grid.readSintefLegacyFormat(slf_grid_prefix);
}
