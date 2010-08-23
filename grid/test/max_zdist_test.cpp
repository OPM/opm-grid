//===========================================================================
//                                                                           
// File: max_zdist_test.cpp                                                  
//                                                                           
// Created: Mon Oct 12 13:08:12 2009                                         
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
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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

#include <dune/common/EclipseGridParser.hpp>
#include <dune/common/EclipseGridInspector.hpp>

using namespace Dune;

int main(int argc, char** argv)
{
    if (argc != 2) {
	std::cerr << "Usage: " << argv[0] << " grdecl-file\n";
	return EXIT_FAILURE;
    }

    EclipseGridParser parser(argv[1]);
    EclipseGridInspector insp(parser);

    boost::array<int, 3> dims = insp.gridSize();

    std::vector<int> actnum;
    if (!parser.hasField("ACTNUM")) {
	actnum.resize(dims[0]* dims[1]*dims[2], 1);
    } else {
	actnum = parser.getIntegerValue("ACTNUM");
    }

    double global_zdiff_max = 1e100;
    int n = 0;
    for (int k = 0; k < dims[2]; ++k) {
	for (int j = 0; j < dims[1]; ++j) {
	    for (int i = 0; i < dims[0]; ++i, ++n) {
		if (actnum[n]) {
		    boost::array<double, 8> cellz = insp.cellZvals(i, j, k);
		    double zdiff, zdiff_max = 0.0;
		    for (int dd = 0; dd < 4; ++dd) {
			zdiff = cellz[dd+4] - cellz[dd];
			if (zdiff > zdiff_max) {
			    zdiff_max = zdiff;
			}
		    }
		    if (zdiff_max != 0.0) {
			if (zdiff_max < global_zdiff_max ) {
			    global_zdiff_max = zdiff_max;
			}
		    }
		}
	    }
	}
    }
    std::cout.precision(17);
    std::cout << "\nMaximal z-tolerance = " << global_zdiff_max << std::endl;

}
