//===========================================================================
//
// File: CpGrid.cpp
//
// Created: Thu Jun  4 12:55:28 2009
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

#include <fstream>
#include "../CpGrid.hpp"

namespace Dune
{




    /// Initialize the grid.
    void CpGrid::init(const parameter::ParameterGroup& param)
    {
	std::string fileformat = param.get<std::string>("fileformat");
	if (fileformat == "sintef_legacy") {
	    std::string grid_prefix = param.get<std::string>("grid_prefix");
	    readSintefLegacyFormat(grid_prefix);
	} else if (fileformat == "eclipse") {
	    std::string filename = param.get<std::string>("filename");
	    double z_tolerance = param.getDefault<double>("z_tolerance", 0.0);
	    readEclipseFormat(filename, z_tolerance);
	} else {
	    THROW("Unknown file format string: " << fileformat);
	}
    }





} // namespace Dune
