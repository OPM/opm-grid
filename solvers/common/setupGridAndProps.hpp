//===========================================================================
//
// File: setupGridAndProps.hpp
//
// Created: Tue Aug 11 14:47:35 2009
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

#ifndef OPENRS_SETUPGRIDANDPROPS_HEADER
#define OPENRS_SETUPGRIDANDPROPS_HEADER

#include <dune/common/param/ParameterGroup.hpp>
#include <dune/grid/CpGrid.hpp>
#include <dune/solvers/common/ReservoirPropertyCapillary.hpp>

namespace Dune
{

    /// @brief
    /// @todo Doc me!
    /// @param
    inline void setupGridAndProps(const parameter::ParameterGroup& param,
				  CpGrid& grid,
				  ReservoirPropertyCapillary<3>& res_prop)
    {
	// Initialize grid and reservoir properties.
	// Parts copied from CpGrid::init().
	std::string fileformat = param.getDefault<std::string>("fileformat", "eclipse");
	if (fileformat == "sintef_legacy") {
	    std::string grid_prefix = param.get<std::string>("grid_prefix");
	    grid.readSintefLegacyFormat(grid_prefix);
	    MESSAGE("Warning: We do not yet read legacy reservoir properties. Using defaults.");
	    res_prop.init(grid.size(0));
	} else if (fileformat == "eclipse") {
	    EclipseGridParser parser(param.get<std::string>("filename"));
	    double z_tolerance = param.getDefault<double>("z_tolerance", 0.0);
	    bool periodic_extension = param.getDefault<bool>("periodic_extension", false);
	    grid.processEclipseFormat(parser, z_tolerance, periodic_extension);
	    std::string rock_list = param.getDefault<std::string>("rock_list", "no_list");
	    std::string* rl_ptr = (rock_list == "no_list") ? 0 : &rock_list;
	    res_prop.init(parser, grid.globalCell(), rl_ptr);
	} else if (fileformat == "cartesian") {
	    array<int, 3> dims = {{ param.get<int>("nx"),
				    param.get<int>("ny"),
				    param.get<int>("nz") }};
	    array<double, 3> cellsz = {{ param.get<double>("dx"),
					 param.get<double>("dy"),
					 param.get<double>("dz") }};
	    grid.createCartesian(dims, cellsz);
	    MESSAGE("Warning: For generated cartesian grids, we use default reservoir properties.");
	    res_prop.init(grid.size(0));
	} else {
	    THROW("Unknown file format string: " << fileformat);
	}
    }

} // namespace Dune


#endif // OPENRS_SETUPGRIDANDPROPS_HEADER
