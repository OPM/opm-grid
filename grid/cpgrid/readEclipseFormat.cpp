//===========================================================================
//
// File: readEclipseFormat.cpp
//
// Created: Fri Jun 12 09:16:59 2009
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
#include "EclipseGridParser.hpp"
#include "EclipseGridInspector.hpp"
#include "../preprocess/grdecl.h"

namespace Dune
{


    /// Read the Sintef legacy grid format ('topogeom').
    void CpGrid::readEclipseFormat(const std::string& filename)
    {
	EclipseGridParser parser(filename);
	EclipseGridInspector inspector(parser);
	grdecl g;
	g.dims[0] = inspector.gridSize()[0];
	g.dims[1] = inspector.gridSize()[1];
	g.dims[2] = inspector.gridSize()[2];
	g.coord = &(parser.getFloatingPointValue("COORD")[0]);
	g.zcorn = &(parser.getFloatingPointValue("ZCORN")[0]);
	g.actnum = &(parser.getIntegerValue("ACTNUM")[0]);
    }

} // namespace Dune
