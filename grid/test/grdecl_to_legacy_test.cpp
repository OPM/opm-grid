//===========================================================================
//
// File: grdecl_to_legacy_test.cpp
//
// Created: Thu Dec 17 11:19:44 2009
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


#if HAVE_CONFIG_H
#include "config.h"
#endif
#include <dune/common/param/ParameterGroup.hpp>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/CpGrid.hpp>
#include <dune/solvers/common/ReservoirPropertyCapillary.hpp>
#include <dune/solvers/common/setupGridAndProps.hpp>

using namespace Dune;

int main(int argc, char** argv)
{
    parameter::ParameterGroup param(argc, argv);
    CpGrid grid;
    ReservoirPropertyCapillary<3> res_prop;
    setupGridAndProps(param, grid, res_prop);

    std::string grid_prefix = param.get<std::string>("grid_prefix");
    grid.writeSintefLegacyFormat(grid_prefix);
    res_prop.writeSintefLegacyFormat(grid_prefix);
}

