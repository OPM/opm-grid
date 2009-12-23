//===========================================================================
//
// File: make_vtk_test.cpp
//
// Created: Mon Aug 17 16:28:33 2009
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
    std::vector<double> porosity(grid.size(0));
    std::vector<double> perm_xx(grid.size(0));
    std::vector<double> perm_yy(grid.size(0));
    std::vector<double> perm_zz(grid.size(0));
    for (int i = 0; i < grid.size(0); ++i) {
	porosity[i] = res_prop.porosity(i);
	perm_xx[i] = res_prop.permeability(i)(0,0);
	perm_yy[i] = res_prop.permeability(i)(1,1);
	perm_zz[i] = res_prop.permeability(i)(2,2);
    }
    VTKWriter<CpGrid::LeafGridView> vtkwriter(grid.leafView());
    vtkwriter.addCellData(porosity, "porosity");
    vtkwriter.addCellData(perm_xx, "perm_xx");
    vtkwriter.addCellData(perm_yy, "perm_yy");
    vtkwriter.addCellData(perm_zz, "perm_zz");
    std::string fname = param.get<std::string>("filename");
    vtkwriter.write(fname.substr(0, fname.find_last_of('.')) + "-cells", VTKOptions::ascii);
}

