//===========================================================================
//
// File: euler_upstream_test.cpp
//
// Created: Tue Jun 16 14:31:39 2009
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
#include "../EulerUpstream.hpp"
#include "../GridInterfaceEuler.hpp"
#include "../ReservoirPropertyInterface.hpp"
#include "../BoundaryConditions.hpp"
#include <dune/grid/CpGrid.hpp>
#include <dune/grid/yaspgrid.hh>
#include <dune/common/mpihelper.hh>
#include <dune/grid/common/SparseVector.hpp>

int main(int argc, char** argv)
{
    Dune::parameter::ParameterGroup param(argc, argv);
    Dune::MPIHelper::instance(argc,argv);

    // Make a grid.
//     const int dim = 3;
//     typedef Dune::YaspGrid<dim> GridType;
//     typedef Dune::FieldVector<int,dim> iTuple;
//     typedef Dune::FieldVector<double,dim> fTuple;
//     typedef Dune::FieldVector<bool,dim> bTuple;
//     fTuple cell_sz(1.0);
//     iTuple dims(3);
//     bTuple periodic(false);
//     int overlap = 1;
//     Dune::YaspGrid<dim> grid(cell_sz, dims, periodic, overlap);
//     grid.globalRefine(2);

    // Make a grid
    typedef Dune::CpGrid GridType;
    Dune::CpGrid grid;
    Dune::EclipseGridParser parser(param.get<std::string>("filename"));
    double z_tolerance = param.getDefault<double>("z_tolerance", 0.0);
    grid.processEclipseFormat(parser, z_tolerance);

    // Make the grid interface
    Dune::GridInterfaceEuler<GridType> g(grid);

    // Reservoir properties.
    Dune::ReservoirPropertyInterface<3> res_prop;
    res_prop.init(parser);

    // Make pressure equation boundary conditions.
    // Pressure 1.0e5 on the left, 0.0 on the right.
    // Recall that the boundary ids range from 1 to 6 for the cartesian edges,
    // and that boundary id 0 means interiour face/intersection.
    typedef Dune::FlowBoundaryCondition BC;
    Dune::FlowBoundaryConditions bcond(7);
    bcond[1] = BC(BC::Dirichlet, 1.0e5);
    bcond[2] = BC(BC::Dirichlet, 0.0);

    // Make a solver.
    Dune::EulerUpstream transport_solver;


    //transport_solver.transportSolve(sat, time, resdata, boundary, face_fluxes);
}

