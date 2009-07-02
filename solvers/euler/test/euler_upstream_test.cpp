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

#define VERBOSE

#include "config.h"
#include "../EulerUpstream.hpp"
#include "../GridInterfaceEuler.hpp"
#include "../ReservoirPropertyInterface.hpp"
#include "../BoundaryConditions.hpp"
#include <dune/grid/CpGrid.hpp>
#include <dune/grid/yaspgrid.hh>
#include <dune/common/mpihelper.hh>
#include <dune/grid/common/SparseVector.hpp>
#include <dune/grid/common/SparseTable.hpp>
#include <dune/grid/common/Volumes.hpp>

namespace Dune
{

    template <class GridInterface>
    class TestSolution
    {
    public:
	TestSolution(const GridInterface& g,
		     const typename GridInterface::CellIterator::Vector& v)
	{
	    // Make a flux consistent with a constant velocity field.
	    typename GridInterface::CellIterator c = g.cellbegin();
	    std::vector<double> cell_fluxes;
	    for (; c != g.cellend(); ++c) {
		cell_fluxes.clear();
		typename GridInterface::CellIterator::FaceIterator f = c->facebegin();
		for (; f != c->faceend(); ++f) {
		    //double flux = Dune::template inner<typename GridInterface::CellIterator::Vector>(v, f->normal())*f->area();
		    double flux = inner(v, f->normal())*f->area();
		    cell_fluxes.push_back(flux);
		}
		halfface_fluxes_.appendRow(cell_fluxes.begin(), cell_fluxes.end());
	    }
	}

	typedef typename GridInterface::CellIterator::FaceIterator FaceIter;
	double outflux(const FaceIter& f) const
	{
	    return halfface_fluxes_[f->cellIndex()][f->localIndex()];
	}

    private:
	SparseTable<double> halfface_fluxes_;
    };

} // namespace Dune


using namespace Dune;


int main(int argc, char** argv)
{
    parameter::ParameterGroup param(argc, argv);
    MPIHelper::instance(argc,argv);

    // Make a grid.
//     const int dim = 3;
//     typedef YaspGrid<dim> GridType;
//     typedef FieldVector<int,dim> iTuple;
//     typedef FieldVector<double,dim> fTuple;
//     typedef FieldVector<bool,dim> bTuple;
//     fTuple cell_sz(1.0);
//     iTuple dims(3);
//     bTuple periodic(false);
//     int overlap = 1;
//     YaspGrid<dim> grid(cell_sz, dims, periodic, overlap);
//     grid.globalRefine(2);

    // Make a grid
    typedef CpGrid GridType;
    CpGrid grid;
    EclipseGridParser parser(param.get<std::string>("filename"));
    double z_tolerance = param.getDefault<double>("z_tolerance", 0.0);
    grid.processEclipseFormat(parser, z_tolerance);

    // Make the grid interface
    typedef GridInterfaceEuler<GridType> GridInterface;
    GridInterface g(grid);

    // Reservoir properties.
    ReservoirPropertyInterface<3> res_prop;
    res_prop.init(parser);

    // Make flow equation boundary conditions.
    // Pressure 1.0e5 on the left, 0.0 on the right.
    // Recall that the boundary ids range from 1 to 6 for the cartesian edges,
    // and that boundary id 0 means interiour face/intersection.
    typedef FlowBoundaryCondition BC;
    FlowBoundaryConditions flow_bcond(7);
    flow_bcond[1] = BC(BC::Dirichlet, 1.0e5);
    flow_bcond[2] = BC(BC::Dirichlet, 0.0);

    // Make transport equation boundary conditions.
    // The default one is fine (sat = 1.0 on inflow).
    SaturationBoundaryConditions sat_bcond(7);

    // No injection or production.
    SparseVector<double> injection_rates(g.numberOfCells());

    // Make a solver.
    typedef EulerUpstream<GridInterface, ReservoirPropertyInterface<3>, SaturationBoundaryConditions> TransportSolver;
    TransportSolver transport_solver(g, res_prop, sat_bcond, injection_rates);

    // Define a flow field with constant velocity
    FieldVector<double, 3> vel(0.0);
    vel[0] = 2.0;
    vel[1] = 1.0;
    TestSolution<GridInterface> flow_solution(g, vel);
    // Solve a step.
    double time = 1.0;
    std::vector<double> sat(g.numberOfCells(), 0.0);
    FieldVector<double, 3> gravity(0.0);
    gravity[2] = -9.81;
    transport_solver.transportSolve(sat, time, gravity, flow_solution);
}

