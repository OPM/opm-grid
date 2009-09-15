//===========================================================================
//
// File: parsolver_test.cpp
//
// Created: Tue Sep  8 10:34:37 2009
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

#ifdef HAVE_MPI // #else clause at bottom of file


#define VERBOSE

#include "../ParIncompFlowSolverHybrid.hpp"
#include <dune/common/mpihelper.hh>
#include <dune/common/param/ParameterGroup.hpp>
#include <dune/grid/CpGrid.hpp>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/solvers/common/ReservoirPropertyCapillary.hpp>
#include <dune/solvers/common/setupGridAndProps.hpp>
#include <dune/solvers/common/GridInterfaceEuler.hpp>
#include <dune/solvers/common/SimulatorUtilities.hpp>
#include <dune/solvers/mimetic/MimeticIPEvaluator.hpp>
#include <dune/grid/common/GridPartitioning.hpp>


template<int dim, class RI>
void assign_permeability(RI& r, int nc, double k)
{
    typedef typename RI::SharedPermTensor Tensor;
    for (int c = 0; c < nc; ++c) {
        Tensor K = r.permeabilityModifiable(c);
        for (int i = 0; i < dim; ++i) {
            K(i,i) = k;
        }
    }
}


template<int dim, class Grid, class RI>
void test_flowsolver(const Grid& grid, const RI& r,
		     const std::vector<int>& partition,
		     int rank,
		     bool output_is_vtk = true)
{
    typedef Dune::GridInterfaceEuler<Grid>          GI;
    typedef typename GI::CellIterator               CI;
    typedef typename CI::FaceIterator               FI;
    typedef Dune::BoundaryConditions<true, false>   FBC;
    typedef Dune::ParIncompFlowSolverHybrid<GI, RI, FBC, Dune::MimeticIPEvaluator> FlowSolver;

    GI g(grid, true);
    FlowSolver solver;

    typedef Dune::FlowBC BC;
    FBC flow_bc(7);
    //flow_bc.flowCond(1) = BC(BC::Dirichlet, 1.0*Dune::unit::barsa);
    //flow_bc.flowCond(2) = BC(BC::Dirichlet, 0.0*Dune::unit::barsa);
    flow_bc.flowCond(5) = BC(BC::Dirichlet, 100.0*Dune::unit::barsa);
    flow_bc.flowCond(6) = BC(BC::Dirichlet, 0.0);

    solver.init(g, r, flow_bc, partition, rank);
    // solver.printStats(std::cout);

    //solver.assembleStatic(g, r);
    //solver.printIP(std::cout);

    std::vector<double> src(g.numberOfCells(), 0.0);
    std::vector<double> sat(g.numberOfCells(), 0.0);
#if 0
    if (g.numberOfCells() > 1) {
        src[0]     = 1.0;
        src.back() = -1.0;
    }
#endif

    typename CI::Vector gravity(0.0);
    //gravity[2] = Dune::unit::gravity;
    solver.solve(r, sat, flow_bc, src, gravity);

    if (output_is_vtk) {
	std::vector<double> cell_velocity;
	estimateCellVelocity(cell_velocity, g, solver.getSolution());
	std::vector<double> cell_pressure;
	getCellPressure(cell_pressure, g, solver.getSolution());
	Dune::VTKWriter<Dune::CpGrid::LeafGridView> vtkwriter(grid.leafView());
	vtkwriter.addCellData(cell_velocity, "velocity");
	vtkwriter.addCellData(cell_pressure, "pressure");
	vtkwriter.write("parsolver_test_output", Dune::VTKOptions::ascii);
    } else {
	solver.printSystem("system");
	typedef typename FlowSolver::SolutionType FlowSolution;
	FlowSolution soln = solver.getSolution();
	std::cout << "Cell Pressure:\n" << std::scientific << std::setprecision(15);
	for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
	    std::cout << '\t' << soln.pressure(c) << '\n';
	}
	std::cout << "Cell (Out) Fluxes:\n";
	std::cout << "flux = [\n";
	for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
	    for (FI f = c->facebegin(); f != c->faceend(); ++f) {
		std::cout << soln.outflux(f) << ' ';
	    }
	    std::cout << "\b\n";
	}
	std::cout << "]\n";
    }
}



int main(int argc , char ** argv)
{
    // Initialize MPI, finalize is done automatically on exit.
    Dune::MPIHelper::instance(argc,argv).rank();
    int mpi_rank = 1;//Dune::MPIHelper::instance(argc,argv).rank();
    int mpi_size = 2;//Dune::MPIHelper::instance(argc,argv).size();
    std::cout << "Hello from rank " << mpi_rank << std::endl;

    Dune::CpGrid grid;
    Dune::ReservoirPropertyCapillary<3> res_prop;
    // Get parameters.
    Dune::parameter::ParameterGroup param(argc, argv);
    // Make grid and reservoir properties.
    Dune::setupGridAndProps(param, grid, res_prop);
    assign_permeability<3>(res_prop, grid.size(0), 0.1*Dune::unit::darcy);
    // Partitioning.
//     boost::array<int, 3> split = {{ param.getDefault("sx", 1), 
// 				    param.getDefault("sy", 1),
// 				    param.getDefault("sz", 1) }};
    boost::array<int, 3> split = {{ 1, 1, mpi_size }};
    int num_part = 0;
    std::vector<int> partition;
    Dune::partition(grid, split, num_part, partition);
    if (num_part != mpi_size) {
	std::cerr << "Right now, parsolver_test() wants as many partitions as mpi tasks.\n"
		  << "Number of partitions: " << num_part << "   Mpi parallells: " << mpi_size << '\n';
	return EXIT_FAILURE;
    }
    test_flowsolver<3>(grid, res_prop, partition, mpi_rank);

    return EXIT_SUCCESS;
}

#else

// We do not have MPI

#include <iostream>

int main()
{
    std::cerr << "This program does nothing if HAVE_MPI is undefined.\n"
	"To enable MPI, pass --enable-parallel to configure (or dunecontrol) when setting up dune.\n";
}


#endif
