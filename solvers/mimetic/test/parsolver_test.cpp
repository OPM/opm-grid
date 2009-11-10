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

#if HAVE_MPI // #else clause at bottom of file


#define VERBOSE

#include "../ParIncompFlowSolverHybrid.hpp"
#include <dune/common/mpihelper.hh>
#include <dune/common/param/ParameterGroup.hpp>
#include <dune/grid/CpGrid.hpp>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/solvers/common/ReservoirPropertyCapillary.hpp>
#include <dune/solvers/common/setupGridAndProps.hpp>
#include <dune/solvers/common/GridInterfaceEuler.hpp>
#include <dune/solvers/common/SimulatorUtilities.hpp>
#include <dune/solvers/mimetic/MimeticIPEvaluator.hpp>
#include <dune/grid/common/GridPartitioning.hpp>



template<class Grid, class RI>
void test_flowsolver(const Grid& grid,
		     const RI& r,
		     const std::vector<int>& partition,
		     int rank,
		     double tol = 1e-8,
		     int verbosity = 1,
		     bool output_vtk = true)
{
    typedef Dune::GridInterfaceEuler<Grid>          GI;
    typedef typename GI::CellIterator               CI;
    typedef typename CI::FaceIterator               FI;
    typedef Dune::BasicBoundaryConditions<true, false>   FBC;
    typedef Dune::ParIncompFlowSolverHybrid<GI, RI, FBC, Dune::MimeticIPEvaluator> FlowSolver;

    GI g(grid, true);
    FlowSolver solver;

    typedef Dune::FlowBC BC;
    FBC flow_bc(2*Grid::dimension + 1);
    //flow_bc.flowCond(1) = BC(BC::Dirichlet, 1.0*Dune::unit::barsa);
    //flow_bc.flowCond(2) = BC(BC::Dirichlet, 0.0*Dune::unit::barsa);
    flow_bc.flowCond(1) = BC(BC::Dirichlet, 1.0);//100.0*Dune::unit::barsa);
    flow_bc.flowCond(2) = BC(BC::Dirichlet, 2.0);//200.0*Dune::unit::barsa);

    solver.init(g, r, flow_bc, partition, rank);
    // solver.printStats(std::cout);
    // solver.assembleStatic(g, r);
    // solver.printIP(std::cout);

    std::vector<double> src(g.numberOfCells(), 0.0);
    std::vector<double> sat(g.numberOfCells(), 1.0);
#if 0
    if (g.numberOfCells() > 1) {
        src[0]     = 1.0;
        src.back() = -1.0;
    }
#endif

    typename CI::Vector gravity(0.0);
    //gravity[2] = Dune::unit::gravity;
    solver.solve(r, sat, flow_bc, src, gravity, tol, verbosity);

    if (output_vtk) {
	std::vector<double> cell_velocity;
	estimateCellVelocity(cell_velocity, g, solver.getSolution(), partition, rank);
	std::vector<double> cell_pressure;
	getCellPressure(cell_pressure, g, solver.getSolution(), partition, rank);
	if (rank == 0) {
	    Dune::VTKWriter<typename Grid::LeafGridView> vtkwriter(grid.leafView());
	    vtkwriter.addCellData(cell_velocity, "velocity");
	    vtkwriter.addCellData(cell_pressure, "pressure");
	    vtkwriter.write("parsolver_test_output", Dune::VTKOptions::ascii);
	}
    }
}



int main(int argc , char ** argv)
{
    // Initialize MPI, finalize is done automatically on exit.
    Dune::MPIHelper::instance(argc,argv).rank();
    int mpi_rank = Dune::MPIHelper::instance(argc,argv).rank();
    int mpi_size = Dune::MPIHelper::instance(argc,argv).size();
    std::cout << "Hello from rank " << mpi_rank << std::endl;

    // Get parameters.
    Dune::parameter::ParameterGroup param(argc, argv);
    mpi_rank = param.getDefault("force_rank", mpi_rank);
    mpi_size = param.getDefault("force_size", mpi_size);
    double tol = param.getDefault("residual_tolerance", 1e-8);
    int verbosity = param.getDefault("linsolver_verbosity", 1);

    // Make grid and reservoir properties.
    Dune::CpGrid grid;
    Dune::ReservoirPropertyCapillary<3> res_prop;
    Dune::setupGridAndProps(param, grid, res_prop);
    // Partitioning.
//     boost::array<int, 3> split = {{ param.getDefault("sx", 1), 
// 				    param.getDefault("sy", 1),
// 				    param.getDefault("sz", 1) }};
    boost::array<int, 3> split = {{ mpi_size, 1, 1 }};
    int num_part = 0;
    std::vector<int> partition;
    Dune::partition(grid, split, num_part, partition);
    if (num_part != mpi_size) {
	std::cerr << "Right now, parsolver_test() wants as many partitions as mpi tasks.\n"
		  << "Number of partitions: " << num_part << "   Mpi parallells: " << mpi_size << '\n';
	return EXIT_FAILURE;
    }
    const int sgrid_dim = 2;
    const int n[3] = { param.getDefault("nx", 1),
		       param.getDefault("ny", 1),
		       param.getDefault("nz", 1) };
    double d[3] =  { param.getDefault<double>("dx", 1.0)*n[0],
		     param.getDefault<double>("dy", 1.0)*n[1],
		     param.getDefault<double>("dz", 1.0)*n[2] };
    Dune::SGrid<sgrid_dim, sgrid_dim> sgrid(n, d);
    Dune::ReservoirPropertyCapillary<sgrid_dim> sres_prop;
    sres_prop.init(sgrid.size(0), 0.2, 1e-3);
    test_flowsolver(sgrid, sres_prop, partition, mpi_rank, tol, verbosity);
    // test_flowsolver(grid, res_prop, partition, mpi_rank, tol, verbosity);

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
