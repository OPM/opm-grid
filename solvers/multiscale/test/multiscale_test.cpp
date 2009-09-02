//===========================================================================
//
// File: multiscale_test.cpp
//
// Created: Wed Sep  2 13:13:12 2009
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
#include <dune/common/mpihelper.hh>
#include <dune/common/param/ParameterGroup.hpp>
#include <dune/grid/CpGrid.hpp>
#include <dune/solvers/common/ReservoirPropertyCapillary.hpp>
#include <dune/solvers/common/setupGridAndProps.hpp>
#include <dune/solvers/common/GridInterfaceEuler.hpp>
#include <dune/solvers/multiscale/MultiscaleFlowSolver.hpp>

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


template<int dim, class GI, class RI>
void test_flowsolver(const GI& g, const RI& r)
{
    typedef typename GI::CellIterator               CI;
    typedef typename CI::FaceIterator               FI;
    typedef Dune::BoundaryConditions<true, false>   FBC;
    typedef Dune::MultiscaleFlowSolver<GI, RI, FBC> FlowSolver;

    FlowSolver solver;

    typedef Dune::FlowBC BC;
    FBC flow_bc(7);
    //flow_bc.flowCond(1) = BC(BC::Dirichlet, 1.0*Dune::unit::barsa);
    //flow_bc.flowCond(2) = BC(BC::Dirichlet, 0.0*Dune::unit::barsa);
    flow_bc.flowCond(5) = BC(BC::Dirichlet, 100.0*Dune::unit::barsa);

    solver.init(g, r, flow_bc);
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

    typename CI::Vector gravity;
    gravity[0] = gravity[1] = 0.0;
    gravity[2] = Dune::unit::gravity;
    solver.solve(r, sat, flow_bc, src, gravity);

#if 0
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
#endif
}



int main(int argc , char ** argv)
{
    // Initialize MPI, finalize is done automatically on exit.
    Dune::MPIHelper::instance(argc,argv);
    // Get parameters.
    Dune::parameter::ParameterGroup param(argc, argv);
    // Make grid and reservoir properties.
    Dune::CpGrid grid;
    Dune::ReservoirPropertyCapillary<3> res_prop;
    Dune::setupGridAndProps(param, grid, res_prop);
    Dune::GridInterfaceEuler<Dune::CpGrid> g(grid);

    assign_permeability<3>(res_prop, g.numberOfCells(), 0.1*Dune::unit::darcy);
    test_flowsolver<3>(g, res_prop);

    return EXIT_SUCCESS;
}
