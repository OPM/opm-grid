//===========================================================================
//
// File: mimetic_periodic_test.cpp
//
// Created: Wed Aug 19 13:42:04 2009
//
// Author(s): Bård Skaflestad     <bard.skaflestad@sintef.no>
//            Atgeirr F Rasmussen <atgeirr@sintef.no>
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

#include <iostream>

#include <boost/array.hpp>

#include <dune/common/Units.hpp>
#include <dune/common/param/ParameterGroup.hpp>

#include <dune/grid/CpGrid.hpp>

#include <dune/solvers/common/PeriodicHelpers.hpp>
#include <dune/solvers/common/BoundaryConditions.hpp>
#include <dune/solvers/common/GridInterfaceEuler.hpp>
#include <dune/solvers/common/ReservoirPropertyCapillary.hpp>

#include <dune/solvers/mimetic/MimeticIPEvaluator.hpp>
#include <dune/solvers/mimetic/IncompFlowSolverHybrid.hpp>

using namespace Dune;

int main(int argc, char** argv)
{
    typedef Dune::GridInterfaceEuler<CpGrid>           GI;
    typedef GI  ::CellIterator                         CI;
    typedef CI  ::FaceIterator                         FI;
    typedef Dune::MimeticIPEvaluator<CI,3,true>        IP;
    typedef Dune::BoundaryConditions<true, false>      BCs;
    typedef Dune::ReservoirPropertyCapillary<3>        RI;
    typedef Dune::IncompFlowSolverHybrid<GI,RI,BCs,IP> FlowSolver;

    parameter::ParameterGroup param(argc, argv);
    CpGrid grid;
    grid.init(param);
    grid.setUniqueBoundaryIds(true);
    GridInterfaceEuler<CpGrid> g(grid);
    typedef FlowBC FBC;
    boost::array<FBC, 6> cond = {{ FBC(FBC::Periodic,  1.0*unit::barsa),
                                   FBC(FBC::Periodic, -1.0*unit::barsa),
                                   FBC(FBC::Neumann,   0.0),
                                   FBC(FBC::Neumann,   0.0),
                                   FBC(FBC::Neumann,   0.0),
                                   FBC(FBC::Neumann,   0.0) }};
    BCs fbc;
    createPeriodic(fbc, g, cond);

    RI r;
    r.init(g.numberOfCells());

    FlowSolver solver;
    solver.init(g, r, fbc);

    std::vector<double> src(g.numberOfCells(), 0.0);
    std::vector<double> sat(g.numberOfCells(), 0.0);

    CI::Vector gravity;
    gravity[0] = gravity[1] = gravity[2] = 0.0;
#if 0
    gravity[2] = Dune::unit::gravity;
#endif
    solver.solve(r, sat, fbc, src, gravity);

    FlowSolver::SolutionType soln = solver.getSolution();
#if 1
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
    
    return 0;
}

