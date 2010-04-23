//===========================================================================
//
// File: mimetic_solver_test.cpp
//
// Created: Tue Jul  7 15:35:21 2009
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

#include <config.h>

#include <algorithm>
#include <iostream>
#include <iomanip>

#include <boost/static_assert.hpp>

#include <dune/common/array.hh>
#include <dune/common/mpihelper.hh>
#include <dune/common/Units.hpp>

#if HAVE_ALUGRID
#include <dune/common/shared_ptr.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/alugrid.hh>
#endif

#include <dune/solvers/common/SimulatorUtilities.hpp>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/CpGrid.hpp>
#include <dune/grid/common/EclipseGridParser.hpp>
#include <dune/grid/common/EclipseGridInspector.hpp>

#include <dune/solvers/common/fortran.hpp>
#include <dune/solvers/common/blas_lapack.hpp>
#include <dune/solvers/common/Matrix.hpp>
#include <dune/solvers/common/GridInterfaceEuler.hpp>
#include <dune/solvers/common/ReservoirPropertyCapillary.hpp>
#include <dune/solvers/common/BoundaryConditions.hpp>

#include <dune/solvers/mimetic/MimeticIPEvaluator.hpp>
#include <dune/solvers/mimetic/IncompFlowSolverHybrid.hpp>
#include <dune/common/param/ParameterGroup.hpp>


void build_grid(const Dune::EclipseGridParser& parser,
                const double z_tol, Dune::CpGrid& grid,
                boost::array<int,3>& cartDims)
{
    Dune::EclipseGridInspector insp(parser);

    grdecl g;
    cartDims[0] = g.dims[0] = insp.gridSize()[0];
    cartDims[1] = g.dims[1] = insp.gridSize()[1];
    cartDims[2] = g.dims[2] = insp.gridSize()[2];

    g.coord = &parser.getFloatingPointValue("COORD")[0];
    g.zcorn = &parser.getFloatingPointValue("ZCORN")[0];

    if (parser.hasField("ACTNUM")) {
        g.actnum = &parser.getIntegerValue("ACTNUM")[0];
        grid.processEclipseFormat(g, z_tol, false, false);
    } else {
        std::vector<int> dflt_actnum(g.dims[0] * g.dims[1] * g.dims[2], 1);

        g.actnum = &dflt_actnum[0];
        grid.processEclipseFormat(g, z_tol, false, false);
    }
}


#if USE_ALUGRID
template<class GType>
Dune::shared_ptr<GType>
make_gmsh(const std::string& msh_file)
{
    return Dune::shared_ptr<GType>(Dune::GmshReader<GType>::read(msh_file));
}
#endif


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
    typedef typename GI::CellIterator                   CI;
    typedef typename CI::FaceIterator                   FI;
    typedef Dune::BasicBoundaryConditions<true, false>  FBC;
    typedef Dune::IncompFlowSolverHybrid<GI, RI, FBC,
                                         Dune::MimeticIPEvaluator> FlowSolver;

    FlowSolver solver;

    typedef Dune::FlowBC BC;
    FBC flow_bc(7);

#if !USE_ALUGRID
    //flow_bc.flowCond(1) = BC(BC::Dirichlet, 1.0*Dune::unit::barsa);
    //flow_bc.flowCond(2) = BC(BC::Dirichlet, 0.0*Dune::unit::barsa);
    flow_bc.flowCond(5) = BC(BC::Dirichlet, 100.0*Dune::unit::barsa);
#endif

    typename CI::Vector gravity;
    gravity[0] = gravity[1] = 0.0;
    gravity[2] = Dune::unit::gravity;

    solver.init(g, r, gravity, flow_bc);

    std::vector<double> src(g.numberOfCells(), 0.0);
    std::vector<double> sat(g.numberOfCells(), 0.0);
#if 1
    if (g.numberOfCells() > 1) {
        src[0]     = 1.0;
        src.back() = -1.0;
    }
#endif

    solver.solve(r, sat, flow_bc, src, 5e-9, 3, 0);

#if 1
    typedef typename FlowSolver::SolutionType FlowSolution;
    FlowSolution soln = solver.getSolution();

    std::vector<typename GI::Vector> cell_velocity;
    estimateCellVelocity(cell_velocity, g, soln);
    // Dune's vtk writer wants multi-component data to be flattened.
    std::vector<double> cell_velocity_flat(&*cell_velocity.front().begin(),
                                           &*cell_velocity.back().end());
    std::vector<double> cell_pressure;
    getCellPressure(cell_pressure, g, soln);

    Dune::VTKWriter<typename GI::GridType::LeafGridView> vtkwriter(g.grid().leafView());
    vtkwriter.addCellData(cell_velocity_flat, "velocity", dim);
    vtkwriter.addCellData(cell_pressure, "pressure");
    vtkwriter.write("testsolution-" + boost::lexical_cast<std::string>(0),
                    Dune::VTKOptions::ascii);
#else    
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


using namespace Dune;

int main(int argc, char** argv)
{
    Dune::parameter::ParameterGroup param(argc, argv);
    Dune::MPIHelper::instance(argc,argv);

    // Make a grid
#if !USE_ALUGRID
    Dune::CpGrid grid;

    Dune::EclipseGridParser parser(param.get<std::string>("filename"));
    double z_tol = param.getDefault<double>("z_tolerance", 0.0);
    boost::array<int,3> cartDims;
    build_grid(parser, z_tol, grid, cartDims);

    // Make the grid interface
    Dune::GridInterfaceEuler<Dune::CpGrid> g(grid);

    // Reservoir properties.
    ReservoirPropertyCapillary<3> res_prop;
    res_prop.init(parser, grid.globalCell());
    assign_permeability<3>(res_prop, g.numberOfCells(), 0.1*Dune::unit::darcy);
#else
    typedef Dune::ALUSimplexGrid<3,3> GType;
    Dune::shared_ptr<GType> pgrid = make_gmsh<GType>(param.get<std::string>("filename"));
    Dune::GridInterfaceEuler<GType> g(*pgrid);

    ReservoirPropertyCapillary<3> res_prop;
    res_prop.init(g.numberOfCells());
#endif

    test_flowsolver<3>(g, res_prop);
}
