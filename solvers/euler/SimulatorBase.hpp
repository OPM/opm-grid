//===========================================================================
//
// File: SimulatorBase.hpp
//
// Created: Tue Aug 11 15:01:48 2009
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

#ifndef OPENRS_SIMULATORBASE_HEADER
#define OPENRS_SIMULATORBASE_HEADER




#include <dune/solvers/euler/EulerUpstream.hpp>
#include <dune/solvers/euler/GridInterfaceEuler.hpp>
#include <dune/solvers/euler/ReservoirPropertyCapillary.hpp>
#include <dune/solvers/euler/BoundaryConditions.hpp>
#include <dune/solvers/euler/setupGridAndProps.hpp>
#include <dune/common/param/ParameterGroup.hpp>
#include <dune/grid/CpGrid.hpp>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/common/SparseVector.hpp>
#include <dune/grid/common/SparseTable.hpp>
#include <dune/grid/common/Volumes.hpp>
#include <dune/grid/common/Units.hpp>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/solvers/mimetic/IncompFlowSolverHybrid.hpp>
#include <dune/solvers/mimetic/MimeticIPEvaluator.hpp>
#include <fstream>
#include <iterator>
#include <boost/lexical_cast.hpp>

namespace Dune
{



    class SimulatorBase
    {
    public:
	typedef FieldVector<double, 3> Vector;

	SimulatorBase()
	    : simulation_steps_(1),
	      stepsize_(1.0*Dune::units::DAYS2SECONDS)
	{
	}

	void init(const parameter::ParameterGroup& param)
	{
	    simulation_steps_ = param.getDefault("simulation_steps", simulation_steps_);
	    stepsize_ = param.getDefault("stepsize", stepsize_)*Dune::units::DAYS2SECONDS;

	    setupGridAndProps(param, grid_, res_prop_);

	    // Make grid interface.
	    ginterf_.init(grid_);
	    // Make flow equation boundary conditions.
	    // Default is pressure 1.0e5 on the left, 0.0 on the right.
	    // Recall that the boundary ids range from 1 to 6 for the cartesian edges,
	    // and that boundary id 0 means interiour face/intersection.
	    flow_bcond_.resize(7);
	    std::string flow_bc_type = param.getDefault<std::string>("flow_bc_type", "dirichlet");
	    FBC::BCType bct = FBC::Dirichlet;
	    double leftval = 1.0e5;
	    double rightval = 0.0;
	    if (flow_bc_type == "neumann") {
		bct = FBC::Neumann;
		leftval = param.get<double>("left_flux");
		rightval = param.getDefault<double>("right_flux", -leftval);
	    } else if (flow_bc_type == "dirichlet") {
		leftval = param.getDefault<double>("left_pressure", leftval);
		rightval = param.getDefault<double>("right_pressure", rightval);
	    } else {
		THROW("Unknown flow boundary condition type " << flow_bc_type);
	    }
	    flow_bcond_[1] = FBC(bct, leftval);
	    flow_bcond_[2] = FBC(bct, rightval);
	    // Make transport equation boundary conditions.
	    // The default ones are fine (sat = 1.0 on inflow).
	    transport_bcond_.resize(7); // Again 7 conditions, see comment above.
	    // Make flow solver.
	    flow_solver_.init(ginterf_);
	    flow_solver_.assembleStatic(ginterf_, res_prop_);
	}

	template <class FlowSol>
	void estimateCellVelocity(std::vector<double>& cell_velocity,
				  const FlowSol& flow_solution)
	{
	    // Algorithm used is same as in halfFaceFluxToCellVelocity.hpp
	    // in the Sintef legacy c++ code.
	    cell_velocity.clear();
	    cell_velocity.resize(ginterf_.numberOfCells());
	    for (CellIter c = ginterf_.cellbegin(); c != ginterf_.cellend(); ++c) {
		int numf = 0;
		Vector cell_v(0.0);
		for (FaceIter f = c->facebegin(); f != c->faceend(); ++f, ++numf) {
		    double flux = flow_solution.outflux(f);
		    Vector v = f->centroid();
		    v -= c->centroid();
		    v *= flux/c->volume();
		    cell_v += v;
		}
		cell_velocity[c->index()] = cell_v.two_norm();
	    }
	}


	template <class FlowSol>
	void getCellPressure(std::vector<double>& cell_pressure,
			     const FlowSol& flow_solution)
	{
	    cell_pressure.clear();
	    cell_pressure.resize(ginterf_.numberOfCells());
	    for (CellIter c = ginterf_.cellbegin(); c != ginterf_.cellend(); ++c) {
		cell_pressure[c->index()] = flow_solution.pressure(c);
	    }
	}


	void run()
	{
	    // No injection or production.
	    SparseVector<double> injection_rates(ginterf_.numberOfCells());
	    std::vector<double> src(ginterf_.numberOfCells());
	    // Make a transport solver.
	    TransportSolver transport_solver(ginterf_, res_prop_, transport_bcond_, injection_rates);
	    // Initial saturation.
	    std::vector<double> sat(ginterf_.numberOfCells(), 0.0);
	    // Gravity.
	    FieldVector<double, 3> gravity(0.0);
	    // gravity[2] = -9.81;
	    // Compute flow field.
	    if (gravity.two_norm() > 0.0) {
		MESSAGE("Warning: Gravity not handled by flow solver.");
	    }

	    // Solve some steps.
	    for (int i = 0; i < simulation_steps_; ++i) {
		std::cout << "================    Simulation step number " << i << "    ===============" << std::endl;
		// Flow.
		flow_solver_.solve(ginterf_, res_prop_, sat, flow_bcond_, src);
// 		if (i == 0) {
// 		    flow_solver_.printSystem("linsys_dump_mimetic");
// 		}
		// Transport.
		transport_solver.transportSolve(sat, stepsize_, gravity, flow_solver_.getSolution());
		// Output.
		std::vector<double> cell_velocity;
		estimateCellVelocity(cell_velocity, flow_solver_.getSolution());
		std::vector<double> cell_pressure;
		getCellPressure(cell_pressure, flow_solver_.getSolution());
		Dune::VTKWriter<GridType::LeafGridView> vtkwriter(grid_.leafView());
		vtkwriter.addCellData(cell_velocity, "velocity");
		vtkwriter.addCellData(sat, "saturation");
		vtkwriter.addCellData(cell_pressure, "pressure");
		vtkwriter.write("testsolution-" + boost::lexical_cast<std::string>(i), Dune::VTKOptions::ascii);
	    }
	}

    protected:
	typedef CpGrid                                      GridType;
	typedef GridInterfaceEuler<GridType>                GridInterface;
	typedef EulerUpstream<GridInterface, ReservoirPropertyCapillary<3>, SaturationBoundaryConditions> TransportSolver;
	typedef GridInterface::CellIterator                 CellIter;
	typedef CellIter::FaceIterator                      FaceIter;
	typedef Dune::MimeticIPEvaluator<CellIter, 3, true> InnerProd;
	typedef Dune::FlowBoundaryCondition                 FBC;
	typedef Dune::FlowBoundaryConditions                FBCs;
	typedef Dune::IncompFlowSolverHybrid<GridInterface, FBCs, InnerProd> FlowSolver;

	int simulation_steps_;
	double stepsize_;
	GridType grid_;
	GridInterface ginterf_;
	ReservoirPropertyCapillary<3> res_prop_;
	FBCs flow_bcond_;
	SaturationBoundaryConditions transport_bcond_;
	FlowSolver flow_solver_;
	//TransportSolver transport_solver_;
    };



} // namespace Dune



#endif // OPENRS_SIMULATORBASE_HEADER
