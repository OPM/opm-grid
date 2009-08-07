//===========================================================================
//
// File: SimulatorTester.hpp
//
// Created: Fri Aug  7 09:21:55 2009
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

#ifndef OPENRS_SIMULATORTESTER_HEADER
#define OPENRS_SIMULATORTESTER_HEADER




#include "../EulerUpstream.hpp"
#include "../GridInterfaceEuler.hpp"
#include "../ReservoirPropertyCapillary.hpp"
#include "../BoundaryConditions.hpp"
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



    class SimulatorTester
    {
    public:
	typedef FieldVector<double, 3> Vector;

	SimulatorTester()
	    : simulation_steps_(1),
	      stepsize_(1.0*Dune::units::DAYS2SECONDS)
	{
	}

	void init(const parameter::ParameterGroup& param)
	{
	    simulation_steps_ = param.getDefault("simulation_steps", simulation_steps_);
	    stepsize_ = param.getDefault("stepsize", stepsize_)*Dune::units::DAYS2SECONDS;
	    // Initialize grid and reservoir properties.
	    // Parts copied from CpGrid::init().
	    std::string fileformat = param.getDefault<std::string>("fileformat", "eclipse");
	    if (fileformat == "sintef_legacy") {
		std::string grid_prefix = param.get<std::string>("grid_prefix");
		grid_.readSintefLegacyFormat(grid_prefix);
		MESSAGE("Warning: We do not yet read legacy reservoir properties. Using defaults.");
		res_prop_.init(grid_.size(0));
	    } else if (fileformat == "eclipse") {
		EclipseGridParser parser(param.get<std::string>("filename"));
		double z_tolerance = param.getDefault<double>("z_tolerance", 0.0);
		std::string rock_list = param.getDefault<std::string>("rock_list", "no_list");
		grid_.processEclipseFormat(parser, z_tolerance);
		std::string* rl_ptr = (rock_list == "no_list") ? 0 : &rock_list;
		res_prop_.init(parser, grid_.globalCell(), rl_ptr);
	    } else if (fileformat == "cartesian") {
		array<int, 3> dims = {{ param.get<int>("nx"),
					param.get<int>("ny"),
					param.get<int>("nz") }};
		array<double, 3> cellsz = {{ param.get<double>("dx"),
					     param.get<double>("dy"),
					     param.get<double>("dz") }};
		grid_.createCartesian(dims, cellsz);
		MESSAGE("Warning: For generated cartesian grids, we use default reservoir properties.");
		res_prop_.init(grid_.size(0));
	    } else {
		THROW("Unknown file format string: " << fileformat);
	    }
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
		// Transport.
		transport_solver.transportSolve(sat, stepsize_, gravity, flow_solver_.getSolution());
		// Output.
		std::vector<double> cell_velocity;
		estimateCellVelocity(cell_velocity, flow_solver_.getSolution());
		Dune::VTKWriter<GridType::LeafGridView> vtkwriter(grid_.leafView());
		vtkwriter.addCellData(cell_velocity, "velocity");
		vtkwriter.addCellData(sat, "saturation");
		vtkwriter.write("testsolution-" + boost::lexical_cast<std::string>(i), Dune::VTKOptions::ascii);
	    }
	}


	template <class CellData>
	void output(const std::string& filename, const std::string& fieldname, const CellData& celldata)
	{
	    // VTK output.
	    Dune::VTKWriter<typename GridType::LeafGridView> vtkwriter(grid_.leafView());
	    vtkwriter.addCellData(celldata, fieldname);
	    vtkwriter.write(filename, Dune::VTKOptions::ascii);
	    // Dumping the saturation.
// 	    std::ofstream os((filename + "_" + fieldname).c_str());
// 	    std::copy(celldata.begin(), celldata.end(), std::ostream_iterator<double>(os, "\n"));
	}




    private:
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
    };



} // namespace Dune


#endif // OPENRS_SIMULATORTESTER_HEADER
