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


#include <fstream>
#include <iterator>

#include <boost/lexical_cast.hpp>

#include <dune/common/param/ParameterGroup.hpp>

#include <dune/common/SparseVector.hpp>
#include <dune/common/SparseTable.hpp>
#include <dune/common/Units.hpp>

#include <dune/grid/common/Volumes.hpp>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/CpGrid.hpp>
#include <dune/grid/yaspgrid.hh>

#include <dune/solvers/common/GridInterfaceEuler.hpp>
#include <dune/solvers/common/ReservoirPropertyCapillary.hpp>
#include <dune/solvers/common/BoundaryConditions.hpp>
#include <dune/solvers/common/setupGridAndProps.hpp>
#include <dune/solvers/common/setupBoundaryConditions.hpp>

#include <dune/solvers/euler/EulerUpstream.hpp>

#include <dune/solvers/mimetic/MimeticIPEvaluator.hpp>
#include <dune/solvers/mimetic/IncompFlowSolverHybrid.hpp>

namespace Dune
{




    /// @brief
    /// @todo Doc me!
    /// @tparam
    class SimulatorBase
    {
    public:

	/// @brief
	/// @todo Doc me!
	/// @tparam
	typedef FieldVector<double, 3> Vector;


	/// @brief
	/// @todo Doc me!
	SimulatorBase()
	    : simulation_steps_(1),
	      stepsize_(1.0),   // init() expects units of days! Yes, but now the meaning of stepsize_ changes
	                        // from days (here) to seconds (after init()). Solution to that?
	      init_saturation_(0.0)
	{
	}

	/// @brief Initialization from parameters.
	/// @param param a parameter object
	void init(const parameter::ParameterGroup& param)
	{
	    initControl(param);
	    initGridAndProps(param);
	    initInitialConditions(param);
	    initBoundaryConditions(param);
	    initSolvers(param);

	    // Write any unused parameters.
	    std::cout << "====================   Unused parameters:   ====================\n";
	    param.displayUsage();
	    std::cout << "================================================================\n";
	}

	/// @brief Estimates a scalar cell velocity from outgoing fluxes.
	/// @tparam FlowSol a flow solution type.
	/// @param[out] cell_velocity the estimated velocities.
	/// @param[in] flow_solution the object containing the fluxes.
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


	/// @brief
	/// @todo Doc me!
	/// @tparam
	/// @param
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

    protected:
	typedef CpGrid                                         GridType;
 	enum { Dimension = GridType::dimension };
 	typedef ReservoirPropertyCapillary<Dimension>          ResProp;
	typedef GridInterfaceEuler<GridType>                   GridInterface;
	typedef GridInterface::CellIterator                    CellIter;
	typedef CellIter::FaceIterator                         FaceIter;
	typedef MimeticIPEvaluator<CellIter, Dimension, true>  InnerProd;
	typedef BoundaryConditions<true, true>                 BCs;
	typedef IncompFlowSolverHybrid<GridInterface,
				       ResProp,
				       BCs,
				       InnerProd>              FlowSolver;
	typedef EulerUpstream<GridInterface,
			      ResProp,
			      BCs>                             TransportSolver;

	int simulation_steps_;
	double stepsize_;
	double init_saturation_;
	GridType grid_;
	GridInterface ginterf_;
	ReservoirPropertyCapillary<3> res_prop_;
	BCs bcond_;
	FlowSolver flow_solver_;
	TransportSolver transport_solver_;


	virtual void initControl(const parameter::ParameterGroup& param)
	{
	    simulation_steps_ = param.getDefault("simulation_steps", simulation_steps_);
	    stepsize_ = Dune::unit::convert::from(param.getDefault("stepsize", stepsize_),
                                                  Dune::unit::day);
	}

	virtual void initGridAndProps(const parameter::ParameterGroup& param)
	{
	    setupGridAndProps(param, grid_, res_prop_);
	    ginterf_.init(grid_);
	}

	virtual void initInitialConditions(const parameter::ParameterGroup& param)
	{
	    init_saturation_ = param.getDefault("init_saturation", init_saturation_);
	}

	virtual void initBoundaryConditions(const parameter::ParameterGroup& param)
	{
	    setupBoundaryConditions(param, ginterf_, bcond_);
	}

	virtual void initSolvers(const parameter::ParameterGroup& param)
	{
	    // Initialize flow solver.
	    flow_solver_.init(ginterf_, res_prop_, bcond_);
	    //flow_solver_.assembleStatic(ginterf_, res_prop_);
	    // Initialize transport solver.
	    transport_solver_.init(param, ginterf_, res_prop_, bcond_);
	}


    };



} // namespace Dune



#endif // OPENRS_SIMULATORBASE_HEADER
