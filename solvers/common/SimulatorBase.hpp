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
#include <dune/grid/CpGrid.hpp>
#include <dune/grid/yaspgrid.hh>

#include <dune/solvers/common/GridInterfaceEuler.hpp>
#include <dune/solvers/common/ReservoirPropertyCapillary.hpp>
#include <dune/solvers/common/BoundaryConditions.hpp>
#include <dune/solvers/common/setupGridAndProps.hpp>
#include <dune/solvers/common/setupBoundaryConditions.hpp>
#include <dune/solvers/common/SimulatorUtilities.hpp>

#include <dune/solvers/euler/EulerUpstream.hpp>

#include <dune/solvers/mimetic/MimeticIPEvaluator.hpp>
#include <dune/solvers/mimetic/IncompFlowSolverHybrid.hpp>


namespace Dune
{




    /// @brief
    /// @todo Doc me!
    /// @tparam
    template <template <int> class ResPropT = ReservoirPropertyCapillary,
	      template <class, int, bool> class InnerProd = MimeticIPEvaluator>
    class SimulatorBase
    {
    public:


	/// @brief
	/// @todo Doc me!
	SimulatorBase()
	    : simulation_steps_(1),
	      stepsize_(1.0),   // init() expects units of days! Yes, but now the meaning of stepsize_ changes
	                        // from days (here) to seconds (after init()). Solution to that?
	      init_saturation_(0.0),
              residual_tolerance_(1e-8),
              linsolver_verbosity_(1),
              linsolver_type_(1)
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

    protected:
	typedef CpGrid                                         GridType;
 	enum { Dimension = GridType::dimension };
	typedef FieldVector<double, Dimension>                 Vector;
 	typedef ResPropT<Dimension>                            ResProp;
	typedef GridInterfaceEuler<GridType>                   GridInterface;
	typedef GridInterface::CellIterator                    CellIter;
	typedef CellIter::FaceIterator                         FaceIter;
	typedef BasicBoundaryConditions<true, true>                 BCs;
	typedef IncompFlowSolverHybrid<GridInterface,
				       ResProp,
				       BCs,
				       InnerProd>     FlowSolver;
	typedef EulerUpstream<GridInterface,
			      ResProp,
			      BCs>                             TransportSolver;

	int simulation_steps_;
	double stepsize_;
	double init_saturation_;
        Vector gravity_;
	double residual_tolerance_;
	int linsolver_verbosity_;
	int linsolver_type_;

	GridType grid_;
	GridInterface ginterf_;
	ResProp res_prop_;
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

            gravity_[0] = param.getDefault("gx", 0.0);
            gravity_[1] = param.getDefault("gy", 0.0);
            gravity_[2] = param.getDefault("gz", 0.0); //Dune::unit::gravity);
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
	    flow_solver_.init(ginterf_, res_prop_, gravity_, bcond_);
            residual_tolerance_ = param.getDefault("residual_tolerance", residual_tolerance_);
            linsolver_verbosity_ = param.getDefault("linsolver_verbosity", linsolver_verbosity_);
            linsolver_type_ = param.getDefault("linsolver_type", linsolver_type_);
	    //flow_solver_.assembleStatic(ginterf_, res_prop_);
	    // Initialize transport solver.
	    transport_solver_.init(param, ginterf_, res_prop_, bcond_);
	}


    };



} // namespace Dune



#endif // OPENRS_SIMULATORBASE_HEADER
