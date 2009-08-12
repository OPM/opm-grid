//===========================================================================
//
// File: Upscaler.hpp
//
// Created: Mon Aug 10 09:29:04 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
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

#ifndef OPENRS_UPSCALER_HEADER
#define OPENRS_UPSCALER_HEADER

#include "config.h"
#include <dune/solvers/euler/SimulatorBase.hpp>
// #include "FluxChecker.hpp"
// #include "UpscalingBoundary.hpp"
// #include "NewFluid.hpp"
// #include "TransportSolverDummy.hpp"

namespace Dune
{
    /**
       @brief A class for doing upscaling.
       @author Atgeirr F. Rasmussen <atgeirr@sintef.no>
       @date Thu Aug 28 14:45:51 2008
    */
//     template <class grid_t,
// 	      class reservoir_properties_t,
// 	      class pressure_solver_t,
// 	      class transport_solver_t>
    class Upscaler : public SimulatorBase
    {
    public:
	// ------- Typedefs and enums -------

	typedef CpGrid                                      GridType;
	enum { Dimension = GridType::dimension };
	typedef GridInterfaceEuler<GridType>                GridInterface;
	typedef GridInterface::CellIterator                 CellIter;
	typedef CellIter::FaceIterator                      FaceIter;
	typedef ReservoirPropertyCapillary<Dimension>       ResProp;
	typedef MimeticIPEvaluator<CellIter,
				   Dimension,
				   true>                    InnerProd;
	typedef FlowBoundaryCondition                       FBC;
	typedef FlowBoundaryConditions                      FBCs;
	typedef IncompFlowSolverHybrid<GridInterface,
				       FBCs,
				       InnerProd>           FlowSolver;
	typedef EulerUpstream<GridInterface,
			      ResProp,
			      SaturationBoundaryConditions> TransportSolver;


	// ------- Methods -------

	/// Default constructor.
	//	Upscaler();

	/// Initializes the upscaler.
	void init(parameter::ParameterGroup& param);

	/// A type for the upscaled permeability.
	typedef ResProp::PermTensor permtensor_t;

	/// Does a single-phase upscaling.
	/// @return an upscaled permeability tensor.
	permtensor_t upscaleSinglePhase();

	/// Does a steady-state upscaling.
	/// @param initial_saturation the initial saturation profiles for the steady-state computations.
	/// There should be one such saturation vector for each cardinal direction, each vector having 
	/// size equal to the number of cells in the grid.
	/// @param boundary_saturation the saturation of fluid flowing in across the boundary,
	/// only needed for nonperiodic upscaling.
	/// @param pressure_drop the pressure drop in Pascal over the domain.
	/// @param upscaled_perm typically the output of upscaleSinglePhase().
	/// @return the upscaled relative permeability matrix of the first phase (usually water).
	/// The relative permeability matrix, call it k, is such that if K_w is the phase
	/// permeability and K the absolute permeability, K_w = k*K.
	permtensor_t upscaleSteadyState(const Dune::array<std::vector<double>, Dimension>& initial_saturations,
					const double boundary_saturation,
					const double pressure_drop,
					const permtensor_t& upscaled_perm);

	/// Accessor for the fine-scale fluid property object. For testing purposes.
	// const fluid_t& fluid() const;

	/// Accessor for the grid object.
	const GridType& grid() const;

	/// Accessor for the steady state saturation fields. This is empty until
	/// upscaleSteadyState() is called, at which point it will
	/// contain the last computed (steady) saturation states (one for each cardinal direction).
	const Dune::array<std::vector<double>, Dimension>& lastSaturations() const;

    private:
	// Methods.
	// std::vector<double> setupInitialSaturation(double target_saturation);
	double computeAverageVelocity(const std::vector<double>& hface_fluxes,
				      const int flow_dir) const;
	double computeAveragePhaseVelocity(const std::vector<double>& hface_fluxes,
					   const std::vector<double>& saturations,
					   const int flow_dir) const;
	double computeDelta(const int flow_dir) const;



	// FluxChecker flux_checker_;
	// typename grid_t::point_t gravity_;
	Dune::array<std::vector<double>, Dimension> last_saturations_;
    };

} // namespace Dune

#include "Upscaler_impl.hpp"


#endif // OPENRS_UPSCALER_HEADER
