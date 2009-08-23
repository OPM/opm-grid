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
#include <dune/solvers/common/SimulatorBase.hpp>
#include <boost/array.hpp>

namespace Dune
{
    /**
       @brief A class for doing upscaling.
       @author Atgeirr F. Rasmussen <atgeirr@sintef.no>
       @date Thu Aug 28 14:45:51 2008
    */
    class Upscaler : public SimulatorBase
    {
    public:
	// ------- Methods -------

	/// Default constructor.
	Upscaler();

	/// Initializes the upscaler.
	// void init(const parameter::ParameterGroup& param); // The inherited one is fine, I hope.

	/// A type for the upscaled permeability.
	typedef ResProp::MutablePermTensor permtensor_t;

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
	permtensor_t upscaleSteadyState(const boost::array<std::vector<double>, Dimension>& initial_saturations,
					const double boundary_saturation,
					const double pressure_drop,
					const permtensor_t& upscaled_perm);

	/// Accessor for the grid object.
	const GridType& grid() const;

	/// Accessor for the steady state saturation fields. This is empty until
	/// upscaleSteadyState() is called, at which point it will
	/// contain the last computed (steady) saturation states (one for each cardinal direction).
	const Dune::array<std::vector<double>, Dimension>& lastSaturations() const;

    protected:
	virtual void initInitialConditions(const parameter::ParameterGroup& param);
	virtual void initBoundaryConditions(const parameter::ParameterGroup& param);
	virtual void initSolvers(const parameter::ParameterGroup& param);

    private:
	// ------- Typedefs and enums -------
	enum BoundaryConditionType { Fixed = 0, Linear = 1, Periodic = 2, PeriodicSingleDirection = 3, Noflow = 4 };

	// ------- Methods -------
	// std::vector<double> setupInitialSaturation(double target_saturation);
	template <class FlowSol>
	double computeAverageVelocity(const FlowSol& flow_solution,
				      const int flow_dir,
				      const int pdrop_dir) const;
	template <class FlowSol>
	double computeAveragePhaseVelocity(const FlowSol& flow_solution,
					   const std::vector<double>& saturations,
					   const int flow_dir,
					   const int pdrop_dir) const;
	double computeDelta(const int flow_dir) const;



	// ------- Data members -------
	// FluxChecker flux_checker_;
	// typename grid_t::point_t gravity_;
	Dune::array<std::vector<double>, Dimension> last_saturations_;
	BoundaryConditionType bctype_;
	bool twodim_hack_;
    };

} // namespace Dune

#include "Upscaler_impl.hpp"


#endif // OPENRS_UPSCALER_HEADER
