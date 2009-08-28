//===========================================================================
//
// File: SteadyStateUpscaler.hpp
//
// Created: Fri Aug 28 14:01:19 2009
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

#ifndef OPENRS_STEADYSTATEUPSCALER_HEADER
#define OPENRS_STEADYSTATEUPSCALER_HEADER


#include "config.h"
#include "SinglePhaseUpscaler.hpp"
#include <dune/solvers/euler/EulerUpstream.hpp>
#include <boost/array.hpp>

namespace Dune
{
    /**
       @brief A class for doing steady state upscaling.
       @author Atgeirr F. Rasmussen <atgeirr@sintef.no>
    */
    class SteadyStateUpscaler : public SinglePhaseUpscaler
    {
    public:
	// ------- Methods -------

	/// Default constructor.
	SteadyStateUpscaler();

	/// Initializes the upscaler.
	virtual void init(const parameter::ParameterGroup& param);

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

	/// Accessor for the steady state saturation fields. This is empty until
	/// upscaleSteadyState() is called, at which point it will
	/// contain the last computed (steady) saturation states (one for each cardinal direction).
	const boost::array<std::vector<double>, Dimension>& lastSaturations() const;

    protected:
	// ------- Typedefs -------
	typedef EulerUpstream<GridInterface, ResProp, BCs> TransportSolver;
	// ------- Methods -------
	template <class FlowSol>
	double computeAveragePhaseVelocity(const FlowSol& flow_solution,
					   const std::vector<double>& saturations,
					   const int flow_dir,
					   const int pdrop_dir) const;
	// ------- Data members -------
	boost::array<std::vector<double>, Dimension> last_saturations_;
	bool output_;
	int simulation_steps_;
	double stepsize_;
	TransportSolver transport_solver_;
    };

} // namespace Dune

#include "SteadyStateUpscaler_impl.hpp"


#endif // OPENRS_STEADYSTATEUPSCALER_HEADER
