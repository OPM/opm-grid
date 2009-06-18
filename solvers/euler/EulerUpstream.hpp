//===========================================================================
//
// File: EulerUpstream.hpp
//
// Created: Tue Jun 16 14:07:53 2009
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

#ifndef OPENRS_EULERUPSTREAM_HEADER
#define OPENRS_EULERUPSTREAM_HEADER


#include "../../common/param/ParameterGroup.hpp"
#include <tr1/unordered_map>


namespace Dune {

    /// Class for doing simple transport by explicit Euler upstream method for general grid.
    //template <class GridInterface, class Problem, class ReservoirProperties>
    class EulerUpstream
    {
    public:
	EulerUpstream();
// 	EulerUpstream(const GridInterface& grid,
// 		      const Problem& problem,
// 		      const ReservoirProperties& resprop);
	void init(const parameter::ParameterGroup& param);
	void display();
	/**
	   Set the courant number. that is dt=dt_cfl*courant_number.
	   For this explicite methode it should be < 1.
	*/
	void setCourantNumber(double cn);

	//typedef typename grid_t::point_t point_t;
	template <class grid_t,class flowsys_t, class rock_t, class well_t, class boundary_t>
	void transportSolve(std::vector<double>& saturation,
			    const double time,
			    const flowsys_t& flow_system,
			    const grid_t& my_grid,
			    const rock_t& rock_data,
			    const well_t& well_data,
			    const boundary_t& boundary,
			    const typename grid_t::point_t& gravity,
			    const std::vector<double>& flux_vector) const;
	/**
	   Routione to do one euler upstream time step
	   /param saturation output saturation
	   /param dt time step
	   /param flowsys class for calculating fluid properties
	   /param velocity input velocity
	   /param gravity input gravity vector
	   /param permeability input vector of matrixes of permability
	   /param porosity input porosity
	   /param source vector of rate for each cell
	*/

    protected:
	template <class grid_t,class flowsys_t, class rock_t, class well_t, class boundary_t>
	void smallTimeStep(std::vector<double>& saturation,
			   const double dt,
			   const flowsys_t& flowsys,
			   const grid_t& grid,
			   const rock_t& rock_data,
			   const well_t& well_data,
			   const boundary_t& boundary,
			   const typename grid_t::point_t& gravity,
			   const std::vector<double>& flux_vector) const;
	
	template <class grid_t,class flowsys_t, class rock_t, class well_t, class boundary_t>
	void computeSatDelta(const std::vector<double>& saturation,
			     const double dt,
			     const flowsys_t& flowsys,
			     const grid_t& grid,
			     const rock_t& rock_data,
			     const well_t& well_data,
			     const boundary_t& boundary,
			     const typename grid_t::point_t& gravity,
			     const std::vector<double>& flux_vector,
			     std::vector<double>& sat_change) const;


	template <class grid_t, class rock_t, class fluid_t>
	typename grid_t::point_t
	estimateCapPressureGradient(const grid_t& grid,
				    const rock_t& rock,
				    const fluid_t& fluid,
				    int cell,
				    int face,
				    const std::vector<double>& sat) const;

	void checkAndPossiblyClampSat(std::vector<double>& s) const;

	template <class boundary_t>
	void initPeriodics(const boundary_t& boundary) const;

	template <class grid_t, class fluid_t, class rock_t>
	double computeCflTime(const grid_t& my_grid, const fluid_t& flow_system,
			      const rock_t rock_data, const std::vector<double>& flux_vector,
			      double time, const typename grid_t::point_t& gravity) const;

	bool method_viscous_;
	bool method_gravity_;
	bool method_capillary_;
	double courant_number_;
	int minimum_small_steps_;
	bool check_sat_;
	bool clamp_sat_;
	/**
	 * courant_number is the used to multiply with the cfl time to get the time step
	 */
	// We store the periodic boundary conditions for fast access while awaiting
	// a rewrite of the boundary objects that does the same...
	mutable std::tr1::unordered_map<int, int> periodic_partner_;
    };
} // namespace Dune

#include "EulerUpstream_impl.hpp"


#endif // OPENRS_EULERUPSTREAM_HEADER
