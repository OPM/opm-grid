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

#include <tr1/unordered_map>

#include <dune/common/param/ParameterGroup.hpp>
#include <dune/common/SparseVector.hpp>


namespace Dune {

    /// Class for doing simple transport by explicit Euler upstream method for general grid.
    /// @tparam
    template <class GridInterface, class ReservoirProperties, class BoundaryConditions>
    class EulerUpstream
    {
    public:
	/// @brief
	/// @todo Doc me
	EulerUpstream();
	/// @brief
	/// @todo Doc me
	/// @param
 	EulerUpstream(const GridInterface& grid,
 		      const ReservoirProperties& resprop,
		      const BoundaryConditions& boundary);
	/// @brief
	/// @todo Doc me
	/// @param
	void init(const parameter::ParameterGroup& param);
	/// @brief
	/// @todo Doc me
	/// @param
	void init(const parameter::ParameterGroup& param,
		  const GridInterface& grid,
		  const ReservoirProperties& resprop,
		  const BoundaryConditions& boundary);
	/// @brief
	/// @todo Doc me
	/// @param
	void display();

	/// \brief Set the Courant number.
	/// That is dt = dt_cfl*courant_number.
	/// For this explicit method it should be < 1.
	void setCourantNumber(double cn);

	/// \brief Solve transport equation, evolving \param saturation
	/// for \param time seconds.
	/// Cfl type conditions may force many explicit timesteps to
	/// be taken, before the function returns.
	/// @tparam
	/// @param
	template <class PressureSolution>
	void transportSolve(std::vector<double>& saturation,
			    const double time,
			    const typename GridInterface::Vector& gravity,
			    const PressureSolution& pressure_sol,
			    const SparseVector<double>& injection_rates) const;

    protected:
	typedef typename GridInterface::CellIterator CIt;
	typedef typename CIt::FaceIterator FIt;

	void initPeriodics();

	template <class PressureSolution>
	double computeCflTime(const std::vector<double>& saturation,
			      const double time,
			      const typename GridInterface::Vector& gravity,
			      const PressureSolution& pressure_sol) const;

	template <class PressureSolution>
	void smallTimeStep(std::vector<double>& saturation,
			   const double time,
			   const typename GridInterface::Vector& gravity,
			   const PressureSolution& pressure_sol) const;
	
	template <class PressureSolution>
	void computeSatDelta(std::vector<double>& saturation,
			     const double time,
			     const typename GridInterface::Vector& gravity,
			     const PressureSolution& pressure_sol,
			     std::vector<double>& sat_change) const;


	typename GridInterface::Vector
	estimateCapPressureGradient(FIt f, const std::vector<double>& sat) const;

	void checkAndPossiblyClampSat(std::vector<double>& s) const;


	const GridInterface* pgrid_;
	const ReservoirProperties* preservoir_properties_;
	const BoundaryConditions* pboundary_;
	bool method_viscous_;
	bool method_gravity_;
	bool method_capillary_;
	// The courant_number is the multiplied with the cfl time to get the time step.
	double courant_number_;
	int minimum_small_steps_;
	bool check_sat_;
	bool clamp_sat_;

	// We store the periodic boundary conditions for fast access while awaiting
	// a rewrite of the boundary objects that does the same...
	typedef std::map<FIt, FIt> PartnerMapType;
	PartnerMapType periodic_partner_;
    };

} // namespace Dune

#include "EulerUpstream_impl.hpp"

#endif // OPENRS_EULERUPSTREAM_HEADER
