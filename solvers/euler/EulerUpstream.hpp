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

#ifndef OPENRS_EULERUPSTREAM_HEADER
#define OPENRS_EULERUPSTREAM_HEADER

#include <dune/solvers/euler/EulerUpstreamResidual.hpp>

#include <tr1/unordered_map>

#include <dune/common/param/ParameterGroup.hpp>
#include <dune/common/SparseVector.hpp>


namespace Dune {

    namespace {
        // Forward declaration for friendship purposes.
        template <class UpstreamSolver, class PressureSolution>
        struct UpdateForCell;
    } // anon namespace


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
	void initObj(const GridInterface& grid,
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
        template <class S, class P>
        friend class UpdateForCell;
	typedef typename GridInterface::CellIterator CIt;
	typedef typename CIt::FaceIterator FIt;
	typedef typename FIt::Vector Vector;
        typedef ReservoirProperties RP;


	template <class PressureSolution>
	double computeCflTime(const std::vector<double>& saturation,
			      const double time,
			      const typename GridInterface::Vector& gravity,
			      const PressureSolution& pressure_sol) const;

	template <class PressureSolution>
	void smallTimeStep(std::vector<double>& saturation,
			   const double time,
			   const typename GridInterface::Vector& gravity,
			   const PressureSolution& pressure_sol,
                           const SparseVector<double>& injection_rates) const;

	void checkAndPossiblyClampSat(std::vector<double>& s) const;


        EulerUpstreamResidual<GridInterface,
                              ReservoirProperties,
                              BoundaryConditions> residual_;

	bool method_viscous_;
	bool method_gravity_;
	bool method_capillary_;
	bool use_cfl_viscous_;
	bool use_cfl_gravity_;
	bool use_cfl_capillary_;
	// The courant_number is the multiplied with the cfl time to get the time step.
	double courant_number_;
	int minimum_small_steps_;
	int maximum_small_steps_;
	bool check_sat_;
	bool clamp_sat_;

	// Storing sat_change_ so that we won't have to reallocate it for every step.
	mutable std::vector<double> sat_change_;
    };

} // namespace Dune

#include "EulerUpstream_impl.hpp"

#endif // OPENRS_EULERUPSTREAM_HEADER
