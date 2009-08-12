//===========================================================================
//
// File: Upscaler_impl.hpp
//
// Created: Mon Aug 10 09:46:00 2009
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

#ifndef OPENRS_UPSCALER_IMPL_HEADER
#define OPENRS_UPSCALER_IMPL_HEADER


//#include "GridFactory.hpp"
//#include "WellsPressureAndRate.hpp"
//#include "MatrixInverse.hpp"
//#include "Units.hpp"

//#include <algorithm>

//#include <dune/solvers/euler/setupGridAndProps.hpp>

namespace Dune
{

	/*
	// Gravity is not used yet.
	std::fill(gravity_.begin(), gravity_.end(), 0.0);

	// Set up flux checker.
	flux_checker_.init(param);


	*/

#if 0

    namespace Upscaler_helper
    {
	class SinglePhaseFluid
	{
	public:
	    typedef double saturation_t;
	    double mixtureDensity(int /*cell*/, const double /*S*/) const
	    {
		return 1.0;
	    }
	    double totalMobility(int /*cell*/, const double /*S*/) const
	    {
		return 1.0;
	    }
	    double calculateDensityOfMixture(const double /*S*/) const
	    {
		return 1.0;
	    }
	};
    }

//     template <class grid_t,
// 	      class rock_t,
// 	      class pressure_solver_t,
// 	      class transport_solver_t>
//     inline typename rock_t::permtensor_t
//     Upscaler<grid_t, rock_t, pressure_solver_t, transport_solver_t>::upscaleSinglePhase()
    inline Upscaler::permtensor_t
    Upscaler::upscaleSinglePhase()
    {
	// Set up (no) wells and variables.
	WellsPressureAndRate wells;
	int num_cells = grid_.template numberOf<grid::CellType>();
	std::vector<double> saturation(num_cells, 1.0);
	std::vector<double> flux, cell_pressure, face_pressure;
	typename rock_t::permtensor_t upscaled_K;
	Upscaler_helper::SinglePhaseFluid single_fluid;
	for (int pdd = 0; pdd < Dimension; ++pdd) {
	    // Set boundary direction for current pressure drop direction (0 = x, 1 = y, 2 = z).
	    boundary_.setPressureDropDirection(static_cast<UpscalingBoundary::PressureDropDirection>(pdd));
	    boundary_.setPressureDrop(1.0);
	    boundary_.makeBoundaryConditions(grid_, gravity_, single_fluid);

	    // Run pressure solver.
	    pressure_solver_.pressureSolveAndUpdateWells(saturation, single_fluid, grid_, rock_,
							 wells, boundary_, gravity_,
							 flux, cell_pressure, face_pressure);

	    // Check and fix fluxes.
	    flux_checker_.checkDivergence(grid_, wells, flux);
	    flux_checker_.fixFlux(grid_, wells, boundary_, flux);

	    // Compute upscaled K.
	    double Q[Dimension];
	    switch (boundary_.getConditionType()) {
	    case UpscalingBoundary::Fixed:
		std::fill(Q, Q+Dimension, 0); // resetting Q
		Q[pdd] = computeAverageVelocity(flux, pdd);
		break;
	    case UpscalingBoundary::Linear:
	    case UpscalingBoundary::Periodic:
		for (int i = 0; i < Dimension; ++i) {
		    Q[i] = computeAverageVelocity(flux, i);
		}
		break;
	    default:
		THROW("Unknown boundary type: " << boundary_.getConditionType());
	    }
	    double delta = computeDelta(pdd);
	    for (int i = 0; i < Dimension; ++i) {
		upscaled_K(i, pdd) = Q[i] * delta;
	    }
	}
	// Commented out, because it is wrong for the single-phase upscaling,
	// since we set the total mobility to 1.0 in SinglePhaseFluid.
	// upscaled_K *= fluid_.getViscosityOne();
	return upscaled_K;
    }




//     template <class grid_t,
// 	      class rock_t,
// 	      class pressure_solver_t,
// 	      class transport_solver_t>
//     inline typename rock_t::permtensor_t
//     Upscaler<grid_t, rock_t, pressure_solver_t, transport_solver_t>::
    inline Upscaler::permtensor_t
    Upscaler::
    upscaleSteadyState(const Go::Array<std::vector<double>, Dimension>& initial_saturations,
		       const double boundary_saturation,
		       const double pressure_drop,
		       const permtensor_t& upscaled_perm)
    {
	// Set up (no) wells and some variables.
	WellsPressureAndRate wells;
	std::vector<double> flux, cell_pressure, face_pressure;

	// Loop over the three pressure drop directions (x, y, z).
	permtensor_t relative_K; // v_w = -relative_K grad p, so relative_K = (k_rw/mu_w)K
	for (int pdd = 0; pdd < Dimension; ++pdd) {
	    // Set boundary direction for current pressure drop direction (0 = x, 1 = y, 2 = z).
	    boundary_.setPressureDropDirection(static_cast<UpscalingBoundary::PressureDropDirection>(pdd));
	    boundary_.setPressureDrop(pressure_drop);
	    boundary_.setBoundarySaturation(boundary_saturation);
	    boundary_.makeBoundaryConditions(grid_, gravity_, fluid_);

	    // Set up initial saturation profile.
	    // std::vector<double> saturation = setupInitialSaturation(target_saturation);
	    std::vector<double> saturation = initial_saturations[pdd];

	    // Do an initial pressure solve.
	    pressure_solver_.pressureSolveAndUpdateWells(saturation, fluid_, grid_, rock_,
							 wells, boundary_, gravity_,
							 flux, cell_pressure, face_pressure);

	    // Do a run till steady state. For now, we just do some pressure and transport steps...
	    for (int iter = 0; iter < num_steps_; ++iter) {
		// Check and fix fluxes.
		flux_checker_.checkDivergence(grid_, wells, flux);
		flux_checker_.fixFlux(grid_, wells, boundary_, flux);

		// Run transport solver.
		transport_solver_.transportSolve(saturation,
						 time_step_,
						 fluid_,
						 grid_,
						 rock_,
						 wells,
						 boundary_,
						 gravity_,
						 flux);

		// Run pressure solver.
		pressure_solver_.pressureSolveAndUpdateWells(saturation, fluid_, grid_, rock_,
							     wells, boundary_, gravity_,
							     flux, cell_pressure, face_pressure);
	    }

	    // A check on the final fluxes.
	    flux_checker_.checkDivergence(grid_, wells, flux);
	    flux_checker_.fixFlux(grid_, wells, boundary_, flux);

	    // Compute upscaled relperm.
	    double Q[Dimension];
	    switch (boundary_.getConditionType()) {
	    case UpscalingBoundary::Fixed:
		std::fill(Q, Q+Dimension, 0); // resetting Q
		Q[pdd] = computeAveragePhaseVelocity(flux, saturation, pdd);
		break;
	    case UpscalingBoundary::Linear:
	    case UpscalingBoundary::Periodic:
		for (int i = 0; i < Dimension; ++i) {
		    Q[i] = computeAveragePhaseVelocity(flux, saturation, i);
		}
		break;
	    default:
		THROW("Unknown boundary type: " << boundary_.getConditionType());
	    }
	    double delta = computeDelta(pdd);
	    for (int i = 0; i < Dimension; ++i) {
		relative_K(i, pdd) = Q[i] * delta / pressure_drop;
	    }

	    // Set the steady state saturation fields for eventual outside access.
	    last_saturations_[pdd].swap(saturation);
	}

	// Compute the relative K tensor.
	relative_K *= fluid_.getViscosityOne();
	permtensor_t relperm_matrix(relative_K*inverse(upscaled_perm));
	// std::cout << relperm_matrix << std::endl;
	return relperm_matrix;
    }




//     template <class grid_t,
// 	      class rock_t,
// 	      class pressure_solver_t,
// 	      class transport_solver_t>
//     inline const typename Upscaler<grid_t, rock_t, pressure_solver_t, transport_solver_t>::fluid_t&
//     Upscaler<grid_t, rock_t, pressure_solver_t, transport_solver_t>::fluid() const
//     {
// 	return fluid_;
//     }



//     template <class grid_t,
// 	      class rock_t,
// 	      class pressure_solver_t,
// 	      class transport_solver_t>
//     inline const grid_t&
//     Upscaler<grid_t, rock_t, pressure_solver_t, transport_solver_t>::grid() const
    inline const Upscaler::GridType&
    Upscaler::grid() const
    {
	return grid_;
    }




//     template <class grid_t,
// 	      class rock_t,
// 	      class pressure_solver_t,
// 	      class transport_solver_t>
//     inline const Go::Array<std::vector<double>, grid_t::Dimension>&
//     Upscaler<grid_t, rock_t, pressure_solver_t, transport_solver_t>::lastSaturations() const
    inline const Dune::array<std::vector<double>, Upscaler::Dimension>&
    Upscaler::lastSaturations() const
    {
	return last_saturations_;
    }




//     template <class grid_t,
// 	      class rock_t,
// 	      class pressure_solver_t,
// 	      class transport_solver_t>
//     inline std::vector<double>
//     Upscaler<grid_t, rock_t, pressure_solver_t, transport_solver_t>::setupInitialSaturation(double target_saturation)
//     {
// 	// For now, a very simple variant.
// 	// Later, we should include more variants, and choose by means of a parameter.
// 	int num_cells = grid_.template numberOf<grid::CellType>();
// 	std::vector<double> s(num_cells, target_saturation);
// 	return s;
//     }



//     template <class grid_t,
// 	      class rock_t,
// 	      class pressure_solver_t,
// 	      class transport_solver_t>
//     inline double
//     Upscaler<grid_t,
// 	     rock_t,
// 	     pressure_solver_t,
// 	     transport_solver_t>::computeAverageVelocity(const std::vector<double>& hface_fluxes,
// 							 const int flow_dir) const
    double Upscaler::computeAverageVelocity(const std::vector<double>& hface_fluxes,
					    const int flow_dir) const
    {
	return boundary_.computeAverageVelocity(grid_, hface_fluxes, flow_dir);

	/*
	double side1_flux = 0.0;
	double side2_flux = 0.0;
	double side1_area = 0.0;
	double side2_area = 0.0;

	const std::vector<int>& bdyfaces = boundary_.boundaryFaces();
	const std::vector<int>& canon_pos = boundary_.canonicalPositions();
	UpscalingBoundary::PressureDropDirection pdd = boundary_.getPressureDropDirection();
	const int num_bdyfaces = bdyfaces.size();
	for (int i = 0; i < num_bdyfaces; ++i) {
	    const int face = bdyfaces[i];
	    if (canon_pos[i]/2 == flow_dir) {
		typename grid_t::range_t hfaces = grid_.template neighbours<grid::FaceType, grid::HalfFaceType>(face);
		ASSERT(hfaces.size() == 1);
		const int hface = hfaces[0];
		double flux = hface_fluxes[hface];
		double area = grid_.getFaceArea(face);
		// @@@ Check that using this is correct:
		double norm_comp = grid_.getHalfFaceNormal(hface)[flow_dir];
		if (canon_pos[i] == 2*flow_dir) {
		    if (flow_dir == pdd && flux > 0.0) {
			std::cerr << "Flow may be in wrong direction at face: " << face
				  << " Magnitude: " << std::fabs(flux) << std::endl;
			// THROW("Detected outflow at entry face: " << face);
		    }
		    side1_flux += flux*norm_comp;
		    side1_area += area;
		} else {
		    if (flow_dir == pdd && flux < 0.0) {
			std::cerr << "Flow may be in wrong direction at face: " << face
				  << " Magnitude: " << std::fabs(flux) << std::endl;
			// THROW("Detected inflow at exit face: " << face);
		    }
		    side2_flux += flux*norm_comp;
		    side2_area += area;
		}
	    }
	}
	// q is the average velocity.
	return 0.5*(side1_flux/side1_area + side2_flux/side2_area);
	*/
    }




    template <class grid_t,
	      class rock_t,
	      class pressure_solver_t,
	      class transport_solver_t>
    inline double
    Upscaler<grid_t,
	     rock_t,
	     pressure_solver_t,
	     transport_solver_t>::computeAveragePhaseVelocity(const std::vector<double>& hface_fluxes,
							      const std::vector<double>& saturations,
							      const int flow_dir) const
    {
	return boundary_.computeAveragePhaseVelocity(grid_, fluid_, hface_fluxes, saturations, flow_dir);
    }




    template <class grid_t,
	      class rock_t,
	      class pressure_solver_t,
	      class transport_solver_t>
    inline double
    Upscaler<grid_t,
	     rock_t,
	     pressure_solver_t,
	     transport_solver_t>::computeDelta(const int flow_dir) const
    {
	return boundary_.computeDelta(grid_, flow_dir);
	/*
	double side1_pos = 0.0;
	double side2_pos = 0.0;
	double side1_area = 0.0;
	double side2_area = 0.0;

	const std::vector<int>& bdyfaces = boundary_.boundaryFaces();
	const std::vector<int>& canon_pos = boundary_.canonicalPositions();
	const int num_bdyfaces = bdyfaces.size();
	for (int i = 0; i < num_bdyfaces; ++i) {
	    const int face = bdyfaces[i];
	    if (canon_pos[i]/2 == flow_dir) {
		typename grid_t::range_t hfaces = grid_.template neighbours<grid::FaceType, grid::HalfFaceType>(face);
		ASSERT(hfaces.size() == 1);
		double area = grid_.getFaceArea(face);
		double pos_comp = grid_.getFaceCentroid(face)[flow_dir];
		if (canon_pos[i] == 2*flow_dir) {
		    // At flow entry face.
		    side1_pos += area*pos_comp;
		    side1_area += area;
		} else {
		    // At flow exit face.
		    side2_pos += area*pos_comp;
		    side2_area += area;
		}
	    }
	}
	// delta is the average length.
	return  side2_pos/side2_area - side1_pos/side1_area;
	*/
    }
#endif

} // namespace Dune

#endif // OPENRS_UPSCALER_IMPL_HEADER
