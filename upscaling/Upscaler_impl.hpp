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


namespace Dune
{

    inline Upscaler::Upscaler()
	: bctype_(Fixed),
	  periodic_dir_(0)
    {
    }



    inline void Upscaler::init(const parameter::ParameterGroup& param)
    {
	SimulatorBase::init(param);
	int bctype = param.get<int>("boundary_condition_type");
	if (bctype < 0 || bctype > 4) {
	    THROW("Illegal boundary condition type (0-4 are legal): " << bctype);
	} else {
	    bctype_ = static_cast<BoundaryConditionType>(bctype);
	}
	if (bctype_ == PeriodicSingleDirection) {
	    periodic_dir_ = param.getDefault("periodic_direction", periodic_dir_);
	}
    }


    inline Upscaler::permtensor_t
    Upscaler::upscaleSinglePhase()
    {
	int num_cells = ginterf_.numberOfCells();
	// No source or sink.
	std::vector<double> src(num_cells, 0.0);
	// Just water.
	std::vector<double> sat(num_cells, 1.0);
	// Gravity.
	FieldVector<double, 3> gravity(0.0);
	// gravity[2] = -Dune::unit::gravity;
	if (gravity.two_norm() > 0.0) {
	    MESSAGE("Warning: Gravity not yet handled by flow solver.");
	}

	permtensor_t upscaled_K;
	boost::array<FlowBC, 2*Dimension> fcond;
	for (int pdd = 0; pdd < Dimension; ++pdd) {
	    // Set boundary direction for current pressure drop direction (0 = x, 1 = y, 2 = z).
// 	    boundary_.setPressureDropDirection(static_cast<UpscalingBoundary::PressureDropDirection>(pdd));
// 	    boundary_.setPressureDrop(1.0);
// 	    boundary_.makeBoundaryConditions(grid_, gravity_, single_fluid);

	    // Run pressure solver.
	    flow_solver_.solve(res_prop_, sat, bcond_, src, gravity);

	    // Check and fix fluxes.
// 	    flux_checker_.checkDivergence(grid_, wells, flux);
// 	    flux_checker_.fixFlux(grid_, wells, boundary_, flux);

	    // Compute upscaled K.
	    double Q[Dimension];
	    switch (bctype_) {
	    case Fixed:
		std::fill(Q, Q+Dimension, 0); // resetting Q
		Q[pdd] = computeAverageVelocity(flow_solver_.getSolution(), pdd);
		break;
	    case Linear:
	    case Periodic:
		for (int i = 0; i < Dimension; ++i) {
		    Q[i] = computeAverageVelocity(flow_solver_.getSolution(), i);
		}
		break;
	    default:
		THROW("Unknown boundary type: " << bctype_);
	    }
	    double delta = computeDelta(pdd);
	    for (int i = 0; i < Dimension; ++i) {
		upscaled_K(i, pdd) = Q[i] * delta;
	    }
	}
	upscaled_K *= res_prop_.viscosityFirstPhase();
	return upscaled_K;
    }


#if 0

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
#endif

    inline const Upscaler::GridType&
    Upscaler::grid() const
    {
	return grid_;
    }


    inline const Dune::array<std::vector<double>, Upscaler::Dimension>&
    Upscaler::lastSaturations() const
    {
	return last_saturations_;
    }


    template <class FlowSol>
    double Upscaler::computeAverageVelocity(const FlowSol& flow_solution,
					    const int flow_dir) const
    {
	return 0.0; //boundary_.computeAverageVelocity(grid_, hface_fluxes, flow_dir);

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




    template <class FlowSol>
    inline double Upscaler::computeAveragePhaseVelocity(const FlowSol& flow_solution,
							const std::vector<double>& saturations,
							const int flow_dir) const
    {
	return 0.0;
	//return boundary_.computeAveragePhaseVelocity(grid_, fluid_, hface_fluxes, saturations, flow_dir);
    }




    inline double Upscaler::computeDelta(const int flow_dir) const
    {
	return 0.0;
	//return boundary_.computeDelta(grid_, flow_dir);
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

} // namespace Dune

#endif // OPENRS_UPSCALER_IMPL_HEADER
