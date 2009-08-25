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


#include <dune/solvers/common/MatrixInverse.hpp>


namespace Dune
{

    inline Upscaler::Upscaler()
	: bctype_(Fixed),
	  twodim_hack_(false),
	  residual_tolerance_(1e-8)
    {
    }

    inline void Upscaler::initControl(const parameter::ParameterGroup& param)
    {
	SimulatorBase::initControl(param);
	residual_tolerance_ = param.getDefault("residual_tolerance", residual_tolerance_);
    }

    inline void Upscaler::initInitialConditions(const parameter::ParameterGroup& param)
    {
	init_saturation_ = 1e100; // Not used. Yet.
    }

    inline void Upscaler::initBoundaryConditions(const parameter::ParameterGroup& param)
    {
        int bct = param.get<int>("boundary_condition_type");
	bctype_ = static_cast<BoundaryConditionType>(bct);
	twodim_hack_ = param.getDefault("2d_hack", false);
    }

    inline void Upscaler::initSolvers(const parameter::ParameterGroup& param)
    {
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

	permtensor_t upscaled_K(3, 3, (double*)0);
	for (int pdd = 0; pdd < Dimension; ++pdd) {
	    setupUpscalingConditions(ginterf_, bctype_, pdd, 1.0, 1.0, twodim_hack_, bcond_);
	    if (pdd == 0) {
		// Only on first iteration, since we do not change the
		// structure of the system, the way the flow solver is
		// implemented.
		flow_solver_.init(ginterf_, res_prop_, bcond_);
	    }

	    // Run pressure solver.
	    flow_solver_.solve(res_prop_, sat, bcond_, src, gravity, residual_tolerance_);

	    // Check and fix fluxes.
// 	    flux_checker_.checkDivergence(grid_, wells, flux);
// 	    flux_checker_.fixFlux(grid_, wells, boundary_, flux);

	    // Compute upscaled K.
	    double Q[Dimension];
	    switch (bctype_) {
	    case Fixed:
		std::fill(Q, Q+Dimension, 0); // resetting Q
		Q[pdd] = computeAverageVelocity(flow_solver_.getSolution(), pdd, pdd);
		break;
	    case Linear:
	    case Periodic:
		for (int i = 0; i < Dimension; ++i) {
		    Q[i] = computeAverageVelocity(flow_solver_.getSolution(), i, pdd);
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


    inline Upscaler::permtensor_t
    Upscaler::
    upscaleSteadyState(const boost::array<std::vector<double>, Dimension>& initial_saturations,
		       const double boundary_saturation,
		       const double pressure_drop,
		       const permtensor_t& upscaled_perm)
    {
	int num_cells = ginterf_.numberOfCells();
	// No source or sink.
	std::vector<double> src(num_cells, 0.0);
	SparseVector<double> injection(num_cells);
	// Just water.
	std::vector<double> sat(num_cells, 1.0);
	// Gravity.
	FieldVector<double, 3> gravity(0.0);
	// gravity[2] = -Dune::unit::gravity;
	if (gravity.two_norm() > 0.0) {
	    MESSAGE("Warning: Gravity not yet handled by flow solver.");
	}

	permtensor_t upscaled_K(3, 3, (double*)0);

	// Loop over the three pressure drop directions (x, y, z).
	permtensor_t relative_K; // v_w = -relative_K grad p, so relative_K = (k_rw/mu_w)K
	for (int pdd = 0; pdd < Dimension; ++pdd) {
	    // Set up initial saturation profile.
	    // std::vector<double> saturation = setupInitialSaturation(target_saturation);
	    std::vector<double> saturation = initial_saturations[pdd];

	    // Set up boundary conditions.
	    setupUpscalingConditions(ginterf_, bctype_, pdd, 1.0, 1.0, twodim_hack_, bcond_);

	    // Set up solvers.
	    if (pdd == 0) {
		flow_solver_.init(ginterf_, res_prop_, bcond_);
	    }
	    transport_solver_.initObj(ginterf_, res_prop_, bcond_);

	    // Run pressure solver.
	    flow_solver_.solve(res_prop_, sat, bcond_, src, gravity, residual_tolerance_);

	    // Do a run till steady state. For now, we just do some pressure and transport steps...
	    for (int iter = 0; iter < simulation_steps_; ++iter) {
		// Check and fix fluxes.
// 		flux_checker_.checkDivergence(grid_, wells, flux);
// 		flux_checker_.fixFlux(grid_, wells, boundary_, flux);

		// Run transport solver.
		transport_solver_.transportSolve(saturation, stepsize_, gravity, flow_solver_.getSolution(), injection);

		// Run pressure solver.
		flow_solver_.solve(res_prop_, sat, bcond_, src, gravity, residual_tolerance_);
	    }

	    // A check on the final fluxes.
// 	    flux_checker_.checkDivergence(grid_, wells, flux);
// 	    flux_checker_.fixFlux(grid_, wells, boundary_, flux);

	    // Compute upscaled relperm.
	    double Q[Dimension];
	    switch (bctype_) {
	    case Fixed:
		std::fill(Q, Q+Dimension, 0); // resetting Q
		Q[pdd] = computeAveragePhaseVelocity(flow_solver_.getSolution(), saturation, pdd, pdd);
		break;
	    case Linear:
	    case Periodic:
		for (int i = 0; i < Dimension; ++i) {
		    Q[i] = computeAveragePhaseVelocity(flow_solver_.getSolution(), saturation, i, pdd);
		}
		break;
	    default:
		THROW("Unknown boundary type: " << bctype_);
	    }
	    double delta = computeDelta(pdd);
	    for (int i = 0; i < Dimension; ++i) {
		relative_K(i, pdd) = Q[i] * delta / pressure_drop;
	    }

	    // Set the steady state saturation fields for eventual outside access.
	    last_saturations_[pdd].swap(saturation);
	}

	// Compute the relative K tensor.
	relative_K *= res_prop_.viscosityFirstPhase();
	permtensor_t relperm_matrix(matprod(relative_K, inverse3x3(upscaled_perm)));
	// std::cout << relperm_matrix << std::endl;
	return relperm_matrix;
    }

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
					    const int flow_dir,
					    const int pdrop_dir) const
    {
	double side1_flux = 0.0;
	double side2_flux = 0.0;
	double side1_area = 0.0;
	double side2_area = 0.0;

	for (CellIter c = ginterf_.cellbegin(); c != ginterf_.cellend(); ++c) {
	    for (FaceIter f = c->facebegin(); f != c->faceend(); ++f) {
		if (f->boundary()) {
		    int canon_bid = bcond_.getCanonicalBoundaryId(f->boundaryId());
		    if ((canon_bid - 1)/2 == flow_dir) {
			double flux = flow_solution.outflux(f);
			double area = f->area();
			double norm_comp = f->normal()[flow_dir];
			if (canon_bid - 1 == 2*flow_dir) {
			    if (flow_dir == pdrop_dir && flux > 0.0) {
				std::cerr << "Flow may be in wrong direction at bid: " << f->boundaryId()
					  << " Magnitude: " << std::fabs(flux) << std::endl;
				// THROW("Detected outflow at entry face: " << face);
			    }
			    side1_flux += flux*norm_comp;
			    side1_area += area;
			} else {
			    if (flow_dir == pdrop_dir && flux < 0.0) {
				std::cerr << "Flow may be in wrong direction at bid: " << f->boundaryId()
					  << " Magnitude: " << std::fabs(flux) << std::endl;
				// THROW("Detected inflow at exit face: " << face);
			    }
			    side2_flux += flux*norm_comp;
			    side2_area += area;
			}
		    }		    
		}
	    }
	}
	// q is the average velocity.
	return 0.5*(side1_flux/side1_area + side2_flux/side2_area);
    }




    template <class FlowSol>
    inline double Upscaler::computeAveragePhaseVelocity(const FlowSol& flow_solution,
							const std::vector<double>& saturations,
							const int flow_dir,
							const int pdrop_dir) const
    {
	// Apart from the two lines defining frac_flow and flux below, this code
	// is identical to computeAverageVelocity().
	// \todo Unify. Also, there is something fishy about using the cell's fractional flow.
	// Should use the periodic partner's, perhaps?
	// Or maybe just do this for outflow?
	double side1_flux = 0.0;
	double side2_flux = 0.0;
	double side1_area = 0.0;
	double side2_area = 0.0;

	for (CellIter c = ginterf_.cellbegin(); c != ginterf_.cellend(); ++c) {
	    for (FaceIter f = c->facebegin(); f != c->faceend(); ++f) {
		if (f->boundary()) {
		    int canon_bid = bcond_.getCanonicalBoundaryId(f->boundaryId());
		    if ((canon_bid - 1)/2 == flow_dir) {
			double frac_flow = res_prop_.fractionalFlow(c->index(), saturations[c->index()]);
			double flux = flow_solution.outflux(f)*frac_flow;
			double area = f->area();
			double norm_comp = f->normal()[flow_dir];
			if (canon_bid - 1 == 2*flow_dir) {
			    if (flow_dir == pdrop_dir && flux > 0.0) {
				std::cerr << "Flow may be in wrong direction at bid: " << f->boundaryId()
					  << " Magnitude: " << std::fabs(flux) << std::endl;
				// THROW("Detected outflow at entry face: " << face);
			    }
			    side1_flux += flux*norm_comp;
			    side1_area += area;
			} else {
			    if (flow_dir == pdrop_dir && flux < 0.0) {
				std::cerr << "Flow may be in wrong direction at bid: " << f->boundaryId()
					  << " Magnitude: " << std::fabs(flux) << std::endl;
				// THROW("Detected inflow at exit face: " << face);
			    }
			    side2_flux += flux*norm_comp;
			    side2_area += area;
			}
		    }		    
		}
	    }
	}
	// q is the average velocity.
	return 0.5*(side1_flux/side1_area + side2_flux/side2_area);
    }




    inline double Upscaler::computeDelta(const int flow_dir) const
    {
	double side1_pos = 0.0;
	double side2_pos = 0.0;
	double side1_area = 0.0;
	double side2_area = 0.0;
	for (CellIter c = ginterf_.cellbegin(); c != ginterf_.cellend(); ++c) {
	    for (FaceIter f = c->facebegin(); f != c->faceend(); ++f) {
		if (f->boundary()) {
		    int canon_bid = bcond_.getCanonicalBoundaryId(f->boundaryId());
		    if ((canon_bid - 1)/2 == flow_dir) {
			double area = f->area();
			double pos_comp = f->centroid()[flow_dir];
			if (canon_bid - 1 == 2*flow_dir) {
			    side1_pos += area*pos_comp;
			    side1_area += area;
			} else {
			    side2_pos += area*pos_comp;
			    side2_area += area;
			}
		    }		    
		}
	    }
	}
	// delta is the average length.
	return  side2_pos/side2_area - side1_pos/side1_area;
    }

} // namespace Dune

#endif // OPENRS_UPSCALER_IMPL_HEADER
