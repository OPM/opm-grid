//===========================================================================
//
// File: SteadyStateUpscaler_impl.hpp
//
// Created: Fri Aug 28 14:07:51 2009
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

#ifndef OPENRS_STEADYSTATEUPSCALER_IMPL_HEADER
#define OPENRS_STEADYSTATEUPSCALER_IMPL_HEADER


#include <boost/lexical_cast.hpp>
#include <dune/solvers/common/MatrixInverse.hpp>
#include <dune/solvers/common/SimulatorUtilities.hpp>
#include <dune/grid/io/file/vtk/vtkwriter.hh>


namespace Dune
{

    inline SteadyStateUpscaler::SteadyStateUpscaler()
	: SinglePhaseUpscaler(),
	  output_(false),
	  simulation_steps_(10),
	  stepsize_(0.1)
    {
    }

    inline void SteadyStateUpscaler::initImpl(const parameter::ParameterGroup& param)
    {
	SinglePhaseUpscaler::initImpl(param);
	output_ = param.getDefault("output", output_);
	simulation_steps_ = param.getDefault("simulation_steps", simulation_steps_);
	stepsize_ = Dune::unit::convert::from(param.getDefault("stepsize", stepsize_),
					      Dune::unit::day);
	transport_solver_.init(param);
    }


    inline SteadyStateUpscaler::permtensor_t
    SteadyStateUpscaler::
    upscaleSteadyState(const boost::array<std::vector<double>, Dimension>& initial_saturations,
		       const double boundary_saturation,
		       const double pressure_drop,
		       const permtensor_t& upscaled_perm)
    {
	static int count = 0;
	++count;
	int num_cells = ginterf_.numberOfCells();
	// No source or sink.
	std::vector<double> src(num_cells, 0.0);
	SparseVector<double> injection(num_cells);
	// Gravity.
	FieldVector<double, 3> gravity(0.0);
	// gravity[2] = -Dune::unit::gravity;
	if (gravity.two_norm() > 0.0) {
	    MESSAGE("Warning: Gravity not yet handled by flow solver.");
	}

	permtensor_t upscaled_K(3, 3, (double*)0);

	// v_w = -relative_K grad p, so relative_K = (k_rw/mu_w)K
	permtensor_t relative_K(3, 3, (double*)0);

	// Loop over the three pressure drop directions (x, y, z).
	for (int pdd = 0; pdd < Dimension; ++pdd) {
	    // Set up initial saturation profile.
	    // std::vector<double> saturation = setupInitialSaturation(target_saturation);
	    std::vector<double> saturation = initial_saturations[pdd];

	    // Set up boundary conditions.
	    setupUpscalingConditions(ginterf_, bctype_, pdd, pressure_drop, boundary_saturation, twodim_hack_, bcond_);

	    // Set up solvers.
	    if (pdd == 0) {
		flow_solver_.init(ginterf_, res_prop_, gravity, bcond_);
	    }
	    transport_solver_.initObj(ginterf_, res_prop_, bcond_);

	    // Run pressure solver.
	    flow_solver_.solve(res_prop_, saturation, bcond_, src, residual_tolerance_, linsolver_verbosity_);

	    // Do a run till steady state. For now, we just do some pressure and transport steps...
	    for (int iter = 0; iter < simulation_steps_; ++iter) {
		// Check and fix fluxes.
// 		flux_checker_.checkDivergence(grid_, wells, flux);
// 		flux_checker_.fixFlux(grid_, wells, boundary_, flux);

		// Run transport solver.
		transport_solver_.transportSolve(saturation, stepsize_, gravity, flow_solver_.getSolution(), injection);

		// Run pressure solver.
		flow_solver_.solve(res_prop_, saturation, bcond_, src, residual_tolerance_, linsolver_verbosity_);

		// Output.
		if (output_) {
		    std::vector<double> cell_velocity;
		    estimateCellVelocity(cell_velocity, ginterf_, flow_solver_.getSolution());
		    std::vector<double> cell_pressure;
		    getCellPressure(cell_pressure, ginterf_, flow_solver_.getSolution());
		    Dune::VTKWriter<GridType::LeafGridView> vtkwriter(grid_.leafView());
		    vtkwriter.addCellData(cell_velocity, "velocity");
		    vtkwriter.addCellData(saturation, "saturation");
		    vtkwriter.addCellData(cell_pressure, "pressure");
		    vtkwriter.write(std::string("output-steadystate")
				    + '-' + boost::lexical_cast<std::string>(count)
				    + '-' + boost::lexical_cast<std::string>(pdd)
				    + '-' + boost::lexical_cast<std::string>(iter),
				    Dune::VTKOptions::ascii);
		}
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



    inline const boost::array<std::vector<double>, SteadyStateUpscaler::Dimension>&
    SteadyStateUpscaler::lastSaturations() const
    {
	return last_saturations_;
    }



    template <class FlowSol>
    inline double SteadyStateUpscaler::computeAveragePhaseVelocity(const FlowSol& flow_solution,
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


} // namespace Dune


#endif // OPENRS_STEADYSTATEUPSCALER_IMPL_HEADER
