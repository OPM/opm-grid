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

#ifndef OPENRS_STEADYSTATEUPSCALER_IMPL_HEADER
#define OPENRS_STEADYSTATEUPSCALER_IMPL_HEADER


#include <boost/lexical_cast.hpp>
#include <dune/solvers/common/MatrixInverse.hpp>
#include <dune/solvers/common/SimulatorUtilities.hpp>
#include <dune/solvers/common/ReservoirPropertyFixedMobility.hpp>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <algorithm>

namespace Dune
{

    inline SteadyStateUpscaler::SteadyStateUpscaler()
	: SinglePhaseUpscaler(),
	  output_(false),
	  simulation_steps_(10),
	  stepsize_(0.1),
	  relperm_threshold_(1.0e-4)
    {
    }

    inline void SteadyStateUpscaler::initImpl(const parameter::ParameterGroup& param)
    {
	SinglePhaseUpscaler::initImpl(param);
	output_ = param.getDefault("output", output_);
	simulation_steps_ = param.getDefault("simulation_steps", simulation_steps_);
	stepsize_ = Dune::unit::convert::from(param.getDefault("stepsize", stepsize_),
					      Dune::unit::day);
	relperm_threshold_ = param.getDefault("relperm_threshold", relperm_threshold_);

	transport_solver_.init(param);
        // Set viscosities and densities if given.
        double v1_default = res_prop_.viscosityFirstPhase();
        double v2_default = res_prop_.viscositySecondPhase();
        res_prop_.setViscosities(param.getDefault("viscosity1", v1_default), param.getDefault("viscosity2", v2_default));
        double d1_default = res_prop_.densityFirstPhase();
        double d2_default = res_prop_.densitySecondPhase();
        res_prop_.setDensities(param.getDefault("density1", d1_default), param.getDefault("density2", d2_default));
    }


    inline std::pair<SteadyStateUpscaler::permtensor_t, SteadyStateUpscaler::permtensor_t>
    SteadyStateUpscaler::
    upscaleSteadyState(const int flow_direction,
                       const std::vector<double>& initial_saturation,
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

        // Set up initial saturation profile.
        std::vector<double> saturation = initial_saturation;

        // Set up boundary conditions.
        setupUpscalingConditions(ginterf_, bctype_, flow_direction, pressure_drop, boundary_saturation, twodim_hack_, bcond_);

        // Set up solvers.
        if (flow_direction == 0) {
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
                std::vector<GridInterface::Vector> cell_velocity;
                estimateCellVelocity(cell_velocity, ginterf_, flow_solver_.getSolution());
                // Dune's vtk writer wants multi-component data to be flattened.
                std::vector<double> cell_velocity_flat(&*cell_velocity.front().begin(),
                                                       &*cell_velocity.back().end());
                std::vector<double> cell_pressure;
                getCellPressure(cell_pressure, ginterf_, flow_solver_.getSolution());
                std::vector<double> cap_pressure;
                computeCapPressure(cap_pressure, res_prop_, saturation);
                Dune::VTKWriter<GridType::LeafGridView> vtkwriter(grid_.leafView());
                vtkwriter.addCellData(cell_velocity_flat, "velocity", GridInterface::Vector::dimension);
                vtkwriter.addCellData(saturation, "saturation");
                vtkwriter.addCellData(cell_pressure, "pressure");
                vtkwriter.addCellData(cap_pressure, "capillary pressure");
                vtkwriter.write(std::string("output-steadystate")
                                + '-' + boost::lexical_cast<std::string>(count)
                                + '-' + boost::lexical_cast<std::string>(flow_direction)
                                + '-' + boost::lexical_cast<std::string>(iter),
                                Dune::VTKOptions::ascii);
            }
        }

        // A check on the final fluxes.
        // 	    flux_checker_.checkDivergence(grid_, wells, flux);
        // 	    flux_checker_.fixFlux(grid_, wells, boundary_, flux);

        // Compute phase mobilities.
        std::vector<double> mob1(num_cells, 0.0);
        std::vector<double> mob2(num_cells, 0.0);
        const double mob1_threshold = relperm_threshold_ / res_prop_.viscosityFirstPhase();
        const double mob2_threshold = relperm_threshold_ / res_prop_.viscositySecondPhase();
        for (int c = 0; c < num_cells; ++c) {
            mob1[c] = std::max(res_prop_.mobilityFirstPhase(c, saturation[c]), mob1_threshold);
            mob2[c] = std::max(res_prop_.mobilitySecondPhase(c, saturation[c]), mob2_threshold);
        }

        // Compute upscaled relperm for each phase.
        ReservoirPropertyFixedMobility fluid_first(mob1);
        permtensor_t eff_Kw = upscaleEffectivePerm(fluid_first);
        ReservoirPropertyFixedMobility fluid_second(mob2);
        permtensor_t eff_Ko = upscaleEffectivePerm(fluid_second);

        // Set the steady state saturation fields for eventual outside access.
        last_saturations_[flow_direction].swap(saturation);

	// Compute the (anisotropic) upscaled mobilities.
        // eff_Kw := lambda_w*K
        //  =>  lambda_w = eff_Kw*inv(K); 
	permtensor_t lambda_w(matprod(eff_Kw, inverse3x3(upscaled_perm)));
	permtensor_t lambda_o(matprod(eff_Ko, inverse3x3(upscaled_perm)));

        // Compute (anisotropic) upscaled relative permeabilities.
        // lambda = k_r/mu
        permtensor_t k_rw(lambda_w);
        k_rw *= res_prop_.viscosityFirstPhase();
        permtensor_t k_ro(lambda_o);
        k_ro *= res_prop_.viscositySecondPhase();
	return std::make_pair(k_rw, k_ro);
    }



    inline const boost::array<std::vector<double>, SteadyStateUpscaler::Dimension>&
    SteadyStateUpscaler::lastSaturations() const
    {
	return last_saturations_;
    }


    double SteadyStateUpscaler::lastSaturationUpscaled(int flow_direction) const
    {
        double pore_vol = 0.0;
        double sat_vol = 0.0;
        for (CellIter c = ginterf_.cellbegin(); c != ginterf_.cellend(); ++c) {
            double cell_pore_vol = c->volume()*res_prop_.porosity(c->index());
            pore_vol += cell_pore_vol;
            sat_vol += cell_pore_vol*last_saturations_[flow_direction][c->index()];
        }
        // Dividing by pore volume gives average saturations.
        return sat_vol/pore_vol;
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
