//===========================================================================
//
// File: SinglePhaseUpscaler_impl.hpp
//
// Created: Fri Aug 28 13:46:07 2009
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

#ifndef OPENRS_SINGLEPHASEUPSCALER_IMPL_HEADER
#define OPENRS_SINGLEPHASEUPSCALER_IMPL_HEADER


#include <dune/solvers/common/setupGridAndProps.hpp>
#include <dune/solvers/common/setupBoundaryConditions.hpp>
#include <dune/solvers/common/ReservoirPropertyTracerFluid.hpp>

namespace Dune
{

    inline SinglePhaseUpscaler::SinglePhaseUpscaler()
	: bctype_(Fixed),
	  twodim_hack_(false),
	  residual_tolerance_(1e-8),
	  linsolver_verbosity_(0),
          linsolver_type_(1)
    {
    }




    inline void SinglePhaseUpscaler::init(const parameter::ParameterGroup& param)
    {
	initImpl(param);
	initFinal(param);
    }



    inline void SinglePhaseUpscaler::initImpl(const parameter::ParameterGroup& param)
    {
        // Request the boundary condition type parameter early since,
        // depending on the actual type, we may have to manufacture
        // and insert other parameters into the ParameterGroup prior
        // to constructing the grid and associated properties.
        //
        int bct = param.get<int>("boundary_condition_type");
	bctype_ = static_cast<BoundaryConditionType>(bct);

        // Import control parameters pertaining to reduced physical
        // dimensionality ("2d_hack = true" precludes periodic
        // boundary conditions in the Z direction), and linear solves.
        //
        twodim_hack_ = param.getDefault("2d_hack", twodim_hack_);
	residual_tolerance_ = param.getDefault("residual_tolerance", residual_tolerance_);
	linsolver_verbosity_ = param.getDefault("linsolver_verbosity", linsolver_verbosity_);
        linsolver_type_ = param.getDefault("linsolver_type", linsolver_type_);

        // Ensure sufficient grid support for requested boundary
        // condition type.
        //
	parameter::ParameterGroup temp_param = param;
        if (bctype_ == Linear || bctype_ == Periodic) {
            if (!temp_param.has("use_unique_boundary_ids")) {
                temp_param.insertParameter("use_unique_boundary_ids", "true");
            }
        }
        if (bctype_ == Periodic) {
            if (!temp_param.has("periodic_extension")) {
                temp_param.insertParameter("periodic_extension", "true");
            }
        }

	setupGridAndProps(temp_param, grid_, res_prop_);
	ginterf_.init(grid_);
    }



    inline void SinglePhaseUpscaler::initFinal(const parameter::ParameterGroup& param)
    {
	// Report any unused parameters.
	std::cout << "====================   Unused parameters:   ====================\n";
	param.displayUsage();
	std::cout << "================================================================\n";
    }




    inline void SinglePhaseUpscaler::init(const EclipseGridParser& parser,
                                          BoundaryConditionType bctype,
                                          double perm_threshold,
                                          double z_tolerance,
                                          double residual_tolerance,
                                          int linsolver_verbosity,
                                          int linsolver_type,
                                          bool twodim_hack)
    {
	bctype_ = bctype;
	residual_tolerance_ = residual_tolerance;
	linsolver_verbosity_ = linsolver_verbosity;
        linsolver_type_ = linsolver_type;
	twodim_hack_ = twodim_hack;

	// Faking some parameters depending on bc type.
        bool periodic_ext = (bctype_ == Periodic);
        bool turn_normals = false;
        bool clip_z = (bctype_ == Periodic);
        bool unique_bids = (bctype_ == Linear || bctype_ == Periodic);
        std::string rock_list("no_list");
	setupGridAndPropsEclipse(parser, z_tolerance,
                                 periodic_ext, turn_normals, clip_z, unique_bids,
                                 perm_threshold, rock_list,
                                 useJ<ResProp>(), 1.0, 0.0,
                                 grid_, res_prop_);
	ginterf_.init(grid_);
    }



    inline const SinglePhaseUpscaler::GridType&
    SinglePhaseUpscaler::grid() const
    {
	return grid_;
    }




    inline void
    SinglePhaseUpscaler::setBoundaryConditionType(BoundaryConditionType type)
    {
        if ((type == Periodic && bctype_ != Periodic)
            || (type != Periodic && bctype_ == Periodic)) {
            THROW("Cannot switch to or from Periodic boundary condition, "
                  "periodic must be set in init() params.");
        } else {
            bctype_ = type;
            if (type == Periodic || type == Linear) {
                grid_.setUniqueBoundaryIds(true);
            } else {
                grid_.setUniqueBoundaryIds(false);
            }
        }
    }



    inline void
    SinglePhaseUpscaler::setPermeability(const int cell_index, const permtensor_t& k)
    {
        res_prop_.permeabilityModifiable(cell_index) = k;
    }




    inline SinglePhaseUpscaler::permtensor_t
    SinglePhaseUpscaler::upscaleSinglePhase()
    {
	int num_cells = ginterf_.numberOfCells();
	// No source or sink.
	std::vector<double> src(num_cells, 0.0);
	// Just water.
	std::vector<double> sat(num_cells, 1.0);
	// Gravity.
	FieldVector<double, 3> gravity(0.0);
	// gravity[2] = -Dune::unit::gravity;

	permtensor_t upscaled_K(3, 3, (double*)0);
	for (int pdd = 0; pdd < Dimension; ++pdd) {
	    setupUpscalingConditions(ginterf_, bctype_, pdd, 1.0, 1.0, twodim_hack_, bcond_);
	    if (pdd == 0) {
		// Only on first iteration, since we do not change the
		// structure of the system, the way the flow solver is
		// implemented.
		flow_solver_.init(ginterf_, res_prop_, gravity, bcond_);
	    }

	    // Run pressure solver.
            ReservoirPropertyTracerFluid fluid;
	    flow_solver_.solve(fluid, sat, bcond_, src, residual_tolerance_, linsolver_verbosity_, linsolver_type_);

	    // Check and fix fluxes.
// 	    flux_checker_.checkDivergence(grid_, wells, flux);
// 	    flux_checker_.fixFlux(grid_, wells, boundary_, flux);

	    // Compute upscaled K.
	    double Q[Dimension] =  { 0 };
	    switch (bctype_) {
	    case Fixed:
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
	return upscaled_K;
    }




    template <class FlowSol>
    double SinglePhaseUpscaler::computeAverageVelocity(const FlowSol& flow_solution,
					    const int flow_dir,
					    const int pdrop_dir) const
    {
	double side1_flux = 0.0;
	double side2_flux = 0.0;
	double side1_area = 0.0;
	double side2_area = 0.0;

	int num_faces = 0;
	int num_bdyfaces = 0;
	int num_side1 = 0;
	int num_side2 = 0;

	for (CellIter c = ginterf_.cellbegin(); c != ginterf_.cellend(); ++c) {
	    for (FaceIter f = c->facebegin(); f != c->faceend(); ++f) {
		++num_faces;
		if (f->boundary()) {
		    ++num_bdyfaces;
		    int canon_bid = bcond_.getCanonicalBoundaryId(f->boundaryId());
		    if ((canon_bid - 1)/2 == flow_dir) {
			double flux = flow_solution.outflux(f);
			double area = f->area();
			double norm_comp = f->normal()[flow_dir];
			// std::cout << "bid " << f->boundaryId() << "   area " << area << "   n " << norm_comp << std::endl;
			if (canon_bid - 1 == 2*flow_dir) {
			    ++num_side1;
			    if (flow_dir == pdrop_dir && flux > 0.0) {
				std::cerr << "Flow may be in wrong direction at bid: " << f->boundaryId()
					  << " Magnitude: " << std::fabs(flux) << std::endl;
				// THROW("Detected outflow at entry face: " << face);
			    }
			    side1_flux += flux*norm_comp;
			    side1_area += area;
			} else {
			    ASSERT(canon_bid - 1 == 2*flow_dir + 1);
			    ++num_side2;
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
// 	std::cout << "Faces: " << num_faces << "   Boundary faces: " << num_bdyfaces
// 		  << "   Side 1 faces: " << num_side1 << "   Side 2 faces: " << num_side2 << std::endl;
	// q is the average velocity.
	return 0.5*(side1_flux/side1_area + side2_flux/side2_area);
    }




    inline double SinglePhaseUpscaler::computeDelta(const int flow_dir) const
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




    double SinglePhaseUpscaler::upscalePorosity() const
    {
        double total_vol = 0.0;
        double total_pore_vol = 0.0;
	for (CellIter c = ginterf_.cellbegin(); c != ginterf_.cellend(); ++c) {
            total_vol += c->volume();
            total_pore_vol += c->volume()*res_prop_.porosity(c->index());
        }
        return total_pore_vol/total_vol;
    }




} // namespace Dune



#endif // OPENRS_SINGLEPHASEUPSCALER_IMPL_HEADER
