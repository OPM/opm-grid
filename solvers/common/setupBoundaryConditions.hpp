//===========================================================================
//
// File: setupBoundaryConditions.hpp
//
// Created: Fri Aug 21 09:07:09 2009
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

#ifndef OPENRS_SETUPBOUNDARYCONDITIONS_HEADER
#define OPENRS_SETUPBOUNDARYCONDITIONS_HEADER

#include <dune/common/param/ParameterGroup.hpp>
#include <dune/solvers/common/BoundaryConditions.hpp>
#include <dune/solvers/common/PeriodicHelpers.hpp>

namespace Dune
{

    /// @brief Setup boundary conditions for a simulation.
    /// It is assumed that the boundary ids are 1-6, similar to cartesian case/Yaspgrid,
    /// unless periodic, in which case we assume unique boundary ids.
    template <class GridType, class BCs>
    inline void setupBoundaryConditions(const parameter::ParameterGroup& param,
					const GridType& g,
					BCs& bcs)
    {
	if (param.getDefault("upscaling", false)) {
	    setupUpscalingConditions(param, g, bcs);
	    return;
	}
	// Make flow equation boundary conditions.
	// Default is pressure 1.0e5 on the left, 0.0 on the right.
	// Recall that the boundary ids range from 1 to 6 for the cartesian edges,
	// and that boundary id 0 means interiour face/intersection.
	std::string flow_bc_type = param.getDefault<std::string>("flow_bc_type", "dirichlet");
	FlowBC::BCType bct = FlowBC::Dirichlet;
	double leftval = 1.0*Dune::unit::barsa;
	double rightval = 0.0;
	if (flow_bc_type == "neumann") {
	    bct = FlowBC::Neumann;
	    leftval = param.get<double>("left_flux");
	    rightval = param.getDefault<double>("right_flux", -leftval);
	} else if (flow_bc_type == "dirichlet") {
	    leftval = param.getDefault<double>("left_pressure", leftval);
	    rightval = param.getDefault<double>("right_pressure", rightval);
	} else if (flow_bc_type == "periodic") {
	    THROW("Periodic conditions not here yet.");
	} else {
	    THROW("Unknown flow boundary condition type " << flow_bc_type);
	}
	bcs.resize(7);
	bcs.flowCond(1) = FlowBC(bct, leftval);
	bcs.flowCond(2) = FlowBC(bct, rightval);

	// Default transport boundary conditions are used.
    }

    /// @brief
    /// @todo Doc me!
    /// @param
    template <class GridType, class BCs>
    inline void setupUpscalingConditions(const parameter::ParameterGroup& param,
					 const GridType& g,
					 BCs& bcs)
    {
	// Caution: This enum is copied from Upscaler.hpp.
	enum BoundaryConditionType { Fixed = 0, Linear = 1, Periodic = 2, PeriodicSingleDirection = 3, Noflow = 4 };
        int bct = param.get<int>("boundary_condition_type");
        if (bct < 0 || bct > 2) {
            THROW("Illegal boundary condition type (0-2 are legal): " << bct); // Later on, we may allow 3 and 4.
        }
	BoundaryConditionType bctype = static_cast<BoundaryConditionType>(bct);
        int pddir = param.getDefault("pressure_drop_direction", 0);
        ASSERT(pddir >=0 && pddir <= 2);
        double boundary_pressuredrop = param.getDefault("boundary_pressuredrop", 1.0e5);

	// Flow conditions.
	switch (bctype) {
	case Fixed:
	    {
		// ASSERT(!g.uniqueBoundaryIds());
		bcs.clear();
		bcs.resize(7);
		bcs.flowCond(2*pddir + 1) = FlowBC(FlowBC::Dirichlet, boundary_pressuredrop);
		bcs.flowCond(2*pddir + 2) = FlowBC(FlowBC::Dirichlet, 0.0);
		double boundary_saturation = param.getDefault("boundary_saturation", 1.0);
		bcs.satCond(2*pddir + 1) = SatBC(SatBC::Dirichlet, boundary_saturation); // The only possible inflow location.
		break;
	    }
	case Linear:
	    {
		THROW("Linear boundary conditions not done yet");
		break;
	    }
	case Periodic:
	    {
		// ASSERT(g.uniqueBoundaryIds());
		FlowBC fb(FlowBC::Periodic, 0.0);
		boost::array<FlowBC, 6> fcond = {{ fb, fb, fb, fb, fb, fb }};
		fcond[2*pddir] = FlowBC(FlowBC::Periodic, boundary_pressuredrop);
		fcond[2*pddir + 1] = FlowBC(FlowBC::Periodic, -boundary_pressuredrop);
		SatBC sb(SatBC::Periodic, 0.0);
		boost::array<SatBC, 6> scond = {{ sb, sb, sb, sb, sb, sb }};
		if (param.getDefault("2d_hack", false)) {
		    fcond[2] = FlowBC(FlowBC::Neumann, 0.0);
		    fcond[3] = FlowBC(FlowBC::Neumann, 0.0);
		    fcond[4] = FlowBC(FlowBC::Neumann, 0.0);
		    fcond[5] = FlowBC(FlowBC::Neumann, 0.0);
		    scond[2] = SatBC(SatBC::Dirichlet, 1.0);
		    scond[3] = SatBC(SatBC::Dirichlet, 1.0);
		    scond[4] = SatBC(SatBC::Dirichlet, 1.0);
		    scond[5] = SatBC(SatBC::Dirichlet, 1.0);
		}
		createPeriodic(bcs, g, fcond, scond);
		break;
	    }
	default:
	    THROW("Error in switch statement, should never be here.");
	}

	// Default transport boundary conditions are used.
    }

} // namespace Dune


#endif // OPENRS_SETUPBOUNDARYCONDITIONS_HEADER
