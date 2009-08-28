//===========================================================================
//
// File: SimulatorUtilities.hpp
//
// Created: Fri Aug 28 15:00:15 2009
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

#ifndef OPENRS_SIMULATORUTILITIES_HEADER
#define OPENRS_SIMULATORUTILITIES_HEADER

namespace Dune
{


    /// @brief Estimates a scalar cell velocity from outgoing fluxes.
    /// @tparam GridInterface a grid interface.
    /// @tparam FlowSol a flow solution type.
    /// @param[out] cell_velocity the estimated velocities.
    /// @param[in] ginterf an interface to the grid.
    /// @param[in] flow_solution the object containing the fluxes.
    template <class GridInterface, class FlowSol>
    void estimateCellVelocity(std::vector<double>& cell_velocity,
			      const GridInterface& ginterf,
			      const FlowSol& flow_solution)
    {
	// Algorithm used is same as in halfFaceFluxToCellVelocity.hpp
	// in the Sintef legacy c++ code.
	cell_velocity.clear();
	cell_velocity.resize(ginterf.numberOfCells());
	for (typename GridInterface::CellIterator c = ginterf.cellbegin(); c != ginterf.cellend(); ++c) {
	    int numf = 0;
	    typename GridInterface::Vector cell_v(0.0);
	    typename GridInterface::CellIterator::FaceIterator f = c->facebegin();
	    for (; f != c->faceend(); ++f, ++numf) {
		double flux = flow_solution.outflux(f);
		typename GridInterface::Vector v = f->centroid();
		v -= c->centroid();
		v *= flux/c->volume();
		cell_v += v;
	    }
	    cell_velocity[c->index()] = cell_v.two_norm();
	}
    }


    /// @brief
    /// @todo Doc me!
    /// @tparam
    /// @param
    template <class GridInterface, class FlowSol>
    void getCellPressure(std::vector<double>& cell_pressure,
			 const GridInterface& ginterf,
			 const FlowSol& flow_solution)
    {
	cell_pressure.clear();
	cell_pressure.resize(ginterf.numberOfCells());
	for (typename GridInterface::CellIterator c = ginterf.cellbegin(); c != ginterf.cellend(); ++c) {
	    cell_pressure[c->index()] = flow_solution.pressure(c);
	}
    }


} // namespace Dune


#endif // OPENRS_SIMULATORUTILITIES_HEADER
