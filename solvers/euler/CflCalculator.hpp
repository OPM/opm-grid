//===========================================================================
//
// File: CflCalculator.hpp
//
// Created: Wed Jul  1 09:35:59 2009
//
// Author(s): Halvor M Nilsen     <hnil@sintef.no>
//            Atgeirr F Rasmussen <atgeirr@sintef.no>
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

#ifndef OPENRS_CFLCALCULATOR_HEADER
#define OPENRS_CFLCALCULATOR_HEADER

#include <dune/common/ErrorMacros.hpp>
#include <dune/common/Average.hpp>

#include <dune/solvers/common/Matrix.hpp>

namespace Dune {
    namespace cfl_calculator {


	/// @brief
	/// @todo Doc me!
	/// @tparam
	/// @param
	template <class Grid, class ReservoirProperties, class PressureSolution>
	double findCFLtimeVelocity(const Grid& grid,
				   const ReservoirProperties& resprop,
				   const PressureSolution& pressure_sol)
	{
	    double dt = 1e100;
	    typename Grid::CellIterator c = grid.cellbegin();
	    for (; c != grid.cellend(); ++c) {
		double flux_p = 0.0;
		double flux_n = 0.0;
		typename Grid::CellIterator::FaceIterator f = c->facebegin();
		for (; f != c->faceend(); ++f) {
		    const double loc_flux = pressure_sol.outflux(f);
		    if (loc_flux > 0) {
			flux_p += loc_flux;
		    } else {
			flux_n -= loc_flux;
		    }
		}
		double flux = std::max(flux_n, flux_p);
		double loc_dt = (resprop.cflFactor()*c->volume()*resprop.porosity(c->index()))/flux;
		if (loc_dt == 0.0) {
		    THROW("Cfl computation gave dt = 0.0");
		}
		if (loc_dt < dt){
		    dt = loc_dt;
		}
	    }
	    return dt;
	}


	/// @brief
	/// @todo Doc me!
	/// @tparam
	/// @param
	template <class Grid, class ReservoirProperties>
	double findCFLtimeGravity(const Grid& grid,
				  const ReservoirProperties& resprop,
				  const typename Grid::Vector& gravity)
	{
	    typedef typename ReservoirProperties::PermTensor PermTensor;
	    typedef typename ReservoirProperties::MutablePermTensor MutablePermTensor;
	    const int dimension = Grid::Vector::dimension;
	    double dt = 1e100;
	    typename Grid::CellIterator c = grid.cellbegin();
	    for (; c != grid.cellend(); ++c) {
		double flux = 0.0;
		typename Grid::CellIterator::FaceIterator f = c->facebegin();
		for (; f != c->faceend(); ++f) {
		    // UGLY WARNING
		    MutablePermTensor loc_perm_aver;
		    const double* permdata = 0;
		    if (!f->boundary()) {
			PermTensor K0 = resprop.permeability(f->cellIndex());
			PermTensor K1 = resprop.permeability(f->neighbourCellIndex());
			loc_perm_aver = utils::arithmeticAverage<PermTensor, MutablePermTensor>(K0, K1);
			permdata = loc_perm_aver.data();
		    } else {
			permdata = resprop.permeability(f->cellIndex()).data();
		    }
		    PermTensor loc_perm(dimension, dimension, permdata);
		    typename Grid::Vector loc_halfface_normal = f->normal();
		    double loc_gravity_flux = 0.0;
		    for (int k = 0; k < dimension; ++k) {
			for (int q = 0; q < dimension; ++q) {
			    loc_gravity_flux += loc_halfface_normal[q]*(loc_perm(q,k)*gravity[k]*resprop.densityDifference());
			}
		    }
		    loc_gravity_flux *= f->area();
		    if (loc_gravity_flux > 0) {
			flux += loc_gravity_flux;
		    }
		}
		double loc_dt = (resprop.cflFactorGravity()*c->volume()*resprop.porosity(c->index()))/flux;
		if (loc_dt < dt){
		    dt = loc_dt;
		}
	    }
	    return dt;
	}

    } // namespace cfl_calculator
} // namespace Dune


#endif // OPENRS_CFLCALCULATOR_HEADER
