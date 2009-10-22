//===========================================================================
//
// File: EulerUpstream_impl.hpp
//
// Created: Tue Jun 16 14:25:24 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
//            Halvor M Nilsen     <hnil@sintef.no>
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

#ifndef OPENRS_EULERUPSTREAM_IMPL_HEADER
#define OPENRS_EULERUPSTREAM_IMPL_HEADER



#include <cassert>
#include <cmath>
#include <algorithm>

#include <dune/common/ErrorMacros.hpp>
#include <dune/common/Average.hpp>
#include <dune/common/Units.hpp>
#include <dune/grid/common/Volumes.hpp>
#include <dune/solvers/euler/CflCalculator.hpp>

namespace Dune
{


    template <class GI, class RP, class BC>
    inline EulerUpstream<GI, RP, BC>::EulerUpstream()
	: pgrid_(0),
	  preservoir_properties_(0),
	  pboundary_(0),
	  method_viscous_(true),
	  method_gravity_(true),
	  method_capillary_(true),
	  courant_number_(0.5),
	  minimum_small_steps_(1),
	  check_sat_(true),
	  clamp_sat_(false)
    {
    }

    template <class GI, class RP, class BC>
    inline EulerUpstream<GI, RP, BC>::EulerUpstream(const GI& g, const RP& r, const BC& b)
	: pgrid_(&g),
	  preservoir_properties_(&r),
	  pboundary_(&b),
	  method_viscous_(true),
	  method_gravity_(true),
	  method_capillary_(true),
	  courant_number_(0.5),
	  minimum_small_steps_(1),
	  check_sat_(true),
	  clamp_sat_(false)
    {
	initFinal();
    }



    template <class GI, class RP, class BC>
    inline void EulerUpstream<GI, RP, BC>::init(const parameter::ParameterGroup& param)
    {
	courant_number_ = param.getDefault("courant_number", courant_number_);
	method_viscous_ = param.getDefault("method_viscous", method_viscous_);
	method_gravity_ = param.getDefault("method_gravity", method_gravity_);
	method_capillary_ = param.getDefault("method_capillary", method_capillary_);
	minimum_small_steps_ = param.getDefault("minimum_small_steps", minimum_small_steps_);
	check_sat_ = param.getDefault("check_sat", check_sat_);
	clamp_sat_ = param.getDefault("clamp_sat", clamp_sat_);
    }

    template <class GI, class RP, class BC>
    inline void EulerUpstream<GI, RP, BC>::init(const parameter::ParameterGroup& param,
						const GI& g, const RP& r, const BC& b)
    {
	init(param);
	initObj(g, r, b);
    }


    template <class GI, class RP, class BC>
    inline void EulerUpstream<GI, RP, BC>::initObj(const GI& g, const RP& r, const BC& b)
    {
	pgrid_ = &g;
	preservoir_properties_ = &r;
	pboundary_ = &b;
	initFinal();
    }



    template <class GI, class RP, class BC>
    inline void EulerUpstream<GI, RP, BC>::display()
    {
	using namespace std;
	cout << endl;
	cout <<"Displaying some members of EulerUpstream" << endl;
	cout << endl;
	cout << "courant_number = " << courant_number_ << endl;
    }



    template <class GI, class RP, class BC>
    inline void EulerUpstream<GI, RP, BC>::setCourantNumber(double cn)
    {
	courant_number_ = cn;
    }



    template <class GI, class RP, class BC>
    template <class PressureSolution>
    void EulerUpstream<GI, RP, BC>::transportSolve(std::vector<double>& saturation,
						   const double time,
						   const typename GI::Vector& gravity,
						   const PressureSolution& pressure_sol,
						   const SparseVector<double>& injection_rates) const
    {
	if (injection_rates.nonzeroSize() != 0) {
	    MESSAGE("Warning: EulerUpstream currently ignores source terms.");
	}
	// Compute the cfl time-step.
	double cfl_dt = computeCflTime(saturation, time, gravity, pressure_sol);

	// Compute the number of small steps to take, and the actual small timestep.
	int nr_transport_steps;
	if (cfl_dt > time){
	    nr_transport_steps = minimum_small_steps_;
	} else {
	    nr_transport_steps = std::max(int(std::ceil(time/cfl_dt)), minimum_small_steps_);
	}
	double dt_transport = time/nr_transport_steps;

	// Do the timestepping. The try-catch blocks are there to handle
	// the situation that smallTimeStep throws, which may happen due
	// to saturation out of bounds (if check_sat_ is true).
	// We cannot guarantee that this does not happen, since we do not
	// (yet) compute a capillary cfl condition.
	// Using exception for "alternate control flow" like this is bad
	// design, should rather use error return values for this.
	std::vector<double> saturation_initial(saturation);
	bool finished = false;
	int repeats = 0;
	const int max_repeats = 10;
	while (!finished) {
	    try {
#ifdef VERBOSE
		std::cout << "Doing " << nr_transport_steps
			  << " steps for saturation equation with stepsize "
			  << dt_transport << " in seconds." << std::endl;
#endif // VERBOSE
		for (int q = 0; q < nr_transport_steps; ++q) {
		    smallTimeStep(saturation,
				  dt_transport,
				  gravity,
				  pressure_sol);
		}
		finished = true;
	    }
	    catch (...) {
		++repeats;
		if (repeats > max_repeats) {
		    throw;
		}
		std::cout << "Warning: Transport failed, retrying with more steps." << std::endl;
		nr_transport_steps *= 2;
		dt_transport = time/nr_transport_steps;
		saturation = saturation_initial;
	    }
	}
    }








    template <class GI, class RP, class BC>
    inline void EulerUpstream<GI, RP, BC>::initFinal()
    {
	// Resize the max timestep per face vectors.
	int num_faces = pgrid_->numberOfFaces();
	visc_maxtimes_.resize(num_faces);
	grav_maxtimes_.resize(num_faces);
	cap_maxtimes_.resize(num_faces);

	// Build bid_to_face_ mapping for handling periodic conditions.
	int maxbid = 0;
	for (typename GI::CellIterator c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c) {
	    for (typename GI::CellIterator::FaceIterator f = c->facebegin(); f != c->faceend(); ++f) {
		int bid = f->boundaryId();
		maxbid = std::max(maxbid, bid);
	    }
	}
	bid_to_face_.clear();
	bid_to_face_.resize(maxbid + 1);
	for (typename GI::CellIterator c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c) {
	    for (typename GI::CellIterator::FaceIterator f = c->facebegin(); f != c->faceend(); ++f) {
		int bid = f->boundaryId();
		if (pboundary_->satCond(bid).isPeriodic()) {
		    bid_to_face_[bid] = f;
		}
	    }
	}
    }





    template <class GI, class RP, class BC>
    template <class PressureSolution>
    inline double EulerUpstream<GI, RP, BC>::computeCflTime(const std::vector<double>& /*saturation*/,
#ifdef VERBOSE
							    const double time,
#else
							    const double,
#endif
							    const typename GI::Vector& gravity,
							    const PressureSolution& pressure_sol) const
    {
	// Deal with cfl stuff, compute the necessary number of small time steps.
	double cfl_dt_v = 1e99;
	double cfl_dt_g = 1e99;
	double cfl_dt_c = 1e99;

	// Viscous cfl.
	if (method_viscous_) {
	    cfl_dt_v = cfl_calculator::findCFLtimeVelocity(*pgrid_,
							   *preservoir_properties_,
							   pressure_sol);
#ifdef VERBOSE
	    std::cout << "CFL dt for velocity is "
                      << Dune::unit::convert::to(cfl_dt_v, Dune::unit::day)
		      << " and total impes time is "
                      << Dune::unit::convert::to(time, Dune::unit::day)
		      << " in days." << std::endl;
#endif // VERBOSE
	}

	// Gravity cfl.
	if (method_gravity_) {
	    cfl_dt_g = cfl_calculator::findCFLtimeGravity(*pgrid_, *preservoir_properties_, gravity);
#ifdef VERBOSE
	    std::cout << "CFL dt for gravity is "
                      << Dune::unit::convert::to(cfl_dt_g, Dune::unit::day)
		      << " and total impes time is "
                      << Dune::unit::convert::to(time, Dune::unit::day)
		      << " in days." << std::endl;
#endif // VERBOSE
	}

	// Capillary cfl. Not done yet.
	if (method_capillary_) {
#ifdef VERBOSE
	    std::cout << "CFL dt for capillarity is not implemented yet." << std::endl;
#endif // VERBOSE
	}

	double cfl_dt = std::min(std::min(cfl_dt_v, cfl_dt_g), cfl_dt_c);
	cfl_dt *= courant_number_;
#ifdef VERBOSE
	std::cout << "Final modified CFL dt is "
                  << Dune::unit::convert::to(cfl_dt, Dune::unit::day)
		  << " and total impes time is "
                  << Dune::unit::convert::to(time, Dune::unit::day)
		  << " in days." << std::endl;
#endif // VERBOSE
	return cfl_dt;
    }







    template <class GI, class RP, class BC>
    inline typename GI::Vector
    EulerUpstream<GI, RP, BC>::estimateCapPressureGradient(const FIt& f, const FIt& nbf, const std::vector<double>& sat) const
    {
	typedef typename GI::CellIterator::FaceIterator Face;
	typedef typename Face::Cell Cell;
	typedef typename GI::Vector Vector;

	// At nonperiodic boundaries, we return a zero gradient.
	// That is (sort of) a trivial Neumann (noflow) condition for the capillary pressure.
	if (f->boundary() && !pboundary_->satCond(f->boundaryId()).isPeriodic()) {
	    return Vector(0.0);
	}
	// Find neighbouring cell and face: nbc and nbf.
	// If we are not on a periodic boundary, nbf is of course equal to f.
	Cell c = f->cell();
	Cell nb = f->boundary() ? (f == nbf ? c : nbf->cell()) : f->neighbourCell();

	// Estimate the gradient like a finite difference between
	// cell centers, except that in order to handle periodic
	// conditions we pass through the face centroid(s).
	Vector cell_c = c.centroid();
	Vector nb_c = nb.centroid();
	Vector f_c = f->centroid();
	Vector nbf_c = nbf->centroid();
	double d0 = (cell_c - f_c).two_norm();
	double d1 = (nb_c - nbf_c).two_norm();
	int cell = c.index();
	int nbcell = nb.index();
	double cp0 = cap_pressures_[cell];
	double cp1 = cap_pressures_[nbcell];
	double val = (cp1 - cp0)/(d0 + d1);
	Vector res = nb_c - nbf_c + f_c - cell_c;
	res /= res.two_norm();
	res *= val;
	return res;
    }


    template <class GI, class RP, class BC>
    inline void EulerUpstream<GI, RP, BC>::computeCapPressures(const std::vector<double>& sat) const
    {
	int num_cells = sat.size();
	cap_pressures_.resize(num_cells);
	for (int cell = 0; cell < num_cells; ++cell) {
	    cap_pressures_[cell] = preservoir_properties_->capillaryPressure(cell, sat[cell]);
	}
    }



    template <class GI, class RP, class BC>
    inline void EulerUpstream<GI, RP, BC>::checkAndPossiblyClampSat(std::vector<double>& s) const
    {
	int num_cells = s.size();
	for (int cell = 0; cell < num_cells; ++cell) {
	    if (s[cell] > 1.0 || s[cell] < 0.0) {
		if (clamp_sat_) {
		    s[cell] = std::max(std::min(s[cell], 1.0), 0.0);
		} else if (s[cell] > 1.001 || s[cell] < -0.001) {
		    THROW("Saturation out of range in EulerUpstream: Cell " << cell << "   sat " << s[cell]);
		}
	    }
	}
    }




	
    template <class GI, class RP, class BC>
    template <class PressureSolution>
    inline void EulerUpstream<GI, RP, BC>::smallTimeStep(std::vector<double>& saturation,
							 const double dt,
							 const typename GI::Vector& gravity,
							 const PressureSolution& pressure_sol) const
    {
	computeCapPressures(saturation);
	computeSatDelta(saturation, gravity, pressure_sol);
	int num_cells = saturation.size();
	for (int i = 0; i < num_cells; ++i) {
	    saturation[i] += dt*sat_change_[i];
	}
	if (check_sat_ || clamp_sat_) {
	    checkAndPossiblyClampSat(saturation);
	}
    }



    /*
      template <class grid_t, class flowsys_t, class rock_t, class well_t>
      inline void wellDelta(const std::vector<double>& saturation,
      const double dt,
      const flowsys_t& flowsys,
      const grid_t& grid,
      const rock_t& rock_data,
      const well_t& well_data,
      std::vector<double>& sat_change)
      {
      // Source terms from wells.
      for (int i = 0; i < well_data.getNumberOfWells(); ++i) {
      typename well_t::well_type well = well_data.getWell(i);
      const std::vector<int>& cell_index = well.cellIndices();
      for (int j = 0; j < well.numberOfCells(); ++j){
      double dS = 0.0;
      double source_rate = well.getRate(j);
      if ( source_rate < 0) {
      dS -= source_rate*(flowsys.mobilityOne(cell_index[j], saturation[cell_index[j]])
      /flowsys.totalMobility(cell_index[j], saturation[cell_index[j]]));
      }
      if (source_rate > 0) {
      dS-= source_rate; 
      }
      sat_change[cell_index[j]] -= (dt/rock_data.getPorosity(cell_index[j]))*dS/grid.getCellVolume(cell_index[j]);
      }
      }	 
      }

    */

    namespace
    {
	template <typename T, template <typename> class StoragePolicy, class OrderingPolicy>
	FullMatrix<T, OwnData, OrderingPolicy>
        arithAver(const FullMatrix<T, StoragePolicy, OrderingPolicy>& m1,
                  const FullMatrix<T, StoragePolicy, OrderingPolicy>& m2)
	{
	    return utils::arithmeticAverage<FullMatrix<T, StoragePolicy, OrderingPolicy>,
		FullMatrix<T, OwnData, OrderingPolicy> >(m1, m2);
	}
    } // anon namespace



    template <class GI, class RP, class BC>
    template <class PressureSolution>
    inline void EulerUpstream<GI, RP, BC>::computeSatDelta(const std::vector<double>& saturation,
							   const typename GI::Vector& gravity,
							   const PressureSolution& pressure_sol) const
    {
	typedef typename GI::Vector Vector;
	typedef typename RP::PermTensor PermTensor;
	typedef typename RP::MutablePermTensor MutablePermTensor;

	// Make sure sat_change is zero, and has the right size.
	sat_change_.clear();
	sat_change_.resize(saturation.size(), 0.0);

	// For every face, we will modify sat_change for adjacent cells.
	// We loop over every cell and intersection, and modify only if
	// this cell has lower index than the neighbour, or we are on the boundary.
	typename GI::CellIterator c = pgrid_->cellbegin();
	for (; c != pgrid_->cellend(); ++c) {
	    for (FIt f = c->facebegin(); f != c->faceend(); ++f) {
		// Neighbour face, will be changed if on a periodic boundary.
		FIt nbface = f;
		double dS = 0.0;
		int cell[2];
		double cell_sat[2];
		double cell_pvol[2];
		cell[0] = f->cellIndex();
		cell_sat[0] = saturation[cell[0]];
		cell_pvol[0] = c->volume()*preservoir_properties_->porosity(cell[0]);
		if (f->boundary()) {
		    int bid = f->boundaryId();
		    if (pboundary_->satCond(bid).isPeriodic()) {
			nbface = bid_to_face_[pboundary_->getPeriodicPartner(bid)];
			ASSERT(nbface != f);
			cell[1] = nbface->cellIndex();
			ASSERT(cell[0] != cell[1]);
			// Periodic faces will be visited twice, but only once
			// should they contribute. We make sure that we skip the
			// periodic faces half the time.
			if (cell[0] > cell[1]) {
			    // We skip this face.
			    continue;
			}
			cell_sat[1] = saturation[cell[1]];
			cell_pvol[1] = nbface->cell().volume()*preservoir_properties_->porosity(cell[1]);
		    } else {
			ASSERT(pboundary_->satCond(bid).isDirichlet());
			cell[1] = cell[0];
			cell_sat[1] = pboundary_->satCond(bid).saturation();
			cell_pvol[1] = cell_pvol[0];
		    }
		} else {
		    cell[1] = f->neighbourCellIndex();
		    ASSERT(cell[0] != cell[1]);
		    if (cell[0] > cell[1]) {
			// We skip this face.
			continue;
		    }
		    cell_sat[1] = saturation[cell[1]];
		    cell_pvol[1] = f->neighbourCellVolume()*preservoir_properties_->porosity(cell[1]);
		}

		// Get some local properties, and compute averages.
		const double loc_area = f->area();
		const double loc_flux = pressure_sol.outflux(f);
		const Vector loc_normal = f->normal();
		const int face_index = f->index();
		// Doing arithmetic averages. Should we consider harmonic or geometric instead?
		using utils::arithmeticAverage;
		const MutablePermTensor aver_perm
		    = arithAver(preservoir_properties_->permeability(cell[0]),
				preservoir_properties_->permeability(cell[1]));
		const double aver_sat
		    = arithmeticAverage<double, double>(cell_sat[0], cell_sat[1]);
		const double aver_lambda_one
		    = arithmeticAverage<double, double>(preservoir_properties_->mobilityFirstPhase(cell[0], aver_sat),
							preservoir_properties_->mobilityFirstPhase(cell[1], aver_sat));
		const double aver_lambda_two
		    = arithmeticAverage<double, double>(preservoir_properties_->mobilitySecondPhase(cell[0], aver_sat), 
							preservoir_properties_->mobilitySecondPhase(cell[1], aver_sat));

		// The local gravity flux is needed for finding the correct phase mobilities.
		double loc_gravity_flux = 0.0;
		if (method_gravity_) {
		    double grav_comp = inner(loc_normal, prod(aver_perm, gravity));
		    loc_gravity_flux = loc_area*preservoir_properties_->densityDifference()*grav_comp;
		}
		// Find the correct phasemobilities to use
		const double flux_one = loc_flux + loc_gravity_flux*aver_lambda_two;
		const double flux_two = loc_flux - loc_gravity_flux*aver_lambda_one;
		double lambda_one;
		double lambda_two;
		// total velocity term
		if (flux_one > 0){
		    lambda_one = preservoir_properties_->mobilityFirstPhase(cell[0], cell_sat[0]);
		} else {
		    lambda_one = preservoir_properties_->mobilityFirstPhase(cell[1], cell_sat[1]); 
		}
		if (flux_two > 0){
		    lambda_two = preservoir_properties_->mobilitySecondPhase(cell[0], cell_sat[0] );
		} else {
		    lambda_two = preservoir_properties_->mobilitySecondPhase(cell[1], cell_sat[1] );
		}

		// Viscous (pressure driven) term.
		if (method_viscous_) {
		    const double visc_change = loc_flux*(lambda_one/(lambda_two+lambda_one));
		    const double satdiff = cell_sat[1] - cell_sat[0];
		    dS += visc_change;
		}

		// Gravity term.
		if (method_gravity_) {
		    if (cell[0] != cell[1]) {
			// We only add gravity flux on internal or periodic faces.
			const double grav_change = loc_gravity_flux*(lambda_one*lambda_two/(lambda_two+lambda_one));
			dS += grav_change;
		    }
		}

		// Capillary term.
		if (method_capillary_) {
		    // J(s_w) = \frac{p_c(s_w)\sqrt{k/\phi}}{\sigma \cos\theta}
		    // p_c = \frac{J \sigma \cos\theta}{\sqrt{k/\phi}}
		    const double cap_vel = inner(loc_normal, prod(aver_perm, estimateCapPressureGradient(f, nbface, saturation)));
		    const double loc_cap_flux = cap_vel*loc_area;
		    const double cap_change = loc_cap_flux*(aver_lambda_two*aver_lambda_one/(aver_lambda_one + aver_lambda_two));
		    dS += cap_change;
		}

		// Modify saturation.
		if (cell[0] != cell[1]){
		    sat_change_[cell[0]] -= dS/cell_pvol[0];
		    sat_change_[cell[1]] += dS/cell_pvol[1];
		} else {
		    ASSERT(cell[0] == cell[1]);
		    sat_change_[cell[0]] -= dS/cell_pvol[0];
		}
	    }
	}
	// wellDelta(saturation, dt, flowsys, grid, rock_data, well_data, sat_change);
    }

} // end namespace Dune


#endif // OPENRS_EULERUPSTREAM_IMPL_HEADER
