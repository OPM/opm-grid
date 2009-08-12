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
#include <dune/grid/common/ErrorMacros.hpp>
#include <dune/grid/common/Average.hpp>
#include <dune/grid/common/Units.hpp>
#include "CflCalculator.hpp"

namespace Dune
{


    template <class GI, class RP, class BC>
    inline EulerUpstream<GI, RP, BC>::EulerUpstream(const GI& g, const RP& r, const BC& b,
						    const SparseVector<double>& injection_rates)
	: grid_(g),
	  reservoir_properties_(r),
	  boundary_(b),
	  injection_rates_(injection_rates),
	  method_viscous_(true),
	  method_gravity_(true),
	  method_capillary_(true),
	  courant_number_(0.5),
	  minimum_small_steps_(1),
	  check_sat_(true),
	  clamp_sat_(false)
    {
	initPeriodics();
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
						   const PressureSolution& pressure_sol) const
    {
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
    inline void EulerUpstream<GI, RP, BC>::initPeriodics()
    {
	// Store the periodic boundaries, if any.
	periodic_partner_.clear();
	MESSAGE("Not yet supporting periodic boundaries, should this be done in grid instead?");
#if 0
	for (int i = 0; i < boundary_.getNumberOfPeriodicConditions(); ++i) {
	    std::pair<int, int> partners = boundary_.getPeriodicConditionFaces(i);
	    periodic_partner_[partners.first] = partners.second;
	    periodic_partner_[partners.second] = partners.first;
	}
#endif
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
	    cfl_dt_v = cfl_calculator::findCFLtimeVelocity(grid_, reservoir_properties_,
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
	    cfl_dt_g = cfl_calculator::findCFLtimeGravity(grid_, reservoir_properties_, gravity);
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
    EulerUpstream<GI, RP, BC>::estimateCapPressureGradient(FIt f, const std::vector<double>& sat) const
    {
	typedef typename GI::CellIterator::FaceIterator Face;
	typedef typename Face::Cell Cell;
	typedef typename GI::Vector Vector;

	// For now, we return zero on boundaries. Obviously we should have
	// capillary pressure boundary conditions.
	if (f->boundary()) {
	    return Vector(0.0);
	}
	// Find neighbouring cell and face: nbc and nbf.
	// If we are not on a periodic boundary, nbf is of course equal to f.
	Cell c = f->cell();
	Cell nb = f->boundary() ? c : f->neighbourCell();
	Face nbf = f; // Must change for periodic boundaries
// 	if (f->boundary()) {
// 	    std::tr1::unordered_map<int,int>::iterator it = periodic_partner_.find(face);
// 	    if (it == periodic_partner_.end()) {
// 		// At (nonperiodic) boundaries, we just return a zero gradient.
// 		Vector g(0.0);
// 		return g;
// 	    } else {
// 		nbface = it->second;
// 		nbcell = grid.template neighbours<grid::FaceType, grid::CellType>(nbface)[0];
// 	    }
// 	}

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
	double cp0 = reservoir_properties_.capillaryPressure(cell, sat[cell]);
	double cp1 = reservoir_properties_.capillaryPressure(nbcell, sat[nbcell]);
	double val = (cp1 - cp0)/(d0 + d1);
	Vector res = nb_c - nbf_c + f_c - cell_c;
	res /= res.two_norm();
	res *= val;
	return res;
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
	std::vector<double> sat_change;
	computeSatDelta(saturation, dt, gravity, pressure_sol, sat_change);
	int num_cells = saturation.size();
	for (int i = 0; i < num_cells; ++i) {
	    saturation[i] += sat_change[i];
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
    inline void EulerUpstream<GI, RP, BC>::computeSatDelta(std::vector<double>& saturation,
							   const double dt,
							   const typename GI::Vector& gravity,
							   const PressureSolution& pressure_sol,
							   std::vector<double>& sat_change) const
    {
	typedef typename GI::Vector Vector;
	typedef typename RP::PermTensor PermTensor;
	typedef typename RP::MutablePermTensor MutablePermTensor;

	// Make sure sat_change is zero, and has the right size.
	sat_change.clear();
	sat_change.resize(saturation.size(), 0.0);

	// For every face, we will modify sat_change for adjacent cells.
	// We loop over every cell and intersection, and modify only if
	// this cell has lower index than the neighbour, or we are on the boundary.
	typename GI::CellIterator c = grid_.cellbegin();
	for (; c != grid_.cellend(); ++c) {
	    for (FIt f = c->facebegin(); f != c->faceend(); ++f) {
		double dS = 0.0;
		int cell[2];
		double cell_sat[2];
		if (f->boundary()) {
		    cell[0] = f->cellIndex();
		    cell_sat[0] = saturation[cell[0]];
		    typename PartnerMapType::const_iterator it = periodic_partner_.find(f);
		    if (it == periodic_partner_.end()) {
			cell[1] = cell[0];
			cell_sat[1] = boundary_[f->boundaryId()].saturation();
		    } else {
			FIt nbface = it->second;
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
		    }
		} else {
		    cell[0] = f->cellIndex();
		    cell[1] = f->neighbourCellIndex();
		    ASSERT(cell[0] != cell[1]);
		    if (cell[0] > cell[1]) {
			// We skip this face.
			continue;
		    }
		    cell_sat[0] = saturation[cell[0]];
		    cell_sat[1] = saturation[cell[1]];
		}

		// Get some local properties, and compute averages.
		const double loc_area = f->area();
		const double loc_flux = pressure_sol.outflux(f);
		const Vector loc_normal = f->normal();
		// Doing arithmetic averages. Should we consider harmonic or geometric instead?
		using utils::arithmeticAverage;
		const MutablePermTensor loc_perm
		    = arithAver(reservoir_properties_.permeability(cell[0]),
				reservoir_properties_.permeability(cell[1]));
		const double average_saturation
		    = arithmeticAverage<double, double>(cell_sat[0], cell_sat[1]);
		const double aver_mob_phase1
		    = arithmeticAverage<double, double>(reservoir_properties_.mobilityFirstPhase(cell[0], average_saturation),
							reservoir_properties_.mobilityFirstPhase(cell[1], average_saturation));
		const double aver_mob_phase2
		    = arithmeticAverage<double, double>(reservoir_properties_.mobilitySecondPhase(cell[0], average_saturation), 
							reservoir_properties_.mobilitySecondPhase(cell[1], average_saturation));

		// The local gravity flux is needed for finding the correct phase mobilities.
		double loc_gravity_flux = 0.0;
		if (method_gravity_) {
		    double grav_comp = inner(loc_normal, prod(loc_perm, gravity));
		    loc_gravity_flux = loc_area*reservoir_properties_.densityDifference()*grav_comp;
		}
		// Find the correct phasemobilities to use
		const double flux_phase1 = loc_flux + loc_gravity_flux*aver_mob_phase2;
		const double flux_phase2 = loc_flux - loc_gravity_flux*aver_mob_phase1;
		double lambda_one;
		double lambda_two;
		// total velocity term
		if (flux_phase1 > 0){
		    lambda_one = reservoir_properties_.mobilityFirstPhase(cell[0], cell_sat[0]);
		} else {
		    lambda_one = reservoir_properties_.mobilityFirstPhase(cell[1], cell_sat[1]); 
		}
		if (flux_phase2 > 0){
		    lambda_two = reservoir_properties_.mobilitySecondPhase(cell[0], cell_sat[0] );
		} else {
		    lambda_two = reservoir_properties_.mobilitySecondPhase(cell[1], cell_sat[1] );
		}

		// Viscous (pressure driven) term.
		if (method_viscous_) {
		    dS+=loc_flux*(lambda_one/(lambda_two+lambda_one));
		}

		// Gravity term.
		if (method_gravity_) {
		    if (cell[0] != cell[1]) {
			// We only add gravity flux on internal or periodic faces.
			dS+=loc_gravity_flux*(lambda_one*lambda_two/(lambda_two+lambda_one));
		    }
		}

		// Capillary term.
		if (method_capillary_) {
		    // J(s_w) = \frac{p_c(s_w)\sqrt{k/\phi}}{\sigma \cos\theta}
		    // p_c = \frac{J \sigma \cos\theta}{\sqrt{k/\phi}}
		    double cap_term = inner(loc_normal, prod(loc_perm, estimateCapPressureGradient(f, saturation)));
		    dS += cap_term*loc_area*aver_mob_phase2*aver_mob_phase1/(aver_mob_phase1 + aver_mob_phase2);
		}

		// Modify saturation.
		if (cell[0] != cell[1]){
		    sat_change[cell[0]] -= (dt/reservoir_properties_.porosity(cell[0]))*dS/c->volume();
		    sat_change[cell[1]] += (dt/reservoir_properties_.porosity(cell[1]))*dS/f->neighbourCellVolume();
		} else {
		    ASSERT(cell[0] == cell[1]);
		    sat_change[cell[0]] -= (dt/reservoir_properties_.porosity(cell[0]))
			*dS/c->volume();
		}
	    }
	}
	// wellDelta(saturation, dt, flowsys, grid, rock_data, well_data, sat_change);
    }

} // end namespace Dune


#endif // OPENRS_EULERUPSTREAM_IMPL_HEADER
