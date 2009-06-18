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
#include "../../grid/common/ErrorMacros.hpp"
// #include "EntityType.hpp"
// #include "FullMatrix.hpp"
// #include "CflCalculator.hpp"
// #include "GridHelper.hpp"
// #include "Average.hpp"


namespace Dune
{


    inline EulerUpstream::EulerUpstream()
	: method_viscous_(true),
	  method_gravity_(true),
	  method_capillary_(true),
	  courant_number_(0.5),
	  minimum_small_steps_(1),
	  check_sat_(true),
	  clamp_sat_(false)
    {
    }



    inline void EulerUpstream::init(const parameter::ParameterGroup& param)
    {
	courant_number_ = param.getDefault("courant_number", courant_number_);
	method_viscous_ = param.getDefault("method_viscous", method_viscous_);
	method_gravity_ = param.getDefault("method_gravity", method_gravity_);
	method_capillary_ = param.getDefault("method_capillary", method_capillary_);
	minimum_small_steps_ = param.getDefault("minimum_small_steps", minimum_small_steps_);
	check_sat_ = param.getDefault("check_sat", check_sat_);
	clamp_sat_ = param.getDefault("clamp_sat", clamp_sat_);
    }



    inline void EulerUpstream::display()
    {
	using namespace std;
	cout << endl;
	cout <<"Displaying some members of EulerUpstream" << endl;
	cout << endl;
	cout << "courant_number = " << courant_number_ << endl;
    }



    inline void EulerUpstream::setCourantNumber(double cn)
    {
	courant_number_ = cn;
    }
    
#if 0

    template <class grid_t,class flowsys_t, class rock_t, class well_t, class boundary_t>
    inline void EulerUpstream::transportSolve(std::vector<double>& saturation,
					      const double time,
					      const flowsys_t& flow_system,
					      const grid_t& my_grid,
					      const rock_t& rock_data,
					      const well_t& well_data,
					      const boundary_t& boundary,
					      const typename grid_t::point_t& gravity,
					      const std::vector<double>& flux_vector) const
    {
	// Initialize periodic boundary data. It would be more optimal to do this once,
	// instead of for every transportSolve(), but that would not fit well with being able
	// to give a new boundary object (and grid etc.) for every call.
	initPeriodics(boundary);

	// Compute the cfl time-step.
	double cfl_dt = computeCflTime(my_grid, flow_system, rock_data, flux_vector, time, gravity);

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
				  flow_system,
				  my_grid,
				  rock_data,
				  well_data,
				  boundary,
				  gravity,
				  flux_vector);
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






    /*
    template <class grid_t, class rock_t, class fluid_t>
    inline typename grid_t::point_t EulerUpstream::
    estimateCapPressureGradient(const grid_t& grid,
				const rock_t& rock,
				const fluid_t& fluid,
				int cell,
				int face,
				const std::vector<double>& sat) const
    {
	// Find neighbouring cell and face: nbcell and nbface.
	// If we are not on a periodic boundary, nbface is of course equal to face.
	typedef typename grid_t::point_t point_t;
	int nbcell = grid_helper::neighbourCell(grid, cell, face);
	int nbface = face;
	if (nbcell == -1) {
	    std::tr1::unordered_map<int,int>::iterator it = periodic_partner_.find(face);
	    if (it == periodic_partner_.end()) {
		// At (nonperiodic) boundaries, we just return a zero gradient.
		point_t g;
		g.zero();
		return g;
	    } else {
		nbface = it->second;
		nbcell = grid.template neighbours<grid::FaceType, grid::CellType>(nbface)[0];
	    }
	}

	// Estimate the gradient like a finite difference between
	// cell centers, except that in order to handle periodic
	// conditions we pass through the face centroid(s).
	point_t cell_c = grid.getCellCentroid(cell);
	point_t nb_c = grid.getCellCentroid(nbcell);
	point_t f_c = grid.getFaceCentroid(face);
	point_t nbf_c = grid.getFaceCentroid(nbface);
	double d0 = cell_c.dist(f_c);
	double d1 = nb_c.dist(nbf_c);
	double cp0 = fluid.capPressure(cell, sat[cell],
				       rock.getPermeability(cell).trace()/rock_t::Dimension,
				       rock.getPorosity(cell));
	double cp1 = fluid.capPressure(nbcell, sat[nbcell],
				       rock.getPermeability(nbcell).trace()/rock_t::Dimension,
				       rock.getPorosity(nbcell));
	double val = (cp1 - cp0)/(d0 + d1);
	point_t dir = nb_c - nbf_c + f_c - cell_c;
	dir.normalize();
	return dir*val;
    }
    */




    template <class boundary_t>
    inline void EulerUpstream::initPeriodics(const boundary_t& boundary) const
    {
	// Store the periodic boundaries, if any.
	periodic_partner_.clear();
	for (int i = 0; i < boundary.getNumberOfPeriodicConditions(); ++i) {
	    std::pair<int, int> partners = boundary.getPeriodicConditionFaces(i);
	    periodic_partner_[partners.first] = partners.second;
	    periodic_partner_[partners.second] = partners.first;
	}

    }





    template <class grid_t, class fluid_t, class rock_t>
    inline double EulerUpstream::computeCflTime(const grid_t& my_grid,
						     const fluid_t& flow_system,
						     const rock_t rock_data,
						     const std::vector<double>& flux_vector,
						     double time,
						     const typename grid_t::point_t& gravity) const
    {
	// Deal with cfl stuff, compute the necessary number of small time steps.
	double cfl_dt_v = 1e99;
	double cfl_dt_g = 1e99;
	double cfl_dt_c = 1e99;

	// Viscous cfl.
	if (method_viscous_) {
	    cfl_dt_v =
		cfl_calculator::findCFLtimeVelocity(my_grid,
						    flow_system,
						    flux_vector,
						    rock_data);
#ifdef VERBOSE
	    std::cout << "CFL dt for velocity is " << cfl_dt_v/samcode::units::DAYS2SECONDS
		      << " and total impes time is " << time/samcode::units::DAYS2SECONDS
		      << " in days." << std::endl;
#endif // VERBOSE
	}

	// Gravity cfl.
	if (method_gravity_) {
	    cfl_dt_g =
		cfl_calculator::findCFLtimeGravity(my_grid,
						   flow_system,
						   gravity,
						   rock_data);
#ifdef VERBOSE
	    std::cout << "CFL dt for gravity is " << cfl_dt_g/samcode::units::DAYS2SECONDS
		      << " and total impes time is " << time/samcode::units::DAYS2SECONDS
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
	std::cout << "Final modified CFL dt is " << cfl_dt/samcode::units::DAYS2SECONDS
		  << " and total impes time is " << time/samcode::units::DAYS2SECONDS
		  << " in days." << std::endl;
#endif // VERBOSE
	return cfl_dt;
    }






    inline void EulerUpstream::checkAndPossiblyClampSat(std::vector<double>& s) const
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




	
    template <class grid_t,class flowsys_t, class rock_t, class well_t, class boundary_t>
    inline void EulerUpstream::smallTimeStep(std::vector<double>& saturation,
						  const double dt,
						  const flowsys_t& flowsys,
						  const grid_t& grid,
						  const rock_t& rock_data,
						  const well_t& well_data,
						  const boundary_t& boundary,
						  const typename grid_t::point_t& gravity,
						  const std::vector<double>& flux_vector) const
    {
	std::vector<double> sat_change;
	computeSatDelta(saturation, dt, flowsys, grid, rock_data, well_data,
			boundary, gravity, flux_vector, sat_change);
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



    template <class grid_t,class flowsys_t, class rock_t, class well_t, class boundary_t>
    inline void EulerUpstream::computeSatDelta(const std::vector<double>& saturation,
						    const double dt,
						    const flowsys_t& flowsys,
						    const grid_t& grid,
						    const rock_t& rock_data,
						    const well_t& well_data,
						    const boundary_t& boundary,
						    const typename grid_t::point_t& gravity,
						    const std::vector<double>& flux_vector,
						    std::vector<double>& sat_change) const
    {
	using namespace grid;
	typedef typename grid_t::range_t range_t;
	typedef typename grid_t::point_t point_t;
	typedef typename rock_t::permtensor_t permtensor_t;

	// Make sure sat_change is zero, and has the right size.
	sat_change.clear();
	sat_change.resize(saturation.size(), 0.0);

	// For every face, we will modify sat_change for adjacent cells.
	int num_faces = grid.template numberOf<FaceType>();
	for (int face = 0; face < num_faces; ++face){
	    // Find the cell indices and saturations of adjacent cells.
	    double dS = 0.0;
	    const range_t face_halfface = grid.template neighbours<FaceType, HalfFaceType>(face);
	    int cell[2];
	    double cell_sat[2];
	    if (face_halfface.size() == 2) {
		cell[0] = grid.template neighbours<HalfFaceType, CellType>(face_halfface[0])[0];
		cell[1] = grid.template neighbours<HalfFaceType, CellType>(face_halfface[1])[0];
		cell_sat[0] = saturation[cell[0]];
		cell_sat[1] = saturation[cell[1]];
	    } else {
		cell[0] = grid.template neighbours<HalfFaceType, CellType>(face_halfface[0])[0];
		cell_sat[0] = saturation[cell[0]];
		std::tr1::unordered_map<int, int>::iterator it = periodic_partner_.find(face);
		if (it == periodic_partner_.end()) {
		    cell_sat[1] = boundary.getSaturationAtBoundaryFace(face);
		    cell[1] = cell[0];
		} else {
		    int nbface = it->second;
		    ASSERT(nbface != face);
		    // Periodic faces will be visited twice, but only once
		    // should they contribute. We make sure that we skip the
		    // periodic faces half the time.
		    if (nbface < face) continue;
		    cell[1] = grid.template neighbours<grid::FaceType, grid::CellType>(nbface)[0];
		    ASSERT(cell[1] != cell[0]);
		    cell_sat[1] = saturation[cell[1]];
		}
	    }

	    // Get some local properties, and compute averages.
	    const int halfface_index = face_halfface[0];
	    const double loc_area = grid.getFaceArea(face);
	    const double loc_flux = flux_vector[halfface_index];
	    const point_t loc_halfface_normal = grid.getHalfFaceNormal(halfface_index);
	    // Doing arithmetic averages. Should we consider harmonic or geometric instead?
	    using utils::arithmeticAverage;
	    const permtensor_t loc_perm = arithmeticAverage(rock_data.getPermeability(cell[0]),
							    rock_data.getPermeability(cell[1]));
	    const double average_saturation = arithmeticAverage(cell_sat[0], cell_sat[1]);
	    const double aver_mob_phase1 = arithmeticAverage(flowsys.mobilityOne(cell[0], average_saturation),
							     flowsys.mobilityOne(cell[1], average_saturation));
	    const double aver_mob_phase2 = arithmeticAverage(flowsys.mobilityTwo(cell[0], 1.0 - average_saturation), 
							     flowsys.mobilityTwo(cell[1], 1.0 - average_saturation));

	    // The local gravity flux is needed for finding the correct phase mobilities.
	    const double loc_gravity_flux = method_gravity_ ?
		loc_area*flowsys.densityDiff()*(loc_halfface_normal*(loc_perm*gravity)) : 0.0;

	    // Find the correct phasemobilities to use
	    const double flux_phase1 = loc_flux + loc_gravity_flux*aver_mob_phase2;
	    const double flux_phase2 = loc_flux - loc_gravity_flux*aver_mob_phase1;
	    double lambda_one;
	    double lambda_two;
	    // total velocity term
	    if (flux_phase1 > 0){
		lambda_one = flowsys.mobilityOne(cell[0], cell_sat[0]);
	    } else {
		lambda_one = flowsys.mobilityOne(cell[1], cell_sat[1]); 
	    }
	    if (flux_phase2 > 0){
		lambda_two = flowsys.mobilityTwo(cell[0], 1.0 - cell_sat[0] );
	    } else {
		lambda_two = flowsys.mobilityTwo(cell[1], 1.0 - cell_sat[1] );
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
	    /*
	    // Capillary term.
	    if (method_capillary_) {
		// J(s_w) = \frac{p_c(s_w)\sqrt{k/\phi}}{\sigma \cos\theta}
		// p_c = \frac{J \sigma \cos\theta}{\sqrt{k/\phi}}
		point_t cap_term = loc_perm
		    *estimateCapPressureGradient(grid, rock_data, flowsys, cell[0], face, saturation)
		    *aver_mob_phase2*aver_mob_phase1/(aver_mob_phase1 + aver_mob_phase2);
		dS += cap_term*loc_halfface_normal*loc_area;
	    }
	    */
	    // Modify saturation.
	    if (cell[0] != cell[1]){
		sat_change[cell[0]] -= (dt/rock_data.getPorosity(cell[0]))*dS/grid.getCellVolume(cell[0]);
		sat_change[cell[1]] += (dt/rock_data.getPorosity(cell[1]))*dS/grid.getCellVolume(cell[1]);
	    } else {
		assert(cell[0] == cell[1]);
		sat_change[cell[0]] -= (dt/rock_data.getPorosity(cell[0]))
		    *dS/grid.getCellVolume(cell[0]);
	    }
	}

	// wellDelta(saturation, dt, flowsys, grid, rock_data, well_data, sat_change);
    }
#endif

} // end namespace Dune



#endif // OPENRS_EULERUPSTREAM_IMPL_HEADER
