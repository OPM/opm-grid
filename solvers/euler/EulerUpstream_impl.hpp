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
// #define USE_TBB
#ifdef USE_TBB
#include <tbb/parallel_for.h>
#endif

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

        // Build cell_iters_.
        const int num_cells_per_iter = std::min(50, pgrid_->numberOfCells());
        int counter = 0;
	for (typename GI::CellIterator c = pgrid_->cellbegin(); c != pgrid_->cellend(); ++c, ++counter) {
            if (counter % num_cells_per_iter == 0) {
                cell_iters_.push_back(c);
            }
        }
        cell_iters_.push_back(pgrid_->cellend());
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
	double max_ok_dt = 1e100;
	const double tol = 1e-10;
	int num_cells = saturation.size();
	for (int i = 0; i < num_cells; ++i) {
	    const double sc = sat_change_[i];
	    saturation[i] += dt*sc;
	    if (sc > tol) {
		max_ok_dt = std::min(max_ok_dt, (1.0 - saturation[i])/sc);
	    } else if (sc < -tol) {
		max_ok_dt = std::min(max_ok_dt, -saturation[i]/sc);
	    }
	}
// 	std::cout << "Maximum nonviolating timestep is " << max_ok_dt << " seconds\n";
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

        template <class UpstreamSolver, class PressureSolution>
        struct UpdateForCell
        {
            typedef typename UpstreamSolver::Vector Vector;
            typedef typename UpstreamSolver::FIt FIt;
            typedef typename UpstreamSolver::RP::PermTensor PermTensor;
            typedef typename UpstreamSolver::RP::MutablePermTensor MutablePermTensor;

            const UpstreamSolver& s;
            const std::vector<double>& saturation;
            const Vector& gravity;
            const PressureSolution& pressure_sol;

            UpdateForCell(const UpstreamSolver& solver,
                          const std::vector<double>& sat,
                          const Vector& grav,
                          const PressureSolution& psol)
                : s(solver), saturation(sat), gravity(grav), pressure_sol(psol)
            {
            }

            template <class CIt>
            void operator()(const CIt& c) const
            {
                // This is constant for the whole run.
                const double delta_rho = s.preservoir_properties_->densityDifference();

                // Loop over all cell faces.
                for (FIt f = c->facebegin(); f != c->faceend(); ++f) {
                    // Neighbour face, will be changed if on a periodic boundary.
                    FIt nbface = f;
                    double dS = 0.0;
                    int cell[2];
                    double cell_sat[2];
                    double cell_pvol[2];
                    cell[0] = f->cellIndex();
                    cell_sat[0] = saturation[cell[0]];
                    cell_pvol[0] = c->volume()*s.preservoir_properties_->porosity(cell[0]);
                    if (f->boundary()) {
                        int bid = f->boundaryId();
                        if (s.pboundary_->satCond(bid).isPeriodic()) {
                            nbface = s.bid_to_face_[s.pboundary_->getPeriodicPartner(bid)];
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
                            cell_pvol[1] = nbface->cell().volume()*s.preservoir_properties_->porosity(cell[1]);
                        } else {
                            ASSERT(s.pboundary_->satCond(bid).isDirichlet());
                            cell[1] = cell[0];
                            cell_sat[1] = s.pboundary_->satCond(bid).saturation();
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
                        cell_pvol[1] = f->neighbourCellVolume()*s.preservoir_properties_->porosity(cell[1]);
                    }

                    // Get some local properties.
                    const double loc_area = f->area();
                    const double loc_flux = pressure_sol.outflux(f);
                    const Vector loc_normal = f->normal();

                    // We will now try to establish the upstream directions for each
                    // phase. They may be the same, or different (due to gravity).
                    // Recall the equation for v_w (water phase velocity):
                    //   v_w  = lambda_w * (lambda_o + lambda_w)^{-1}
                    //          * (v + lambda_o * K * grad p_{cow} + lambda_o * K * (rho_w - rho_o) * g)
                    //             ^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    //     viscous term       capillary term                    gravity term
                    //
                    // For the purpose of upstream weighting, we only consider the viscous and gravity terms.
                    // The question is, in which direction does v_w and v_o point? That is, what is the sign
                    // of v_w*loc_normal and v_o*loc_normal?
                    //
                    // For the case when the mobilities are scalar, the following analysis applies:
                    // The viscous contribution to v_w is loc_area*loc_normal*f_w*v == f_w*loc_flux.
                    // Then the phase fluxes become
                    //     flux_w = f_w*(loc_flux + loc_area*loc_normal*lambda_o*K*(rho_w - rho_o)*g)
                    //                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                    //                                           := lambda_o*G (only scalar case)
                    //     flux_o = f_o*(loc_flux - lambda_w*G)
                    // In the above, we must decide where to evaluate K, and for this purpose (deciding
                    // upstream directions) we use a K averaged between the two cells.
                    // Since all mobilities and fractional flow functions are positive, the sign
                    // of one of these cases is trivial. If G >= 0, flux_w is in the same direction as
                    // loc_flux, if G <= 0, flux_o is in the same direction as loc_flux.
                    // The phase k for which flux_k and loc_flux are of same sign, is called the trivial
                    // phase in the code below.
                    //
                    // Assuming for the moment that G >=0, we know the direction of the water flux
                    // (same as loc_flux) and evaluate lambda_w in the upstream cell. Then we may use
                    // that lambda_w to evaluate flux_o using the above formula. Knowing flux_o, we know
                    // the direction of the oil flux, and can evaluate lambda_o in the corresponding
                    // upstream cell. Finally, we can use the equation for flux_w to compute that flux.
                    // The opposite case is similar.
                    //
                    // What about tensorial mobilities? In the following code, we make the assumption
                    // that the directions of vectors are not so changed by the multiplication with
                    // mobility tensors that upstream directions change. In other words, we let all
                    // the upstream logic stand as it is. This assumption may need to be revisited.
                    // A worse problem is that
                    // 1) we do not have v, just loc_area*loc_normal*v,
                    // 2) we cannot define G, since the lambdas do not commute with the dot product.

                    typedef typename UpstreamSolver::RP::Mobility Mob;
                    using utils::arithmeticAverage;
                    // Doing arithmetic averages. Should we consider harmonic or geometric instead?
                    const MutablePermTensor aver_perm
                        = arithAver(s.preservoir_properties_->permeability(cell[0]),
                                    s.preservoir_properties_->permeability(cell[1]));
                    // Computing the raw gravity influence vector = (rho_w - rho_o)Kg
                    Vector grav_influence = prod(aver_perm, gravity);
                    grav_influence *= delta_rho;
                    // Computing G. Note that we do not multiply with the mobility,
                    // so this G is wrong in case of anisotropic relperm.
                    const double G = s.method_gravity_ ?
                        loc_area*inner(loc_normal, grav_influence) 
                        : 0.0;
                    const int triv_phase = G >= 0.0 ? 0 : 1;
                    const int ups_cell = loc_flux >= 0.0 ? 0 : 1;
                    // Compute mobility of the trivial phase.
                    Mob m_ups[2];
                    s.preservoir_properties_->phaseMobility(triv_phase, cell[ups_cell],
                                                          cell_sat[ups_cell], m_ups[triv_phase].mob);
                    // Compute gravity flow of the nontrivial phase.
                    double sign_G[2] = { -1.0, 1.0 };
                    double grav_flux_nontriv = sign_G[triv_phase]*loc_area
                        *inner(loc_normal, m_ups[triv_phase].multiply(grav_influence));
                    // Find flow direction of nontrivial phase.
                    const int ups_cell_nontriv = (loc_flux + grav_flux_nontriv >= 0.0) ? 0 : 1;
                    const int nontriv_phase = (triv_phase + 1) % 2;
                    s.preservoir_properties_->phaseMobility(nontriv_phase, cell[ups_cell_nontriv],
                                                          cell_sat[ups_cell_nontriv], m_ups[nontriv_phase].mob);
                    // Now we have the upstream phase mobilities in m_ups[].
                    Mob m_tot;
                    m_tot.setToSum(m_ups[0], m_ups[1]);
                    Mob m_totinv;
                    m_totinv.setToInverse(m_tot);


                    const double aver_sat
                        = arithmeticAverage<double, double>(cell_sat[0], cell_sat[1]);

                    Mob m1c0, m1c1, m2c0, m2c1;
                    s.preservoir_properties_->phaseMobility(0, cell[0], aver_sat, m1c0.mob);
                    s.preservoir_properties_->phaseMobility(0, cell[1], aver_sat, m1c1.mob);
                    s.preservoir_properties_->phaseMobility(1, cell[0], aver_sat, m2c0.mob);
                    s.preservoir_properties_->phaseMobility(1, cell[1], aver_sat, m2c1.mob);
                    Mob m_aver[2];
                    m_aver[0].setToAverage(m1c0, m1c1);
                    m_aver[1].setToAverage(m2c0, m2c1);
                    Mob m_aver_tot;
                    m_aver_tot.setToSum(m_aver[0], m_aver[1]);
                    Mob m_aver_totinv;
                    m_aver_totinv.setToInverse(m_aver_tot);

                    /*
                      const double aver_lambda_one
                      = arithmeticAverage<double, double>(s.preservoir_properties_->mobilityFirstPhase(cell[0], aver_sat),
                      s.preservoir_properties_->mobilityFirstPhase(cell[1], aver_sat));
                      const double aver_lambda_two
                      = arithmeticAverage<double, double>(s.preservoir_properties_->mobilitySecondPhase(cell[0], aver_sat), 
                      s.preservoir_properties_->mobilitySecondPhase(cell[1], aver_sat));
                    */

                    /*
                    // The local gravity flux is needed for finding the correct phase mobilities.
                    double loc_gravity_flux = 0.0;
                    if (method_gravity_) {
		    double grav_comp = inner(loc_normal, prod(aver_perm, gravity));
		    loc_gravity_flux = loc_area*delta_rho*grav_comp;
                    }
                    // Find the correct phasemobilities to use
                    const double flux_one = loc_flux + loc_gravity_flux*aver_lambda_two.mob;
                    const double flux_two = loc_flux - loc_gravity_flux*aver_lambda_one.mob;
                    // const double flux_one = loc_flux + loc_gravity_flux*aver_lambda_two;
                    // const double flux_two = loc_flux - loc_gravity_flux*aver_lambda_one;
                    double lambda_one;
                    double lambda_two;
                    // total velocity term
                    if (flux_one > 0){
		    lambda_one = s.preservoir_properties_->mobilityFirstPhase(cell[0], cell_sat[0]);
                    } else {
		    lambda_one = s.preservoir_properties_->mobilityFirstPhase(cell[1], cell_sat[1]); 
                    }
                    if (flux_two > 0){
		    lambda_two = s.preservoir_properties_->mobilitySecondPhase(cell[0], cell_sat[0] );
                    } else {
		    lambda_two = s.preservoir_properties_->mobilitySecondPhase(cell[1], cell_sat[1] );
                    }
                    */

                    // Viscous (pressure driven) term.
                    if (s.method_viscous_) {
                        // v is not correct for anisotropic relperm.
                        Vector v(loc_normal);
                        v *= loc_flux;
                        const double visc_change = inner(loc_normal, m_ups[0].multiply(m_totinv.multiply(v)));
                        // 		    const double visc_change = (m_ups[0].mob/(m_ups[1].mob + m_ups[0].mob))*loc_flux;
                        // 		    std::cout << "New: " << visc_change_2 << "   old: " << visc_change << '\n';
                        dS += visc_change;
                    }

                    // Gravity term.
                    if (s.method_gravity_) {
                        if (cell[0] != cell[1]) {
                            // We only add gravity flux on internal or periodic faces.
                            const double grav_change = loc_area
                                *inner(loc_normal, m_ups[0].multiply(m_totinv.multiply(m_ups[1].multiply(grav_influence))));
                            // const double grav_change = (lambda_one*lambda_two/(lambda_two+lambda_one))*G;
                            // const double grav_change = (lambda_one*lambda_two/(lambda_two+lambda_one))*loc_gravity_flux;
                            dS += grav_change;
                        }
                    }

                    // Capillary term.
                    if (s.method_capillary_) {
                        // J(s_w) = \frac{p_c(s_w)\sqrt{k/\phi}}{\sigma \cos\theta}
                        // p_c = \frac{J \sigma \cos\theta}{\sqrt{k/\phi}}
                        Vector cap_influence = prod(aver_perm, s.estimateCapPressureGradient(f, nbface, saturation));
                        const double cap_change = loc_area
			    *inner(loc_normal, m_aver[0].multiply(m_aver_totinv.multiply(m_aver[1].multiply(cap_influence))));
                        // 		    const double cap_vel = inner(loc_normal, prod(aver_perm, estimateCapPressureGradient(f, nbface, saturation)));
                        // 		    const double loc_cap_flux = cap_vel*loc_area;
                        // //   		    const double cap_change = loc_cap_flux*(m_aver[1].mob*m_aver[0].mob
                        // //   							    /(m_aver[0].mob + m_aver[1].mob));
                        //  		    const double cap_change = loc_cap_flux*(aver_lambda_two*aver_lambda_one
                        //  							    /(aver_lambda_one + aver_lambda_two));
                        dS += cap_change;
                    }

                    // Modify saturation.
                    if (cell[0] != cell[1]){
                        s.sat_change_[cell[0]] -= dS/cell_pvol[0];
                        s.sat_change_[cell[1]] += dS/cell_pvol[1];
                    } else {
                        ASSERT(cell[0] == cell[1]);
                        s.sat_change_[cell[0]] -= dS/cell_pvol[0];
                    }
                }
            }
        };

        template <typename Iter>
        struct IndirectRange
        {
            typedef Iter Iterator;
            IndirectRange(const std::vector<Iter>& iters)
                : iters_(iters), beg_(0), end_(iters_.size() - 1)
            {
                ASSERT(iters_.size() >= 2);
            }
#ifdef USE_TBB
            IndirectRange(IndirectRange& r, tbb::split)
                : iters_(r.iters_)
            {
                int m = (r.beg_ + r.end_)/2;
                beg_ = m;
                end_ = r.end_;
                r.end_ = m;
            }
#endif
            bool empty() const
            {
                return beg_ == end_;
            }
            bool is_divisible() const
            {
                return end_ - beg_ > 1;
            }
            Iter begin() const
            {
                return iters_[beg_];
            }
            Iter end() const
            {
                return iters_[end_];
            }
        private:
            const std::vector<Iter>& iters_;
            int beg_;
            int end_;
        };

        template <class Updater>
        struct UpdateLoopBody
        {
            UpdateLoopBody(const Updater& upd)
                : updater(upd)
            {
            }
            const Updater& updater;
            template <class Range>
            void operator()(const Range& r) const
            {
                typename Range::Iterator c = r.begin();
                typename Range::Iterator cend = r.end();
                for (; c != cend; ++c) {
                    updater(c);
                }
            }
        };

    } // anon namespace



    template <class GI, class RP, class BC>
    template <class PressureSolution>
    inline void EulerUpstream<GI, RP, BC>::computeSatDelta(const std::vector<double>& saturation,
							   const typename GI::Vector& gravity,
							   const PressureSolution& pressure_sol) const
    {
	// Make sure sat_change is zero, and has the right size.
	sat_change_.clear();
	sat_change_.resize(saturation.size(), 0.0);

	// For every face, we will modify sat_change for adjacent cells.
	// We loop over every cell and intersection, and modify only if
	// this cell has lower index than the neighbour, or we are on the boundary.
        typedef UpdateForCell<EulerUpstream<GI,RP,BC>, PressureSolution> CellUpdater;
        CellUpdater update_cell(*this, saturation, gravity, pressure_sol);
        UpdateLoopBody<CellUpdater> body(update_cell);
        IndirectRange<CIt> r(cell_iters_);
        // #define USE_TBB
#ifdef USE_TBB
        tbb::parallel_for(r, body);
#else
        body(r);
#endif
	// wellDelta(saturation, dt, flowsys, grid, rock_data, well_data, sat_change);
    }

} // end namespace Dune


#endif // OPENRS_EULERUPSTREAM_IMPL_HEADER
