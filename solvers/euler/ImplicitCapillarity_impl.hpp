//===========================================================================
//
// File: ImplicitCapillarity_impl.hpp
//
// Created: Thu May  6 15:36:07 2010
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Jostein R Natvig    <jostein.r.natvig@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2010 Statoil ASA.

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

#ifndef OPENRS_IMPLICITCAPILLARITY_IMPL_HEADER
#define OPENRS_IMPLICITCAPILLARITY_IMPL_HEADER



#include <cassert>
#include <cmath>
#include <algorithm>

#include <dune/common/ErrorMacros.hpp>
#include <dune/common/Average.hpp>
#include <dune/common/Units.hpp>
#include <dune/common/RootFinders.hpp>
#include <dune/grid/common/Volumes.hpp>
#include <dune/solvers/common/ReservoirPropertyFixedMobility.hpp>
#include <dune/solvers/euler/MatchSaturatedVolumeFunctor.hpp>

namespace Dune
{


    template <class GI, class RP, class BC, template <class, class> class IP>
    inline ImplicitCapillarity<GI, RP, BC, IP>::ImplicitCapillarity()
	: method_viscous_(true),
	  method_gravity_(true),
	  check_sat_(true),
	  clamp_sat_(false),
          residual_tolerance_(1e-8),
          linsolver_verbosity_(1),
          linsolver_type_(1),
          update_relaxation_(1.0)
    {
    }

    template <class GI, class RP, class BC, template <class, class> class IP>
    inline ImplicitCapillarity<GI, RP, BC, IP>::ImplicitCapillarity(const GI& g, const RP& r, const BC& b)
	: residual_(g, r, b),
          method_viscous_(true),
	  method_gravity_(true),
	  check_sat_(true),
	  clamp_sat_(false),
          residual_tolerance_(1e-8),
          linsolver_verbosity_(1),
          linsolver_type_(1),
          update_relaxation_(1.0)
    {
        FieldVector<double, GI::Dimension> grav(0.0);
        psolver_.init(g, r, grav, b);
    }



    template <class GI, class RP, class BC, template <class, class> class IP>
    inline void ImplicitCapillarity<GI, RP, BC, IP>::init(const parameter::ParameterGroup& param)
    {
	method_viscous_ = param.getDefault("method_viscous", method_viscous_);
	method_gravity_ = param.getDefault("method_gravity", method_gravity_);
	check_sat_ = param.getDefault("check_sat", check_sat_);
	clamp_sat_ = param.getDefault("clamp_sat", clamp_sat_);
	residual_tolerance_ = param.getDefault("residual_tolerance", residual_tolerance_);
	linsolver_verbosity_ = param.getDefault("linsolver_verbosity", linsolver_verbosity_);
	linsolver_type_ = param.getDefault("linsolver_type", linsolver_type_);
        update_relaxation_ = param.getDefault("update_relaxation", update_relaxation_);
    }

    template <class GI, class RP, class BC, template <class, class> class IP>
    inline void ImplicitCapillarity<GI, RP, BC, IP>::init(const parameter::ParameterGroup& param,
                                                          const GI& g, const RP& r, const BC& b)
    {
	init(param);
	initObj(g, r, b);
    }


    template <class GI, class RP, class BC, template <class, class> class IP>
    inline void ImplicitCapillarity<GI, RP, BC, IP>::initObj(const GI& g, const RP& r, const BC& b)
    {
        residual_.initObj(g, r, b);
        FieldVector<double, GI::Dimension> grav(0.0);
        psolver_.init(g, r, grav, b);

    }


    template <class GI, class RP, class BC, template <class, class> class IP>
    template <class PressureSolution>
    void ImplicitCapillarity<GI, RP, BC, IP>::transportSolve(std::vector<double>& saturation,
                                                             const double /*time*/,
                                                             const typename GI::Vector& gravity,
                                                             const PressureSolution& pressure_sol,
                                                             const SparseVector<double>& injection_rates) const
    {
        // Start a timer.
        time::StopWatch clock;
        clock.start();

        // Compute capillary mobilities.
        typedef typename RP::Mobility Mob;
        int num_cells = saturation.size();
        std::vector<Mob> cap_mob(num_cells);
        for (int c = 0; c < num_cells; ++c) {
            Mob& m = cap_mob[c];
            residual_.reservoirProperties().phaseMobility(0, c, saturation[c], m.mob);
            Mob mob2;
            residual_.reservoirProperties().phaseMobility(1, c, saturation[c], mob2.mob);
            Mob mob_tot;
            mob_tot.setToSum(m, mob2);
            Mob mob_totinv;
            mob_totinv.setToInverse(mob_tot);
            m *= mob_totinv;
            m *= mob2;
        }
        ReservoirPropertyFixedMobility<Mob> capillary_mobilities(cap_mob);

        // Set up boundary conditions.
        BC cap_press_bcs(residual_.boundaryConditions());
        for (int i = 0; i < cap_press_bcs.size(); ++i) {
            if (cap_press_bcs.flowCond(i).isPeriodic()) {
                cap_press_bcs.flowCond(i) = FlowBC(FlowBC::Periodic, 0.0);
            }
        }

        // Compute injection rates from residual.
        std::vector<double> injection_rates_residual(num_cells);
	residual_.computeResidual(saturation, gravity, pressure_sol, injection_rates,
                                  method_viscous_, method_gravity_, false,
                                  injection_rates_residual);
        for (int i = 0; i < num_cells; ++i) {
            injection_rates_residual[i] = -injection_rates_residual[i];
        }

        // Compute capillary pressure.
        // Note that the saturation is just a dummy for this call, since the mobilities are fixed.
        psolver_.solve(capillary_mobilities, saturation, cap_press_bcs, injection_rates_residual,
                        residual_tolerance_, linsolver_verbosity_, linsolver_type_);

        // Solve for constant to change capillary pressure solution by.
        std::vector<double> cap_press(num_cells);
        const PressureSolution& pcapsol = psolver_.getSolution();
        for (CIt c = residual_.grid().cellbegin(); c != residual_.grid().cellend(); ++c) {
            cap_press[c->index()] = pcapsol.pressure(c);
        }
        MatchSaturatedVolumeFunctor<GI, RP> functor(residual_.grid(),
                                                    residual_.reservoirProperties(),
                                                    saturation,
                                                    cap_press);
        double min_cap_press = *std::min_element(cap_press.begin(), cap_press.end());
        double max_cap_press = *std::max_element(cap_press.begin(), cap_press.end());
        double cap_press_range = max_cap_press - min_cap_press;
        double mod_low = 1e100;
        double mod_high = -1e100;
        bracketZero(functor, 0.0, cap_press_range, mod_low, mod_high);
        const int max_iter = 40;
        const double nonlinear_tolerance = 1e-12;
        int iterations_used = -1;
        double mod_correct = modifiedRegulaFalsi(functor, mod_low, mod_high, max_iter, nonlinear_tolerance, iterations_used);
        std::cout << "Moved capillary pressure solution by " << mod_correct << " after "
                  << iterations_used << " iterations." << std::endl;
        // saturation = functor.lastSaturations();
        const std::vector<double>& sat_new = functor.lastSaturations();
        for (int i = 0; i < num_cells; ++i) {
            saturation[i] = (1.0 - update_relaxation_)*saturation[i] + update_relaxation_*sat_new[i];
        }

        // Optionally check and/or clamp results.
	if (check_sat_ || clamp_sat_) {
	    checkAndPossiblyClampSat(saturation);
	}

        // Stop timer and optionally print seconds taken.
        clock.stop();
#ifdef VERBOSE
        std::cout << "Seconds taken by transport solver: " << clock.secsSinceStart() << std::endl;
#endif // VERBOSE
    }



    template <class GI, class RP, class BC, template <class, class> class IP>
    inline void ImplicitCapillarity<GI, RP, BC, IP>::checkAndPossiblyClampSat(std::vector<double>& s) const
    {
	int num_cells = s.size();
	for (int cell = 0; cell < num_cells; ++cell) {
	    if (s[cell] > 1.0 || s[cell] < 0.0) {
		if (clamp_sat_) {
		    s[cell] = std::max(std::min(s[cell], 1.0), 0.0);
		} else if (s[cell] > 1.001 || s[cell] < -0.001) {
		    THROW("Saturation out of range in ImplicitCapillarity: Cell " << cell << "   sat " << s[cell]);
		}
	    }
	}
    }



} // end namespace Dune



#endif // OPENRS_IMPLICITCAPILLARITY_IMPL_HEADER
