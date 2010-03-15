//===========================================================================
//
// File: ReservoirPropertyCapillaryAnisotropicRelperm_impl.hpp
//
// Created: Mon Oct 26 10:12:15 2009
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

#ifndef OPENRS_RESERVOIRPROPERTYCAPILLARYANISOTROPICRELPERM_IMPL_HEADER
#define OPENRS_RESERVOIRPROPERTYCAPILLARYANISOTROPICRELPERM_IMPL_HEADER

#include <boost/lambda/lambda.hpp>

namespace Dune
{

    template <int dim>
    template <class MatrixType>
    void ReservoirPropertyCapillaryAnisotropicRelperm<dim>::phaseMobility(int phase_index,
									  int cell_index,
									  double saturation,
									  MatrixType& phase_mob) const
    {
	const int region = Super::rock_.size() > 0 ? Super::cell_to_rock_[cell_index] : -1;
	phaseMobilityByRock(phase_index, region, saturation, phase_mob);
    }


    template <int dim>
    double ReservoirPropertyCapillaryAnisotropicRelperm<dim>::fractionalFlow(int cell_index, double saturation) const
    {
        // This method is a hack.
	// Assumes that the relperm is diagonal.
	Mobility m;
        double ff_first = 0.0;
        for (int direction = 0; direction < dim; ++direction) {
            phaseMobility(0, cell_index, saturation, m.mob);
            double l1 = m.mob(direction, direction);
            phaseMobility(1, cell_index, saturation, m.mob);
            double l2 = m.mob(direction, direction);
            ff_first += l1/(l1 + l2);
        }
        ff_first /= double(dim);
        return ff_first;
    }


    template <int dim>
    template <class MatrixType>
    void ReservoirPropertyCapillaryAnisotropicRelperm<dim>::phaseMobilityByRock(int phase_index,
										int rock_index,
										double saturation,
										MatrixType& phase_mob) const
    {
	ASSERT ((0 <= phase_index) && (Super::NumberOfPhases < 2));
	BOOST_STATIC_ASSERT(Super::NumberOfPhases == 2);

	double visc = phase_index == 0 ? Super::viscosity1_ : Super::viscosity2_;
	if (rock_index != -1) {
	    ASSERT (rock_index < int(Super::rock_.size()));
	    Super::rock_[rock_index].kr(phase_index, saturation, phase_mob);
	    //using namespace boost::lambda;
	    std::transform(phase_mob.data(), phase_mob.data() + dim*dim,
			   phase_mob.data(), boost::lambda::_1/visc);
	} else {
	    // HACK: With no rocks, we use quadratic relperm functions.
	    double kr = phase_index == 0 ? saturation*saturation : (1.0 - saturation)*(1.0 - saturation);
	    double mob = kr/visc;
	    zero(phase_mob);
            for (int j = 0; j < dim; ++j) {
                phase_mob(j,j) = mob;
            }
	}
    }


    template <int dim>
    void ReservoirPropertyCapillaryAnisotropicRelperm<dim>::cflFracFlows(int rock, int direction,
									 double s, double& ff_first, double& ff_gravity) const
    {
	// Assumes that the relperm is diagonal.
	Mobility m;
	phaseMobilityByRock(0, rock, s, m.mob);
	double l1 = m.mob(direction, direction);
	phaseMobilityByRock(1, rock, s, m.mob);
	double l2 = m.mob(direction, direction);
	ff_first = l1/(l1 + l2);
	ff_gravity = l1*l2/(l1 + l2);
    }




    template <int dim>
    std::pair<double, double> ReservoirPropertyCapillaryAnisotropicRelperm<dim>::computeSingleRockCflFactors(int rock) const
    {
        const int N = 257;
        double delta = 1.0/double(N - 1);
        double last_ff1, last_ffg;
        double max_der1 = -1e100;
        double max_derg = -1e100;
	for (int direction = 0; direction < dim; ++direction) {
	    cflFracFlows(rock, direction, 0.0, last_ff1, last_ffg);
	    for (int i = 1; i < N; ++i) {
		double s = double(i)*delta;
		double ff1, ffg;
		cflFracFlows(rock, direction, s, ff1, ffg);
		double est_deriv_ff1 = std::fabs(ff1 - last_ff1)/delta;
		double est_deriv_ffg = std::fabs(ffg - last_ffg)/delta;
		max_der1 = std::max(max_der1, est_deriv_ff1);
		max_derg = std::max(max_derg, est_deriv_ffg);
		last_ff1 = ff1;
		last_ffg = ffg;
	    }
        }
        return std::make_pair(1.0/max_der1, 1.0/max_derg);
    }




    template <int dim>
    void ReservoirPropertyCapillaryAnisotropicRelperm<dim>::computeCflFactors()
    {
        if (Super::rock_.empty()) {
            std::pair<double, double> fac = computeSingleRockCflFactors(-1);
            Super::cfl_factor_ = fac.first;
            Super::cfl_factor_gravity_ = fac.second;
        } else {
            Super::cfl_factor_ = 1e100;
            Super::cfl_factor_gravity_ = 1e100;
            for (int r = 0; r < int(Super::rock_.size()); ++r) {
                std::pair<double, double> fac = computeSingleRockCflFactors(r);
                Super::cfl_factor_ = std::min(Super::cfl_factor_, fac.first);
                Super::cfl_factor_gravity_ = std::min(Super::cfl_factor_gravity_, fac.second);
            }
        }
    }






} // namespace Dune


#endif // OPENRS_RESERVOIRPROPERTYCAPILLARYANISOTROPICRELPERM_IMPL_HEADER
