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
        ASSERT ((0 <= phase_index) && (Super::NumberOfPhases < 2));
	BOOST_STATIC_ASSERT(Super::NumberOfPhases == 2);

	double visc = phase_index == 0 ? Super::viscosity1_ : Super::viscosity2_;
        if (Super::rock_.size() > 0) {
            const int region = Super::cell_to_rock_[cell_index];
            ASSERT (region < int(Super::rock_.size()));
	    Super::rock_[region].kr(phase_index, saturation, phase_mob);
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
    void ReservoirPropertyCapillaryAnisotropicRelperm<dim>::computeCflFactors()
    {
	MESSAGE("Warning: Not currently computing cfl factors properly.");
	Super::cfl_factor_ = 1.0;
	Super::cfl_factor_gravity_ = 1.0;
    }

} // namespace Dune


#endif // OPENRS_RESERVOIRPROPERTYCAPILLARYANISOTROPICRELPERM_IMPL_HEADER
