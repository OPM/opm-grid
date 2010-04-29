//===========================================================================
//
// File: ReservoirPropertyFixedMobility.hpp
//
// Created: Thu Apr 15 14:50:22 2010
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

#ifndef OPENRS_RESERVOIRPROPERTYFIXEDMOBILITY_HEADER
#define OPENRS_RESERVOIRPROPERTYFIXEDMOBILITY_HEADER


#include <vector>


namespace Dune {

    template <class Mobility>
    class ReservoirPropertyFixedMobility {
    public:
        ReservoirPropertyFixedMobility(const std::vector<Mobility>& mobs)
            : mobs_(mobs)
        {
        }

        enum { NumberOfPhases = 2 };

        template<class Vector>
        void phaseMobilities(int cell_index, double /*saturation*/, Vector& mobility) const
        {
            mobility[0] = mobs_[cell_index].mob;
            mobility[1] = 0.0;
        }

        template <class ActualMobType>
        void phaseMobility(int phase, int cell_index, double /*saturation*/, ActualMobType& mobility) const
        {
            if (phase == 0) {
                mobility = mobs_[cell_index].mob;
            }
        }

        template<class Vector>
        void phaseDensities(int /*cell_index*/, Vector& rho) const
        {
            rho[0] = 1.0;
            rho[1] = 0.0;
        }
    private:
        const std::vector<Mobility>& mobs_;
    };
}


#endif // OPENRS_RESERVOIRPROPERTYFIXEDMOBILITY_HEADER
