//===========================================================================
//
// File: ReservoirPropertyTracerFluid.hpp
//
// Created: Thu Apr 15 10:45:00 2010
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

#ifndef OPENRS_RESERVOIRPROPERTYTRACERFLUID_HEADER
#define OPENRS_RESERVOIRPROPERTYTRACERFLUID_HEADER

namespace Dune {
    class ReservoirPropertyTracerFluid {
    public:
        enum { NumberOfPhases = 2 };

        template<class Vector>
        void phaseMobilities(int /*cell_index*/, double saturation, Vector& mobility) const
        {
            mobility[0] = saturation;
            mobility[1] = 1 - saturation;
        }

        template<class Vector>
        void phaseDensities(int /*cell_index*/, Vector& rho) const
        {
            rho[0] = 1.0;
            rho[1] = 1.0;
        }
    };
}

#endif /* OPENRS_RESERVOIRPROPERTYTRACERFLUID_HEADER */
