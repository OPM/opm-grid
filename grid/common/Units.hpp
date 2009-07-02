//===========================================================================
//
// File: Units.hpp
//
// Created: Thu Jul  2 09:19:08 2009
//
// Author(s): Halvor M Nilsen <hnil@sintef.no>
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

#ifndef OPENRS_UNITS_HEADER
#define OPENRS_UNITS_HEADER

namespace Dune
{

    namespace units
    {
	//     const double MILLIDARCY = 1.0;//9.86923e-16;
	//     const double VISCOSITY_UNIT = 1.0;//1e-3;
	//     const double DAYS2SECONDS = 1.0;//86400;
	const double MILLIDARCY = 9.86923e-16;
	const double VISCOSITY_UNIT = 1e-3;
	const double DAYS2SECONDS = 86400;
	const double FEET = 0.30479999798832;
	const double WELL_INDEX_UNIT = VISCOSITY_UNIT/(DAYS2SECONDS*1e5);
    }


}



#endif // OPENRS_UNITS_HEADER
