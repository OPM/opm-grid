//===========================================================================
//
// File: upscaling_test.cpp
//
// Created: Mon Aug 10 08:21:44 2009
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


#include <dune/upscaling/Upscaler.hpp>
#include <dune/common/Units.hpp>

using namespace Dune;
using namespace Dune::prefix;
using namespace Dune::unit;

int main(int argc, char** argv)
{
    parameter::ParameterGroup param(argc, argv);
    // MPIHelper::instance(argc,argv);
    Upscaler upscaler;
    upscaler.init(param);
    Upscaler::permtensor_t upscaled_K = upscaler.upscaleSinglePhase();
    upscaled_K *= (1.0/(milli*darcy));
    std::cout.precision(15);
    std::cout << "Upscaled K in millidarcy:\n" << upscaled_K << std::endl;
}
