//===========================================================================
//
// File: steadystate_test.cpp
//
// Created: Fri Aug 28 14:11:03 2009
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



#include <dune/upscaling/SteadyStateUpscaler.hpp>
#include <dune/common/Units.hpp>

using namespace Dune;
using namespace Dune::prefix;
using namespace Dune::unit;


int main(int argc, char** argv)
{
    // Initialize.
    parameter::ParameterGroup param(argc, argv);
    // MPIHelper::instance(argc,argv);
    SteadyStateUpscaler upscaler;
    upscaler.init(param);

    // First, compute an upscaled permeability.
    SteadyStateUpscaler::permtensor_t upscaled_K = upscaler.upscaleSinglePhase();
    upscaled_K *= (1.0/(milli*darcy));
    std::cout.precision(15);
    std::cout << "Upscaled K in millidarcy:\n" << upscaled_K << std::endl;

    // Then, compute some upscaled relative permeabilities.
    int num_cells = upscaler.grid().numberOfCells();
    const int num_sats = 7;
    double saturations[num_sats] = { 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 };
    // const int num_sats = 11;
    // double saturations[num_sats] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    const int num_pdrops = 5;
    double pdrops[num_pdrops] = { 1e1, 1e2, 1e3, 1e4, 1e5};
    for (int i = 0; i < num_sats; ++i) {
	// Starting every computation with a trio of uniform profiles.
	std::vector<double> init_sat(num_cells, saturations[i]);
	boost::array<std::vector<double>, 3> init_sats = {{ init_sat, init_sat, init_sat }};
	for (int j = 0; j < num_pdrops; ++j) {
	    double pdrop = pdrops[j];
	    SteadyStateUpscaler::permtensor_t upscaled_relperm
		= upscaler.upscaleSteadyState(init_sats, saturations[i], pdrop, upscaled_K);
	    std::cout << "Tensor of upscaled relperms for saturation " << saturations[i]
		      << " and pressure drop " << pdrop << ":\n" << upscaled_relperm << std::endl;
	    // Changing initial saturations for next pressure drop to equal the steady state of the last
	    init_sats = upscaler.lastSaturations();
	}
    }
}
