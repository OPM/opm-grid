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


//#define VERBOSE

#include <dune/upscaling/SteadyStateUpscalerManager.hpp>
#include <dune/upscaling/UpscalingTraits.hpp>

using namespace Dune;

int main(int argc, char** argv)
{
    // Initialize.
    parameter::ParameterGroup param(argc, argv);
    // MPIHelper::instance(argc,argv);

    SteadyStateUpscalerManager<UpscalingTraitsBasic> mgr;
    mgr.upscale(param);
}
