//===========================================================================
//
// File: simulator_test.cpp
//
// Created: Fri Aug  7 10:08:17 2009
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


#define VERBOSE
//#define USE_TBB

#include "config.h"
#include "../SimulatorTester.hpp"
#include <dune/common/mpihelper.hh>

#ifdef USE_TBB
#include <tbb/task_scheduler_init.h>
#endif

using namespace Dune;


int main(int argc, char** argv)
{
    parameter::ParameterGroup param(argc, argv);
    MPIHelper::instance(argc,argv);
#ifdef USE_TBB
    int num_threads = param.getDefault("num_threads", tbb::task_scheduler_init::default_num_threads());
    tbb::task_scheduler_init init(num_threads);
#endif
    SimulatorTester<> tester;
    tester.init(param);
    tester.run();
}

