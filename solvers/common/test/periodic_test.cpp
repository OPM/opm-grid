//===========================================================================
//
// File: periodic_test.cpp
//
// Created: Fri Aug 14 13:03:17 2009
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

#include "../PeriodicHelpers.hpp"
#include <dune/solvers/common/BoundaryConditions.hpp>
#include <dune/solvers/common/GridInterfaceEuler.hpp>
#include <dune/grid/CpGrid.hpp>
#include <boost/array.hpp>
#include <dune/common/param/ParameterGroup.hpp>

using namespace Dune;

int main(int argc, char** argv)
{
    parameter::ParameterGroup param(argc, argv);
    CpGrid grid;
    grid.init(param);
    grid.setUniqueBoundaryIds(true);
    GridInterfaceEuler<CpGrid> gi(grid);
    typedef FlowBC FBC;
    boost::array<FBC, 6> fcond = {{ FBC(FBC::Periodic, 1.0e5),
				    FBC(FBC::Periodic, -1.0e5),
				    FBC(FBC::Periodic, 0.0),
				    FBC(FBC::Periodic, 0.0),
				    FBC(FBC::Neumann, 0.0),
				    FBC(FBC::Neumann, 0.0) }};
    typedef SatBC SBC;
    boost::array<SBC, 6> scond = {{ SBC(SBC::Periodic, 0.0),
				    SBC(SBC::Periodic, 0.0),
				    SBC(SBC::Periodic, 0.0),
				    SBC(SBC::Periodic, 0.0),
				    SBC(SBC::Dirichlet, 0.0),
				    SBC(SBC::Dirichlet, 0.0) }};
    BasicBoundaryConditions<true, true> bcs;
    createPeriodic(bcs, gi, fcond, scond);
    std::cout << bcs;
}

