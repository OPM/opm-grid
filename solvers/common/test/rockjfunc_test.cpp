//===========================================================================
//
// File: rockjfunc_test.cpp
//
// Created: Tue Jul 20 15:28:39 2010
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

#include "../RockJfunc.hpp"
#include <dune/solvers/common/Matrix.hpp>
#include <dune/common/Units.hpp>

using namespace Dune;

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " file\n";
        return 1;
    }
    RockJfunc r;
    r.read("", argv[1]);


    OwnCMatrix perm;
    eye(perm);
    perm *= 100*prefix::milli*unit::darcy;
    double poro = 0.2;

    double sat;
    std::cout << "sat = " << std::flush;
    std::cin >> sat;
    double cp = r.capPress(perm, poro, sat);
    std::cout << "cap press = " << cp << std::endl;
    double sat2 = r.satFromCapPress(perm, poro, cp);
    std::cout << "sat2 = " << sat2 << std::endl;
}
