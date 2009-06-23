//===========================================================================
//
// File: cpgrid_test.cpp
//
// Created: Fri May 29 14:07:20 2009
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

#include "../CpGrid.hpp"

//#include <config.h>
#include <iostream>

#include <dune/grid/test/gridcheck.cc>
#include <dune/grid/test/checkcommunicate.cc>
#include <dune/grid/test/checkgeometryinfather.cc>
#include <dune/grid/test/checkintersectionit.cc>


void check_cpgrid()
{
    const int dim = 3;
    std::cout << "\nCpGrid\n" << std::endl;

    //typedef Dune::FieldVector<int,dim> iTuple;
    //typedef Dune::FieldVector<double,dim> fTuple;
    //typedef Dune::FieldVector<bool,dim> bTuple;
    //    fTuple cell_sz(1.0);
    //    iTuple dims(3);
    //     bTuple p(false);
    //     p[0] = p0;
    //     int overlap = 1;
    // Dune::YaspGrid<dim> grid(Len,s,p,overlap);
    // grid.globalRefine(2);

    Dune::CpGrid grid;
    Dune::array<int, dim> dims = {{ 1, 1, 1 }};
    Dune::array<double, dim> cell_sz = {{ 1.0, 1.0, 1.0 }};
    // grid.createCartesian(dims, cell_sz);

    gridcheck(grid);

    // check communication interface
    checkCommunication(grid,-1,Dune::dvverb);
    for (int l=0; l<=grid.maxLevel(); ++l)
        checkCommunication(grid,l,Dune::dvverb);

    // check the method geometryInFather()
    //checkGeometryInFather(grid);
    // check the intersection iterator and the geometries it returns
    checkIntersectionIterator(grid);
}


int main(int /*argc*/ , char** /*argv*/)
{
    try {
        check_cpgrid();
    } catch (Dune::Exception &e) {
        std::cerr << e << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Generic exception!" << std::endl;
        return 2;
    }
    return 0;
}
