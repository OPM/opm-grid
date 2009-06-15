//===========================================================================
//
// File: gie_test.cpp
//
// Created: Mon Jun 15 13:19:31 2009
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


#include <config.h>
#include <iostream>
#include <dune/grid/yaspgrid.hh>

#include "../GridInterfaceEuler.hpp"


template <class Interface>
void test_interface(const Interface& g)
{
    std::cout << "Called test_interface()" << std::endl;
    int count = 0;
    typename Interface::CellIterator c = g.cellbegin();
    for (; c != g.cellend(); ++c) {
	std::cout << "\nCell number: " << count++
		  << "\n    Cell volume   = " << c->volume()
		  << "\n    Cell centroid = " << c->centroid() << '\n';
	typename Interface::CellIterator::FaceIterator f = c->facebegin();
	int fcount = 0;
	for (; f != c->faceend(); ++f) {
	    std::cout << "        Face number: " << fcount++
		      << "\n        Face area     = " << f->area()
		      << "\n        Face centroid = " << f->centroid()
		      << "\n        Face normal   = " << f->normal() << '\n';
	}
	std::cout << "    Total face count for cell: " << fcount << '\n';
    }
    std::cout << "\nTotal cell count = " << count << std::endl;
}


template <int dim>
void check_yasp(bool p0=false) {
    typedef Dune::FieldVector<int,dim> iTupel;
    typedef Dune::FieldVector<double,dim> fTupel;
    typedef Dune::FieldVector<bool,dim> bTupel;

    std::cout << std::endl << "YaspGrid<" << dim << ">";
    if (p0) std::cout << " periodic\n";
    std::cout << std::endl << std::endl;

    fTupel Len; Len = 1.0;
    iTupel s; s = 1;
    bTupel p; p = false;
    p[0] = p0;
    int overlap = 1;

#if HAVE_MPI
    Dune::YaspGrid<dim> grid(MPI_COMM_WORLD,Len,s,p,overlap);
#else
    Dune::YaspGrid<dim> grid(Len,s,p,overlap);
#endif
    grid.globalRefine(1);

    // Test the interface
    Dune::GridInterfaceEuler<Dune::YaspGrid<dim> > gie(grid);
    test_interface(gie);
};

int main (int argc , char **argv) {
    try {
#if HAVE_MPI
	// initialize MPI
	MPI_Init(&argc,&argv);
	// get own rank
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif    
	check_yasp<1>();
	check_yasp<2>();
	check_yasp<3>();
    } catch (Dune::Exception &e) {
	std::cerr << e << std::endl;
	return 1;
    } catch (...) {
	std::cerr << "Generic exception!" << std::endl;
	return 2;
    }
  
#if HAVE_MPI
    // Terminate MPI
    MPI_Finalize();
#endif

    return 0;
};

