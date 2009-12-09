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

#include <boost/static_assert.hpp>

#include <dune/common/array.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/CpGrid.hpp>


#include "../GridInterfaceEuler.hpp"


template <class Interface>
void test_interface(const Interface& g, bool verbose = true)
{
    std::cout << "Called test_interface()" << std::endl;
    typename Interface::CellIterator c = g.cellbegin();
    int count = 0;
    for (; c != g.cellend(); ++c, ++count) {
	if (verbose) {
	    std::cout << "\nCell number: " << count
		      << "\nCell index : " << c->index()
		      << "\n    Cell volume   = " << c->volume()
		      << "\n    Cell centroid = " << c->centroid() << '\n';
	} else {
	    std::cout << "\nCell index: " << c->index()
		      << "  Cell centroid: " << c->centroid()
		      << "\n   Neighbour cells:";
	}
	typename Interface::CellIterator::FaceIterator f = c->facebegin();
        int fcount = 0;
	for (; f != c->faceend(); ++f, ++fcount) {
	    if (verbose) {
		std::cout << "        Face number: " << fcount
			  << "\n            Local index    = " << f->localIndex()
			  << "\n            Boundary       = " << f->boundary()
			  << "\n            Boundary Id    = " << f->boundaryId()
			  << "\n            My cell        = " << f->cellIndex()
			  << "\n            Neighbour cell = " << f->neighbourCellIndex()
			  << "\n            Face area      = " << f->area()
			  << "\n            Face centroid  = " << f->centroid()
			  << "\n            Face normal    = " << f->normal() << '\n';
	    } else {
		std::cout << ' ' << f->neighbourCellIndex();
	    }
	    typename Interface::CellIterator::FaceIterator::Cell this_c = f->cell();
	}
	std::cout << "    Total face count for cell: " << fcount << '\n';
    }
    std::cout << "\nTotal cell count = " << count << std::endl;
}


template <int dim, int refinement>
void check_yasp(bool p0 = false) {
    typedef Dune::FieldVector<int,dim> iTupel;
    typedef Dune::FieldVector<double,dim> fTupel;
    typedef Dune::FieldVector<bool,dim> bTupel;

    std::cout << std::endl << "YaspGrid, dim = " << dim << ", refinement = " << refinement;
    if (p0) std::cout << " periodic\n";
    std::cout << std::endl << std::endl;

    fTupel Len; Len = 1.0;
    iTupel s; s = 2;
    bTupel p; p = false;
    p[0] = p0;
    int overlap = 0;

#if HAVE_MPI
    Dune::YaspGrid<dim> grid(MPI_COMM_WORLD,Len,s,p,overlap);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
    Dune::YaspGrid<dim> grid(Len,s,p,overlap);
#endif
    grid.globalRefine(refinement);

    // Test the interface
    Dune::GridInterfaceEuler<Dune::YaspGrid<dim> > gie(grid);
#if HAVE_MPI
    //    if (rank == 0) {
    test_interface(gie, false);
    //    }
#else
    test_interface(gie);
#endif
}

template <int dim, int refinement>
void check_sgrid()
{
    typedef Dune::FieldVector<int,dim> iTupel;
    typedef Dune::FieldVector<double,dim> fTupel;

    std::cout << std::endl << "SGrid, dim = " << dim << ", refinement = " << refinement;
    std::cout << std::endl << std::endl;

    fTupel len(1.0);
    iTupel num(1);
    fTupel pos(0.0);

    Dune::SGrid<dim, dim, double> grid(num, pos, len);
    grid.globalRefine(refinement);

    // Test the interface
    Dune::GridInterfaceEuler<Dune::SGrid<dim, dim, double> > gie(grid);
    test_interface(gie);
}


//-----------------------------------------------------------------------------
template <int refinement>
void check_cpgrid()
//-----------------------------------------------------------------------------
{
    std::cout << '\n' << "CpGrid<" << refinement << ">\n" << std::endl;

    Dune::CpGrid grid;

    Dune::array<int   , 3> dims;
    std::fill(dims.begin(), dims.end(), 1 << refinement);

    Dune::array<double, 3> cell_sz;
    std::fill(cell_sz.begin(), cell_sz.end(), 1.0 / (1 << refinement));

    grid.createCartesian(dims, cell_sz);

    // Test the interface
    Dune::GridInterfaceEuler<Dune::CpGrid> gie(grid);
    test_interface(gie);
}


int main (int argc , char **argv) {
    try {
#if HAVE_MPI
	// initialize MPI
	MPI_Init(&argc,&argv);
	// get own rank
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
	check_yasp<1,1>();
	//check_yasp<2,1>();

	//check_yasp<3,0>();  // 3D, 1 x 1 x 1 cell
	//check_cpgrid<0>();

#ifdef REFINE
#undef REFINE
#endif

#define REFINE 1

#if REFINE > 0
	check_yasp<3,1>();  // 3D, 2 x 2 x 2 cells
	check_sgrid<3,1>();  // 3D, 2 x 2 x 2 cells
	check_cpgrid<1>();
#endif

#if REFINE > 1
	check_yasp<3,2>();  // 3D, 4 x 4 x 4 cells
	check_sgrid<3,2>();  // 3D, 4 x 4 x 4 cells
	check_cpgrid<2>();
#endif

#if REFINE > 2
	check_yasp<3,3>();  // 3D, 8 x 8 x 8 cells
	check_sgrid<3,3>();  // 3D, 8 x 8 x 8 cells
	check_cpgrid<3>();
#endif

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
}
