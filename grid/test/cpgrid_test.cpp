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


template <class Grid>
void mygridcheck (Grid &g)
{
  /*
   * first do the compile-test: this will not produce any code but
   * fails if an interface-component is missing
   */
  GridInterface<Grid>();

  enum { dim      = Grid :: dimension };
  enum { dimworld = Grid :: dimensionworld };
  typedef typename Grid  :: ctype ctype;
  typedef typename Grid  :: GridFamily GridFamily;

  // type of GridInterface == GridDefaultImplementation 
  typedef Dune::GridDefaultImplementation<dim,dimworld,ctype,GridFamily> GridIF;
  const GridIF & gridIF = g;
  // check functionality when grid is interpreted as reference to interface
  GridInterface<GridIF>::check(gridIF);
  /*
   * now the runtime-tests
   */
  const Grid & cg = g;
  iteratorEquals(g);
  iteratorEquals(cg);
  iterate<true>(g);
  iterate<false>(cg);
  //zeroEntityConsistency(g);
  //zeroEntityConsistency(cg);
  assertNeighbor(g);
  assertNeighbor(cg);
  // note that for some grid this might fail
  // then un comment this test 
  Dune :: checkIndexSet( g, g.leafView(), Dune :: dvverb );
  for( int level = 0; level <= g.maxLevel(); ++level )
    Dune :: checkIndexSet( g, g.levelView( level ), Dune :: dvverb, true );
}


void check_cpgrid()
{
    const int dim = 3;
    typedef Dune::FieldVector<int,dim> iTupel;
    typedef Dune::FieldVector<double,dim> fTupel;
    typedef Dune::FieldVector<bool,dim> bTupel;

    std::cout << "\nCpGrid\n" << std::endl;

//     fTupel Len;
//     Len = 1.0;
//     iTupel s;
//     s = 3;
//     bTupel p;
//     p = false;
//     p[0] = p0;
//     int overlap = 1;

    // Dune::YaspGrid<dim> grid(Len,s,p,overlap);
    // grid.globalRefine(2);

    Dune::CpGrid grid;

#if 0
    mygridcheck(grid);

    // check communication interface
    checkCommunication(grid,-1,Dune::dvverb);
    for (int l=0; l<=grid.maxLevel(); ++l)
        checkCommunication(grid,l,Dune::dvverb);

    // check the method geometryInFather()
    //checkGeometryInFather(grid);
    // check the intersection iterator and the geometries it returns
    checkIntersectionIterator(grid);
#endif
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
