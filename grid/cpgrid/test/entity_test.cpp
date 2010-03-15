//===========================================================================
//
// File: entity_test.cpp
//
// Created: Fri May 29 14:04:50 2009
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

#define BOOST_TEST_DYN_LINK
#define NVERBOSE // to suppress our messages when throwing


#define BOOST_TEST_MODULE EntityTests
#include <boost/test/unit_test.hpp>
#include <sstream>

#include "../Entity.hpp"
#include "../../CpGrid.hpp"

using namespace Dune;



BOOST_AUTO_TEST_CASE(entity)
{
    CpGrid g;
    cpgrid::Entity<0, CpGrid> e1(g, 0);
    cpgrid::Entity<0, CpGrid> e2(g, ~0);
    cpgrid::Entity<0, CpGrid> e3(g, 1);
    cpgrid::Entity<0, CpGrid> e4(g, ~1);
    BOOST_CHECK(e1 != e2);
    BOOST_CHECK(e1 != e3);
    BOOST_CHECK(e1 != e4);
    BOOST_CHECK(e2 != e3);
    BOOST_CHECK(e2 != e4);
    BOOST_CHECK(e3 != e4);
    BOOST_CHECK_EQUAL(e1.level(), 0);
    // BOOST_CHECK(e1.type().isSingular()); // Our new type
    BOOST_CHECK(e1.type().isCube());
    BOOST_CHECK_EQUAL(e1.partitionType(), InteriorEntity);
    cpgrid::Entity<3, CpGrid> e5(g, 0);
    BOOST_CHECK(e5.type().isCube());

    // Cannot check other members without a real grid.
    // Put in more checks when it is possible to construct
    // test grids easily.
    //   geometry()
    //   count()
    //   ileafbegin()
    //   ileafend()
}


BOOST_AUTO_TEST_CASE(entity_ptr)
{
    CpGrid g;
    cpgrid::EntityPointer<0, CpGrid> p1(g, ~5);
    const cpgrid::EntityPointer<0, CpGrid> p2(g, 42);
//     cpgrid::Entity<0, CpGrid>& e1 = *p1;
//     const cpgrid::Entity<0, CpGrid>& e2 = *p2;
//     cpgrid::Entity<0, CpGrid> ee1(g, ~5);
//     cpgrid::Entity<0, CpGrid> ee2(g, 42);
//     BOOST_CHECK(e1 == ee1);
//     BOOST_CHECK(e2 == ee2);
}

