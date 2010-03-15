//===========================================================================
//
// File: entityrep_test.cpp
//
// Created: Wed Aug 26 11:15:48 2009
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


#define BOOST_TEST_MODULE EntityRepTests
#include <boost/test/unit_test.hpp>

#include "../EntityRep.hpp"

using namespace Dune;


BOOST_AUTO_TEST_CASE(entity_rep)
{
    cpgrid::EntityRep<0> e1(0);
    cpgrid::EntityRep<0> e2(~0);
    cpgrid::EntityRep<0> e3(1);
    cpgrid::EntityRep<0> e4(~1);
    BOOST_CHECK(e1.orientation());
    BOOST_CHECK(!e2.orientation());
    BOOST_CHECK(e3.orientation());
    BOOST_CHECK(!e4.orientation());
    BOOST_CHECK_EQUAL(e1.index(), 0);
    BOOST_CHECK_EQUAL(e2.index(), 0);
    BOOST_CHECK_EQUAL(e3.index(), 1);
    BOOST_CHECK_EQUAL(e4.index(), 1);
    BOOST_CHECK(e1 < e2);
    BOOST_CHECK(e1 < e3);
    BOOST_CHECK(e1 < e4);
    BOOST_CHECK(e2 < e3);
    BOOST_CHECK(e2 < e4);
    BOOST_CHECK(e3 < e4);
    BOOST_CHECK(e1 == e1);
    BOOST_CHECK(!(e1 == e2));
    BOOST_CHECK(e1 != e2);
    BOOST_CHECK(!(e1 != e1));
    BOOST_CHECK_EQUAL(sizeof e1, sizeof(int));
}

BOOST_AUTO_TEST_CASE(entity_variable)
{
    // EntityVariableBase tests
    const int sz = 2;
    const double array[sz] = { 2.71828, 3.1415};
    cpgrid::EntityVariableBase<double> base;
    BOOST_CHECK(base.empty());
    base.assign(array, array + sz);
    BOOST_CHECK_EQUAL(base.size(), sz);

    // Needing some entityreps for the rest of the checks.
    cpgrid::EntityRep<0> e1(0);
    cpgrid::EntityRep<0> e2(~0);
    cpgrid::EntityRep<0> e3(1);
    cpgrid::EntityRep<0> e4(~1);

    // EntityVariable tests
    cpgrid::EntityVariable<double, 0> var;
    BOOST_CHECK(var.empty());
    var.assign(array, array + sz);
    BOOST_CHECK_EQUAL(var.size(), sz);
    BOOST_CHECK_EQUAL(var[e1], array[0]);
    BOOST_CHECK_EQUAL(var[e2], array[0]);
    BOOST_CHECK_EQUAL(var[e3], array[1]);
    BOOST_CHECK_EQUAL(var[e4], array[1]);

    // SignedEntityVariable tests
    cpgrid::SignedEntityVariable<double, 0> svar;
    BOOST_CHECK(svar.empty());
    svar.assign(array, array + sz);
    BOOST_CHECK_EQUAL(svar.size(), sz);
    BOOST_CHECK_EQUAL(svar[e1], array[0]);
    BOOST_CHECK_EQUAL(svar[e2], -array[0]);
    BOOST_CHECK_EQUAL(svar[e3], array[1]);
    BOOST_CHECK_EQUAL(svar[e4], -array[1]);
}
