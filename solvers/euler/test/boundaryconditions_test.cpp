//===========================================================================
//
// File: boundaryconditions_test.cpp
//
// Created: Mon Jun 29 21:47:17 2009
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


#define BOOST_TEST_DYN_LINK
#define NVERBOSE // to suppress our messages when throwing

#define BOOST_TEST_MODULE BoundaryConditionsTest
#include <boost/test/unit_test.hpp>

#include "../BoundaryConditions.hpp"

using namespace Dune;

BOOST_AUTO_TEST_CASE(boundarycondition)
{
    BoundaryCondition bc1;
    BOOST_CHECK(bc1.isNeumann());
    BOOST_CHECK(bc1.outflux() == 0.0);
    BOOST_CHECK(!bc1.isDirichlet());
    BoundaryCondition bc2(BoundaryCondition::Neumann, 5.0);
    BOOST_CHECK(bc2.isNeumann());
    BOOST_CHECK(bc2.outflux() == 5.0);
    BOOST_CHECK(!bc2.isDirichlet());
    BoundaryCondition bc3(BoundaryCondition::Dirichlet, 10.0);
    BOOST_CHECK(bc3.isDirichlet());
    BOOST_CHECK(bc3.pressure() == 10.0);
    BOOST_CHECK(!bc3.isNeumann());
    // Tests that only run in debug mode.
#ifndef NDEBUG
    BOOST_CHECK_THROW(bc1.pressure(), std::exception);
    BOOST_CHECK_THROW(bc2.pressure(), std::exception);
    BOOST_CHECK_THROW(bc3.outflux(), std::exception);
#endif
}

BOOST_AUTO_TEST_CASE(boundaryconditions)
{
    BoundaryConditions bc1;
    BOOST_CHECK(bc1.empty());
    BoundaryConditions bc2(2);
    BOOST_CHECK(bc2[0].isNeumann());
    BOOST_CHECK(bc2[0].outflux() == 0.0);
    BOOST_CHECK(!bc2[0].isDirichlet());
    BOOST_CHECK(bc2[1].isNeumann());
    BOOST_CHECK(bc2[1].outflux() == 0.0);
    BOOST_CHECK(!bc2[1].isDirichlet());
    bc2[1] = BoundaryCondition(BoundaryCondition::Dirichlet, 10.0);
    BOOST_CHECK(bc2[1].isDirichlet());
    BOOST_CHECK(bc2[1].pressure() == 10.0);
    BOOST_CHECK(!bc2[1].isNeumann());
    // Tests that only run in debug mode.
#ifndef NDEBUG
    BOOST_CHECK_THROW(bc2[0].pressure(), std::exception);
    BOOST_CHECK_THROW(bc2[1].outflux(), std::exception);
#endif
}
