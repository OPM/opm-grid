//===========================================================================
//
// File: orientedentitytable_test.cpp
//
// Created: Wed Aug 26 11:17:20 2009
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


#define BOOST_TEST_MODULE OrientedEntityTableTests
#include <boost/test/unit_test.hpp>
#include <sstream>

#include "../OrientedEntityTable.hpp"

using namespace Dune;


BOOST_AUTO_TEST_CASE(oriented_entity_table)
{
    // The data below corresponds to this case:
    //
    //    --------4--------------6-------
    //    |               |             |
    //    |               |             |
    //    0    cell 0     1   cell 1    2
    //    |               |             |
    //    |               |             |
    //    --------3--------------5-------
    //
    // All face normals are considered to be pointing towards increasing X or Y coordinates.
    // A face and cell has positive mutual orientation if the normal points out of the cell.
    //
    // Cell to face data.
    const int num_data = 8;
    const int data[num_data] = { ~0, 1, ~3, 4, ~1, 2, ~5, 6 };
    const int num_rows = 2;
    const int row_sizes[num_rows] = { 4, 4 };
    // Face to cell data.
    const int num_data2 = 8;
    const int data2[num_data2] = { ~0, 0, ~1, 1, ~0, 0, ~1, 1 };
    const int num_rows2 = 7;
    const int row_sizes2[num_rows2] = { 1, 2, 1, 1, 1, 1, 1 };

    // Needing some entityreps for the rest of the checks.
    const cpgrid::EntityRep<0> e1(0);
    const cpgrid::EntityRep<0> e2(~0);
    const cpgrid::EntityRep<0> e3(1);
    const cpgrid::EntityRep<0> e4(~1);

    // OrientedEntityTable tests.
    const cpgrid::OrientedEntityTable<0, 1> none;
    BOOST_CHECK(none.empty());
    const cpgrid::OrientedEntityTable<0, 1> cell2face(data, data + num_data, row_sizes, row_sizes + num_rows);
    BOOST_CHECK(!cell2face.empty());
    BOOST_CHECK_EQUAL(cell2face.size(), num_rows);
    const cpgrid::OrientedEntityTable<0, 1>::row_type c1 = cell2face[e1];
    const cpgrid::OrientedEntityTable<0, 1>::row_type c2 = cell2face[e2];
    const cpgrid::OrientedEntityTable<0, 1>::row_type c3 = cell2face[e3];
    const cpgrid::OrientedEntityTable<0, 1>::row_type c4 = cell2face[e4];
    BOOST_CHECK_EQUAL(c1.size(), 4);
    BOOST_CHECK_EQUAL(c2.size(), 4);
    BOOST_CHECK_EQUAL(c3.size(), 4);
    BOOST_CHECK_EQUAL(c4.size(), 4);
    // BOOST_CHECK_EQUAL(c1[0], cpgrid::EntityRep<1>(~0)); // Why doesn't this compile?
    BOOST_CHECK(c1[0] == cpgrid::EntityRep<1>(~0));
    BOOST_CHECK(c2[0] == cpgrid::EntityRep<1>(0));
    BOOST_CHECK(c3[0] == cpgrid::EntityRep<1>(~1));
    BOOST_CHECK(c4[0] == cpgrid::EntityRep<1>(1));

    // Testing makeInverseRelation() method.
    cpgrid::OrientedEntityTable<1, 0> face2cell_byinv;
    cell2face.makeInverseRelation(face2cell_byinv);
    const cpgrid::OrientedEntityTable<1, 0> face2cell(data2, data2 + num_data2, row_sizes2, row_sizes2 + num_rows2);
    BOOST_CHECK(face2cell == face2cell_byinv);
    cpgrid::OrientedEntityTable<0, 1> cell2face_byinv;
    face2cell.makeInverseRelation(cell2face_byinv);
    BOOST_CHECK(cell2face == cell2face_byinv);

    // Printing test
    std::string expect1(" -1  1  0 -1  1  0  0\n"
			"  0 -1  1  0  0 -1  1\n");
    std::ostringstream s1;
    cell2face.printRelationMatrix(s1);
    BOOST_CHECK(expect1 == s1.str());
    std::string expect2(" -1  0\n"
			"  1 -1\n"
			"  0  1\n"
			" -1  0\n"
			"  1  0\n"
			"  0 -1\n"
			"  0  1\n");
    std::ostringstream s2;
    face2cell.printRelationMatrix(s2);
    BOOST_CHECK(expect2 == s2.str());
}

