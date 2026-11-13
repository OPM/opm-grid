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

  This file is part of The Open Porous Media project  (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <config.h>

#define NVERBOSE // to suppress our messages when throwing


#define BOOST_TEST_MODULE OrientedEntityTableTests
#include <boost/test/unit_test.hpp>
#include <sstream>

#include <opm/grid/cpgrid/OrientedEntityTable.hpp>

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
    typedef cpgrid::EntityRep<0> E0;
    typedef cpgrid::EntityRep<1> E1;
    const int num_data = 8;
    const E1 data[num_data] = { E1(0, false), E1(1, true), E1(3, false), E1(4, true),
                                E1(1, false), E1(2, true), E1(5, false), E1(6, true) };
    const int num_rows = 2;
    const int row_sizes[num_rows] = { 4, 4 };
    // Face to cell data.
    const int num_data2 = 8;
    const E0 data2[num_data2] = { E0(0, false), E0(0, true), E0(1, false), E0(1, true),
                                  E0(0, false), E0(0, true), E0(1, false), E0(1, true) };
    const int num_rows2 = 7;
    const int row_sizes2[num_rows2] = { 1, 2, 1, 1, 1, 1, 1 };

    // Needing some entityreps for the rest of the checks.
    const cpgrid::EntityRep<0> e1(0, true);
    const cpgrid::EntityRep<0> e2(0, false);
    const cpgrid::EntityRep<0> e3(1, true);
    const cpgrid::EntityRep<0> e4(1, false);

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
    // BOOST_CHECK_EQUAL(c1[0], cpgrid::EntityRep<1>(0, false)); // Why doesn't this compile?
    BOOST_CHECK(c1[0] == cpgrid::EntityRep<1>(0, false));
    BOOST_CHECK(c2[0] == cpgrid::EntityRep<1>(0, true));
    BOOST_CHECK(c3[0] == cpgrid::EntityRep<1>(1, false));
    BOOST_CHECK(c4[0] == cpgrid::EntityRep<1>(1, true));

    // Testing makeInverseRelation() method.
    const auto face2cell_byinv = cell2face.makeInverseRelation();
    const cpgrid::OrientedEntityTable<1, 0> face2cell(data2, data2 + num_data2, row_sizes2, row_sizes2 + num_rows2);
    BOOST_CHECK(face2cell == face2cell_byinv);
    const auto cell2face_byinv = face2cell.makeInverseRelation();
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

BOOST_AUTO_TEST_CASE(inverse_relation_skips_invalid_entities)
{
    using Cell = cpgrid::EntityRep<0>;
    using Face = cpgrid::EntityRep<1>;
    const Cell data[] = {
        Cell(Cell::InvalidIndex, false), Cell(2, true),
        Cell(0, false), Cell(Cell::InvalidIndex, true),
        Cell(2, false), Cell(0, true)
    };
    const int row_sizes[] = {2, 2, 2, 0};
    const cpgrid::OrientedEntityTable<1, 0> faces(data, data + 6, row_sizes, row_sizes + 4);

    // Inversion preserves orientations and orders incidences by face index,
    // even when source rows have different cell ordering and remote neighbours.
    const Face expected_data[] = {Face(1, false), Face(2, true), Face(0, true), Face(2, false)};
    const int expected_sizes[] = {2, 0, 2};
    const cpgrid::OrientedEntityTable<0, 1> expected(
        expected_data, expected_data + 4, expected_sizes, expected_sizes + 3);
    BOOST_CHECK(faces.makeInverseRelation() == expected);
    BOOST_CHECK(faces.makeInverseRelation(3) == expected);

    const int extended_sizes[] = {2, 0, 2, 0, 0};
    const cpgrid::OrientedEntityTable<0, 1> extended(
        expected_data, expected_data + 4, extended_sizes, extended_sizes + 5);
    BOOST_CHECK(faces.makeInverseRelation(5) == extended);
}

BOOST_AUTO_TEST_CASE(inverse_relation_empty_and_invalid_only)
{
    using Cell = cpgrid::EntityRep<0>;
    const cpgrid::OrientedEntityTable<1, 0> empty;
    BOOST_CHECK(empty.makeInverseRelation().empty());
    const auto empty_inverse = empty.makeInverseRelation(3);
    BOOST_REQUIRE_EQUAL(empty_inverse.size(), 3);
    BOOST_CHECK_EQUAL(empty_inverse.dataSize(), 0);
    for (int cell = 0; cell < empty_inverse.size(); ++cell) {
        BOOST_CHECK(empty_inverse[Cell(cell, true)].empty());
    }

    const Cell data[] = {Cell(Cell::InvalidIndex, true), Cell(Cell::InvalidIndex, false)};
    const int row_sizes[] = {1, 1};
    const cpgrid::OrientedEntityTable<1, 0> invalid_only(data, data + 2, row_sizes, row_sizes + 2);
    BOOST_CHECK(invalid_only.makeInverseRelation().empty());
    BOOST_CHECK(invalid_only.makeInverseRelation(0).empty());
    BOOST_CHECK(invalid_only.makeInverseRelation(3) == empty_inverse);

    // The same sentinel rule applies when the destination entities are faces.
    using Face = cpgrid::EntityRep<1>;
    const Face face_data[] = {Face(Face::InvalidIndex, false), Face(0, true)};
    const int face_row_sizes[] = {2};
    const cpgrid::OrientedEntityTable<0, 1> cells(face_data, face_data + 2, face_row_sizes, face_row_sizes + 1);
    const auto inverse = cells.makeInverseRelation();
    BOOST_REQUIRE_EQUAL(inverse.size(), 1);
    const auto row = inverse[Face(0, true)];
    BOOST_REQUIRE_EQUAL(row.size(), 1);
    BOOST_CHECK(row[0] == Cell(0, true));
}

BOOST_AUTO_TEST_CASE(inverse_relation_validates_destination_size)
{
    using Cell = cpgrid::EntityRep<0>;
    const Cell data[] = {Cell(2, true)};
    const int row_sizes[] = {1};
    const cpgrid::OrientedEntityTable<1, 0> faces(data, data + 1, row_sizes, row_sizes + 1);
    BOOST_CHECK_THROW(faces.makeInverseRelation(-1), std::invalid_argument);
    BOOST_CHECK_THROW(faces.makeInverseRelation(0), std::out_of_range);
    BOOST_CHECK_THROW(faces.makeInverseRelation(2), std::out_of_range);
    BOOST_CHECK_NO_THROW(faces.makeInverseRelation(3));
}
