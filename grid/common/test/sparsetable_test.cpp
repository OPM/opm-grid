//===========================================================================
//
// File: test_sparsetable.cpp
//
// Created: Thu May 28 10:01:46 2009
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

#define BOOST_TEST_MODULE SparseTableTest
#include <boost/test/unit_test.hpp>

#include "../SparseTable.hpp"

using namespace Dune;

BOOST_AUTO_TEST_CASE(constructing)
{
    const SparseTable<int> st1;
    BOOST_CHECK(st1.empty());

    // This should be getting us a table like this:
    // ----------------
    // 0
    // <empty row>
    // 1 2
    // 3 4 5 6
    // 7 8 9
    // ----------------
    const int num_elem = 10;
    const int elem[num_elem] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    const int num_rows = 5;
    const int rowsizes[num_rows] = { 1, 0, 2, 4, 3 };
    const SparseTable<int> st2(elem, elem + num_elem, rowsizes, rowsizes + num_rows);
    BOOST_CHECK(!st2.empty());
    BOOST_CHECK_EQUAL(st2.size(), num_rows);
    BOOST_CHECK_EQUAL(st2[0][0], 0);
    BOOST_CHECK(st2[1].empty());
    BOOST_CHECK_EQUAL(st2[3][1], 4);
    BOOST_CHECK_EQUAL(st2[4][2], 9);
    BOOST_CHECK(st2[4].size() == rowsizes[4]);

    // One element too few.
    BOOST_CHECK_THROW(const SparseTable<int> st3(elem, elem + num_elem - 1, rowsizes, rowsizes + num_rows), std::exception);

    // A few elements too many.
    BOOST_CHECK_THROW(const SparseTable<int> st4(elem, elem + num_elem, rowsizes, rowsizes + num_rows - 1), std::exception);

    // Need at least one row.
    BOOST_CHECK_THROW(const SparseTable<int> st5(elem, elem + num_elem, rowsizes, rowsizes), std::exception);

    // No negative row sizes. This is checked by SparseTable only in debug mode.
#ifndef NDEBUG
    const int err_rs[num_rows] = { 1, 0, -1, 7, 3 };
    BOOST_CHECK_THROW(const SparseTable<int> st5(elem, elem + num_elem, err_rs, err_rs + num_rows), std::exception);
#endif
}


// int main()
