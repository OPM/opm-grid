/*
  Copyright 2018 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

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

#define BOOST_TEST_MODULE mappingTest

#include <boost/test/unit_test.hpp>

#include <opm/grid/utility/compressedToCartesian.hpp>
#include <opm/grid/utility/cartesianToCompressed.hpp>

BOOST_AUTO_TEST_CASE(mapping)
{
    const std::vector<int> global_cell{0, 1, 2, 3, 5, 6, 8, 9};
    const int num_cells = global_cell.size();

    const std::unordered_map<int, int> cartesian_to_compressed =
        Opm::cartesianToCompressed(num_cells, global_cell.data());

    BOOST_CHECK_EQUAL(cartesian_to_compressed.at(0), 0);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.at(1), 1);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.at(2), 2);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.at(3), 3);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.at(5), 4);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.at(6), 5);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.at(8), 6);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.at(9), 7);

    const int non_existing_index = 1829;
    BOOST_CHECK_THROW(cartesian_to_compressed.at(non_existing_index), std::out_of_range);

    const std::vector<int> compressed_to_cartesian = Opm::compressedToCartesian(num_cells,  global_cell.data());

    BOOST_CHECK_EQUAL_COLLECTIONS(compressed_to_cartesian.begin(), compressed_to_cartesian.end(),
                                  global_cell.begin(),             global_cell.end());
}

BOOST_AUTO_TEST_CASE(nullmapping)
{
    const int num_cells = 30;

    const std::unordered_map<int, int> cartesian_to_compressed =
        Opm::cartesianToCompressed(num_cells, nullptr);

    for (int i = 0; i < num_cells; ++i) {
        BOOST_CHECK_EQUAL(cartesian_to_compressed.at(i), i);
    }

    const std::vector<int> compressed_to_cartesian = Opm::compressedToCartesian(num_cells,  nullptr);

    for (int i = 0; i < num_cells; ++i) {
        BOOST_CHECK_EQUAL(compressed_to_cartesian[i], i);
    }
}
