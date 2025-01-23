// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 SINTEF Digital

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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>

#define BOOST_TEST_MODULE ElementChunksTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/utility/ElementChunks.hpp>

#include <opm/grid/CpGrid.hpp>

struct Fixture
{
    Fixture()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(m_argc, m_argv);
        Opm::OpmLog::setupSimpleDefaultLogging();
    }
};

using GV = Dune::CpGrid::LeafGridView;
using EC = Opm::ElementChunks<GV>;

std::vector<int> countChunks(const EC& chunks)
{
    std::vector<int> counts(chunks.size(), 0);
    int chunk_num = 0;
    for (const auto& chunk : chunks) {
        for (const auto& elem : chunk) {
            static_cast<void>(elem); // silence unused variable warning
            ++counts[chunk_num];
        }
        ++chunk_num;
    }
    return counts;
}

void testCase(const GV& gv, const std::size_t num_chunks, const std::vector<int>& expected)
{
    EC chunks(gv, num_chunks);
    auto counts = countChunks(chunks);
    BOOST_CHECK_EQUAL_COLLECTIONS(counts.begin(), counts.end(), expected.begin(), expected.end());
}

BOOST_FIXTURE_TEST_CASE(ElementChunksTests, Fixture)
{
    std::array<int, 3> dims = { 8, 1, 1 };
    std::array<double, 3> cellsz = { 1.0, 1.0, 1.0 };
    Dune::CpGrid grid;
    grid.createCartesian(dims, cellsz);
    const auto& gv = grid.leafGridView();
    // num_chunks == 0
    BOOST_CHECK_THROW(EC(gv, 0), std::logic_error);
    // num_chunks == 1
    testCase(gv, 1, { 8 });
    // num_chunks < num_elements, num_elements = K * num_chunks
    testCase(gv, 2, { 4, 4 });
    // num_chunks < num_elements, general case
    testCase(gv, 3, { 2, 2, 4 });
    // num_chunks == num_elements
    testCase(gv, 8, { 1, 1, 1, 1, 1, 1, 1, 1 });
    // num_chunks > num_elements
    testCase(gv, 11, { 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0 });
}
