/*
  Copyright 2023, 2025 Equinor ASA.

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
#include "config.h"

#define BOOST_TEST_MODULE DisjointPatchesTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/cpgrid/LgrHelpers.hpp>

#include <array>
#include <vector>

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

BOOST_GLOBAL_FIXTURE(Fixture);

BOOST_AUTO_TEST_CASE(lgrs_disjointPatches)
{
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,2}, {3,2,2}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,1,1}, {1,1,3}, {4,3,3}};

    BOOST_CHECK( Opm::disjointPatches(startIJK_vec, endIJK_vec) );
}

BOOST_AUTO_TEST_CASE(patches_share_corner)
{
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {1,1,1}, {3,2,2}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,1}, {2,2,2}, {4,3,3}};

    BOOST_CHECK( !Opm::disjointPatches(startIJK_vec, endIJK_vec) );
}

BOOST_AUTO_TEST_CASE(patches_share_corners)
{
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {2,0,1}, {3,2,2}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,1,1}, {3,1,2}, {4,3,3}};

    BOOST_CHECK( !Opm::disjointPatches(startIJK_vec, endIJK_vec) );
}

BOOST_AUTO_TEST_CASE(pathces_share_face)
{
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {2,0,0}, {3,2,2}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,1,1}, {3,1,1}, {4,3,3}};

    BOOST_CHECK( !Opm::disjointPatches(startIJK_vec, endIJK_vec) );
}

BOOST_AUTO_TEST_CASE(invalid_argument_sizes)
{

    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,1}, {1,1,2}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,2,1}, {3,2,1}};

    BOOST_CHECK_THROW(Opm::disjointPatches(startIJK_vec, endIJK_vec), std::logic_error);
}
