/*
  Copyright 2022 Equinor ASA.

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

#define BOOST_TEST_MODULE GeometryTests
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif
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

    static int rank()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        return Dune::MPIHelper::instance(m_argc, m_argv).rank();
    }
};

BOOST_GLOBAL_FIXTURE(Fixture);
void createAndTestShiftCartGrid(std::array<int, 3> shift)
{
    Dune::CpGrid grid;
    std::array<double,3> cell_size = {1.0, 0.1, .5};
    std::array<int,3> cell_dims = {5, 8, 3};
    grid.createCartesian(cell_dims, cell_size, shift);

    for(const auto& cell: elements(grid.leafGridView()))
    {
        std::array<int, 3> ijk, low{0, 0, 0}, high{cell_dims[0] - 1, cell_dims[1] - 1, cell_dims[2] -1};
        grid.getIJK(cell.index(), ijk);
        if (ijk == low)
        {
            Dune::FieldVector<double, 3> center =
                { (shift[0] + .5) * cell_size[0],
                  (shift[1] + .5) * cell_size[1],
                  (shift[2] + .5) * cell_size[2]};
            BOOST_TEST( center == cell.geometry().center(),
                        boost::test_tools::per_element() );
            continue;
        }
        if (ijk == high)
        {
            Dune::FieldVector<double, 3> center =
                { (shift[0] + cell_dims[0] - .5) * cell_size[0],
                  (shift[1] + cell_dims[1] - .5) * cell_size[1],
                  (shift[2] + cell_dims[2] - .5) * cell_size[2]};
            BOOST_TEST( center == cell.geometry().center(),
                        boost::test_tools::per_element() );
            continue;
        }

    }
}


BOOST_AUTO_TEST_CASE(posShiftedCartGrid, *boost::unit_test::tolerance(1e-12))
{
    createAndTestShiftCartGrid({10, 4 , 9});
}

BOOST_AUTO_TEST_CASE(negShiftedCartGrid, *boost::unit_test::tolerance(1e-12))
{
    createAndTestShiftCartGrid({-10, -4 , -9});
}
BOOST_AUTO_TEST_CASE(nonShiftedCartGrid, *boost::unit_test::tolerance(1e-12))
{
    createAndTestShiftCartGrid({0, 0, 0});
}
