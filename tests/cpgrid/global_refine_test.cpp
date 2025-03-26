//===========================================================================
/*
  Copyright 2024 Equinor ASA.

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

#define BOOST_TEST_MODULE GlobalRefineTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <tests/cpgrid/LgrChecks.hpp>

#include <numeric>
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


BOOST_AUTO_TEST_CASE(globalRefineWithParamZeroDoesNothingToTheGrid)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    grid.globalRefine(0);

    BOOST_CHECK_EQUAL( grid.maxLevel(), 0);
    BOOST_CHECK_EQUAL( grid.leafGridView().size(0), grid.levelGridView(0).size(0));
}

BOOST_AUTO_TEST_CASE(globalRefineWithNegativeParamThrows)
{
    Dune::CpGrid grid;
    BOOST_CHECK_THROW(grid.globalRefine(-5), std::logic_error);
}

BOOST_AUTO_TEST_CASE(globalRefineAgridWithCoarseAndRefinedCellsThrows)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}},
                               /* startIJK_vec = */ {{2,0,0}},
                               /* endIJK_vec = */ {{4,1,1}},
                               /* lgr_name_vec = */ {"LGR1"});

    BOOST_CHECK_THROW(grid.globalRefine(1), std::logic_error);
}

BOOST_AUTO_TEST_CASE(globalRefineWithPositiveParamIsSupported)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {2,2,1}, /* cell_sizes = */ {8.0, 8.0, 4.0});
    grid.globalRefine(3);
    // Note: default subdivisions per cell is 2x2x2 in globalRefine().
    // Starting grid dimensions {2,2,1} -1stGR-> {4,4,2} -2ndGR-> {8,8,4} -3rdGR-> {16,16,8}.

    // Create other grid for comparison, equivalent to "grid.globalRefine(2)". Mark all elements
    // and adapt.
    Dune::CpGrid equivalent_grid;
    equivalent_grid.createCartesian(/* grid_dim = */ {8,8,4}, /* cell_sizes = */ {2.0, 2.0, 1.0});
    std::vector<int> markedCells(256); // 256 = 8x8x4
    std::iota(markedCells.begin(), markedCells.end(), 0);
    Opm::adaptGrid(equivalent_grid, markedCells);

    // We set isBlockShape as false, even though global-refinement implies refinement of a block of cells.
    Opm::compareGrids(grid,
                      equivalent_grid,
                      /* lgrsHaveBlockShape = */ false,
                      /* gridHasBeenGlobalRefined = */ true);
}
