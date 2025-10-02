/*
  Copyright 2025 Equinor ASA.

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

#define BOOST_TEST_MODULE AutoRefineTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <tests/cpgrid/lgr/LgrChecks.hpp>

#include <array>
#include <string>
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

void checkGridAfterAutoRefinement(const Dune::CpGrid& grid,
                                  const std::array<int,3>& nxnynz)
{
    // Extract the refined level grids name, excluding level zero grid name ("GLOBAL").
    // Note: in this case there is only one refined level grid, storing the global
    // refinement.
    std::vector<std::string> lgrNames(grid.maxLevel());
    for (const auto& [name, level] : grid.getLgrNameToLevel()) {
        if (level==0) { // skip level zero grid name for the checks
            continue;
        }
        lgrNames[level-1] = name; // Shift the index since level zero has been removed.
    }
    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {nxnynz},
                           /* lgr_name_vec = */ lgrNames,
                           /* gridHasBeenGlobalRefined = */ true);
}

BOOST_AUTO_TEST_CASE(evenRefinementFactorThrows)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    // nxnynz represents the refinement factors in x-,y-,and z-direction.
    BOOST_CHECK_THROW(grid.autoRefine(/* nxnynz = */ {4,3,5}), std::invalid_argument);
    BOOST_CHECK_THROW(grid.autoRefine(/* nxnynz = */ {3,4,5}), std::invalid_argument);
    BOOST_CHECK_THROW(grid.autoRefine(/* nxnynz = */ {3,5,4}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(nonPositiveRefinementFactorThrows)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    // nxnynz represents the refinement factors in x-,y-,and z-direction.
    BOOST_CHECK_THROW(grid.autoRefine(/* nxnynz = */ {0,3,5}), std::invalid_argument);
    BOOST_CHECK_THROW(grid.autoRefine(/* nxnynz = */ {3,-1,5}), std::invalid_argument);
    BOOST_CHECK_THROW(grid.autoRefine(/* nxnynz = */ {3,5,-3}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(autoRefine)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    if (grid.comm().size()>1) { // distribute level zero grid, in parallel
        grid.loadBalance();
    }

    // nxnynz represents the refinement factors in x-,y-,and z-direction.
    grid.autoRefine(/* nxnynz = */ {3,5,7});

    checkGridAfterAutoRefinement(grid, {3,5,7});
}

BOOST_AUTO_TEST_CASE(callGlobalRefineAfterAutoRefine_serial) {

    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    if (grid.comm().size() == 1 ) { // serial

        // nxnynz represents the refinement factors in x-,y-,and z-direction.
        grid.autoRefine(/* nxnynz = */ {3,3,1});

        checkGridAfterAutoRefinement(grid, /* nxnynz = */ {3,3,1});

        grid.globalRefine(2);

        Opm::checkGridWithLgrs(grid,
                               /* cells_per_dim_vec = */ {{2,2,2}, {2,2,2}},
                               /* lgr_name_vec = */ {"GR2", "GR3"},
                               /* gridHasBeenGlobalRefined = */ true,
                               /* preRefineMaxLevel = */ 1,
                               /* isNested = */ true);
    }
}


BOOST_AUTO_TEST_CASE(callAdaptAfterAutoRefine_serial) {

    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    if (grid.comm().size() == 1 ) { // serial

        // nxnynz represents the refinement factors in x-,y-,and z-direction.
        grid.autoRefine(/* nxnynz = */ {3,3,1});

        checkGridAfterAutoRefinement(grid, /* nxnynz = */ {3,3,1});

        for (const auto& element : Dune::elements(grid.leafGridView())) {
            grid.mark(1, element);
        }
        grid.preAdapt();
        grid.adapt();
        grid.postAdapt();

         Opm::checkGridWithLgrs(grid,
                               /* cells_per_dim_vec = */ {{2,2,2}},
                               /* lgr_name_vec = */ {"GR2"},
                               /* gridHasBeenGlobalRefined = */ true,
                               /* preRefineMaxLevel = */ 1,
                               /* isNested = */ true);
    }
}

BOOST_AUTO_TEST_CASE(callAdaptNotAllElementsMarkedAfterAutoRefine_serial) {

    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    if (grid.comm().size() == 1 ) { // serial

        // nxnynz represents the refinement factors in x-,y-,and z-direction.
        grid.autoRefine(/* nxnynz = */ {3,3,1});

        checkGridAfterAutoRefinement(grid, /* nxnynz = */ {3,3,1});

        const int maxCellId =  grid.currentData().back()->globalIdSet().getMaxCodimGlobalId<0>();

        for (const auto& element : Dune::elements(grid.leafGridView())) {
            const int id = grid.globalIdSet().id(element);
            // Mark two elements
            if ( (id == maxCellId - 2) || (id == maxCellId -1)) {
                grid.mark(1, element);
            }
        }
        grid.preAdapt();
        grid.adapt();
        grid.postAdapt();

        Opm::checkGridWithLgrs(grid,
                               /* cells_per_dim_vec = */ {{2,2,2}},
                               /* lgr_name_vec = */ {"LGR2"},
                               /* gridHasBeenGlobalRefined = */ false,
                               /* preRefineMaxLevel = */ 1,
                               /* isNested = */ true);
    }
}
