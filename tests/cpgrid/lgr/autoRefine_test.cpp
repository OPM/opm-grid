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

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/AutoRefinement.hpp>
#include <opm/input/eclipse/EclipseState/Grid/AutoRefManager.hpp>
#include <opm/input/eclipse/EclipseState/Grid/readKeywordAutoRef.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/grid/CpGrid.hpp>
#include <tests/cpgrid/lgr/LgrChecks.hpp>

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
                           /* cells_per_dim_vec = */ {{3,5,7}},
                           /* lgr_name_vec = */ lgrNames,
                           /* gridHasBeenGlobalRefined = */ true);
}

BOOST_AUTO_TEST_CASE(readAutoref) {

    const std::string deck_string = R"(
RUNSPEC
AUTOREF
3 3 1 0. /
DIMENS
 10 10 3 /
GRID
DX
300*1 /
DY
300*1 /
DZ
300*1 /
TOPS
100*1 /
PORO
300*0.15 /
)";

    Opm::Parser parser;
    Opm::Deck deck = parser.parseString(deck_string);

    const auto& autoref_keyword = deck["AUTOREF"][0];

    Opm::AutoRefManager autoRefManager{};

    Opm::readKeywordAutoRef(autoref_keyword.getRecord(0), autoRefManager);
    const auto autoRef = autoRefManager.getAutoRef();

    BOOST_CHECK_EQUAL( autoRef.NX(), 3);
    BOOST_CHECK_EQUAL( autoRef.NY(), 3);
    BOOST_CHECK_EQUAL( autoRef.NZ(), 1);
    BOOST_CHECK_EQUAL( autoRef.OPTION_TRANS_MULT(), 0.);

    Opm::EclipseState ecl_state(deck);
    Opm::EclipseGrid ecl_grid = ecl_state.getInputGrid();

    Dune::CpGrid grid;
    grid.processEclipseFormat(&ecl_grid, &ecl_state, false, false, false);

    if (grid.comm().size()>1) { // distribute level zero grid, in parallel
        grid.loadBalance();
    }

    // nxnynz represents the refinement factors in x-,y-,and z-direction.
    grid.autoRefine(/* nxnynz = */ { autoRef.NX(), autoRef.NY(), autoRef.NZ()});

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
                           /* cells_per_dim_vec = */ {{autoRef.NX(), autoRef.NY(), autoRef.NZ()}},
                           /* lgr_name_vec = */ lgrNames,
                           /* gridHasBeenGlobalRefined = */ true);
}
