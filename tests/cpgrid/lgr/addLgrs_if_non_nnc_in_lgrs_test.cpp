/*
  Copyright 2023 Equinor ASA.

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

#define BOOST_TEST_MODULE CpGridWithNNC

#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/grid/CpGrid.hpp>

#include <tests/cpgrid/lgr/LgrChecks.hpp>

#include <stdexcept>
#include <string>

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

void testCase(Dune::CpGrid& grid,
              const std::string& deckString,
              const Opm::NNC& nnc)
{
    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);
    Opm::EclipseState es(deck);
    es.appendInputNNC(nnc.input());

    grid.processEclipseFormat(&es.getInputGrid(), &es, false, false, false);
}

BOOST_GLOBAL_FIXTURE(Fixture);

const std::string deckString =
    R"( RUNSPEC
        DIMENS
        1  1  5 /
        GRID
        COORD
        0 0 0
        0 0 1
        1 0 0
        1 0 1
        0 1 0
        0 1 1
        1 1 0
        1 1 1
        /
        ZCORN
        4*0
        8*1
        8*2
        8*3
        8*4
        4*5
        /
        ACTNUM
        5*1
        /
        PORO
        5*0.15
        /)";

BOOST_AUTO_TEST_CASE(addLgrsWithNNCatAnLgrThrows)
{
    Opm::NNC nnc;
    nnc.addNNC(2, 4, 1.0); // connect cell 2 (does not belong to any LGR) and cell 4 (belongs to LGR2)
    // LGR1 parent cell indices = {0,1}.
    // LGR2 parent cell indices = {4}.
    Dune::CpGrid grid;
    testCase(grid, deckString, nnc);

    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}, {3,3,3}},
                                                 /* startIJK_vec = */ {{0,0,0}, {0,0,4}},
                                                 /* endIJK_vec = */ {{2,1,1}, {1,1,5}},
                                                 /* lgr_name_vec = */ {"LGR1", "LGR2"}), std::logic_error);
}

BOOST_AUTO_TEST_CASE(addLgrsWithNNCAtSeveralLgrsThrows)
{
    Opm::NNC nnc;
    nnc.addNNC(0, 1, 1.0); // connect cell 0 and cell 1 (both belong to LGR1)
    nnc.addNNC(2, 4, 1.0); // connect cell 2 (belongs to LGR2) and cell 4 (belongs to LGR3)
    // LGR1 parent cell index = {0}.
    // LGR2 parent cell index = {2}.
    // LGR3 parent cell index = {4}.
    Dune::CpGrid grid;
    testCase(grid, deckString, nnc);

    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}, {3,3,3}, {4,4,4}},
                                                 /* startIJK_vec = */ {{0,0,0}, {0,0,2}, {0,0,4}},
                                                 /* endIJK_vec = */  {{1,1,1}, {1,1,3}, {1,1,5}},
                                                 /* lgr_name_vec = */ {"LGR1", "LGR2", "LGR3"}), std::logic_error);
}

BOOST_AUTO_TEST_CASE(LgrWithNNC_and_lgrsWithoutNNC_addLgrsThrows)
{
    Opm::NNC nnc;
    nnc.addNNC(0, 2, 1.0); // connect cell 0 and cell 2 (both belong to LGR1). LGR2 does not have NNCs.
    // LGR1 parent cell indices = {0,1,2}.
    // LGR2 parent cell index = {4}.
    Dune::CpGrid grid;
    testCase(grid, deckString, nnc);

    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}, {4,4,4}},
                                                 /* startIJK_vec = */ {{0,0,0}, {0,0,4}},
                                                 /* endIJK_vec = */  {{1,1,3}, {1,1,5}},
                                                 /* lgr_name_vec = */ {"LGR1", "LGR2"}), std::logic_error);
}

BOOST_AUTO_TEST_CASE(addLgrsIfNNCoutsideAllLgrsIsSupported)
{
    Opm::NNC nnc;
    nnc.addNNC(1, 3, 1.0); // connect cell 1 and cell 3 (both do NOT belong to any LGR)
    // LGR1 parent cell index = {0}.
    // LGR2 parent cell index = {2}.
    // LGR3 parent cell index = {4}.
    Dune::CpGrid grid;
    testCase(grid, deckString, nnc);

    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}, {3,3,3}, {4,4,4}},
                               /* startIJK_vec = */ {{0,0,0}, {0,0,2}, {0,0,4}},
                               /* endIJK_vec = */  {{1,1,1}, {1,1,3}, {1,1,5}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2", "LGR3"});

    Opm::checkGridBasicHiearchyInfo(grid, /* cells_per_dim_vec = */ {{2,2,2}, {3,3,3}, {4,4,4}});
}
