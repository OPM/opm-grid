/*
  Copyright 2026 NORCE Norwegian Research Centre AS.

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

// A dual-continuum grid cannot be processed on more than one process, and every
// process has to refuse it, also those that are not given the input grid.

#include <config.h>

#define BOOST_TEST_MODULE DualContinuumRefusalTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <stdexcept>

namespace {

struct Fixture
{
    Fixture()
    {
        int argc = boost::unit_test::framework::master_test_suite().argc;
        char** argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(argc, argv);
    }
};

BOOST_GLOBAL_FIXTURE(Fixture);

} // anonymous namespace

BOOST_AUTO_TEST_CASE(EveryProcessRefusesADualContinuumGrid)
{
    Dune::CpGrid grid;
    if (grid.comm().size() == 1) {
        return; // Nothing to refuse on a single process.
    }

    // Two matrix cells and two fracture cells.
    const auto deck = Opm::Parser{}.parseString(
        "RUNSPEC\nOIL\nWATER\nDIMENS\n 2 1 2 /\nDUALPORO\nGRID\n"
        "DX\n 4*100 /\nDY\n 4*100 /\nDZ\n 4*10 /\nTOPS\n 2*2000 /\n"
        "PORO\n 4*0.2 /\nPERMX\n 4*1 /\nSIGMA\n 0.1 /\n");
    Opm::EclipseState eclState(deck);

    // The input grid is only given on the first process, the EclipseState on all.
    const bool root = grid.comm().rank() == 0;
    BOOST_CHECK_THROW(grid.processEclipseFormat(root ? &eclState.getInputGrid() : nullptr, &eclState,
                                                false, false, false),
                      std::invalid_argument);
}
