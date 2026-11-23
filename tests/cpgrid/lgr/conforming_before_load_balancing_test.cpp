/*
  Copyright 2026 Equinor ASA.

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

#define BOOST_TEST_MODULE ConformingBeforeLoadBalancingTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/refinement/conforming/ConformingBlockBuilder.hpp>
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>

#include <tests/cpgrid/lgr/ConformingTestHelpers.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <array>
#include <exception>
#include <memory>
#include <string>

// The conforming refinement of a grid that is not yet load balanced. Run on one rank and
// on several: before load balancing the whole grid is on rank 0, and the other ranks hold
// no cells.

using namespace Opm::ConformingTest;

BOOST_GLOBAL_FIXTURE(MPIFixture);

// A box with an inactive host cell is refined on rank 0 like any other box.
BOOST_AUTO_TEST_CASE(boxWithInactiveHostRefinesBeforeLoadBalancing)
{
    // 12x12x4 unit cells; cell (2,2,2) (I,J,K from 1), inside the box, is inactive.
    const std::string deckString = R"(RUNSPEC
DIMENS
 12 12 4 /
GRID
DX
 576*1 /
DY
 576*1 /
DZ
 576*1 /
BOX
 1 12 1 12 1 1 /
TOPS
 144*0 /
ENDBOX
ACTNUM
 157*1 0 418*1 /
PORO
 576*0.2 /
)";
    const auto deck = Opm::Parser{}.parseString(deckString);
    Opm::EclipseState state(deck);
    auto eclGrid = state.getInputGrid();

    Dune::CpGrid grid;
    grid.setLgrBackend(Opm::Refinement::Backend::Conforming);
    grid.processEclipseFormat(&eclGrid, &state, false, false, false);
    const auto comm = grid.comm();
    BOOST_REQUIRE_EQUAL(grid.size(0), (comm.rank() == 0) ? 575 : 0);

    // Every rank refines from the same corner-point input.
    grid.setRefinementBuilder(std::make_shared<Opm::Refinement::ConformingBlockBuilder>(
        std::array<int,3>{12, 12, 4}, eclGrid.getCOORD(), eclGrid.getZCORN(), eclGrid.getACTNUM()));

    // The box is I = 1..6 (all J and K). A rank that refuses it throws; the ranks agree
    // on the outcome before the grid is checked.
    int refused = 0;
    try {
        grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2, 2, 2}},
                                   /* startIJK_vec = */      {{0, 0, 0}},
                                   /* endIJK_vec = */        {{6, 12, 4}},
                                   /* lgr_name_vec = */      {"LGR1"});
    }
    catch (const std::exception& e) {
        refused = 1;
        BOOST_ERROR("Rank " << comm.rank() << " refused the box: " << e.what());
    }
    BOOST_REQUIRE_EQUAL(comm.max(refused), 0);

    // Rank 0 holds the refined grid: the 288 cells outside the box, and 8 children
    // of each of the 287 active cells inside it.
    BOOST_CHECK_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_EQUAL(grid.size(0), (comm.rank() == 0) ? 288 + 287*8 : 0);
    BOOST_CHECK_CLOSE(comm.sum(totalVolume(grid)), 575.0, 1e-10);
    // The refined grid keeps level zero's communicator, so load balancing sees every rank.
    BOOST_CHECK_EQUAL(grid.comm().size(), comm.size());
}
