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

#define BOOST_TEST_MODULE LgrBackendSwitchTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

struct MPIFixture
{
    MPIFixture()
    {
        auto& argc = boost::unit_test::framework::master_test_suite().argc;
        auto& argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(argc, argv);
    }
};
BOOST_GLOBAL_FIXTURE(MPIFixture);

namespace
{

using Opm::Refinement::Backend;
using Opm::Refinement::BlockRefinement;

void makeGrid(Dune::CpGrid& grid)
{
    grid.createCartesian({{4, 3, 3}}, {{1.0, 1.0, 1.0}});
}

BlockRefinement box(const std::string& name, std::array<int,3> start, std::array<int,3> end)
{
    BlockRefinement request;
    request.name = name;
    request.cellsPerDim = {2, 2, 2};
    request.startIJK = start;
    request.endIJK = end;
    return request;
}

struct RecordingBuilder : Opm::Refinement::Builder
{
    void build(Dune::CpGrid&, const std::vector<BlockRefinement>& requests) override
    {
        recorded = requests;
    }
    std::vector<BlockRefinement> recorded{};
};

void processDeck(Dune::CpGrid& grid)
{
    const std::string deckString = R"(RUNSPEC
DIMENS
 4 3 3 /
GRID
CARFIN
'LGR1' 2 3 1 1 1 1 4 2 2 /
ENDFIN
DX
 36*1 /
DY
 36*1 /
DZ
 36*1 /
TOPS
 12*0 /
PORO
 36*0.2 /
)";
    const auto deck = Opm::Parser{}.parseString(deckString);
    Opm::EclipseState state(deck);
    auto eclGrid = state.getInputGrid();
    grid.processEclipseFormat(&eclGrid, &state, false, false, false);
}
} // anonymous namespace

// Under the default backend the request form refines exactly as upstream's vector form.
BOOST_AUTO_TEST_CASE(trilinearRequestsMatchVectors)
{
    Dune::CpGrid byVectors, byRequests;
    makeGrid(byVectors);
    makeGrid(byRequests);
    BOOST_CHECK(byRequests.lgrBackend() == Backend::Trilinear);

    byVectors.addLgrsUpdateLeafView({{2,2,2}}, {{1,1,1}}, {{3,2,2}}, {"LGR1"});
    byRequests.addLgrsUpdateLeafView(std::vector<BlockRefinement>{ box("LGR1", {1,1,1}, {3,2,2}) });

    BOOST_CHECK_EQUAL(byRequests.maxLevel(), 1);
    BOOST_CHECK_EQUAL(byRequests.maxLevel(), byVectors.maxLevel());
    BOOST_CHECK_EQUAL(byRequests.size(0), byVectors.size(0));
    BOOST_CHECK_EQUAL(byRequests.size(1), byVectors.size(1));
    BOOST_CHECK_EQUAL(byRequests.size(3), byVectors.size(3));
    BOOST_CHECK(byRequests.getLgrNameToLevel() == byVectors.getLgrNameToLevel());
}

// What only the conforming builder can do is refused, naming the backend.
BOOST_AUTO_TEST_CASE(trilinearRefusesConformingOnlyRequests)
{
    auto graded = box("LGR1", {1,1,1}, {2,2,2});
    graded.subdivision[0] = { {0, 0, 0}, {0.0, 0.2, 0.5}, {0.2, 0.5, 1.0} };
    auto pillars = box("LGR1", {1,1,1}, {2,2,2});
    pillars.pillarsFromBoxLayer = true;
    auto minpv = box("LGR1", {1,1,1}, {2,2,2});
    minpv.minpvRemoved.assign(8, 0);

    for (const auto& request : { graded, pillars, minpv }) {
        Dune::CpGrid grid;
        makeGrid(grid);
        BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView(std::vector<BlockRefinement>{ request }),
                          std::invalid_argument);
        BOOST_CHECK_EQUAL(grid.maxLevel(), 0);
    }
}

// The two backends never meet on one grid.
BOOST_AUTO_TEST_CASE(backendFixedOnceRefined)
{
    Dune::CpGrid grid;
    makeGrid(grid);
    grid.addLgrsUpdateLeafView({{2,2,2}}, {{1,1,1}}, {{2,2,2}}, {"LGR1"});
    BOOST_CHECK_THROW(grid.setLgrBackend(Backend::Conforming), std::logic_error);
    BOOST_CHECK_NO_THROW(grid.setLgrBackend(Backend::Trilinear));
}

// Under the Conforming backend both entry points reach the grid's builder,
// after validation.
BOOST_AUTO_TEST_CASE(conformingForwardsToItsBuilder)
{
    Dune::CpGrid grid;
    makeGrid(grid);
    grid.setLgrBackend(Backend::Conforming);
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView({{2,2,2}}, {{1,1,1}}, {{2,2,2}}, {"LGR1"}),
                      std::logic_error);

    auto builder = std::make_shared<RecordingBuilder>();
    grid.setRefinementBuilder(builder);
    grid.addLgrsUpdateLeafView({{2,2,2}, {3,3,3}}, {{0,0,0}, {2,1,1}}, {{1,1,1}, {4,3,3}}, {"A", "B"});
    BOOST_REQUIRE_EQUAL(builder->recorded.size(), 2u);
    BOOST_CHECK_EQUAL(builder->recorded[0].name, "A");
    BOOST_CHECK_EQUAL(builder->recorded[1].name, "B");
    BOOST_CHECK(builder->recorded[1].cellsPerDim == (std::array<int,3>{3, 3, 3}));
    BOOST_CHECK_EQUAL(builder->recorded[1].parentGridName, "GLOBAL");

    builder->recorded.clear();
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView(std::vector<BlockRefinement>{
                          box("C", {0,0,0}, {2,2,2}), box("D", {1,1,1}, {3,3,3}) }),
                      std::invalid_argument);
    BOOST_CHECK(builder->recorded.empty());
}

// A Conforming grid built from a deck with LGRs refines without an explicit builder.
BOOST_AUTO_TEST_CASE(conformingUsesRetainedDeckInput)
{
    Dune::CpGrid grid;
    grid.setLgrBackend(Backend::Conforming);
    processDeck(grid);
    grid.addLgrsUpdateLeafView({{2, 1, 1}}, {{1, 0, 0}}, {{3, 1, 1}}, {"LGR1"});
    BOOST_CHECK_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_EQUAL(grid.size(0), 36 - 2 + 2*2);
}

// Input is retained only when the backend is chosen before processEclipseFormat().
BOOST_AUTO_TEST_CASE(conformingChosenAfterProcessingThrows)
{
    Dune::CpGrid grid;
    processDeck(grid);
    grid.setLgrBackend(Backend::Conforming);
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView({{2, 1, 1}}, {{1, 0, 0}}, {{3, 1, 1}}, {"LGR1"}),
                      std::logic_error);
}
