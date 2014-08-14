/*
  Copyright 2014 Andreas Lauser

  This file is part of The Open Porous Media project  (OPM).

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

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE IntersectionMapperTests
#include <boost/test/unit_test.hpp>

#include <dune/grid/CpGrid.hpp>
#include <dune/grid/cpgrid/IntersectionMapper.hpp>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>

BOOST_AUTO_TEST_CASE(IntersectionMapper)
{

    int argc = boost::unit_test::framework::master_test_suite().argc;
    char** argv = boost::unit_test::framework::master_test_suite().argv;
    Dune::MPIHelper::instance(argc, argv);

    const char* deckString =
        "RUNSPEC\n"
        "DIMENS\n"
        "2 3 4 /\n"
        "GRID\n"
        "DXV\n"
        "1 1 /\n"
        "DYV\n"
        "1 1 1 /\n"
        "DZV\n"
        "1 1 1 1/\n"
        "DEPTHZ\n"
        "12*0.0/\n"
        "SCHEDULE\n";

    Opm::ParserPtr parser(new Opm::Parser);
    Opm::DeckConstPtr deck = parser->parseString(deckString);
    Opm::EclipseStateConstPtr eclipseState(new Opm::EclipseState(deck));

    Dune::CpGrid grid;
    grid.processEclipseFormat(eclipseState->getEclipseGrid(),
                              /*zTolerance=*/0.0,
                              /*isPeriodic=*/false);

    typedef Dune::CpGrid::LeafGridView GridView;
#if DUNE_VERSION_NEWER(DUNE_GRID, 2,3)
    const auto gridView = grid.leafGridView();
#else
    const auto gridView = grid.leafView();
#endif

    typedef Dune::IntersectionMapper<GridView> IntersectionMapper;
    IntersectionMapper isMap(gridView);
    const auto& cartCells = grid.logicalCartesianSize();
    int numIs = 0;
    for (int i = 0; i < 3; ++i) {
        int n = 0;
        if (i == 0)
            n = cartCells[1]*cartCells[2];
        else if (i == 1)
            n = cartCells[0]*cartCells[2];
        else if (i == 2)
            n = cartCells[0]*cartCells[1];

        numIs += (cartCells[i] + 1)*n;
    }

    BOOST_CHECK_EQUAL(numIs, isMap.size());

    std::vector<bool> isSeen(numIs, false);

    auto elemIt = gridView.begin<0>();
    const auto& elemEndIt = gridView.end<0>();
    for (; elemIt != elemEndIt; ++elemIt) {
        auto isIt = gridView.ibegin(*elemIt);
        const auto& isEndIt = gridView.iend(*elemIt);
        for (; isIt != isEndIt; ++isIt) {
            BOOST_CHECK(0 <= isMap.map(*isIt));
            BOOST_CHECK(isMap.map(*isIt) < numIs);
            isSeen[isMap.map(*isIt)] = true;
        }
    }

    for (int isIdx = 0; isIdx < numIs; ++ isIdx) {
        BOOST_CHECK(isSeen[isIdx]);
    }
}
