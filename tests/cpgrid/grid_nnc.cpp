/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#define NVERBOSE

#define BOOST_TEST_MODULE CpGridWithNNC

#include <boost/test/unit_test.hpp>

#include <dune/common/version.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/grid/CpGrid.hpp>
#include <vector>
#include <utility>

namespace std
{
    ostream& operator<<(ostream& os, const pair<int, int>& x)
    {
        os << "{ " << x.first << ", " << x.second << " }";
        return os;
    }
}

using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LeafGridView>;

struct Fixture
{
    Fixture()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(m_argc, m_argv);
        Opm::OpmLog::setupSimpleDefaultLogging();
    }

    Opm::Parser parser;

    void testCase(const std::string& filename,
                  const Opm::NNC& nnc,
                  const int ex_elemcount,
                  const int ex_intercount,
                  const int ex_bdycount,
                  const std::vector<std::pair<int, int>>& ex_nb,
                  const bool use_deck_porv = false)
    {
        Opm::Deck deck = parser.parseFile(filename);
        Opm::EclipseState es(deck);
        es.appendInputNNC(nnc.input());
        std::vector<double> porv;

        if (use_deck_porv) {
            porv = es.fieldProps().porv(true);
        }

        Dune::CpGrid grid;
        grid.processEclipseFormat(&es.getInputGrid(), &es, false, false, false);
        const auto& gv = grid.leafGridView();
        ElementMapper elmap(gv, Dune::mcmgElementLayout());
        int elemcount = 0;
        int intercount = 0;
        int bdycount = 0;
        std::vector<std::pair<int, int>> nb;
        for (const auto& elem : elements(gv)) {
            ++elemcount;
            for (const auto& inter : intersections(gv, elem)) {
                ++intercount;
                if (inter.boundary()) {
                    ++bdycount;
                } else {
                    // Internal connection, store in vector of pairs, if c1 < c2.
                    const int c1 = elmap.index(elem);
                    const int c2 = elmap.index(inter.outside());
                    if (c1 < c2) {
                        nb.push_back(std::make_pair(c1, c2));
                    }
                }
            }
        }
        BOOST_CHECK_EQUAL(elemcount, ex_elemcount);
        BOOST_CHECK_EQUAL(intercount, ex_intercount);
        BOOST_CHECK_EQUAL(bdycount, ex_bdycount);
        std::ranges::sort(nb);
        //BOOST_TEST(nb == ex_nb, boost::test_tools::per_element());
	BOOST_CHECK_EQUAL_COLLECTIONS(nb.begin(), nb.end(),
                                      ex_nb.begin(), ex_nb.end());
    }

    // A grid of 3x1x2 cells given by DX, DY, DZ and TOPS, the cells are 25 m high.
    // Cell numbers: i + 3*k.
    static std::string blockCentredDeck(const std::string& tops,
                                        const std::string& actnum,
                                        const std::string& extra)
    {
        return "RUNSPEC\n"
               "DIMENS\n 3 1 2 /\n"
               "GRID\n"
               + extra +
               "DX\n 6*100 /\n"
               "DY\n 6*100 /\n"
               "DZ\n 6*25 /\n"
               "BOX\n 1 3 1 1 1 1 /\n"
               "TOPS\n " + tops + " /\n"
               "ENDBOX\n"
               "ACTNUM\n " + actnum + " /\n"
               "PORO\n 6*0.15 /\n"
               "EDIT\n";
    }

    // The sorted connected cell pairs of the grid. With useState false the grid is made
    // without an EclipseState.
    std::vector<std::pair<int, int>> neighbours(const std::string& deckString, const bool useState)
    {
        Opm::Deck deck = parser.parseString(deckString);
        Opm::EclipseState es(deck);
        Dune::CpGrid grid;
        grid.processEclipseFormat(&es.getInputGrid(), useState ? &es : nullptr, false, false, false);
        const auto& gv = grid.leafGridView();
        ElementMapper elmap(gv, Dune::mcmgElementLayout());
        std::vector<std::pair<int, int>> nb;
        for (const auto& elem : elements(gv)) {
            for (const auto& inter : intersections(gv, elem)) {
                if (!inter.boundary()) {
                    const int c1 = elmap.index(elem);
                    const int c2 = elmap.index(inter.outside());
                    if (c1 < c2) {
                        nb.emplace_back(c1, c2);
                    }
                }
            }
        }
        std::ranges::sort(nb);
        return nb;
    }
};

BOOST_AUTO_TEST_SUITE(ConstructingWithNNC)


BOOST_FIXTURE_TEST_CASE(NoNNC, Fixture)
{
    Opm::NNC nnc;
    testCase("FIVE.DATA", nnc, 5, 30, 22, { {0,1}, {1,2}, {2,3}, {3,4} });
}

BOOST_FIXTURE_TEST_CASE(NNCAtExistingFace, Fixture)
{
    Opm::NNC nnc;
    // Should not change grid since cells are already neighbours.
    // Repeated NNCs do not change the grid further.
    nnc.addNNC(2, 3, 1.0);
    nnc.addNNC(3, 2, 1.0);
    nnc.addNNC(2, 3, 1.0);
    testCase("FIVE.DATA", nnc, 5, 30, 22, { {0,1}, {1,2}, {2,3}, {3,4} });
}

BOOST_FIXTURE_TEST_CASE(NNCAtNewFace, Fixture)
{
    Opm::NNC nnc;
    nnc.addNNC(2, 4, 1.0);
    testCase("FIVE.DATA", nnc, 5, 30 + 2, 22, { {0,1}, {1,2}, {2,3}, {2,4}, {3,4} });
}

BOOST_FIXTURE_TEST_CASE(NNCAtSeveralFaces, Fixture)
{
    Opm::NNC nnc;
    nnc.addNNC(2, 4, 1.0);   // new connection
    nnc.addNNC(3, 4, 1.0);
    nnc.addNNC(2, 1, 1.0);
    nnc.addNNC(1, 2, 1.0);
    nnc.addNNC(1, 4, 1.0);   // new connection
    nnc.addNNC(999, 2, 1.0); // invalid
    nnc.addNNC(-1, 2, 1.0);  // invalid
    testCase("FIVE.DATA", nnc, 5, 30 + 4, 22, { {0,1}, {1,2}, {1,4}, {2,3}, {2,4}, {3,4} });
}

BOOST_FIXTURE_TEST_CASE(ActnumNoNNC, Fixture)
{
    Opm::NNC nnc;
    testCase("FIVE_ACTNUM.DATA", nnc, 4, 24, 20, { {0,1}, {2,3} });
}

BOOST_FIXTURE_TEST_CASE(ActnumNNCAtExistingFace, Fixture)
{
    Opm::NNC nnc;
    // Should not change grid since cells are already neighbours.
    // Note that cells 3 and 4 in the NNC spec are logical cartesian (global) faces,
    // since global cell 2 is inactive, this becomes cells 2 and 3 in active numbering.
    nnc.addNNC(3, 4, 1.0);
    testCase("FIVE_ACTNUM.DATA", nnc, 4, 24, 20, { {0,1}, {2,3} });
}

BOOST_FIXTURE_TEST_CASE(ActnumNNCAtNewFace, Fixture)
{
    Opm::NNC nnc;
    nnc.addNNC(2, 4, 1.0); // inactive cell 2, no new connection
    nnc.addNNC(1, 4, 1.0); // new connection
    testCase("FIVE_ACTNUM.DATA", nnc, 4, 24 + 2, 20, { {0,1}, {1,3}, {2,3} });
}

BOOST_FIXTURE_TEST_CASE(ActnumNNCAtSeveralFaces, Fixture)
{
    Opm::NNC nnc;
    nnc.addNNC(2, 4, 1.0);   // inactive cell 2, no new connection
    nnc.addNNC(3, 4, 1.0);
    nnc.addNNC(2, 1, 1.0);
    nnc.addNNC(1, 2, 1.0);
    nnc.addNNC(1, 4, 1.0);   // new connection
    nnc.addNNC(1, 3, 1.0);   // new connection
    nnc.addNNC(999, 2, 1.0); // invalid
    nnc.addNNC(-1, 2, 1.0);  // invalid
    testCase("FIVE_ACTNUM.DATA", nnc, 4, 24 + 4, 20, { {0,1}, {1,2}, {1,3}, {2,3} });
}

BOOST_FIXTURE_TEST_CASE(NNCWithPINCH, Fixture)
{
    Opm::NNC nnc;
    testCase("FIVE_PINCH.DATA", nnc, 4, 24 + 1, 18 + 1, { {0,1}, {1,2}, {2,3} }, true);
}

BOOST_FIXTURE_TEST_CASE(NNCWithPINCHNOGAP, Fixture)
{
    Opm::NNC nnc;
    // Cell with cartindex 1 is isolated because of PINCH NOGAP
    testCase("FIVE_PINCH_NOGAP.DATA", nnc, 4, 4 * 6, 2 * (4 + 5) + 2, { {1,2}, {2,3} }, true);
}

BOOST_FIXTURE_TEST_CASE(NNCWithPINCHNOGAP2, Fixture)
{
    Opm::NNC nnc;
    testCase("FIVE_PINCH_NOGAP2.DATA", nnc, 3, 3*6, 3*6-2, { {1,2} }, true);
}

BOOST_FIXTURE_TEST_CASE(NNCWithPINCHNOGAP3, Fixture)
{
    Opm::NNC nnc;
    testCase("FIVE_PINCH_NOGAP3.DATA", nnc, 2, 2*6, 2*6, { }, true);
}

BOOST_FIXTURE_TEST_CASE(NNCWithPINCHAndMore, Fixture)
{
    Opm::NNC nnc;
    nnc.addNNC(2, 4, 1.0);   // new connection
    nnc.addNNC(3, 4, 1.0);
    // The next causes trouble, and I have disabled it because I do not know
    // what should be the correct behaviour for that case.
    // nnc.addNNC(0, 2, 1.0);   // connection added from pinchout already
    nnc.addNNC(999, 2, 1.0); // invalid
    nnc.addNNC(-1, 2, 1.0);  // invalid
    testCase("FIVE_PINCH.DATA", nnc, 4, 24 + 2 + 1, 18 + 1, { {0,1}, {1,2}, {1,3}, {2,3} }, true);
}


// OLDTRAN (the default for grids given by DX, DY, DZ and TOPS): the cells in a layer are
// connected also when they do not overlap. The columns are stepped by 30 m and the
// cells are 25 m high, so the cells of a layer do not overlap, and the cells of the
// second layer overlap the first layer of the neighbouring column.
BOOST_FIXTURE_TEST_CASE(OldTranConnectsLateralNeighbours, Fixture)
{
    const auto deck = blockCentredDeck("2500 2530 2560", "6*1", "");
    const std::vector<std::pair<int, int>> expected{
        {0,1}, {0,3}, {1,2}, {1,3}, {1,4}, {2,4}, {2,5}, {3,4}, {4,5} };
    for (const bool useState : {true, false}) {
        const auto nb = neighbours(deck, useState);
        // the grid-only construction has no state, and gets the same topology
        BOOST_CHECK_EQUAL_COLLECTIONS(nb.begin(), nb.end(), expected.begin(), expected.end());
    }
}

// NEWTRAN: only the cells that overlap are connected
BOOST_FIXTURE_TEST_CASE(NewTranDoesNotConnectLateralNeighbours, Fixture)
{
    const auto deck = blockCentredDeck("2500 2530 2560", "6*1", "NEWTRAN\n");
    const std::vector<std::pair<int, int>> expected{ {0,3}, {1,3}, {1,4}, {2,4}, {2,5} };
    const auto nb = neighbours(deck, true);
    BOOST_CHECK_EQUAL_COLLECTIONS(nb.begin(), nb.end(), expected.begin(), expected.end());
}

// Inactive cells are not connected
BOOST_FIXTURE_TEST_CASE(OldTranInactiveCells, Fixture)
{
    const auto deck = blockCentredDeck("2500 2530 2560", "1 0 1 1 1 1", "");
    // cell 1 is inactive, the cells are numbered 0, 1, 2, 3, 4 (cell 2 -> 1, ...)
    // active cells 0, 2, 3, 4, 5 become 0, 1, 2, 3, 4
    const std::vector<std::pair<int, int>> expected{
        {0,2}, {1,3}, {1,4}, {2,3}, {3,4} };
    const auto nb = neighbours(deck, true);
    BOOST_CHECK_EQUAL_COLLECTIONS(nb.begin(), nb.end(), expected.begin(), expected.end());
}

// Cells in a flat grid are connected by faces, no connection is added twice
BOOST_FIXTURE_TEST_CASE(OldTranFlatGridSameAsNewTran, Fixture)
{
    const auto oldtran = neighbours(blockCentredDeck("6*2500", "6*1", ""), true);
    const auto newtran = neighbours(blockCentredDeck("6*2500", "6*1", "NEWTRAN\n"), true);
    const std::vector<std::pair<int, int>> expected{
        {0,1}, {0,3}, {1,2}, {1,4}, {2,5}, {3,4}, {4,5} };
    BOOST_CHECK_EQUAL_COLLECTIONS(oldtran.begin(), oldtran.end(), expected.begin(), expected.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(newtran.begin(), newtran.end(), expected.begin(), expected.end());
}

BOOST_AUTO_TEST_SUITE_END()
