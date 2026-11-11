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

#define BOOST_TEST_MODULE FaultedBoundaryTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/cpgrid/LgrHelpers.hpp>
#include <opm/grid/cpgrid/refinement/conforming/FaultedBoundaryFaces.hpp>

#include <tests/cpgrid/lgr/ConformingTestHelpers.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <vector>

using namespace Opm::ConformingTest;

BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_AUTO_TEST_CASE(faultedIBoundarySplitsConnections)
{
    // Fault between i=1 and i=2 (throw 0.6). Box i2-3, refined 2x2x2; its
    // i-minus boundary sits on the fault. Boundary cells that straddle the
    // throw must connect to two coarse neighbours. Indices are zero-based.
    //
    // Plan view (5 x 2 x 2; box = i2-3, all j, all k, drawn as its 4 x 4
    // refined cells). The tested face is '#':
    //        +-----+-----#--+--+--+--+-----+
    //        |     |     #  |  |  |  |     |
    //    j=1 |     |     #--+--+--+--+     |
    //        |     |     #  |  |  |  |     |
    //        +-----+-----#--+--+--+--+-----+
    //        |     |     #  |  |  |  |     |
    //    j=0 |     |     #--+--+--+--+     |
    //        |     |     #  |  |  |  |     |
    //        +-----+-----#--+--+--+--+-----+
    //                    ^ fault at x=2 (axis=0, side=-1)
    //
    // Section through the fault at any j (depth down). The coarse column i=1
    // is refined like the box (2 per cell) to compute the faces, so each box
    // face is also cut at depths 1.0 and 3.0:
    //
    //   depth   coarse i=1    box, i-layer 0
    //    0.0   +----------+
    //    0.6   |   k=0    |  +----------+
    //    1.0   |- - - - - |  |   k=0    |  -> (i=1,k=0) twice
    //    1.6   |          |  +----------+
    //    2.0   +----------+  |   k=1    |  -> (i=1,k=0) and (i=1,k=1): split
    //    2.6   |          |  +----------+
    //    3.0   |- - - - - |  |   k=2    |  -> (i=1,k=1) twice
    //    3.6   |   k=1    |  +----------+
    //    4.0   +----------+  |   k=3    |  -> (i=1,k=1) and the domain (-1)
    //    4.6                 +----------+
    //
    // Each box face is 0.5 wide in j and faces one coarse j, so the 4 (j) x 4 (k)
    // boundary cells get 2 connections each: 32 in all.
    auto depth = [](int i_, int, int k_) {
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;
    };
    auto g = makeGrid({5, 2, 2}, depth);

    const auto conns = Opm::Refinement::boxBoundaryConnections(
        g.dims, g.coord.data(), g.zcorn.data(), /*actnum=*/nullptr,
        Opm::Refinement::BlockRefinement{ .name = "BOX", .cellsPerDim = {2,2,2},
                                          .startIJK = {2,0,0}, .endIJK = {4,2,2} },
        /*axis=*/0, /*side=*/-1, /*edgeConformal=*/true);

    BOOST_REQUIRE(!conns.empty());
    // Each of the 16 refined cells at the i-minus box boundary gets 2 connections.
    BOOST_CHECK_EQUAL(conns.size(), 32u);

    std::map<std::array<int,3>, std::set<int>> neighborsOf; // box cell -> coarse carts
    std::map<std::array<int,3>, int> connectionsOf;         // box cell -> its connections
    int domainBoundaryParts = 0;
    for (const auto& c : conns) {
        BOOST_CHECK_EQUAL(c.boxCell[0], 0);               // all faces lie in the box's i-minus layer (plane x = 2)
        BOOST_CHECK_EQUAL(c.faceNodes.size(), 4u);        // a quadrilateral: no hanging nodes here
        ++connectionsOf[c.boxCell];
        if (c.coarseNeighborCart < 0) {
            ++domainBoundaryParts;                        // below column i=1 (faces the domain)
            BOOST_CHECK_EQUAL(c.boxCell[2], 3);           // only the box's bottom layer reaches there
        }
        else {
            // grid dims {5,2,2}: coarse cart = i + 5*j + 10*k
            BOOST_CHECK_EQUAL(Opm::Lgr::getIJK(c.coarseNeighborCart, g.dims)[0], 1); // coarse neighbour at i=1
            neighborsOf[c.boxCell].insert(c.coarseNeighborCart);
        }
    }
    BOOST_CHECK_EQUAL(connectionsOf.size(), 16u);         // 4 (j) x 4 (k) sub-cells
    for (const auto& [cell, n] : connectionsOf) {
        BOOST_CHECK_EQUAL(n, 2);
    }

    int split = 0;
    for (const auto& [cell, nb] : neighborsOf) {
        if (nb.size() >= 2) {
            ++split;
            BOOST_CHECK_EQUAL(cell[2], 1);                // the box layer across depth 2.0
        }
    }
    BOOST_CHECK_EQUAL(split, 4);                          // the fault splits box layer k=1, all 4 j
    BOOST_CHECK_EQUAL(domainBoundaryParts, 4);            // the throw puts box layer k=3 below column i=1

    // The piece of a box cell's face toward cart: its four corners at x = 2,
    // y in [0, 0.5] (j = 0), depth top to bottom; in any order, up to round-off.
    const auto checkPiece = [&conns](const std::array<int,3>& cell, int cart, double top, double bottom) {
        const auto piece = std::find_if(conns.begin(), conns.end(), [&](const auto& c) {
            return (c.boxCell == cell) && (c.coarseNeighborCart == cart);
        });
        BOOST_REQUIRE(piece != conns.end());
        for (const auto& corner : { std::array{2.0, 0.0, top}, std::array{2.0, 0.5, top},
                                    std::array{2.0, 0.5, bottom}, std::array{2.0, 0.0, bottom} }) {
            BOOST_CHECK(std::any_of(piece->faceNodes.begin(), piece->faceNodes.end(),
                                    [&corner](const auto& node) {
                                        return std::abs(node[0] - corner[0]) < 1e-12
                                            && std::abs(node[1] - corner[1]) < 1e-12
                                            && std::abs(node[2] - corner[2]) < 1e-12;
                                    }));
        }
    };
    // Box cell (0,0,1) straddles the fault: coarse (1,0,0) above depth 2.0, coarse (1,0,1) below.
    checkPiece({0,0,1},  1, 1.6, 2.0);
    checkPiece({0,0,1}, 11, 2.0, 2.6);
    // Box cell (0,0,3) reaches below column i=1: coarse (1,0,1) above depth 4.0, the domain below.
    checkPiece({0,0,3}, 11, 3.6, 4.0);
    checkPiece({0,0,3}, -1, 4.0, 4.6);
}

BOOST_AUTO_TEST_CASE(unfaultedBoundaryGivesOneToOne)
{
    // Same grid; the box's i-plus boundary (toward i=4) is unfaulted (i>=2 all
    // thrown equally), so it is a clean 1:1 refined-to-coarse mosaic.
    auto depth = [](int i_, int, int k_) {
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;
    };
    auto g = makeGrid({5, 2, 2}, depth);

    const auto conns = Opm::Refinement::boxBoundaryConnections(
        g.dims, g.coord.data(), g.zcorn.data(), /*actnum=*/nullptr,
        Opm::Refinement::BlockRefinement{ .name = "BOX", .cellsPerDim = {2,2,2},
                                          .startIJK = {2,0,0}, .endIJK = {4,2,2} },
        /*axis=*/0, /*side=*/+1, /*edgeConformal=*/true);

    BOOST_REQUIRE(!conns.empty());

    std::map<std::array<int,3>, std::set<int>> neighborsOf;
    for (const auto& c : conns) {
        BOOST_CHECK_EQUAL(c.boxCell[0], 3);               // box's last i-layer
        BOOST_CHECK_EQUAL(Opm::Lgr::getIJK(c.coarseNeighborCart, g.dims)[0], 4); // coarse neighbour at i=4
        neighborsOf[c.boxCell].insert(c.coarseNeighborCart);
    }
    BOOST_CHECK_EQUAL(neighborsOf.size(), 16u);
    for (const auto& [cell, nb] : neighborsOf) {
        BOOST_CHECK_EQUAL(nb.size(), 1u);                 // exactly one neighbour each
    }
}

BOOST_AUTO_TEST_CASE(domainBoundarySideFacesTheDomain)
{
    // The box's j-minus side is the domain edge (box starts at j=0): no shell,
    // so each boundary cell gets one face toward the domain.
    auto g = makeGrid({5, 2, 2}, flat);

    const auto conns = Opm::Refinement::boxBoundaryConnections(
        g.dims, g.coord.data(), g.zcorn.data(), /*actnum=*/nullptr,
        Opm::Refinement::BlockRefinement{ .name = "BOX", .cellsPerDim = {2,2,2},
                                          .startIJK = {2,0,0}, .endIJK = {4,2,2} },
        /*axis=*/1, /*side=*/-1, /*edgeConformal=*/true);

    std::set<std::array<int,3>> cells;
    for (const auto& c : conns) {
        BOOST_CHECK_EQUAL(c.coarseNeighborCart, -1);
        BOOST_CHECK_EQUAL(c.boxCell[1], 0);
        cells.insert(c.boxCell);
    }
    BOOST_CHECK_EQUAL(conns.size(), 16u);
    BOOST_CHECK_EQUAL(cells.size(), 16u);
}
