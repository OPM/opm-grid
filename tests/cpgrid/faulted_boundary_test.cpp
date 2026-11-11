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

namespace
{

// The piece of box cell `cell`'s face toward the coarse cell `cart` (-1: the domain) has
// exactly the given corners, in any order, up to round-off.
void checkFaceNodes(const std::vector<Opm::Refinement::BoundaryConnection>& conns,
                    const std::array<int,3>& cell, int cart,
                    const std::vector<std::array<double,3>>& corners)
{
    const auto piece = std::find_if(conns.begin(), conns.end(), [&](const auto& c) {
        return (c.boxCell == cell) && (c.coarseNeighborCart == cart);
    });
    BOOST_REQUIRE(piece != conns.end());
    const auto samePoint = [](const auto& a, const auto& b) {
        return std::abs(a[0] - b[0]) < 1e-12 && std::abs(a[1] - b[1]) < 1e-12
            && std::abs(a[2] - b[2]) < 1e-12;
    };
    BOOST_CHECK(std::is_permutation(piece->faceNodes.begin(), piece->faceNodes.end(),
                                    corners.begin(), corners.end(), samePoint));
}

} // anonymous namespace

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
        checkFaceNodes(conns, cell, cart, { {2.0, 0.0, top}, {2.0, 0.5, top},
                                            {2.0, 0.5, bottom}, {2.0, 0.0, bottom} });
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
    // The same grid; the box's i-plus side (toward i=4) is not faulted: i>=2 is
    // all thrown by 0.6, so each box cell there faces one coarse cell. Indices
    // are zero-based.
    //
    // Plan view (5 x 2 x 2; box = i2-3, all j, all k). The tested face is '#':
    //
    //          i=0   i=1   i=2   i=3   i=4
    //        +-----+-----+-----+-----#-----+
    //    j=1 |     |     |    box    #     |
    //        +-----+-----+           #-----+
    //    j=0 |     |     |           #     |
    //        +-----+-----+-----+-----#-----+
    //                                ^ x=4, no fault (axis=0, side=+1)
    //
    // Section at x=4 (depth down). Both sides are thrown by 0.6, so every
    // coarse layer boundary meets a box sub-layer boundary:
    //
    //   depth   box, i-layer 3   coarse i=4
    //    0.6   +----------+----------+
    //          |   k=0    |          |
    //    1.6   +----------+   k=0    |
    //          |   k=1    |          |
    //    2.6   +----------+----------+
    //          |   k=2    |          |
    //    3.6   +----------+   k=1    |
    //          |   k=3    |          |
    //    4.6   +----------+----------+
    //
    // Each box face is 0.5 wide in j and lies inside one coarse face, so the
    // 4 (j) x 4 (k) boundary cells get 1 connection each: 16 in all.
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

    // Each of the 16 refined cells at the i-plus box boundary gets 1 connection.
    BOOST_CHECK_EQUAL(conns.size(), 16u);

    std::map<std::array<int,3>, std::set<int>> neighborsOf; // box cell -> coarse carts
    std::map<std::array<int,3>, int> connectionsOf;         // box cell -> its connections
    int domainBoundaryParts = 0;
    for (const auto& c : conns) {
        BOOST_CHECK_EQUAL(c.boxCell[0], 3);               // all faces lie in the box's i-plus layer (plane x = 4)
        BOOST_CHECK_EQUAL(c.faceNodes.size(), 4u);        // a quadrilateral: no hanging nodes here
        ++connectionsOf[c.boxCell];
        if (c.coarseNeighborCart < 0) {
            ++domainBoundaryParts;
        }
        else {
            // grid dims {5,2,2}: coarse cart = i + 5*j + 10*k
            BOOST_CHECK_EQUAL(Opm::Lgr::getIJK(c.coarseNeighborCart, g.dims)[0], 4); // coarse neighbour at i=4
            neighborsOf[c.boxCell].insert(c.coarseNeighborCart);
        }
    }
    BOOST_CHECK_EQUAL(connectionsOf.size(), 16u);         // 4 (j) x 4 (k) sub-cells
    for (const auto& [cell, n] : connectionsOf) {
        BOOST_CHECK_EQUAL(n, 1);
    }

    int split = 0;
    for (const auto& [cell, nb] : neighborsOf) {
        if (nb.size() >= 2) {
            ++split;
        }
    }
    BOOST_CHECK_EQUAL(split, 0);                          // no box cell faces two coarse cells
    BOOST_CHECK_EQUAL(domainBoundaryParts, 0);            // column i=4 covers the whole side

    // Box cell (3,0,1) faces coarse (4,0,0), box cell (3,0,2) coarse (4,0,1); each
    // face is x = 4, y in [0, 0.5] (j = 0), one box sub-layer deep.
    checkFaceNodes(conns, {3,0,1},  4, { {4.0, 0.0, 1.6}, {4.0, 0.5, 1.6},
                                         {4.0, 0.5, 2.6}, {4.0, 0.0, 2.6} });
    checkFaceNodes(conns, {3,0,2}, 14, { {4.0, 0.0, 2.6}, {4.0, 0.5, 2.6},
                                         {4.0, 0.5, 3.6}, {4.0, 0.0, 3.6} });
}

BOOST_AUTO_TEST_CASE(domainBoundarySideFacesTheDomain)
{
    // The box's j-minus side is the domain edge (box starts at j=0): no shell,
    // so each boundary cell gets one face toward the domain. Indices are
    // zero-based.
    //
    // Plan view (5 x 2 x 2, no fault; box = i2-3, all j, all k). The tested
    // face is '#':
    //
    //          i=0   i=1   i=2   i=3   i=4
    //        +-----+-----+-----+-----+-----+
    //    j=1 |     |     |    box    |     |
    //        +-----+-----+           +-----+
    //    j=0 |     |     |           |     |
    //        +-----+-----+#####+#####+-----+
    //                    ^ y=0, domain edge (axis=1, side=-1)
    //
    // The face y=0 seen from outside (x right, depth down). Nothing lies
    // beyond it, so all 4 (i) x 4 (k) sub-faces face the domain (-1):
    //
    //      x:  2    2.5    3    3.5    4
    //   depth
    //    0.0   +-----+-----+-----+-----+
    //          | -1  | -1  | -1  | -1  |
    //    1.0   +-----+-----+-----+-----+
    //          | -1  | -1  | -1  | -1  |
    //    2.0   +-----+-----+-----+-----+
    //          | -1  | -1  | -1  | -1  |
    //    3.0   +-----+-----+-----+-----+
    //          | -1  | -1  | -1  | -1  |
    //    4.0   +-----+-----+-----+-----+
    auto g = makeGrid({5, 2, 2}, flat);

    const auto conns = Opm::Refinement::boxBoundaryConnections(
        g.dims, g.coord.data(), g.zcorn.data(), /*actnum=*/nullptr,
        Opm::Refinement::BlockRefinement{ .name = "BOX", .cellsPerDim = {2,2,2},
                                          .startIJK = {2,0,0}, .endIJK = {4,2,2} },
        /*axis=*/1, /*side=*/-1, /*edgeConformal=*/true);

    // Each of the 16 refined cells at the j-minus box boundary gets 1 connection.
    BOOST_CHECK_EQUAL(conns.size(), 16u);

    std::map<std::array<int,3>, int> connectionsOf;         // box cell -> its connections
    int domainBoundaryParts = 0;
    for (const auto& c : conns) {
        BOOST_CHECK_EQUAL(c.boxCell[1], 0);               // all faces lie in the box's j-minus layer (plane y = 0)
        BOOST_CHECK_EQUAL(c.faceNodes.size(), 4u);        // a quadrilateral: no hanging nodes here
        BOOST_CHECK_EQUAL(c.coarseNeighborCart, -1);      // nothing across: the domain
        ++connectionsOf[c.boxCell];
        if (c.coarseNeighborCart < 0) {
            ++domainBoundaryParts;
        }
    }
    BOOST_CHECK_EQUAL(connectionsOf.size(), 16u);         // 4 (i) x 4 (k) sub-cells
    for (const auto& [cell, n] : connectionsOf) {
        BOOST_CHECK_EQUAL(n, 1);
    }
    BOOST_CHECK_EQUAL(domainBoundaryParts, 16);           // every sub-face faces the domain

    // Box cell (0,0,3): x in [2, 2.5], y = 0, depth 3 to 4; box cell (3,0,0): x in
    // [3.5, 4], depth 0 to 1.
    checkFaceNodes(conns, {0,0,3}, -1, { {2.0, 0.0, 3.0}, {2.5, 0.0, 3.0},
                                         {2.5, 0.0, 4.0}, {2.0, 0.0, 4.0} });
    checkFaceNodes(conns, {3,0,0}, -1, { {3.5, 0.0, 0.0}, {4.0, 0.0, 0.0},
                                         {4.0, 0.0, 1.0}, {3.5, 0.0, 1.0} });
}
