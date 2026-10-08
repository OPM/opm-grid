/*
  Copyright 2026 SINTEF Digital.

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

#include <opm/grid/cpgrid/refinement/conforming/FaultedBoundaryFaces.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <array>
#include <functional>
#include <map>
#include <set>
#include <vector>

namespace
{

struct MPIFixture
{
    MPIFixture()
    {
        auto& argv = boost::unit_test::framework::master_test_suite().argv;
        auto& argc = boost::unit_test::framework::master_test_suite().argc;
        Dune::MPIHelper::instance(argc, argv);
    }
};

int cellOf(int doubled) { return doubled / 2; }
int sideOf(int doubled) { return doubled % 2; }

struct Grdecl
{
    std::array<int,3> dims;
    std::vector<double> coord;
    std::vector<double> zcorn;
};

// Vertical-pillar grdecl with per-corner depths from depth(i_, j_, k_) on the
// doubled lattice (same convention as the other refinement tests).
Grdecl makeGrid(const std::array<int,3>& dims,
                const std::function<double(int,int,int)>& depth)
{
    Grdecl g;
    g.dims = dims;
    const auto& [nx, ny, nz] = dims;
    g.coord.resize(6*(nx + 1)*(ny + 1));
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            double* p = &g.coord[6*(static_cast<std::size_t>(j)*(nx + 1) + i)];
            p[0] = i; p[1] = j; p[2] = 0.0;
            p[3] = i; p[4] = j; p[5] = 100.0;
        }
    }
    g.zcorn.resize(8*static_cast<std::size_t>(nx)*ny*nz);
    for (int k_ = 0; k_ < 2*nz; ++k_) {
        for (int j_ = 0; j_ < 2*ny; ++j_) {
            for (int i_ = 0; i_ < 2*nx; ++i_) {
                g.zcorn[static_cast<std::size_t>(i_) + 2*static_cast<std::size_t>(nx)*j_
                        + 4*static_cast<std::size_t>(nx)*ny*k_] = depth(i_, j_, k_);
            }
        }
    }
    return g;
}

} // namespace

BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_AUTO_TEST_CASE(faultedIBoundarySplitsConnections)
{
    // Fault between i=1 and i=2 (throw 0.6). Box i2-3, refined 2x2x2; its
    // i-minus boundary sits on the fault. Boundary cells that straddle the
    // throw must connect to two coarse neighbours.
    auto depth = [](int i_, int, int k_) {
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;
    };
    auto g = makeGrid({5, 2, 2}, depth);

    const auto conns = Opm::Refinement::faultedBoundaryConnections(
        g.dims, g.coord.data(), g.zcorn.data(), nullptr,
        Opm::Refinement::BlockRefinement{ "BOX", "GLOBAL", {2,2,2}, {2,0,0}, {4,2,2} }, /*axis=*/0, /*side=*/-1, /*edgeConformal=*/false);

    BOOST_REQUIRE(!conns.empty());

    std::map<std::array<int,3>, std::set<int>> neighborsOf; // box cell -> coarse carts
    std::set<std::array<int,3>> boxCells;
    int domainParts = 0;
    for (const auto& c : conns) {
        BOOST_CHECK_EQUAL(c.boxCell[0], 0);               // box's i-minus layer
        BOOST_CHECK_GE(c.faceNodes.size(), 3u);           // a real polygon
        boxCells.insert(c.boxCell);
        if (c.coarseNeighborCart < 0) {
            ++domainParts;                                // throw scarp (faces domain)
        }
        else {
            BOOST_CHECK_EQUAL(c.coarseNeighborCart % 5, 1); // coarse neighbour at i=1
            neighborsOf[c.boxCell].insert(c.coarseNeighborCart);
        }
    }
    BOOST_CHECK_EQUAL(boxCells.size(), 16u);              // 4 (j) x 4 (k) sub-cells

    int split = 0;
    for (const auto& [cell, nb] : neighborsOf) {
        if (nb.size() >= 2) ++split;
    }
    BOOST_CHECK_GT(split, 0);                             // the fault splits faces
    BOOST_CHECK_GT(domainParts, 0);                      // the throw lifts the top
                                                          // sub-faces above the grid
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

    const auto conns = Opm::Refinement::faultedBoundaryConnections(
        g.dims, g.coord.data(), g.zcorn.data(), nullptr,
        Opm::Refinement::BlockRefinement{ "BOX", "GLOBAL", {2,2,2}, {2,0,0}, {4,2,2} }, /*axis=*/0, /*side=*/+1, /*edgeConformal=*/false);

    BOOST_REQUIRE(!conns.empty());

    std::map<std::array<int,3>, std::set<int>> neighborsOf;
    for (const auto& c : conns) {
        BOOST_CHECK_EQUAL(c.boxCell[0], 3);               // box's last i-layer
        BOOST_CHECK_EQUAL(c.coarseNeighborCart % 5, 4);   // coarse neighbour at i=4
        neighborsOf[c.boxCell].insert(c.coarseNeighborCart);
    }
    BOOST_CHECK_EQUAL(neighborsOf.size(), 16u);
    for (const auto& [cell, nb] : neighborsOf) {
        BOOST_CHECK_EQUAL(nb.size(), 1u);                 // exactly one neighbour each
    }
}

BOOST_AUTO_TEST_CASE(domainBoundarySideHasNoConnections)
{
    // The box's j-minus side is the domain edge (box starts at j=0): no shell,
    // so no connections.
    auto depth = [](int, int, int k_) { return 2.0*(cellOf(k_) + sideOf(k_)); };
    auto g = makeGrid({5, 2, 2}, depth);

    const auto conns = Opm::Refinement::faultedBoundaryConnections(
        g.dims, g.coord.data(), g.zcorn.data(), nullptr,
        Opm::Refinement::BlockRefinement{ "BOX", "GLOBAL", {2,2,2}, {2,0,0}, {4,2,2} }, /*axis=*/1, /*side=*/-1, /*edgeConformal=*/false);

    BOOST_CHECK(conns.empty());
}
