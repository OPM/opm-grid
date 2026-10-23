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

  Local grid refinement on edge-conformal grids. With edge_conformal=true a
  fault inserts the neighbour columns' nodes onto the shared pillar edges,
  so faces carry more than four nodes. These tests probe whether the
  refinement builder copes (review §3 flagged this as untested).
*/
#include <config.h>

#define BOOST_TEST_MODULE EdgeConformalRefinementTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/refinement/conforming/ConformingBlockBuilder.hpp>
#include <opm/grid/cpgrid/refinement/GridStateWriter.hpp>
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <dune/common/parallel/mpihelper.hh>

#include <array>
#include <functional>
#include <memory>
#include <vector>
#include <set>

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

struct TestGrdecl
{
    std::array<int,3> dims;
    std::vector<double> coord;
    std::vector<double> zcorn;
    std::vector<int> actnum;
    grdecl raw() const
    {
        grdecl g;
        g.dims[0] = dims[0]; g.dims[1] = dims[1]; g.dims[2] = dims[2];
        g.coord = coord.data(); g.zcorn = zcorn.data();
        g.actnum = actnum.empty() ? nullptr : actnum.data();
        return g;
    }
};

int cellOf(int doubled) { return doubled / 2; }
int sideOf(int doubled) { return doubled % 2; }

// Vertical-pillar grid; per-corner depth via depth(i_,j_,k_) on the doubled
// lattice, so faults (z jumps between columns) are expressible.
TestGrdecl makeGrid(const std::array<int,3>& dims,
                    const std::function<double(int,int,int)>& depth)
{
    TestGrdecl g;
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

void useConformingBuilder(Dune::CpGrid& grid, std::unique_ptr<Opm::Refinement::Builder> builder)
{
    grid.setLgrBackend(Opm::Refinement::Backend::Conforming);
    grid.setRefinementBuilder(std::move(builder));
}

double totalVolume(const Dune::CpGrid& grid)
{
    double v = 0.0;
    for (const auto& e : Dune::elements(grid.leafGridView())) {
        v += e.geometry().volume();
    }
    return v;
}

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(MPIFixture);

// Baseline: an edge-conformal grid with a fault, refining a box that does
// not touch the fault. The refined region is away from the edge-conformal
// (multi-node) faces.
BOOST_AUTO_TEST_CASE(refineBoxAwayFromFaultOnEdgeConformalGrid)
{
    // Fault between i-columns 1 and 2, throw 0.6 (< layer thickness 2).
    auto depth = [](int i_, int, int k_) {
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;
    };
    auto g = makeGrid({4, 2, 2}, depth);

    Dune::CpGrid grid;
    auto raw = g.raw();
    grid.processEclipseFormat(raw, /*remove_ij_boundary=*/false,
                              /*turn_normals=*/false, /*edge_conformal=*/true);
    const double v0 = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        g.dims, g.coord, g.zcorn, g.actnum));

    // Box = the single cell (0,0,0)..(0,1,1) on the unfaulted left side,
    // i in [0,1): its boundary faces toward the fault are at i=1, which is
    // not the fault plane (i=2).
    grid.addLgrsUpdateLeafView({{2,2,2}}, {{0,0,0}}, {{1,2,2}}, {"LGR1"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_CLOSE(totalVolume(grid), v0, 1e-8);
}

// The box *contains* the fault (fault plane interior, boundaries
// unfaulted), on an edge-conformal grid. Faults inside a block are
// supported; edge-conformality of level zero must not break that.
BOOST_AUTO_TEST_CASE(refineBoxContainingFaultOnEdgeConformalGrid)
{
    auto depth = [](int i_, int, int k_) {
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;
    };
    auto g = makeGrid({4, 2, 2}, depth);

    Dune::CpGrid grid;
    auto raw = g.raw();
    grid.processEclipseFormat(raw, false, false, /*edge_conformal=*/true);
    const double v0 = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        g.dims, g.coord, g.zcorn, g.actnum));

    // Box i in [1,3): fault plane i=2 is inside, boundaries i=1 and i=3 are
    // unfaulted.
    grid.addLgrsUpdateLeafView({{2,2,2}}, {{1,0,0}}, {{3,2,2}}, {"LGR1"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_CLOSE(totalVolume(grid), v0, 1e-8);
}

// Torture case: a pillar with all four surrounding columns faulted to
// different depths. On an edge-conformal grid that pillar accumulates the
// nodes of all four columns. We build it and refine a box in the flat
// corner, away from the pillar - the complex pillar's edge-conformal faces
// must be preserved through refinement.
BOOST_AUTO_TEST_CASE(fourColumnsFaultedAroundAPillarEdgeConformal)
{
    // 4x4x2 grid. The four cells around the central pillar (2,2) -
    // (1,1),(2,1),(1,2),(2,2) - each get a distinct vertical throw; the rest
    // of the grid is flat.
    auto depth = [](int i_, int j_, int k_) {
        const int ci = cellOf(i_);
        const int cj = cellOf(j_);
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        double th = 0.0;
        if (ci == 1 && cj == 1) th = 0.2;
        else if (ci == 2 && cj == 1) th = 0.4;
        else if (ci == 1 && cj == 2) th = 0.6;
        else if (ci == 2 && cj == 2) th = 0.8;
        return base + th;
    };
    auto g = makeGrid({4, 4, 2}, depth);

    // Confirm edge_conformal actually changes this grid: it inserts the
    // neighbour columns' nodes onto the shared faces' edges, so the total
    // face-node count grows (the face and vertex counts may stay the same).
    const auto totalFaceNodes = [](Dune::CpGrid& cg) {
        return Opm::Refinement::GridStateWriter::faceToPoint(*cg.currentData().back()).dataSize();
    };
    Dune::CpGrid plain;
    auto rawPlain = g.raw();
    plain.processEclipseFormat(rawPlain, false, false, /*edge_conformal=*/false);

    Dune::CpGrid grid;
    auto raw = g.raw();
    grid.processEclipseFormat(raw, false, false, /*edge_conformal=*/true);
    BOOST_CHECK_GT(totalFaceNodes(grid), totalFaceNodes(plain));

    const double v0 = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        g.dims, g.coord, g.zcorn, g.actnum));

    // Refine the flat corner cell (0,0): its boundaries (i=1, j=1) are
    // unfaulted, away from the complex pillar.
    grid.addLgrsUpdateLeafView({{2,2,2}}, {{0,0,0}}, {{1,1,2}}, {"LGR1"});

    BOOST_REQUIRE_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_CLOSE(totalVolume(grid), v0, 1e-8);
}

// The edge-conformal toggle: with it on, refinement makes the leaf
// edge-conformal by inserting the LGR-boundary nodes into the coarse faces
// that share those edges. Off keeps the current (face-conformal) behaviour.
// Geometry (volume) is identical either way - only face node lists change.
BOOST_AUTO_TEST_CASE(edgeConformalToggleInsertsBoundaryNodes)
{
    // Flat 3x3x2 grid; refine the central cell (1,1) so its boundary pillars
    // and lateral edges acquire subdivision nodes that the diagonal coarse
    // neighbours would otherwise leave hanging.
    auto depth = [](int, int, int k_) {
        return static_cast<double>(cellOf(k_) + sideOf(k_));
    };
    auto g = makeGrid({3, 3, 2}, depth);

    const auto faceNodesAfterRefine = [&](bool edgeConformal) {
        Dune::CpGrid grid;
        auto raw = g.raw();
        grid.processEclipseFormat(raw, false);
        useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
            g.dims, g.coord, g.zcorn, g.actnum, edgeConformal));
        grid.addLgrsUpdateLeafView({{3,3,3}}, {{1,1,0}}, {{2,2,2}}, {"LGR1"});
        return std::make_pair(
            Opm::Refinement::GridStateWriter::faceToPoint(*grid.currentData().back()).dataSize(),
            totalVolume(grid));
    };

    const auto [nodesOff, volOff] = faceNodesAfterRefine(false);
    const auto [nodesOn, volOn] = faceNodesAfterRefine(true);

    // On inserts hanging nodes into coarse faces -> strictly more face nodes.
    BOOST_CHECK_GT(nodesOn, nodesOff);
    // Geometry is unchanged by the post-pass.
    BOOST_CHECK_CLOSE(volOn, volOff, 1e-10);
    BOOST_CHECK_CLOSE(volOn, 3.0*3.0*2.0, 1e-8);
}

// The global edge-conformal property: after the pass, NO leaf node lies in
// the interior of ANY face's edge without being listed in that face - no
// matter which cell defined the node. (Off, hanging nodes exist; on, none.)
BOOST_AUTO_TEST_CASE(edgeConformalLeavesNoHangingNode)
{
    auto depth = [](int, int, int k_) {
        return static_cast<double>(cellOf(k_) + sideOf(k_));
    };
    auto g = makeGrid({3, 3, 2}, depth);

    // Count nodes that lie strictly on some face edge but are not listed in
    // that face. Edge-conformal <=> this count is zero.
    const auto countHangingNodes = [&](bool edgeConformal) {
        Dune::CpGrid grid;
        auto raw = g.raw();
        grid.processEclipseFormat(raw, false);
        useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
            g.dims, g.coord, g.zcorn, g.actnum, edgeConformal));
        grid.addLgrsUpdateLeafView({{3,3,3}}, {{1,1,0}}, {{2,2,2}}, {"LGR1"});

        auto& leaf = *grid.currentData().back();
        const auto& cg = Opm::Refinement::GridStateWriter::geometry(leaf)
                               .geomVector(std::integral_constant<int,3>());
        const int nc = cg.size();
        std::vector<Dune::FieldVector<double,3>> P(nc);
        for (int i = 0; i < nc; ++i) P[i] = cg.get(i).center();

        auto& f2p = Opm::Refinement::GridStateWriter::faceToPoint(leaf);
        int hanging = 0;
        for (int face = 0; face < f2p.size(); ++face) {
            auto row = f2p[face];
            const int n = row.size();
            std::set<int> inFace(row.begin(), row.end());
            for (int e = 0; e < n; ++e) {
                const int a = row[e], b = row[(e + 1) % n];
                auto d = P[b]; d -= P[a];
                const double len2 = d.two_norm2();
                if (len2 <= 0) continue;
                for (int c = 0; c < nc; ++c) {
                    if (inFace.count(c)) continue;
                    auto ac = P[c]; ac -= P[a];
                    const double t = (ac * d) / len2;
                    if (t <= 1e-9 || t >= 1.0 - 1e-9) continue;
                    auto proj = d; proj *= t; proj += P[a];
                    auto diff = P[c]; diff -= proj;
                    if (diff.two_norm() < 1e-9) ++hanging;
                }
            }
        }
        return hanging;
    };

    BOOST_CHECK_GT(countHangingNodes(false), 0); // refinement opens hanging nodes
    BOOST_CHECK_EQUAL(countHangingNodes(true), 0); // the pass removes them all
}

// A box boundary on the fault plane: the split faces are rebuilt from the
// corner-point processor and must survive the edge-conformal pass (build,
// volume conserved).
BOOST_AUTO_TEST_CASE(refineBoxBoundaryOnFaultBuilds)
{
    auto depth = [](int i_, int, int k_) {
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;
    };
    auto g = makeGrid({4, 2, 2}, depth);

    Dune::CpGrid grid;
    auto raw = g.raw();
    grid.processEclipseFormat(raw, false, false, /*edge_conformal=*/true);
    const double v0 = totalVolume(grid);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        g.dims, g.coord, g.zcorn, g.actnum));

    // Box i in [0,2): right boundary is the fault plane i=2.
    BOOST_REQUIRE_NO_THROW(grid.addLgrsUpdateLeafView({{2,2,2}}, {{0,0,0}}, {{2,2,2}}, {"LGR1"}));
    BOOST_CHECK_EQUAL(grid.maxLevel(), 1);
    BOOST_CHECK_CLOSE(totalVolume(grid), v0, 1e-8);
}
