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
#include "config.h"

#define BOOST_TEST_MODULE GeometricCheckTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>
#include <opm/grid/cpgrid/GeometricCheck.hpp>
#include <opm/grid/cpgrid/coarsening/CornerPointCoarsening.hpp>

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_OPM_COMMON
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#endif

#include <cmath>
#include <vector>

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

namespace
{
using Opm::Coarsening::CoarsenRequest;
using Opm::Coarsening::Grdecl;

constexpr int n = 4;
constexpr double dx = 100.0, dz = 10.0;

Grdecl uniformGrid()
{
    Grdecl g;
    g.dims = {n, n, n};
    for (int j = 0; j <= n; ++j) {
        for (int i = 0; i <= n; ++i) {
            g.coord.insert(g.coord.end(), {i*dx, j*dx, 0.0, i*dx, j*dx, n*dz});
        }
    }
    g.zcorn.resize(8ull*n*n*n);
    for (int k = 0; k < 2*n; ++k) {
        for (int j = 0; j < 2*n; ++j) {
            for (int i = 0; i < 2*n; ++i) {
                g.zcorn[i + 2ull*n*(j + 2ull*n*k)] = ((k + 1)/2)*dz;
            }
        }
    }
    g.actnum.assign(1ull*n*n*n, 1);
    return g;
}

// ZCORN entry of corner (di, dj, dk) of cell (i, j, k).
double& z(Grdecl& g, int i, int j, int k, int di, int dj, int dk)
{
    return g.zcorn[(2*i + di) + 2ull*n*((2*j + dj) + 2ull*n*(2*k + dk))];
}

grdecl view(const Grdecl& g)
{
    grdecl in{};
    for (int d = 0; d < 3; ++d) {
        in.dims[d] = g.dims[d];
    }
    in.coord = g.coord.data();
    in.zcorn = g.zcorn.data();
    in.actnum = g.actnum.data();
    return in;
}

// Corner-point processing without pinch, as a deck without PINCH.
Dune::cpgrid::GeometricCheck processed(const Grdecl& g, bool edgeConformal)
{
    Dune::CpGrid grid;
    grid.processEclipseFormat(view(g), false, false, edgeConformal);
    return Dune::cpgrid::checkGeometric(grid);
}

// Merge route: edge-conformal with pinch, as the mechanics grid is built.
Dune::cpgrid::GeometricCheck merged(const Grdecl& g, const std::vector<CoarsenRequest>& requests)
{
    const auto layout = Opm::Coarsening::blockLayout(g.dims, requests);
    Dune::CpGrid grid;
    grid.processEclipseFormatCoarsened(view(g), layout.blockOfCartesian, layout.boxes,
                                       /*edge_conformal*/ true, /*collapse_coarse_faces*/ true);
    return Dune::cpgrid::checkGeometric(grid);
}

// Columns (1..2) x (1..2): layer 1 removed to zero thickness at its top.
Grdecl pinchedLayer(bool cellBelowTakesTheVolume)
{
    auto g = uniformGrid();
    for (int j = 1; j < 3; ++j) {
        for (int i = 1; i < 3; ++i) {
            g.actnum[i + n*(j + n*1)] = 0;
            for (int dj = 0; dj < 2; ++dj) {
                for (int di = 0; di < 2; ++di) {
                    z(g, i, j, 1, di, dj, 1) = z(g, i, j, 1, di, dj, 0);
                    if (cellBelowTakesTheVolume) {
                        z(g, i, j, 2, di, dj, 0) = z(g, i, j, 1, di, dj, 0);
                    }
                }
            }
        }
    }
    return g;
}
} // namespace

BOOST_AUTO_TEST_CASE(UniformGridIsGeometric)
{
    const auto check = processed(uniformGrid(), true);
    BOOST_TEST_MESSAGE(check.summary());
    BOOST_CHECK(check.ok());
    BOOST_REQUIRE_EQUAL(check.boundaries.size(), 1u);
    BOOST_CHECK_CLOSE(check.boundaries[0].volume, n*dx*n*dx*n*dz, 1e-9);
}

// Layer 1 lifted off layer 2 in four columns: the gap is a closed surface of its own.
BOOST_AUTO_TEST_CASE(VoidBetweenLayers)
{
    auto g = uniformGrid();
    for (int j = 1; j < 3; ++j) {
        for (int i = 1; i < 3; ++i) {
            for (int dj = 0; dj < 2; ++dj) {
                for (int di = 0; di < 2; ++di) {
                    z(g, i, j, 1, di, dj, 1) -= 4.0;   // bottom of layer 1 up by 4 m
                }
            }
        }
    }
    for (const bool ec : {false, true}) {
        const auto check = processed(g, ec);
        BOOST_TEST_MESSAGE(check.summary());
        BOOST_CHECK(!check.ok());
        BOOST_REQUIRE_EQUAL(check.boundaries.size(), 2u);
        BOOST_CHECK_CLOSE(-check.boundaries[1].volume, 2*dx*2*dx*4.0, 1e-6);
    }
}

// An inactive cell inside the body is a cavity.
BOOST_AUTO_TEST_CASE(InactiveCellIsACavity)
{
    auto g = uniformGrid();
    g.actnum[1 + n*(1 + n*1)] = 0;
    const auto check = processed(g, true);
    BOOST_TEST_MESSAGE(check.summary());
    BOOST_CHECK(!check.ok());
    BOOST_REQUIRE_EQUAL(check.boundaries.size(), 2u);
    BOOST_CHECK_CLOSE(-check.boundaries[1].volume, dx*dx*dz, 1e-6);
}

// A zero-thickness layer, its volume moved to the cell below (as MINPV merging
// does): without pinch the layers on either side do not touch, a crack of zero
// volume; with pinch, as the mechanics grid is processed, the body is whole.
BOOST_AUTO_TEST_CASE(PinchedLayerCracksWithoutPinch)
{
    const auto g = pinchedLayer(/*cellBelowTakesTheVolume*/ true);
    for (const bool ec : {false, true}) {
        const auto check = processed(g, ec);
        BOOST_TEST_MESSAGE("no pinch, edge-conformal " << ec << ": " << check.summary());
        BOOST_CHECK(!check.ok());
        BOOST_REQUIRE_GE(check.boundaries.size(), 2u);
        BOOST_CHECK_SMALL(check.boundaries[1].volume, 1e-6);
        BOOST_CHECK_CLOSE(check.boundaries[1].area, 2*2*dx*dx*2, 1e-6);   // both sides
    }
    std::vector<CoarsenRequest> none;
    const auto check = merged(g, none);
    BOOST_TEST_MESSAGE("pinch: " << check.summary());
    BOOST_CHECK_EQUAL(check.boundaries.size(), 1u);
    BOOST_CHECK_EQUAL(check.unpairedBoundaryEdges, 0);
}

// The same layer removed without giving its volume away: a void, which pinch
// does not close.
BOOST_AUTO_TEST_CASE(PinchedLayerWithoutMergeLeavesAVoid)
{
    auto g = pinchedLayer(/*cellBelowTakesTheVolume*/ false);
    const auto check = merged(g, {});
    BOOST_TEST_MESSAGE(check.summary());
    BOOST_CHECK(!check.ok());
    BOOST_REQUIRE_EQUAL(check.boundaries.size(), 2u);
    BOOST_CHECK_CLOSE(-check.boundaries[1].volume, 2*dx*2*dx*dz, 1e-6);
}

// Coarsened by the merge with collapsed faces, and as a corner-point grid processed
// edge-conformal: both geometric, also where the columns group their layers
// differently. Without edge-conformal processing the thin cells' corners hang on
// the thick cells' faces.
BOOST_AUTO_TEST_CASE(CoarsenedGrids)
{
    const auto g = uniformGrid();
    CoarsenRequest south, north;
    south.startIJK = {0, 0, 0};
    south.endIJK = {4, 2, 4};
    south.cellsPerDim = {2, 2, 2};
    north.startIJK = {0, 2, 0};
    north.endIJK = {4, 4, 4};
    north.cellsPerDim = {2, 2, 4};

    const auto collapsed = merged(g, {south, north});
    BOOST_TEST_MESSAGE("collapse: " << collapsed.summary());
    BOOST_CHECK(collapsed.ok());

    const auto coarse = Opm::Coarsening::coarsenCornerPoint(g, {south, north}).grid;
    const auto cornerPoint = processed(coarse, true);
    BOOST_TEST_MESSAGE("corner-point, edge-conformal: " << cornerPoint.summary());
    BOOST_CHECK(cornerPoint.ok());

    const auto plain = processed(coarse, false);
    BOOST_TEST_MESSAGE("corner-point, not edge-conformal: " << plain.summary());
    BOOST_CHECK_GT(plain.nonConformingCells, 0);
    BOOST_CHECK_EQUAL(plain.boundaries.size(), 1u);
}

// One column's layer boundary moved 1e-7 m on a pillar it shares with three others.
// Without a tolerance that is a second node, which edge-conformal processing lists in
// the neighbours' faces (no hanging node, but a 1e-7 m edge); with a tolerance the
// points merge and the faces keep four corners.
BOOST_AUTO_TEST_CASE(NearCoincidentPillarPoints)
{
    auto g = uniformGrid();
    z(g, 0, 0, 1, 1, 1, 1) += 1e-7;   // bottom of (0,0,1) at the pillar (100, 100)
    z(g, 0, 0, 2, 1, 1, 0) += 1e-7;   // top of (0,0,2) there

    const auto facesAndNodes = [&g](double tolerance) {
        processed_grid out{};
        const auto in = view(g);
        process_grdecl(0, /*edge_conformal*/ 1, tolerance, &in, nullptr, &out);
        // the face between (1,0,2) and (1,1,2): at y = 100, x in [100, 200], z in [20, 30]
        int corners = -1;
        for (int f = 0; f < out.number_of_faces; ++f) {
            const int c0 = out.face_neighbors[2*f], c1 = out.face_neighbors[2*f + 1];
            if (std::min(c0, c1) == 1 + n*(0 + n*2) && std::max(c0, c1) == 1 + n*(1 + n*2)) {
                corners = static_cast<int>(out.face_node_ptr[f + 1] - out.face_node_ptr[f]);
            }
        }
        const std::pair<int,int> result{out.number_of_nodes, corners};
        free_processed_grid(&out);
        return result;
    };
    const auto [nodesExact, cornersExact] = facesAndNodes(0.0);
    const auto [nodesMerged, cornersMerged] = facesAndNodes(1e-6);
    BOOST_CHECK_EQUAL(nodesExact, nodesMerged + 1);
    BOOST_CHECK_EQUAL(cornersExact, 5);
    BOOST_CHECK_EQUAL(cornersMerged, 4);

    const auto check = processed(g, true);
    BOOST_TEST_MESSAGE(check.summary());
    BOOST_CHECK(check.ok());
}

// A fault between i = 1 and i = 2 whose throw changes sign along it, so layer
// boundaries of the two sides cross: the crossings are nodes of the fault faces,
// and edge-conformal processing lists them in the top and bottom faces too.
BOOST_AUTO_TEST_CASE(FaultWithCrossingLayers)
{
    auto g = uniformGrid();
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            for (int i = 2; i < n; ++i) {
                for (int dk = 0; dk < 2; ++dk) {
                    for (int dj = 0; dj < 2; ++dj) {
                        for (int di = 0; di < 2; ++di) {
                            const double y = (j + dj)*dx;                 // 0 .. 400
                            z(g, i, j, k, di, dj, dk) += 6.0*(y/(n*dx)) - 2.5;   // -2.5 .. +3.5 m
                        }
                    }
                }
            }
        }
    }
    {
        processed_grid out{};
        const auto in = view(g);
        process_grdecl(0, 1, 0.0, &in, nullptr, &out);
        BOOST_CHECK_GT(out.number_of_nodes, out.number_of_nodes_on_pillars);   // crossings
        free_processed_grid(&out);
    }
    const auto plain = processed(g, false);
    BOOST_TEST_MESSAGE("not edge-conformal: " << plain.summary());
    BOOST_CHECK_GT(plain.nonConformingCells, 0);

    const auto conformal = processed(g, true);
    BOOST_TEST_MESSAGE("edge-conformal: " << conformal.summary());
    BOOST_CHECK(conformal.ok());
}

// Field coordinates (~1e5 m) must not make small cells look open.
BOOST_AUTO_TEST_CASE(TranslatedSmallCellsAreGeometric)
{
    for (const double offset : {0.0, 1.0e5}) {
        auto g = uniformGrid();
        for (std::size_t p = 0; p < g.coord.size(); p += 3) {
            g.coord[p] = offset + 0.1*g.coord[p];
            g.coord[p + 1] = offset + 0.1*g.coord[p + 1];
            g.coord[p + 2] = 2000.0 + 0.1*g.coord[p + 2];
        }
        for (auto& zc : g.zcorn) {
            zc = 2000.0 + 0.1*zc;
        }
        const auto check = processed(g, true);
        BOOST_TEST_INFO("offset " << offset << ": " << check.summary());
        BOOST_CHECK(check.ok());
    }
}

#if HAVE_OPM_COMMON
// Layer 2 pinches out to the west: 0.1 m thick on average in column 0 (removed
// by PINCH 0.5) but 2.6 m in column 1, whose west corners are only 0.2 m apart.
// Without merging, those corners give a 0.2 m edge on the shared pillar.
BOOST_AUTO_TEST_CASE(PinchOutLeavesNoSliverEdge)
{
    const auto deck = [](bool pinch) {
        return std::string(R"(
RUNSPEC
DIMENS
2 1 3 /
OIL
WATER
GRID
COORD
  0   0 1000    0   0 1030
100   0 1000  100   0 1030
200   0 1000  200   0 1030
  0 100 1000    0 100 1030
100 100 1000  100 100 1030
200 100 1000  200 100 1030 /
ZCORN
8*1000 8*1010
8*1010 1010 1010.2 1010.2 1015 1010 1010.2 1010.2 1015
1010 1010.2 1010.2 1015 1010 1010.2 1010.2 1015 8*1030 /
PORO
6*0.2 /
PERMX
6*100 /
PERMY
6*100 /
PERMZ
6*100 /
)") + (pinch ? "PINCH\n0.5 /\n" : "");
    };
    const auto build = [](const std::string& text, Dune::CpGrid& grid, double mergeTolerance = 0.0) {
        Opm::EclipseState es(Opm::Parser{}.parseString(text));
        grid.processEclipseFormat(&es.getInputGrid(), &es, false, false, false,
                                  /*edge_conformal*/ true, mergeTolerance);
    };
    const auto shortEdges = [](const Dune::CpGrid& grid, double length) {
        int count = 0;
        for (int f = 0; f < grid.numFaces(); ++f) {
            const int nv = grid.numFaceVertices(f);
            for (int v = 0; v < nv; ++v) {
                auto d = grid.vertexPosition(grid.faceVertex(f, v));
                d -= grid.vertexPosition(grid.faceVertex(f, (v + 1) % nv));
                count += d.two_norm() < length ? 1 : 0;
            }
        }
        return count;
    };
    const auto volume = [](const Dune::CpGrid& grid) {
        double sum = 0.0;
        for (int c = 0; c < grid.numCells(); ++c) {
            sum += grid.cellVolume(c);
        }
        return sum;
    };

    Dune::CpGrid input, pinched;
    build(deck(false), input);
    build(deck(true), pinched);
    BOOST_CHECK_GT(shortEdges(input, 0.5), 0);
    BOOST_CHECK_EQUAL(input.numCells(), 6);

    BOOST_CHECK_EQUAL(pinched.numCells(), 5);
    BOOST_CHECK_EQUAL(shortEdges(pinched, 0.5), 0);
    BOOST_CHECK_CLOSE(volume(pinched), volume(input), 1e-10);
    const auto check = Dune::cpgrid::checkGeometric(pinched);
    BOOST_TEST_MESSAGE("pinched, edge-conformal: " << check.summary());
    BOOST_CHECK(check.ok());

    // Without PINCH, the merge tolerance alone does the same.
    Dune::CpGrid merged;
    build(deck(false), merged, 0.5);
    BOOST_CHECK_EQUAL(merged.numCells(), 5);
    BOOST_CHECK_EQUAL(shortEdges(merged, 0.5), 0);
    BOOST_CHECK_CLOSE(volume(merged), volume(input), 1e-10);
    BOOST_CHECK(Dune::cpgrid::checkGeometric(merged).ok());

    // A tolerance below the sliver leaves it.
    Dune::CpGrid fine;
    build(deck(false), fine, 0.05);
    BOOST_CHECK_EQUAL(fine.numCells(), 6);
    BOOST_CHECK_GT(shortEdges(fine, 0.5), 0);
}
#endif

// Scissor fault along pillar line i = 2. Left of it the top cells are wedges with the i+ side
// collapsed, so their boundary top face shares the fault edge with their bottom face.
BOOST_AUTO_TEST_CASE(PinchedWedgeAtAFault)
{
    auto g = uniformGrid();
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            for (int i = 2; i < n; ++i) {
                for (int c = 0; c < 8; ++c) {
                    const int dj = (c >> 1) & 1;
                    z(g, i, j, k, c & 1, dj, c >> 2) += ((j + dj) % 2 == 0) ? -4.0 : 4.0;
                }
            }
        }
    }
    for (int j = 0; j < n; ++j) {
        for (int dj = 0; dj < 2; ++dj) {
            z(g, 1, j, 0, 1, dj, 0) = z(g, 1, j, 0, 1, dj, 1);
        }
    }
    const auto check = processed(g, true);
    BOOST_TEST_MESSAGE(check.summary());
    BOOST_CHECK(check.ok());
}
