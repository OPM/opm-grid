/*
  Copyright 2026 Equinor ASA

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

#define BOOST_TEST_MODULE UnstructuredGridPartitioningTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/UnstructuredGrid.h>
#include <opm/grid/cart_grid.h>
#include <opm/grid/common/UnstructuredGridPartitioning.hpp>

#include <algorithm>
#include <memory>
#include <numeric>
#include <set>
#include <vector>

using namespace Dune::unstructured;

namespace
{

struct GridPtr
{
    explicit GridPtr(UnstructuredGrid* g) : grid(g, &destroy_grid) {}
    const UnstructuredGrid& operator*() const { return *grid; }
    std::unique_ptr<UnstructuredGrid, void (*)(UnstructuredGrid*)> grid;
};

} // anonymous namespace

BOOST_AUTO_TEST_CASE(OverlapLayersOnALine)
{
    GridPtr g(create_grid_cart3d(4, 1, 1));
    const std::vector<int> part{0, 0, 1, 1};

    auto none = computeOverlapCells(*g, part, 2, 0);
    BOOST_CHECK(none[0].empty() && none[1].empty());

    auto one = computeOverlapCells(*g, part, 2, 1);
    BOOST_CHECK((one[0] == std::vector<int>{2}));
    BOOST_CHECK((one[1] == std::vector<int>{1}));

    auto two = computeOverlapCells(*g, part, 2, 2);
    BOOST_CHECK((two[0] == std::vector<int>{2, 3}));
    BOOST_CHECK((two[1] == std::vector<int>{0, 1}));

    // More layers than cells must not change anything or loop forever.
    auto many = computeOverlapCells(*g, part, 2, 10);
    BOOST_CHECK(many == two);
}

BOOST_AUTO_TEST_CASE(InvalidPartitionThrows)
{
    GridPtr g(create_grid_cart3d(2, 1, 1));
    BOOST_CHECK_THROW(computeOverlapCells(*g, {0}, 2, 1), std::invalid_argument);
    BOOST_CHECK_THROW(computeOverlapCells(*g, {0, 2}, 2, 1), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(ExtractLocalGridKeepsCellOrderAndTopology)
{
    GridPtr g(create_grid_cart3d(4, 1, 1));
    // Deliberately not sorted: the local order must follow the argument.
    const std::vector<int> cells{2, 1};
    auto local = extractLocalGrid(*g, cells);
    const auto& lg = *local.grid;

    BOOST_CHECK_EQUAL(lg.number_of_cells, 2);
    // two cells sharing a face: 6 + 6 - 1 faces and 12 nodes
    BOOST_CHECK_EQUAL(lg.number_of_faces, 11);
    BOOST_CHECK_EQUAL(lg.number_of_nodes, 12);
    BOOST_CHECK_EQUAL(lg.global_cell[0], 2);
    BOOST_CHECK_EQUAL(lg.global_cell[1], 1);
    BOOST_CHECK_EQUAL(lg.cartdims[0], 4);

    // Geometry follows the cells.
    BOOST_CHECK_CLOSE(lg.cell_centroids[0], g.grid->cell_centroids[3 * 2], 1e-12);
    BOOST_CHECK_CLOSE(lg.cell_centroids[3], g.grid->cell_centroids[3 * 1], 1e-12);

    // Exactly one interior face; all others are boundary faces.
    int interior = 0;
    for (int f = 0; f < lg.number_of_faces; ++f) {
        interior += (lg.face_cells[2 * f] >= 0 && lg.face_cells[2 * f + 1] >= 0);
    }
    BOOST_CHECK_EQUAL(interior, 1);
}

BOOST_AUTO_TEST_CASE(ExtractLocalGridConsistency)
{
    GridPtr g(create_grid_cart3d(4, 3, 2));
    const auto& G = *g;
    // Partition into two slabs in x.
    std::vector<int> part(G.number_of_cells);
    for (int c = 0; c < G.number_of_cells; ++c) {
        part[c] = ((c % 4) < 2) ? 0 : 1;
    }

    const int layers = 1;
    const auto overlap = computeOverlapCells(G, part, 2, layers);

    std::set<int> owned_union;
    for (int p = 0; p < 2; ++p) {
        std::vector<int> cells;
        for (int c = 0; c < G.number_of_cells; ++c) {
            if (part[c] == p) {
                cells.push_back(c);
            }
        }
        const auto numOwned = cells.size();
        owned_union.insert(cells.begin(), cells.end());
        cells.insert(cells.end(), overlap[p].begin(), overlap[p].end());

        auto local = extractLocalGrid(G, cells);
        const auto& lg = *local.grid;
        BOOST_REQUIRE_EQUAL(lg.number_of_cells, static_cast<int>(cells.size()));
        BOOST_CHECK(numOwned + overlap[p].size() == cells.size());

        // Every cell/face incidence is mirrored in face_cells, and every
        // local face is used by at least one local cell.
        std::vector<int> used(lg.number_of_faces, 0);
        for (int c = 0; c < lg.number_of_cells; ++c) {
            for (auto pos = lg.cell_facepos[c]; pos < lg.cell_facepos[c + 1]; ++pos) {
                const int f = lg.cell_faces[pos];
                BOOST_REQUIRE(f >= 0 && f < lg.number_of_faces);
                BOOST_CHECK(lg.face_cells[2 * f] == c || lg.face_cells[2 * f + 1] == c);
                ++used[f];
            }
        }
        for (int f = 0; f < lg.number_of_faces; ++f) {
            BOOST_CHECK(used[f] >= 1 && used[f] <= 2);
            // Global/local face mapping preserves area and normal.
            BOOST_CHECK_CLOSE(lg.face_areas[f], G.face_areas[local.globalFace[f]], 1e-12);
        }

        // Volumes are preserved.
        double vol = 0.0, expected = 0.0;
        for (int c = 0; c < lg.number_of_cells; ++c) {
            vol += lg.cell_volumes[c];
            expected += G.cell_volumes[cells[c]];
        }
        BOOST_CHECK_CLOSE(vol, expected, 1e-12);
    }
    BOOST_CHECK_EQUAL(static_cast<int>(owned_union.size()), G.number_of_cells);
}

BOOST_AUTO_TEST_CASE(ExtractLocalGridRejectsBadInput)
{
    GridPtr g(create_grid_cart3d(2, 1, 1));
    BOOST_CHECK_THROW(extractLocalGrid(*g, {0, 0}), std::invalid_argument);
    BOOST_CHECK_THROW(extractLocalGrid(*g, {5}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(PackUnpackRoundtrip)
{
    GridPtr g(create_grid_cart3d(3, 2, 2));
    const std::vector<int> cells{7, 3, 4, 8, 0};
    auto local = extractLocalGrid(*g, cells);
    const std::vector<char> owner{1, 1, 0, 1, 0};

    const auto packed = packGrid(*local.grid, owner, cells);
    const auto unpacked = unpackGrid(packed.ints, packed.doubles);
    const auto& a = *local.grid;
    const auto& b = *unpacked.grid;

    BOOST_CHECK(unpacked.isOwner == owner);
    BOOST_CHECK(unpacked.globalCells == cells);
    BOOST_REQUIRE_EQUAL(a.number_of_cells, b.number_of_cells);
    BOOST_REQUIRE_EQUAL(a.number_of_faces, b.number_of_faces);
    BOOST_REQUIRE_EQUAL(a.number_of_nodes, b.number_of_nodes);
    BOOST_CHECK(std::equal(a.cartdims, a.cartdims + 3, b.cartdims));
    BOOST_CHECK(std::equal(a.global_cell, a.global_cell + a.number_of_cells, b.global_cell));
    BOOST_CHECK(std::equal(a.face_cells, a.face_cells + 2 * a.number_of_faces, b.face_cells));
    BOOST_CHECK(std::equal(a.cell_facepos, a.cell_facepos + a.number_of_cells + 1, b.cell_facepos));
    BOOST_CHECK(std::equal(a.cell_faces, a.cell_faces + a.cell_facepos[a.number_of_cells], b.cell_faces));
    BOOST_CHECK(std::equal(a.cell_facetag, a.cell_facetag + a.cell_facepos[a.number_of_cells], b.cell_facetag));
    BOOST_CHECK(std::equal(a.face_nodepos, a.face_nodepos + a.number_of_faces + 1, b.face_nodepos));
    BOOST_CHECK(std::equal(a.face_nodes, a.face_nodes + a.face_nodepos[a.number_of_faces], b.face_nodes));
    BOOST_CHECK(std::equal(a.node_coordinates, a.node_coordinates + 3 * a.number_of_nodes, b.node_coordinates));
    BOOST_CHECK(std::equal(a.face_areas, a.face_areas + a.number_of_faces, b.face_areas));
    BOOST_CHECK(std::equal(a.cell_volumes, a.cell_volumes + a.number_of_cells, b.cell_volumes));
    BOOST_CHECK(std::equal(a.cell_centroids, a.cell_centroids + 3 * a.number_of_cells, b.cell_centroids));

    // Corrupt data is rejected.
    auto broken = packed.ints;
    broken.pop_back();
    BOOST_CHECK_THROW(unpackGrid(broken, packed.doubles), std::invalid_argument);
}
