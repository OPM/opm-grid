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

#define BOOST_TEST_MODULE LevelGridAssemblerTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/cpgrid/refinement/conforming/LevelGridAssembler.hpp>
#include <opm/grid/cpgrid/refinement/GridStateWriter.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <tests/cpgrid/lgr/ConformingTestHelpers.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <utility>
#include <vector>

using namespace Opm::ConformingTest;

namespace
{

// Number of K faces between each pair of cells (k of the cell above, k of the cell
// below) in each column (i, j) of a grid, by the grid's own Cartesian indices.
std::map<std::array<int,4>, int> verticalConnections(const Dune::cpgrid::CpGridData& grid)
{
    const auto& faceToCell = Opm::Refinement::GridStateWriter::faceToCell(grid);
    const auto& faceTag = Opm::Refinement::GridStateWriter::faceTag(grid);
    std::map<std::array<int,4>, int> connections;
    for (int face = 0; face < faceToCell.size(); ++face) {
        const auto row = faceToCell[Dune::cpgrid::EntityRep<1>(face, true)];
        if ((row.size() != 2) || (faceTag.get(face) != K_FACE)) {
            continue;
        }
        std::array<int,3> a, b;
        grid.getIJK(row[0].index(), a);
        grid.getIJK(row[1].index(), b);
        BOOST_REQUIRE(a[0] == b[0] && a[1] == b[1]);
        ++connections[{a[0], a[1], std::min(a[2], b[2]), std::max(a[2], b[2])}];
    }
    return connections;
}

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_AUTO_TEST_CASE(levelGridWithParentRelations)
{
    // Distorted (but vertical-pillar) parent grid.
    auto depth = [](int i_, int j_, int k_) {
        const double x = cellOf(i_) + sideOf(i_);
        const double y = cellOf(j_) + sideOf(j_);
        return 2.0*(cellOf(k_) + sideOf(k_)) + 0.25*std::sin(x) + 0.15*std::cos(0.7*y);
    };
    auto parent = makeGrid({4, 3, 2}, depth);

    Dune::CpGrid parentGrid;
    auto rawParent = parent.raw();
    parentGrid.processEclipseFormat(rawParent, false);

    const auto req = box("LGR1", {2, 3, 2}, {1, 0, 0}, {3, 2, 2});

    auto& storage = parentGrid.currentData();
    auto level = Opm::Refinement::assembleBlockLevelGrid(*storage[0], parent.dims,
                                                         parent.coord.data(), parent.zcorn.data(),
                                                         nullptr, req, /* levelIndex = */ 1,
                                                         /* parentLevel = */ 0,
                                                         storage);
    storage.push_back(level);

    const int childrenPerParent = 2*3*2;
    BOOST_REQUIRE_EQUAL(level->size(0), 2*2*2*childrenPerParent);
    BOOST_CHECK_EQUAL(level->getGridIdx(), 1);
    BOOST_CHECK(level->logicalCartesianSize() == (std::array<int,3>{4, 6, 4}));

    // Every refined cell has a father in level zero; children of one father
    // partition its volume; geometryInFather() has reference volume 1/12.
    std::map<int, double> volumePerParent;
    std::map<int, int> countPerParent;
    for (int cell = 0; cell < level->size(0); ++cell) {
        const Dune::cpgrid::Entity<0> child(*level, cell, true);
        BOOST_REQUIRE(child.hasFather());
        const auto father = child.father();
        BOOST_CHECK_EQUAL(father.level(), 0);
        volumePerParent[father.index()] += child.geometry().volume();
        countPerParent[father.index()] += 1;

        const auto geomInFather = child.geometryInFather();
        BOOST_CHECK_CLOSE(geomInFather.volume(), 1.0/childrenPerParent, 1e-9);
    }

    BOOST_CHECK_EQUAL(volumePerParent.size(), 8u); // 2x2x2 parents in the box
    for (const auto& [parentIdx, childVolume] : volumePerParent) {
        const Dune::cpgrid::Entity<0> parentCell(*storage[0], parentIdx, true);
        BOOST_CHECK_CLOSE(childVolume, parentCell.geometry().volume(), 1e-8);
        BOOST_CHECK_EQUAL(countPerParent[parentIdx], childrenPerParent);
    }
}

BOOST_AUTO_TEST_CASE(inactiveParentsHaveNoChildren)
{
    auto parent = makeGrid({3, 2, 2}, [](int, int, int k_) {
        return static_cast<double>(cellOf(k_) + sideOf(k_));
    });
    parent.actnum.assign(3*2*2, 1);
    parent.actnum[4] = 0; // cell (1,1,0)

    Dune::CpGrid parentGrid;
    auto rawParent = parent.raw();
    parentGrid.processEclipseFormat(rawParent, false);

    const auto req = box("LGR1", {2, 2, 2}, {0, 0, 0}, {3, 2, 1});

    auto& storage = parentGrid.currentData();
    auto level = Opm::Refinement::assembleBlockLevelGrid(*storage[0], parent.dims,
                                                         parent.coord.data(), parent.zcorn.data(),
                                                         parent.actnum.data(), req, 1,
                                                         /* parentLevel = */ 0,
                                                         storage);
    storage.push_back(level);

    BOOST_CHECK_EQUAL(level->size(0), (3*2*1 - 1)*8);
    for (int cell = 0; cell < level->size(0); ++cell) {
        const Dune::cpgrid::Entity<0> child(*level, cell, true);
        BOOST_REQUIRE(child.hasFather());
    }
}

BOOST_AUTO_TEST_CASE(faultInsideBlockLevelGrid)
{
    // Fault between i-columns 1 and 2 (throw 0.6, layer thickness 2).
    auto depth = [](int i_, int j_, int k_) {
        (void)j_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + 0.6 : base;
    };
    auto parent = makeGrid({4, 2, 2}, depth);

    Dune::CpGrid parentGrid;
    auto rawParent = parent.raw();
    parentGrid.processEclipseFormat(rawParent, false);

    const auto req = box("LGR1", {2, 2, 2}, {1, 0, 0}, {3, 2, 2});

    auto& storage = parentGrid.currentData();
    auto level = Opm::Refinement::assembleBlockLevelGrid(*storage[0], parent.dims,
                                                         parent.coord.data(), parent.zcorn.data(),
                                                         nullptr, req, 1,
                                                         /* parentLevel = */ 0,
                                                         storage);
    storage.push_back(level);

    BOOST_REQUIRE_EQUAL(level->size(0), 2*2*2*8);

    // Volume is still partitioned correctly across the faulted block.
    std::map<int, double> volumePerParent;
    for (int cell = 0; cell < level->size(0); ++cell) {
        const Dune::cpgrid::Entity<0> child(*level, cell, true);
        volumePerParent[child.father().index()] += child.geometry().volume();
    }
    for (const auto& [parentIdx, childVolume] : volumePerParent) {
        const Dune::cpgrid::Entity<0> parentCell(*storage[0], parentIdx, true);
        BOOST_CHECK_CLOSE(childVolume, parentCell.geometry().volume(), 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(pinchConnectionCarriedToChildren)
{
    // 3x1x4 parent, unit columns; layers 2 thick, except layer k = 2, which is
    // collapsed (zero thickness) at depth 4.
    const std::array<double,5> nodeDepth = {0.0, 2.0, 4.0, 4.0, 6.0};
    auto parent = makeGrid({3, 1, 4}, [&nodeDepth](int, int, int k_) {
        return nodeDepth[cellOf(k_) + sideOf(k_)];
    });

    // Processed with pinch-outs active, as for a deck with PINCH: the collapsed
    // layer is removed, and layers k = 1 and k = 3 are joined by one face.
    Dune::CpGrid parentGrid;
    auto& storage = parentGrid.currentData();
    std::array<std::set<std::pair<int,int>>, 2> parentNnc;
    storage[0]->processEclipseFormat(parent.raw(), nullptr, parentNnc,
                                     /* remove_ij_boundary = */ false,
                                     /* turn_normals = */ false,
                                     /* pinchActive = */ true,
                                     /* tolerance_unique_points = */ 0.0,
                                     /* edge_conformal = */ false);
    BOOST_REQUIRE_EQUAL(storage[0]->size(0), 3*3);
    const auto parentConnections = verticalConnections(*storage[0]);
    for (int i = 0; i < 3; ++i) {
        BOOST_REQUIRE_EQUAL(parentConnections.count({i, 0, 1, 3}), 1u);
    }

    // The box covers columns i = 0, 1 and every layer; each parent becomes 2x1x2
    // children, so parent layer k holds refined layers 2k and 2k+1.
    const auto req = box("LGR1", {2, 1, 2}, {0, 0, 0}, {2, 1, 4});
    auto level = Opm::Refinement::assembleBlockLevelGrid(*storage[0], parent.dims,
                                                         parent.coord.data(), parent.zcorn.data(),
                                                         nullptr, req, /* levelIndex = */ 1,
                                                         /* parentLevel = */ 0,
                                                         storage);
    storage.push_back(level);

    // Refined layers 4 and 5 (children of the collapsed layer) are removed.
    BOOST_REQUIRE_EQUAL(level->size(0), 4*1*6);
    // In every refined column: ordinary neighbours in layers 0-3 and 6-7, and the
    // carried pinch connection between layer 3, just above the collapsed layer,
    // and layer 6, just below it. Each pair has one face, and cells one layer
    // apart (such as layers 1 and 3) have none.
    const std::map<std::array<int,2>, int> expectedPerColumn = {
        {{0, 1}, 1}, {{1, 2}, 1}, {{2, 3}, 1}, {{3, 6}, 1}, {{6, 7}, 1}};
    std::map<std::array<int,4>, int> expected;
    for (int ir = 0; ir < 4; ++ir) {
        for (const auto& [layers, count] : expectedPerColumn) {
            expected[{ir, 0, layers[0], layers[1]}] = count;
        }
    }
    const auto connections = verticalConnections(*level);
    BOOST_CHECK(connections == expected);
    for (int ir = 0; ir < 4; ++ir) {
        BOOST_CHECK_EQUAL(connections.count({ir, 0, 3, 6}), 1u);
        BOOST_CHECK_EQUAL(connections.count({ir, 0, 1, 3}), 0u);
    }
}
