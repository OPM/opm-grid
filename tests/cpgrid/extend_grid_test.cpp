/*
  Copyright 2026 NORCE Norwegian Research Centre AS.

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

// Tests for CpGrid::extendGrid(): appending a second corner-point grid to a
// first one as a disjoint union, as needed for the two continua of a
// dual-porosity model.

#include <config.h>

#define BOOST_TEST_MODULE ExtendGridTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>

#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>

#include <dune/grid/common/mcmgmapper.hh>

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

namespace {

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

// Description of a corner-point grid built from explicit COORD and ZCORN.
struct GridSpec
{
    std::array<int, 3> dims{3, 2, 2};
    double dx{10.0};
    double dy{10.0};
    double dz{2.0};
    double x0{0.0};     // x position of the first pillar
    double top{0.0};    // depth of the top layer
    double shear{0.0};  // tilt of the pillars: dx per unit of depth
    double fault_throw{0.0}; // vertical offset of the cells with i >= NX/2
    std::vector<int> actnum{}; // all active when empty
};

Opm::EclipseGrid makeEclipseGrid(const GridSpec& spec)
{
    const auto [nx, ny, nz] = spec.dims;

    // Pillars: 6 values (top xyz, bottom xyz) for each of (nx+1)*(ny+1) lines.
    const double z_top = spec.top - 1.0;
    const double z_bot = spec.top + nz * spec.dz + spec.fault_throw + 1.0;
    std::vector<double> coord;
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            const double x = spec.x0 + i * spec.dx;
            const double y = j * spec.dy;
            coord.insert(coord.end(), {x + spec.shear * (z_top - spec.top), y, z_top,
                                       x + spec.shear * (z_bot - spec.top), y, z_bot});
        }
    }

    // Corner depths in Eclipse order: i fastest, then j, then top/bottom
    // face of the layer, then k.
    std::vector<double> zcorn;
    zcorn.reserve(8 * nx * ny * nz);
    for (int k = 0; k < nz; ++k) {
        for (int face = 0; face < 2; ++face) {
            for (int j = 0; j < ny; ++j) {
                for (int jj = 0; jj < 2; ++jj) {
                    for (int i = 0; i < nx; ++i) {
                        const double shift = (i >= nx / 2) ? spec.fault_throw : 0.0;
                        for (int ii = 0; ii < 2; ++ii) {
                            zcorn.push_back(spec.top + (k + face) * spec.dz + shift);
                        }
                    }
                }
            }
        }
    }

    return Opm::EclipseGrid(spec.dims, coord, zcorn,
                            spec.actnum.empty() ? nullptr : spec.actnum.data());
}

std::unique_ptr<Dune::CpGrid> makeCpGrid(const GridSpec& spec)
{
    auto grid = std::make_unique<Dune::CpGrid>();
    const auto eclGrid = makeEclipseGrid(spec);
    grid->processEclipseFormat(&eclGrid, nullptr, false, false, false);
    return grid;
}

// Cell (0,0,0) and a few others inactive.
std::vector<int> someInactive(const std::array<int, 3>& dims, const std::set<int>& inactive)
{
    std::vector<int> actnum(dims[0] * dims[1] * dims[2], 1);
    for (const int c : inactive) {
        actnum[c] = 0;
    }
    return actnum;
}

constexpr double tol = 1e-12;

void checkVector(const Dune::FieldVector<double, 3>& a, const Dune::FieldVector<double, 3>& b)
{
    for (int d = 0; d < 3; ++d) {
        BOOST_CHECK_CLOSE_FRACTION(a[d] + 1.0, b[d] + 1.0, tol);
    }
}

// The two parts used in most tests: different pillars, different activity,
// and a fault in each.
GridSpec firstSpec()
{
    GridSpec s;
    s.dims = {3, 2, 2};
    s.shear = 0.1;
    s.fault_throw = 0.5;
    s.actnum = someInactive(s.dims, {1, 7});
    return s;
}

GridSpec secondSpec()
{
    GridSpec s;
    s.dims = {3, 2, 3};
    s.dz = 3.0;
    s.top = 7.0;
    s.shear = -0.2;
    s.fault_throw = 1.0;
    s.actnum = someInactive(s.dims, {0, 4, 11, 12});
    return s;
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(CountsAndCartesianSize)
{
    auto grid = makeCpGrid(firstSpec());
    const auto second = makeCpGrid(secondSpec());
    const auto first = makeCpGrid(firstSpec());

    grid->extendGrid(*second);

    BOOST_CHECK_EQUAL(grid->numCells(), first->numCells() + second->numCells());
    BOOST_CHECK_EQUAL(grid->numFaces(), first->numFaces() + second->numFaces());
    BOOST_CHECK_EQUAL(grid->numVertices(), first->numVertices() + second->numVertices());

    const std::array<int, 3> expected_size{3, 2, 5};
    BOOST_CHECK(grid->logicalCartesianSize() == expected_size);

    // Appending must not modify the appended grid.
    const auto second_ref = makeCpGrid(secondSpec());
    BOOST_CHECK_EQUAL(second->numCells(), second_ref->numCells());
    BOOST_CHECK_EQUAL(second->numFaces(), second_ref->numFaces());
    BOOST_CHECK(second->logicalCartesianSize() == second_ref->logicalCartesianSize());
}

#if HAVE_MPI
BOOST_AUTO_TEST_CASE(CellIndexSetHoldsTheCellsOfBothGrids)
{
    // The index set of the cells is what the interfaces used to distribute a grid
    // are made from, so it has to know every cell of the joined grid.
    auto grid = makeCpGrid(firstSpec());
    const auto second = makeCpGrid(secondSpec());
    const int numFirst = grid->numCells();

    grid->extendGrid(*second);

    const auto& indexSet = grid->getCellIndexSet();
    BOOST_CHECK_EQUAL(indexSet.size(), static_cast<std::size_t>(grid->numCells()));
    BOOST_CHECK_EQUAL(indexSet.size(), static_cast<std::size_t>(numFirst + second->numCells()));
}
#endif

BOOST_AUTO_TEST_CASE(CartesianIndexOfAppendedCellsUsesLogicalSize)
{
    auto grid = makeCpGrid(firstSpec());
    const auto first = makeCpGrid(firstSpec());
    const auto second = makeCpGrid(secondSpec());
    grid->extendGrid(*second);

    const int n1 = first->numCells();
    const int logical_size_first = 3 * 2 * 2;
    // Fewer active cells than logical cells, so shifting by the active
    // count would give the wrong index.
    BOOST_REQUIRE(n1 < logical_size_first);

    const auto& global = grid->globalCell();
    BOOST_REQUIRE_EQUAL(static_cast<int>(global.size()), grid->numCells());
    for (int c = 0; c < n1; ++c) {
        BOOST_CHECK_EQUAL(global[c], first->globalCell()[c]);
    }
    for (int c = 0; c < second->numCells(); ++c) {
        BOOST_CHECK_EQUAL(global[n1 + c], second->globalCell()[c] + logical_size_first);
    }
    BOOST_CHECK(std::is_sorted(global.begin(), global.end()));
    BOOST_CHECK(std::adjacent_find(global.begin(), global.end()) == global.end());
}

BOOST_AUTO_TEST_CASE(TwinCellsHaveCartesianIndexOffsetByHalf)
{
    // Equal parts with the same activity: the usual dual-continuum layout.
    GridSpec spec = firstSpec();
    auto grid = makeCpGrid(spec);
    const auto second = makeCpGrid(spec);
    const int n = grid->numCells();
    grid->extendGrid(*second);

    const int half = grid->logicalCartesianSize()[0] * grid->logicalCartesianSize()[1]
        * grid->logicalCartesianSize()[2] / 2;
    for (int c = 0; c < n; ++c) {
        BOOST_CHECK_EQUAL(grid->globalCell()[n + c], grid->globalCell()[c] + half);
    }
}

BOOST_AUTO_TEST_CASE(AppendedPartMatchesStandaloneGrid)
{
    auto grid = makeCpGrid(firstSpec());
    const auto first = makeCpGrid(firstSpec());
    const auto second = makeCpGrid(secondSpec());
    auto copy_of_second = makeCpGrid(secondSpec());
    grid->extendGrid(*second);

    const int nc1 = first->numCells();
    const int nf1 = first->numFaces();
    const int nv1 = first->numVertices();

    // Compare the parts of the combined grid with the standalone grids.
    const auto check_part = [&](const Dune::CpGrid& ref, const int cell_offset,
                                const int face_offset, const int vertex_offset)
    {
        for (int c = 0; c < ref.numCells(); ++c) {
            const int cc = cell_offset + c;
            BOOST_CHECK_CLOSE_FRACTION(grid->cellVolume(cc), ref.cellVolume(c), tol);
            checkVector(grid->cellCentroid(cc), ref.cellCentroid(c));

            BOOST_REQUIRE_EQUAL(grid->numCellFaces(cc), ref.numCellFaces(c));
            for (int k = 0; k < ref.numCellFaces(c); ++k) {
                BOOST_CHECK_EQUAL(grid->cellFace(cc, k), face_offset + ref.cellFace(c, k));
            }
        }

        for (int f = 0; f < ref.numFaces(); ++f) {
            const int ff = face_offset + f;
            BOOST_CHECK_CLOSE_FRACTION(grid->faceArea(ff), ref.faceArea(f), tol);
            checkVector(grid->faceCentroid(ff), ref.faceCentroid(f));
            checkVector(grid->faceNormal(ff), ref.faceNormal(f));

            for (int k = 0; k < 2; ++k) {
                const int ref_cell = ref.faceCell(f, k);
                const int cell = grid->faceCell(ff, k);
                BOOST_CHECK_EQUAL(cell, ref_cell < 0 ? -1 : cell_offset + ref_cell);
            }

            BOOST_REQUIRE_EQUAL(grid->numFaceVertices(ff), ref.numFaceVertices(f));
            for (int k = 0; k < ref.numFaceVertices(f); ++k) {
                BOOST_CHECK_EQUAL(grid->faceVertex(ff, k), vertex_offset + ref.faceVertex(f, k));
            }

            // Boundary ids are derived from the face tags (1-6 for i-, i+, ...).
            BOOST_CHECK_EQUAL(grid->boundaryId(ff), ref.boundaryId(f));
        }

        for (int v = 0; v < ref.numVertices(); ++v) {
            checkVector(grid->vertexPosition(vertex_offset + v), ref.vertexPosition(v));
        }
    };

    check_part(*first, 0, 0, 0);
    check_part(*copy_of_second, nc1, nf1, nv1);
}

BOOST_AUTO_TEST_CASE(CellCornersMatchStandaloneGrid)
{
    // The cell geometries refer to the corners by pointer; check that they
    // still give the right corners after the tables have been extended.
    auto grid = makeCpGrid(firstSpec());
    const auto first = makeCpGrid(firstSpec());
    const auto second = makeCpGrid(secondSpec());
    grid->extendGrid(*second);

    const auto corners = [](const Dune::CpGrid& g)
    {
        std::vector<std::array<Dune::FieldVector<double, 3>, 8>> all;
        for (const auto& element : Dune::elements(g.leafGridView())) {
            std::array<Dune::FieldVector<double, 3>, 8> c;
            for (int i = 0; i < 8; ++i) {
                c[i] = element.geometry().corner(i);
            }
            all.push_back(c);
        }
        return all;
    };

    const auto combined = corners(*grid);
    const auto ref1 = corners(*first);
    const auto ref2 = corners(*second);
    BOOST_REQUIRE_EQUAL(combined.size(), ref1.size() + ref2.size());

    for (std::size_t c = 0; c < combined.size(); ++c) {
        const auto& ref = (c < ref1.size()) ? ref1[c] : ref2[c - ref1.size()];
        for (int i = 0; i < 8; ++i) {
            checkVector(combined[c][i], ref[i]);
        }
    }

    // Centre of the local coordinates must agree with the stored centroid
    // up to the usual difference between the two definitions.
    for (const auto& element : Dune::elements(grid->leafGridView())) {
        const auto& geometry = element.geometry();
        BOOST_CHECK_CLOSE_FRACTION(geometry.volume(), grid->cellVolume(element.index()), tol);
        checkVector(geometry.center(), grid->cellCentroid(element.index()));
    }
}

BOOST_AUTO_TEST_CASE(NoFaceConnectsTheTwoParts)
{
    auto grid = makeCpGrid(firstSpec());
    const int nc1 = grid->numCells();
    const auto second = makeCpGrid(secondSpec());
    grid->extendGrid(*second);

    for (int f = 0; f < grid->numFaces(); ++f) {
        const int a = grid->faceCell(f, 0);
        const int b = grid->faceCell(f, 1);
        if (a >= 0 && b >= 0) {
            BOOST_CHECK_EQUAL(a < nc1, b < nc1);
        }
    }

    // The same through the Dune interface.
    for (const auto& element : Dune::elements(grid->leafGridView())) {
        for (const auto& is : Dune::intersections(grid->leafGridView(), element)) {
            if (is.neighbor()) {
                BOOST_CHECK_EQUAL(element.index() < nc1, is.outside().index() < nc1);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(IntersectionsMatchStandaloneGrid)
{
    auto grid = makeCpGrid(firstSpec());
    const int nc1 = grid->numCells();
    const auto second = makeCpGrid(secondSpec());
    const auto ref = makeCpGrid(secondSpec());
    grid->extendGrid(*second);

    const auto elements = [](const Dune::CpGrid& g)
    {
        std::vector<Dune::CpGrid::LeafGridView::Codim<0>::Entity> v;
        for (const auto& e : Dune::elements(g.leafGridView())) {
            v.push_back(e);
        }
        return v;
    };
    const auto combined_elements = elements(*grid);
    const auto ref_elements = elements(*ref);
    BOOST_REQUIRE_EQUAL(combined_elements.size(), static_cast<std::size_t>(nc1 + ref->numCells()));

    for (int c = 0; c < ref->numCells(); ++c) {
        const auto& e = combined_elements[nc1 + c];
        const auto& r = ref_elements[c];

        std::vector<std::pair<int, std::pair<double, int>>> got, expected;
        for (const auto& is : Dune::intersections(grid->leafGridView(), e)) {
            got.push_back({is.indexInInside(),
                           {is.geometry().volume(), is.neighbor() ? is.outside().index() - nc1 : -1}});
        }
        for (const auto& is : Dune::intersections(ref->leafGridView(), r)) {
            expected.push_back({is.indexInInside(),
                                {is.geometry().volume(), is.neighbor() ? is.outside().index() : -1}});
        }
        BOOST_REQUIRE_EQUAL(got.size(), expected.size());
        for (std::size_t i = 0; i < got.size(); ++i) {
            BOOST_CHECK_EQUAL(got[i].first, expected[i].first);
            BOOST_CHECK_CLOSE_FRACTION(got[i].second.first, expected[i].second.first, tol);
            BOOST_CHECK_EQUAL(got[i].second.second, expected[i].second.second);
        }
    }
}

BOOST_AUTO_TEST_CASE(IdsAreUniqueAndBoundaryIdsConsistent)
{
    auto grid = makeCpGrid(firstSpec());
    const auto second = makeCpGrid(secondSpec());
    grid->extendGrid(*second);

    const auto view = grid->leafGridView();
    std::set<Dune::CpGrid::GlobalIdSet::IdType> cell_ids, vertex_ids;
    for (const auto& element : Dune::elements(view)) {
        cell_ids.insert(grid->globalIdSet().id(element));
        for (int i = 0; i < 8; ++i) {
            vertex_ids.insert(grid->globalIdSet().id(element.template subEntity<3>(i)));
        }
    }
    BOOST_CHECK_EQUAL(static_cast<int>(cell_ids.size()), grid->numCells());
    BOOST_CHECK_EQUAL(static_cast<int>(vertex_ids.size()), grid->numVertices());

    // With unique boundary ids each boundary face has its own, positive, id.
    grid->setUniqueBoundaryIds(true);
    std::set<int> boundary_ids;
    int num_boundary_faces = 0;
    for (int f = 0; f < grid->numFaces(); ++f) {
        if (grid->faceCell(f, 0) < 0 || grid->faceCell(f, 1) < 0) {
            ++num_boundary_faces;
            BOOST_CHECK(grid->boundaryId(f) > 0);
            boundary_ids.insert(grid->boundaryId(f));
        } else {
            BOOST_CHECK_EQUAL(grid->boundaryId(f), 0);
        }
    }
    BOOST_CHECK_EQUAL(static_cast<int>(boundary_ids.size()), num_boundary_faces);
}

BOOST_AUTO_TEST_CASE(TotalVolumeIsSumOfParts)
{
    auto grid = makeCpGrid(firstSpec());
    const auto first = makeCpGrid(firstSpec());
    const auto second = makeCpGrid(secondSpec());
    grid->extendGrid(*second);

    const auto total = [](const Dune::CpGrid& g)
    {
        double sum = 0.0;
        for (int c = 0; c < g.numCells(); ++c) {
            sum += g.cellVolume(c);
        }
        return sum;
    };
    BOOST_CHECK_CLOSE_FRACTION(total(*grid), total(*first) + total(*second), 1e-10);
}

BOOST_AUTO_TEST_CASE(AppendTwice)
{
    auto grid = makeCpGrid(firstSpec());
    const auto second = makeCpGrid(secondSpec());
    const auto third = makeCpGrid(firstSpec());
    const int expected_cells = grid->numCells() + second->numCells() + third->numCells();

    grid->extendGrid(*second);
    grid->extendGrid(*third);

    BOOST_CHECK_EQUAL(grid->numCells(), expected_cells);
    const std::array<int, 3> expected_size{3, 2, 7};
    BOOST_CHECK(grid->logicalCartesianSize() == expected_size);

    const int last_first_part_cell = grid->numCells() - third->numCells();
    BOOST_CHECK_EQUAL(grid->globalCell()[last_first_part_cell],
                      third->globalCell()[0] + 3 * 2 * (2 + 3));
}

BOOST_AUTO_TEST_CASE(RejectsIncompatibleGridsAndLeavesGridUnchanged)
{
    auto grid = makeCpGrid(firstSpec());
    const int cells = grid->numCells();
    const int faces = grid->numFaces();
    const int vertices = grid->numVertices();
    const auto size = grid->logicalCartesianSize();

    GridSpec wrong = firstSpec();
    wrong.dims = {4, 2, 2};
    wrong.actnum.clear();
    const auto other = makeCpGrid(wrong);

    BOOST_CHECK_THROW(grid->extendGrid(*other), std::invalid_argument);
    BOOST_CHECK_THROW(grid->extendGrid(*grid), std::invalid_argument);

    BOOST_CHECK_EQUAL(grid->numCells(), cells);
    BOOST_CHECK_EQUAL(grid->numFaces(), faces);
    BOOST_CHECK_EQUAL(grid->numVertices(), vertices);
    BOOST_CHECK(grid->logicalCartesianSize() == size);
}

BOOST_AUTO_TEST_CASE(AddNncConnectsCellsAcrossTheParts)
{
    // Two equal parts, cells 1 and 7 of each inactive.
    const GridSpec spec = firstSpec();
    auto grid = makeCpGrid(spec);
    grid->extendGrid(*makeCpGrid(spec));
    auto reference = makeCpGrid(spec);
    reference->extendGrid(*makeCpGrid(spec));

    const int n = reference->numCells() / 2;
    const int logical_half = 3 * 2 * 2;

    // The twins, with duplicates, a reversed pair, a pair of cells that share a
    // face, and pairs with an inactive cell: all but the first ten twins are ignored.
    std::vector<std::pair<int, int>> nnc;
    for (int g = 0; g < logical_half; ++g) {
        nnc.emplace_back(g, g + logical_half);
    }
    nnc.emplace_back(logical_half, 0);
    nnc.emplace_back(0, 0 + logical_half);
    nnc.emplace_back(0, 3);
    grid->addNnc(nnc);

    const int num_nnc = 10; // 12 pairs, minus those with the inactive cells 1 and 7
    BOOST_REQUIRE_EQUAL(grid->numFaces(), reference->numFaces() + num_nnc);
    BOOST_CHECK_EQUAL(grid->numCells(), reference->numCells());

    // The new faces come first, and connect each cell to its twin.
    for (int f = 0; f < num_nnc; ++f) {
        const int a = grid->faceCell(f, 0);
        const int b = grid->faceCell(f, 1);
        BOOST_REQUIRE(a >= 0 && a < n);
        BOOST_REQUIRE(b >= n);
        BOOST_CHECK_EQUAL(grid->globalCell()[b], grid->globalCell()[a] + logical_half);
        BOOST_CHECK_EQUAL(grid->numFaceVertices(f), 0);
        BOOST_CHECK_EQUAL(grid->faceArea(f), 1.0);
        BOOST_CHECK_EQUAL(grid->boundaryId(f), 0);
    }

    // The other faces are unchanged but for the new numbering.
    for (int f = 0; f < reference->numFaces(); ++f) {
        const int ff = f + num_nnc;
        BOOST_CHECK_CLOSE_FRACTION(grid->faceArea(ff), reference->faceArea(f), tol);
        checkVector(grid->faceCentroid(ff), reference->faceCentroid(f));
        checkVector(grid->faceNormal(ff), reference->faceNormal(f));
        BOOST_CHECK_EQUAL(grid->faceCell(ff, 0), reference->faceCell(f, 0));
        BOOST_CHECK_EQUAL(grid->faceCell(ff, 1), reference->faceCell(f, 1));
        BOOST_CHECK_EQUAL(grid->boundaryId(ff), reference->boundaryId(f));
        BOOST_REQUIRE_EQUAL(grid->numFaceVertices(ff), reference->numFaceVertices(f));
        for (int k = 0; k < reference->numFaceVertices(f); ++k) {
            BOOST_CHECK_EQUAL(grid->faceVertex(ff, k), reference->faceVertex(f, k));
        }
    }

    int cells_with_nnc = 0;
    for (int c = 0; c < grid->numCells(); ++c) {
        const int faces = grid->numCellFaces(c);
        const bool has_nnc = faces == reference->numCellFaces(c) + 1;
        BOOST_REQUIRE(has_nnc || faces == reference->numCellFaces(c));
        cells_with_nnc += has_nnc;

        // The NNC face is first, and the faces after it are the old ones in the old order.
        const int skip = has_nnc ? 1 : 0;
        if (has_nnc) {
            BOOST_CHECK(grid->cellFace(c, 0) < num_nnc);
        }
        for (int k = 0; k < reference->numCellFaces(c); ++k) {
            BOOST_CHECK_EQUAL(grid->cellFace(c, k + skip), reference->cellFace(c, k) + num_nnc);
        }
    }
    BOOST_CHECK_EQUAL(cells_with_nnc, 2 * num_nnc);

    // Through the Dune interface the NNC is an intersection with index -1.
    int nnc_intersections = 0;
    for (const auto& element : Dune::elements(grid->leafGridView())) {
        for (const auto& is : Dune::intersections(grid->leafGridView(), element)) {
            if (is.indexInInside() == -1) {
                ++nnc_intersections;
                BOOST_CHECK(is.neighbor());
                BOOST_CHECK_EQUAL((element.index() < n), (is.outside().index() >= n));
            }
        }
    }
    BOOST_CHECK_EQUAL(nnc_intersections, 2 * num_nnc);

    // Connections that exist are not added again.
    grid->addNnc(nnc);
    BOOST_CHECK_EQUAL(grid->numFaces(), reference->numFaces() + num_nnc);
}
