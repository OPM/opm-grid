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

#define BOOST_TEST_MODULE LgrBackendComparisonTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>
#include <opm/grid/cpgrid/refinement/conforming/ConformingBlockBuilder.hpp>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

struct MPIFixture
{
    MPIFixture()
    {
        auto& argc = boost::unit_test::framework::master_test_suite().argc;
        auto& argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(argc, argv);
    }
};

BOOST_GLOBAL_FIXTURE(MPIFixture);

namespace
{

using Vec = Dune::FieldVector<double, 3>;

struct Grdecl
{
    std::array<int,3> dims;
    std::vector<double> coord;
    std::vector<double> zcorn;
    std::vector<int> actnum;
};

// Unit cells on vertical pillars; depth(i, j, k) gives the corner-point depth of node (i, j, k).
Grdecl verticalPillarGrid(const std::array<int,3>& dims, const std::function<double(int,int,int)>& depth)
{
    const auto [nx, ny, nz] = dims;
    Grdecl g{dims, {}, std::vector<double>(8 * nx * ny * nz), std::vector<int>(nx * ny * nz, 1)};
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            g.coord.insert(g.coord.end(), {double(i), double(j), 0.0, double(i), double(j), 1.0});
        }
    }
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                for (int c = 0; c < 8; ++c) {
                    const int di = c & 1, dj = (c >> 1) & 1, dk = (c >> 2) & 1;
                    const std::size_t idx = (2*i + di) + 2*nx*(2*j + dj) + 4*nx*ny*(2*k + dk);
                    g.zcorn[idx] = depth(i + di, j + dj, k + dk);
                }
            }
        }
    }
    return g;
}

struct Request
{
    std::array<int,3> cellsPerDim;
    std::array<int,3> start;
    std::array<int,3> end;
};

void refine(Dune::CpGrid& grid, const Grdecl& g, const std::vector<Request>& requests, bool conforming)
{
    grdecl raw;
    raw.dims[0] = g.dims[0]; raw.dims[1] = g.dims[1]; raw.dims[2] = g.dims[2];
    raw.coord = g.coord.data();
    raw.zcorn = const_cast<double*>(g.zcorn.data());
    raw.actnum = g.actnum.data();
    grid.processEclipseFormat(raw, false);
    if (conforming) {
        grid.setLgrBackend(Opm::Refinement::Backend::Conforming);
        grid.setRefinementBuilder(std::make_shared<Opm::Refinement::ConformingBlockBuilder>(
            g.dims, g.coord, g.zcorn, g.actnum));
    }
    std::vector<std::array<int,3>> cellsPerDim, start, end;
    std::vector<std::string> names;
    for (const auto& r : requests) {
        cellsPerDim.push_back(r.cellsPerDim);
        start.push_back(r.start);
        end.push_back(r.end);
        names.push_back("LGR" + std::to_string(names.size() + 1));
    }
    grid.addLgrsUpdateLeafView(cellsPerDim, start, end, names);
}

// A leaf cell as both backends name it: its level-zero Cartesian cell and its place in that parent.
using Key = std::pair<int,int>;

struct Face
{
    double area = 0.0;
    Vec weightedCentre = Vec(0.0);
};

struct Cell
{
    std::array<Vec,8> corners{};
    double volume = 0.0;
    Vec centre = Vec(0.0);
    double boundaryArea = 0.0;
    std::map<Key, Face> neighbours;
};

std::map<Key, Cell> describe(const Dune::CpGrid& grid)
{
    const auto& level0Cartesian = grid.currentData().front()->globalCell();
    const auto key = [&](const auto& element) {
        return Key{level0Cartesian[element.getOrigin().index()],
                   element.hasFather() ? element.getIdxInParentCell() : -1};
    };
    std::map<Key, Cell> cells;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        auto& cell = cells[key(element)];
        for (int c = 0; c < 8; ++c) {
            cell.corners[c] = element.geometry().corner(c);
        }
        cell.volume = element.geometry().volume();
        cell.centre = element.geometry().center();
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            const double area = is.geometry().volume();
            if (!is.neighbor()) {
                cell.boundaryArea += area;
                continue;
            }
            auto& face = cell.neighbours[key(is.outside())];
            face.area += area;
            auto centre = is.geometry().center();
            centre *= area;
            face.weightedCentre += centre;
        }
    }
    return cells;
}

void checkClose(double a, double b, double scale, const std::string& what)
{
    BOOST_CHECK_MESSAGE(std::abs(a - b) <= 1e-10 * scale, what << ": " << a << " vs " << b);
}

// Both backends place the same nodes. A cell with non-planar faces gets its volume and centroids
// from different formulas, so planarFaces = false compares those to a tolerance only.
void compareBackends(const Grdecl& g, const std::vector<Request>& requests, bool planarFaces = true)
{
    Dune::CpGrid trilinear, conforming;
    refine(trilinear, g, requests, false);
    refine(conforming, g, requests, true);
    BOOST_REQUIRE_EQUAL(trilinear.size(0), conforming.size(0));
    BOOST_CHECK_EQUAL(trilinear.maxLevel(), conforming.maxLevel());

    const auto a = describe(trilinear);
    const auto b = describe(conforming);
    BOOST_REQUIRE_EQUAL(a.size(), b.size());
    double volumeA = 0.0, volumeB = 0.0;
    for (const auto& [key, ca] : a) {
        const auto it = b.find(key);
        BOOST_REQUIRE_MESSAGE(it != b.end(), "cell " << key.first << "/" << key.second << " missing");
        const auto& cb = it->second;
        const std::string name = "cell " + std::to_string(key.first) + "/" + std::to_string(key.second);
        for (int c = 0; c < 8; ++c) {
            for (int d = 0; d < 3; ++d) {
                checkClose(ca.corners[c][d], cb.corners[c][d], 1.0, name + " corner");
            }
        }
        volumeA += ca.volume;
        volumeB += cb.volume;
        if (planarFaces) {
            checkClose(ca.volume, cb.volume, 1.0, name + " volume");
            for (int d = 0; d < 3; ++d) {
                checkClose(ca.centre[d], cb.centre[d], 1.0, name + " centre");
            }
        }
        else {
            BOOST_CHECK_CLOSE(ca.volume, cb.volume, 0.2);
        }
        checkClose(ca.boundaryArea, cb.boundaryArea, 1.0, name + " boundary area");
        BOOST_REQUIRE_EQUAL(ca.neighbours.size(), cb.neighbours.size());
        for (const auto& [nb, fa] : ca.neighbours) {
            const auto fb = cb.neighbours.find(nb);
            BOOST_REQUIRE_MESSAGE(fb != cb.neighbours.end(), name << " lacks neighbour " << nb.first << "/" << nb.second);
            checkClose(fa.area, fb->second.area, 1.0, name + " face area");
            for (int d = 0; d < 3 && planarFaces; ++d) {
                checkClose(fa.weightedCentre[d], fb->second.weightedCentre[d], 1.0, name + " face centre");
            }
        }
    }
    checkClose(volumeA, volumeB, volumeA, "total volume");
}

} // anonymous namespace

// On vertical pillars and unfaulted parents both backends must build the same refined grid.
BOOST_AUTO_TEST_CASE(singleBoxFlatLayers)
{
    const auto g = verticalPillarGrid({4, 3, 3}, [](int, int, int k) { return double(k); });
    compareBackends(g, {{{2, 2, 2}, {1, 1, 1}, {3, 2, 2}}});
}

BOOST_AUTO_TEST_CASE(anisotropicBoxTiltedLayers)
{
    const auto g = verticalPillarGrid({5, 4, 3}, [](int i, int j, int k) { return 0.1*i + 0.05*j + k; });
    compareBackends(g, {{{3, 2, 1}, {1, 1, 0}, {4, 3, 2}}});
}

BOOST_AUTO_TEST_CASE(boxOnGridBoundaryCurvedLayers)
{
    const auto g = verticalPillarGrid({4, 4, 3}, [](int i, int j, int k) {
        return k * (1.0 + 0.02*i*j) + 0.03*i*i;
    });
    compareBackends(g, {{{2, 3, 2}, {0, 0, 0}, {2, 2, 3}}}, /* planarFaces = */ false);
}

BOOST_AUTO_TEST_CASE(twoSeparatedBoxes)
{
    const auto g = verticalPillarGrid({7, 4, 2}, [](int i, int, int k) { return 0.05*i + k; });
    compareBackends(g, {{{2, 2, 2}, {0, 0, 0}, {2, 2, 1}},
                        {{3, 3, 1}, {4, 1, 0}, {6, 3, 2}}});
}
