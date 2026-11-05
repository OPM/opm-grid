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

#define BOOST_TEST_MODULE PolyhedralGridFaceGeometryTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/utility/platform_dependent/disable_warnings.h>
#include <dune/common/parallel/mpihelper.hh>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/polyhedralgrid.hh>
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>

#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>

#include "CornerPointTestGrid.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

// The face centers and the area normals that the transmissibilities of the CpGrid are computed
// from, are compared with the ones of the PolyhedralGrid, for the same grid.

using PolyGrid = Dune::PolyhedralGrid<3, 3>;
using Opm::TestGrids::CornerPointGrid;
using Opm::TestGrids::layersWithShift;

namespace
{

struct MpiFixture
{
    MpiFixture()
    {
        int argc = boost::unit_test::framework::master_test_suite().argc;
        char** argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(argc, argv);
    }
};

using Vec = std::array<double, 3>;

// What is known about one side of a cell, or of a face of the cell.
struct SideInfo
{
    int side;
    int neighbor;       // the cell on the other side of the face, or -1
    Vec faceCenter;     // the center of the face, to find the same face in the other grid
    Vec sideCenter;     // the center of the side of the cell
    Vec areaNormal;     // the area normal of the face
};

template <class Vector>
Vec toArray(const Vector& v)
{
    return {v[0], v[1], v[2]};
}

double distance(const Vec& a, const Vec& b)
{
    return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
}

Vec negated(const Vec& a)
{
    return {-a[0], -a[1], -a[2]};
}

Vec collectEclCenter(const Dune::CpGrid& grid, const int cell, const int side, const Dune::cpgrid::Intersection& is)
{
    return toArray(grid.faceCenterEcl(cell, side, is));
}

Vec collectEclCenter(const PolyGrid& grid, const int cell, const int side, const PolyGrid::LeafIntersection& is)
{
    return toArray(grid.faceCenterEcl(cell, side, is));
}

// The index of the face of an intersection.
int faceIndex(const Dune::cpgrid::Intersection& is)
{
    return is.id();
}

int faceIndex(const PolyGrid::LeafIntersection& is)
{
    return is.impl().id();
}

template <class Grid>
std::vector<std::vector<SideInfo>> collect(const Grid& grid)
{
    const auto gridView = grid.leafGridView();
    std::vector<std::vector<SideInfo>> info(gridView.size(0));
    for (const auto& element : elements(gridView)) {
        const int cell = gridView.indexSet().index(element);
        for (const auto& is : intersections(gridView, element)) {
            info[cell].push_back({is.indexInInside(),
                                  is.neighbor() ? static_cast<int>(gridView.indexSet().index(is.outside())) : -1,
                                  toArray(is.geometry().center()),
                                  collectEclCenter(grid, cell, is.indexInInside(), is),
                                  toArray(grid.faceAreaNormalEcl(faceIndex(is)))});
        }
    }
    return info;
}

// The same grid as a CpGrid and as a PolyhedralGrid.
void compareWithCpGrid(const CornerPointGrid& cp)
{
    const Opm::EclipseGrid eclipseGrid(cp.dims, cp.coord, cp.zcorn);

    Dune::CpGrid cpGrid;
    cpGrid.processEclipseFormat(&eclipseGrid, nullptr, false, false, false);
    const PolyGrid polyGrid(eclipseGrid);

    const auto fromCp = collect(cpGrid);
    const auto fromPoly = collect(polyGrid);
    BOOST_REQUIRE_EQUAL(fromCp.size(), fromPoly.size());

    constexpr double tolerance = 1e-9;
    for (std::size_t cell = 0; cell < fromCp.size(); ++cell) {
        BOOST_REQUIRE_EQUAL(fromCp[cell].size(), fromPoly[cell].size());

        // Find each face of the polyhedral grid in the other grid, by its side, its neighbor and its center.
        std::vector<bool> used(fromCp[cell].size(), false);
        for (const auto& poly : fromPoly[cell]) {
            int best = -1;
            double bestDistance = 1e300;
            for (std::size_t i = 0; i < fromCp[cell].size(); ++i) {
                if (used[i] || fromCp[cell][i].side != poly.side || fromCp[cell][i].neighbor != poly.neighbor) {
                    continue;
                }
                const double d = distance(fromCp[cell][i].faceCenter, poly.faceCenter);
                if (d < bestDistance) {
                    bestDistance = d;
                    best = static_cast<int>(i);
                }
            }
            BOOST_REQUIRE_MESSAGE(best >= 0,
                                  "cell " << cell << ": no face on the side " << poly.side << " towards the cell "
                                  << poly.neighbor << " in the CpGrid");
            used[best] = true;
            const auto& ref = fromCp[cell][best];

            BOOST_CHECK_MESSAGE(distance(ref.sideCenter, poly.sideCenter) < tolerance,
                                "cell " << cell << " side " << poly.side << ": the center of the side is ("
                                << poly.sideCenter[0] << ", " << poly.sideCenter[1] << ", " << poly.sideCenter[2]
                                << ") but the CpGrid has (" << ref.sideCenter[0] << ", " << ref.sideCenter[1] << ", "
                                << ref.sideCenter[2] << ")");

            // The sign of the area normal does not matter.
            const double scale = 1.0 + distance(ref.areaNormal, {0.0, 0.0, 0.0});
            const double mismatch = std::min(distance(ref.areaNormal, poly.areaNormal),
                                             distance(negated(ref.areaNormal), poly.areaNormal));
            BOOST_CHECK_MESSAGE(mismatch < tolerance * scale,
                                "cell " << cell << " side " << poly.side << ": the area normal is ("
                                << poly.areaNormal[0] << ", " << poly.areaNormal[1] << ", " << poly.areaNormal[2]
                                << ") but the CpGrid has (" << ref.areaNormal[0] << ", " << ref.areaNormal[1] << ", "
                                << ref.areaNormal[2] << ")");
        }
    }
}

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(MpiFixture);

BOOST_AUTO_TEST_CASE(RegularGrid)
{
    compareWithCpGrid(CornerPointGrid(2, 2, 2, layersWithShift([](int, int, int, int, int) { return 0.0; })));
}

// The faces are not plane and the cells are skewed, but there are no faults.
BOOST_AUTO_TEST_CASE(WarpedGridWithPlaneAndNonPlaneFaces)
{
    compareWithCpGrid(CornerPointGrid(3, 3, 2, [](int i, int j, int k, int zs, int ix, int jy) {
        return k + zs + 0.2 * std::sin(1.3 * (i + ix) + 0.7 * (j + jy) + 0.9 * (k + zs));
    }));
}

BOOST_AUTO_TEST_CASE(FaultWithDividingLinesInTheMiddleOfTheEdges)
{
    compareWithCpGrid(CornerPointGrid(2, 1, 2, layersWithShift([](int i, int, int, int, int) { return i == 1 ? 0.5 : 0.0; })));
}

BOOST_AUTO_TEST_CASE(FaultWithDividingLinesEndingInCorners)
{
    compareWithCpGrid(CornerPointGrid(2, 1, 2, layersWithShift([](int i, int, int, int, int jy) { return i == 1 ? 0.5 * (1 - jy) : 0.0; })));
}

BOOST_AUTO_TEST_CASE(FaultsInBothDirections)
{
    compareWithCpGrid(CornerPointGrid(3, 3, 2, layersWithShift([](int i, int j, int, int ix, int jy) {
        return (i == 1 ? 0.35 * (1 + ix) : 0.0) + (j == 2 ? 0.45 * (1 + jy) : 0.0);
    })));
}

// Faults and skewed cells together.
BOOST_AUTO_TEST_CASE(WarpedGridWithFaults)
{
    compareWithCpGrid(CornerPointGrid(3, 3, 2, [](int i, int j, int k, int zs, int ix, int jy) {
        const double warp = 0.15 * std::sin(1.1 * (i + ix) + 0.8 * (j + jy) + 0.6 * (k + zs));
        return k + zs + warp + (i == 1 ? 0.35 * (1 + ix) : 0.0) + (j == 2 ? 0.45 * (1 + jy) : 0.0);
    }));
}

BOOST_AUTO_TEST_CASE(CollapsedEdge)
{
    compareWithCpGrid(CornerPointGrid(1, 1, 3, [](int, int, int k, int zs, int ix, int jy) {
        return (ix == 1 && jy == 0 && k + zs == 2) ? 1.0 : double(k + zs);
    }));
}

BOOST_AUTO_TEST_CASE(CollapsedSide)
{
    compareWithCpGrid(CornerPointGrid(1, 1, 1, [](int, int, int, int zs, int ix, int) {
        return ix == 1 ? 0.5 : double(zs);
    }));
}
