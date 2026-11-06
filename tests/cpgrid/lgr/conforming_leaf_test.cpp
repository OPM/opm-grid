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

#define BOOST_TEST_MODULE ConformingLeafTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/refinement/conforming/ConformingBlockBuilder.hpp>
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>
#include <opm/grid/cpgrid/refinement/RefinementRequest.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <dune/common/parallel/mpihelper.hh>

#include <algorithm>
#include <array>
#include <functional>
#include <memory>
#include <string>
#include <vector>

// The leaf of the conforming refinement: cells whose children a box's MINPV removes,
// the corner order of a box on a left-handed grid, and the level-zero face an LGR
// boundary face comes from.

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

using Pillar = std::function<std::array<double,6>(int, int)>;

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
        g.coord = coord.data();
        g.zcorn = zcorn.data();
        g.actnum = actnum.empty() ? nullptr : actnum.data();
        return g;
    }
};

int cellOf(int doubled) { return doubled / 2; }
int sideOf(int doubled) { return doubled % 2; }
int nodeOf(int doubled) { return cellOf(doubled) + sideOf(doubled); }

std::array<double,6> vertical(int i, int j)
{
    return { double(i), double(j), 0.0, double(i), double(j), 100.0 };
}

std::array<double,6> leftHanded(int i, int j)
{
    return { double(i), -double(j), 0.0, double(i), -double(j), 100.0 };
}

TestGrdecl makeGrid(const std::array<int,3>& dims, const Pillar& pillar,
                    const std::function<double(int,int,int)>& depth)
{
    TestGrdecl g;
    g.dims = dims;
    const auto& [nx, ny, nz] = dims;
    g.coord.resize(6*(nx + 1)*(ny + 1));
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            const auto p = pillar(i, j);
            std::copy(p.begin(), p.end(), &g.coord[6*(static_cast<std::size_t>(j)*(nx + 1) + i)]);
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

Opm::Refinement::BlockRefinement box(std::array<int,3> cellsPerDim, std::array<int,3> start,
                                     std::array<int,3> end)
{
    Opm::Refinement::BlockRefinement request;
    request.name = "LGR1";
    request.cellsPerDim = cellsPerDim;
    request.startIJK = start;
    request.endIJK = end;
    return request;
}

void refine(Dune::CpGrid& grid, const TestGrdecl& g, Opm::Refinement::BlockRefinement request)
{
    grid.processEclipseFormat(g.raw(), false);
    grid.setLgrBackend(Opm::Refinement::Backend::Conforming);
    grid.setRefinementBuilder(std::make_shared<Opm::Refinement::ConformingBlockBuilder>(
        g.dims, g.coord, g.zcorn, g.actnum));
    grid.addLgrsUpdateLeafView({std::move(request)});
}

double worstClosure(const Dune::CpGrid& grid)
{
    double worst = 0.0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        Dune::FieldVector<double,3> sum(0.0);
        double area = 0.0;
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            auto n = is.centerUnitOuterNormal();
            n *= is.geometry().volume();
            sum += n;
            area += is.geometry().volume();
        }
        worst = std::max(worst, sum.two_norm() / area);
    }
    return worst;
}

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(MPIFixture);

// The box's MINPV removes both children of one host: the host is refined into nothing,
// and its coarse neighbour keeps its face toward it as a boundary face.
BOOST_AUTO_TEST_CASE(hostWithNoChildLeft)
{
    const auto g = makeGrid({4, 1, 1}, vertical, [](int, int, int k_) { return 2.0*nodeOf(k_); });
    auto request = box({2,1,1}, {1,0,0}, {3,1,1});
    request.minpvRemoved = {1, 1, 0, 0};                 // children of host i = 1
    Dune::CpGrid grid;
    refine(grid, g, request);

    BOOST_CHECK_EQUAL(grid.size(0), 1 + 2 + 1);
    BOOST_CHECK_SMALL(worstClosure(grid), 1e-12);
    for (const auto& element : Dune::elements(grid.levelGridView(0))) {
        const int i = grid.currentData().front()->globalCell()[element.index()];
        BOOST_CHECK_EQUAL(element.isLeaf(), (i == 0) || (i == 3));
    }
}

// A box one parent wide in J on a left-handed grid has a level grid one cell thick in J;
// its corners still take the level-zero corner order, so each leaf cell's corners
// surround it.
BOOST_AUTO_TEST_CASE(thinBoxOnLeftHandedGrid)
{
    const auto g = makeGrid({3, 2, 2}, leftHanded, [](int, int, int k_) { return 2.0*nodeOf(k_); });
    Dune::CpGrid grid;
    refine(grid, g, box({2,1,2}, {1,0,0}, {2,1,2}));

    for (const auto& element : Dune::elements(grid.leafGridView())) {
        const auto& geometry = element.geometry();
        Dune::FieldVector<double,3> mean(0.0);
        for (int corner = 0; corner < 8; ++corner) {
            mean += geometry.corner(corner);
        }
        mean /= 8.0;
        BOOST_CHECK_SMALL((mean - geometry.center()).two_norm(), 1e-12);
    }
    BOOST_CHECK_SMALL(worstClosure(grid), 1e-12);
}

// Each LGR boundary face comes from the level-zero face between the two cells' origins
// that holds it.
BOOST_AUTO_TEST_CASE(parentFaceHoldsTheLgrBoundaryFace)
{
    const auto g = makeGrid({4, 2, 2}, vertical, [](int i_, int, int k_) {
        return 2.0*nodeOf(k_) + ((cellOf(i_) >= 2) ? 0.6 : 0.0);
    });
    Dune::CpGrid grid;
    refine(grid, g, box({2,2,2}, {2,0,0}, {4,2,2}));

    int checked = 0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (!is.neighbor() || (element.level() == is.outside().level())) {
                continue;
            }
            const auto parent = grid.getParentIntersectionFromLgrBoundaryFace(is);
            BOOST_CHECK_EQUAL(parent.inside().index(), element.getOrigin().index());
            BOOST_CHECK_EQUAL(parent.outside().index(), is.outside().getOrigin().index());
            // The leaf face lies on the parent face: its centre is within the parent's area
            // around the parent's centre along the fault plane.
            const auto offset = is.geometry().center() - parent.geometry().center();
            BOOST_CHECK_SMALL(std::abs(offset * parent.centerUnitOuterNormal()), 1e-12);
            ++checked;
        }
    }
    BOOST_CHECK_GT(checked, 0);
}
