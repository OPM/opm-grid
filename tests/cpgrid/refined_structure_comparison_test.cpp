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

  Structural-equivalence tests for the refinement builder: the refined
  result must match a grid built *directly* from the refined corner-point
  description (the upstream LgrChecks methodology). Faulted parents are
  the interesting case: the directly processed grid is the oracle for how
  the fault must look after refinement.
*/
#include <config.h>

#define BOOST_TEST_MODULE RefinedStructureComparisonTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/refinement/conforming/ConformingBlockBuilder.hpp>
#include <opm/grid/cpgrid/refinement/conforming/GrdeclRefinement.hpp>
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <dune/common/parallel/mpihelper.hh>

#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
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

TestGrdecl makeVerticalPillarGrid(const std::array<int,3>& dims,
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

// Faulted, distorted, vertical-pillar test grid: fault plane between
// i-columns (faultI-1 | faultI) with the given throw.
TestGrdecl makeFaultedGrid(const std::array<int,3>& dims, int faultI, double throwZ)
{
    return makeVerticalPillarGrid(dims, [faultI, throwZ](int i_, int j_, int k_) {
        const double y = cellOf(j_) + sideOf(j_);
        const double base = 2.0*(cellOf(k_) + sideOf(k_)) + 0.15*std::cos(0.7*y);
        return (cellOf(i_) >= faultI) ? base + throwZ : base;
    });
}

void useConformingBuilder(Dune::CpGrid& grid, std::unique_ptr<Opm::Refinement::Builder> builder)
{
    grid.setLgrBackend(Opm::Refinement::Backend::Conforming);
    grid.setRefinementBuilder(std::move(builder));
}

// Map from refined-local Cartesian index to (volume, centroid) of a grid
// built directly from a refined description.
std::map<int, std::pair<double, Dune::FieldVector<double,3>>>
cellGeometryByCartesian(const Dune::CpGrid& grid)
{
    std::map<int, std::pair<double, Dune::FieldVector<double,3>>> result;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        result[grid.globalCell()[element.index()]] =
            { element.geometry().volume(), element.geometry().center() };
    }
    return result;
}

// Count interior (two-sided) faces of a grid view via intersections.
template <class GridView>
int countInteriorFacePairs(const GridView& view)
{
    int twice = 0;
    for (const auto& element : Dune::elements(view)) {
        for (const auto& intersection : Dune::intersections(view, element)) {
            if (intersection.neighbor()) {
                ++twice;
            }
        }
    }
    BOOST_REQUIRE_EQUAL(twice % 2, 0);
    return twice / 2;
}

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(MPIFixture);

// Whole-grid box on a faulted grid: the leaf must be structurally
// identical to the directly processed refined description.
BOOST_AUTO_TEST_CASE(globalBoxMatchesDirectlyRefinedFaultedGrid)
{
    auto parent = makeFaultedGrid({4, 2, 2}, 2, 0.6);

    // A: refinement through the builder, one box covering everything.
    Dune::CpGrid gridA;
    auto rawParent = parent.raw();
    gridA.processEclipseFormat(rawParent, false);

    useConformingBuilder(gridA, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    Opm::Refinement::BlockRefinement req;
    req.name = "GR";
    req.cellsPerDim = {2, 2, 2};
    req.startIJK = {0, 0, 0};
    req.endIJK = parent.dims;
    gridA.addLgrsUpdateLeafView({req.cellsPerDim}, {req.startIJK}, {req.endIJK}, {req.name});

    // B: the oracle - process the refined description directly.
    const auto refined = Opm::Refinement::refineBlock(parent.dims, parent.coord.data(),
                                                      parent.zcorn.data(), nullptr, req);
    Dune::CpGrid gridB;
    const auto rawRefined = grdecl{ {refined.dims[0], refined.dims[1], refined.dims[2]},
                                    refined.coord.data(), refined.zcorn.data(),
                                    refined.actnum.data() };
    gridB.processEclipseFormat(rawRefined, false);

    BOOST_REQUIRE_EQUAL(gridA.size(0), gridB.size(0));
    BOOST_CHECK_EQUAL(gridA.size(3), gridB.size(3));
    BOOST_CHECK_EQUAL(countInteriorFacePairs(gridA.leafGridView()),
                      countInteriorFacePairs(gridB.leafGridView()));

    // Per-cell match through the refined-local Cartesian index: every leaf
    // cell of A corresponds to exactly one cell of B with the same volume
    // and centroid - including the cells along the refined fault.
    const auto oracle = cellGeometryByCartesian(gridB);
    int matched = 0;
    for (const auto& element : Dune::elements(gridA.leafGridView())) {
        BOOST_REQUIRE(element.hasFather());
        const int levelIdx = element.getLevelElem().index();
        const int refinedCart = gridA.currentData()[1]->globalCell()[levelIdx];
        const auto it = oracle.find(refinedCart);
        BOOST_REQUIRE(it != oracle.end());
        BOOST_CHECK_CLOSE(element.geometry().volume(), it->second.first, 1e-10);
        for (int c = 0; c < 3; ++c) {
            BOOST_CHECK_SMALL(element.geometry().center()[c] - it->second.second[c], 1e-10);
        }
        ++matched;
    }
    BOOST_CHECK_EQUAL(matched, gridB.size(0));
}

// Local box containing a fault: the refined region of the leaf must match
// the directly processed block description cell by cell, and the interior
// connectivity of the refined region (including the refined fault
// connections) must agree with the oracle's.
BOOST_AUTO_TEST_CASE(localFaultedBoxMatchesDirectBlockGrid)
{
    auto parent = makeFaultedGrid({5, 3, 2}, 2, 0.6);

    Dune::CpGrid gridA;
    auto rawParent = parent.raw();
    gridA.processEclipseFormat(rawParent, false);

    useConformingBuilder(gridA, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        parent.dims, parent.coord, parent.zcorn, parent.actnum));

    // Box {1..4} x {0..3} x {0..2}: fault plane i=2 strictly inside.
    Opm::Refinement::BlockRefinement req;
    req.name = "LGR1";
    req.cellsPerDim = {3, 2, 2};
    req.startIJK = {1, 0, 0};
    req.endIJK = {4, 3, 2};
    gridA.addLgrsUpdateLeafView({req.cellsPerDim}, {req.startIJK}, {req.endIJK}, {req.name});

    const auto refined = Opm::Refinement::refineBlock(parent.dims, parent.coord.data(),
                                                      parent.zcorn.data(), nullptr, req);
    Dune::CpGrid gridB;
    const auto rawRefined = grdecl{ {refined.dims[0], refined.dims[1], refined.dims[2]},
                                    refined.coord.data(), refined.zcorn.data(),
                                    refined.actnum.data() };
    gridB.processEclipseFormat(rawRefined, false);

    // Cell-by-cell match of the refined region.
    const auto oracle = cellGeometryByCartesian(gridB);
    int refinedCells = 0;
    for (const auto& element : Dune::elements(gridA.leafGridView())) {
        if (!element.hasFather()) {
            continue;
        }
        ++refinedCells;
        const int levelIdx = element.getLevelElem().index();
        const int refinedCart = gridA.currentData()[1]->globalCell()[levelIdx];
        const auto it = oracle.find(refinedCart);
        BOOST_REQUIRE(it != oracle.end());
        BOOST_CHECK_CLOSE(element.geometry().volume(), it->second.first, 1e-10);
        for (int c = 0; c < 3; ++c) {
            BOOST_CHECK_SMALL(element.geometry().center()[c] - it->second.second[c], 1e-10);
        }
    }
    BOOST_REQUIRE_EQUAL(refinedCells, gridB.size(0));

    // Interior connectivity of the refined region: refined-refined leaf
    // face pairs (including refined fault connections) must equal the
    // oracle's interior face pairs.
    int refinedPairsTwice = 0;
    for (const auto& element : Dune::elements(gridA.leafGridView())) {
        if (!element.hasFather()) {
            continue;
        }
        for (const auto& intersection : Dune::intersections(gridA.leafGridView(), element)) {
            if (intersection.neighbor() && intersection.outside().hasFather()) {
                ++refinedPairsTwice;
            }
        }
    }
    BOOST_REQUIRE_EQUAL(refinedPairsTwice % 2, 0);
    BOOST_CHECK_EQUAL(refinedPairsTwice / 2, countInteriorFacePairs(gridB.leafGridView()));

    // The refined fault is actually there: more interior connections than
    // an unfaulted twin of the same box would have.
    auto parentNoFault = makeFaultedGrid({5, 3, 2}, 2, 0.0);
    const auto refinedNoFault = Opm::Refinement::refineBlock(parentNoFault.dims,
                                                             parentNoFault.coord.data(),
                                                             parentNoFault.zcorn.data(),
                                                             nullptr, req);
    Dune::CpGrid gridC;
    const auto rawNoFault = grdecl{ {refinedNoFault.dims[0], refinedNoFault.dims[1], refinedNoFault.dims[2]},
                                    refinedNoFault.coord.data(), refinedNoFault.zcorn.data(),
                                    refinedNoFault.actnum.data() };
    gridC.processEclipseFormat(rawNoFault, false);
    BOOST_CHECK_GT(countInteriorFacePairs(gridB.leafGridView()),
                   countInteriorFacePairs(gridC.leafGridView()));
}

// processEclipseFormat orders a cell's corners differently on a left-handed
// grid; the leaf's children must still carry their own corners.
BOOST_AUTO_TEST_CASE(childCornersMatchOnBothHandednesses)
{
    for (const double ySign : {1.0, -1.0}) {
        auto parent = makeFaultedGrid({5, 3, 2}, 2, 0.6);
        for (std::size_t p = 0; p < parent.coord.size(); p += 3) {
            parent.coord[p + 1] *= ySign;
        }

        Dune::CpGrid gridA;
        auto rawParent = parent.raw();
        gridA.processEclipseFormat(rawParent, false);

        useConformingBuilder(gridA, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
            parent.dims, parent.coord, parent.zcorn, parent.actnum));

        Opm::Refinement::BlockRefinement req;
        req.name = "LGR1";
        req.cellsPerDim = {3, 2, 2};
        req.startIJK = {1, 0, 0};
        req.endIJK = {4, 3, 2};
        gridA.addLgrsUpdateLeafView({req.cellsPerDim}, {req.startIJK}, {req.endIJK}, {req.name});

        const auto refined = Opm::Refinement::refineBlock(parent.dims, parent.coord.data(),
                                                          parent.zcorn.data(), nullptr, req);
        Dune::CpGrid gridB;
        const auto rawRefined = grdecl{ {refined.dims[0], refined.dims[1], refined.dims[2]},
                                        refined.coord.data(), refined.zcorn.data(),
                                        refined.actnum.data() };
        gridB.processEclipseFormat(rawRefined, false);
        std::map<int, std::array<double,3>> oracle;
        for (const auto& element : Dune::elements(gridB.leafGridView())) {
            oracle[gridB.globalCell()[element.index()]] = gridB.getEclCentroid(element.index());
        }

        for (const auto& element : Dune::elements(gridA.leafGridView())) {
            if (!element.hasFather()) {
                continue;
            }
            const int levelIdx = element.getLevelElem().index();
            const auto it = oracle.find(gridA.currentData()[1]->globalCell()[levelIdx]);
            BOOST_REQUIRE(it != oracle.end());
            const auto centroid = gridA.getEclCentroid(element.index());
            for (int c = 0; c < 3; ++c) {
                BOOST_CHECK_SMALL(centroid[c] - it->second[c], 1e-10);
            }
        }
    }
}
