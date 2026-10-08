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
//! \file
//! \brief Demonstration that the Dune adapt-framework entry point
//!        CpGrid::globalRefine() is served by the new refinement builder
//!        (the engine AdaptiveCpGrid drives), without any of the removed Dune
//!        mark/adapt machinery. globalRefine(1) refines every cell 2x2x2 and is
//!        bit-equivalent to a whole-grid addLgrsUpdateLeafView.
#include <config.h>

#define BOOST_TEST_MODULE GlobalRefineViaBuilderTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/refinement/conforming/ConformingBlockBuilder.hpp>
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <dune/common/parallel/mpihelper.hh>

#include <algorithm>
#include <array>
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

struct Grdecl
{
    std::array<int,3> dims;
    std::vector<double> coord;
    std::vector<double> zcorn;

    grdecl raw() const
    {
        grdecl g{};
        g.dims[0] = dims[0]; g.dims[1] = dims[1]; g.dims[2] = dims[2];
        g.coord = const_cast<double*>(coord.data());
        g.zcorn = const_cast<double*>(zcorn.data());
        g.actnum = nullptr;
        return g;
    }
};

// Uniform Cartesian unit-cell grid as COORD/ZCORN, so the grid carries a real
// corner-point description the conforming builder can refine.
Grdecl makeCartesian(int nx, int ny, int nz)
{
    Grdecl g;
    g.dims = {nx, ny, nz};
    g.coord.resize(6 * (nx + 1) * (ny + 1));
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            double* p = &g.coord[6 * (static_cast<std::size_t>(j) * (nx + 1) + i)];
            p[0] = i; p[1] = j; p[2] = 0.0;
            p[3] = i; p[4] = j; p[5] = nz;
        }
    g.zcorn.resize(8 * static_cast<std::size_t>(nx) * ny * nz);
    for (int k = 0; k < 2 * nz; ++k)
        for (int j = 0; j < 2 * ny; ++j)
            for (int i = 0; i < 2 * nx; ++i)
                g.zcorn[static_cast<std::size_t>(i) + 2 * static_cast<std::size_t>(nx) * j
                        + 4 * static_cast<std::size_t>(nx) * ny * k] = 0.5 * (k + 1);
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

// Every leaf cell is closed (necessary for a conformal leaf).
void checkClosedLeaf(const Dune::CpGrid& grid)
{
    for (const auto& e : Dune::elements(grid.leafGridView())) {
        Dune::FieldVector<double,3> closure(0.0);
        for (const auto& is : Dune::intersections(grid.leafGridView(), e)) {
            auto n = is.centerUnitOuterNormal();
            n *= is.geometry().volume();
            closure += n;
        }
        BOOST_CHECK_SMALL(closure.two_norm(), 1e-9);
    }
}

std::unique_ptr<Opm::Refinement::Builder> makeBuilder(const Grdecl& g)
{
    return std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        g.dims, g.coord, g.zcorn, std::vector<int>{});
}

} // namespace

BOOST_GLOBAL_FIXTURE(MPIFixture);

// globalRefine(1) refines every cell 2x2x2 through the conforming builder and is
// the same leaf as a whole-grid addLgrsUpdateLeafView({2,2,2}).
BOOST_AUTO_TEST_CASE(globalRefineOneMatchesWholeGridLgr)
{
    auto g = makeCartesian(3, 2, 2);

    Dune::CpGrid viaGlobal;
    { auto raw = g.raw(); viaGlobal.processEclipseFormat(raw, false); }
    const int coarse = viaGlobal.size(0);
    useConformingBuilder(viaGlobal, makeBuilder(g));
    viaGlobal.globalRefine(1);

    BOOST_REQUIRE_EQUAL(viaGlobal.maxLevel(), 1);
    BOOST_CHECK_EQUAL(viaGlobal.size(0), coarse * 8);   // 2x2x2 per cell
    checkClosedLeaf(viaGlobal);

    // The reserved level-zero name map entry survives; the refinement gets its
    // own name ("GR"), never "GLOBAL".
    BOOST_CHECK_EQUAL(viaGlobal.getLgrNameToLevel().at("GLOBAL"), 0);
    BOOST_CHECK_EQUAL(viaGlobal.getLgrNameToLevel().at("GR1"), 1);

    // Reference: the same refinement requested directly as one whole-grid LGR.
    Dune::CpGrid viaLgr;
    { auto raw = g.raw(); viaLgr.processEclipseFormat(raw, false); }
    useConformingBuilder(viaLgr, makeBuilder(g));
    viaLgr.addLgrsUpdateLeafView({{2,2,2}}, {{0,0,0}}, {{3,2,2}}, {"WHOLE"});

    BOOST_CHECK_EQUAL(viaGlobal.size(0), viaLgr.size(0));
    BOOST_CHECK_EQUAL(viaGlobal.maxLevel(), viaLgr.maxLevel());
    BOOST_CHECK_CLOSE(totalVolume(viaGlobal), totalVolume(viaLgr), 1e-9);
}

// globalRefine(0) is a no-op.
BOOST_AUTO_TEST_CASE(globalRefineZeroIsNoop)
{
    auto g = makeCartesian(2, 2, 1);
    Dune::CpGrid grid;
    { auto raw = g.raw(); grid.processEclipseFormat(raw, false); }
    grid.setLgrBackend(Opm::Refinement::Backend::Conforming);
    const int coarse = grid.size(0);
    grid.globalRefine(0);
    BOOST_CHECK_EQUAL(grid.maxLevel(), 0);
    BOOST_CHECK_EQUAL(grid.size(0), coarse);
}

// One builder pass gives one level, so globalRefine(2) cannot keep Dune's n-levels contract.
BOOST_AUTO_TEST_CASE(globalRefineTwoThrows)
{
    auto g = makeCartesian(3, 2, 2);
    Dune::CpGrid grid;
    { auto raw = g.raw(); grid.processEclipseFormat(raw, false); }
    useConformingBuilder(grid, makeBuilder(g));
    BOOST_CHECK_THROW(grid.globalRefine(2), std::logic_error);
    BOOST_CHECK_EQUAL(grid.maxLevel(), 0);
}

// Negative count and an already-refined grid throw a clear error.
BOOST_AUTO_TEST_CASE(globalRefineUnsupportedCasesThrow)
{
    auto g = makeCartesian(2, 2, 1);
    Dune::CpGrid grid;
    { auto raw = g.raw(); grid.processEclipseFormat(raw, false); }

    BOOST_CHECK_THROW(grid.globalRefine(-1), std::logic_error);

    useConformingBuilder(grid, makeBuilder(g));
    grid.globalRefine(1);                                        // now refined
    BOOST_CHECK_THROW(grid.globalRefine(1), std::logic_error);   // already refined
}

// autoRefine(nxnynz) is arbitrary anisotropic global refinement = one full-grid
// LGR with cells_per_dim = nxnynz. Demonstrates that "any (odd) division of the
// whole grid in one LGR" is fully supported -- it just needs the per-direction
// interface, not Dune's single-int globalRefine.
BOOST_AUTO_TEST_CASE(autoRefineAnyOddDivisionMatchesWholeGridLgr)
{
    auto g = makeCartesian(3, 2, 2);

    Dune::CpGrid viaAuto;
    { auto raw = g.raw(); viaAuto.processEclipseFormat(raw, false); }
    const int coarse = viaAuto.size(0);
    useConformingBuilder(viaAuto, makeBuilder(g));
    viaAuto.autoRefine({3, 5, 7});                       // anisotropic, odd

    BOOST_REQUIRE_EQUAL(viaAuto.maxLevel(), 1);
    BOOST_CHECK_EQUAL(viaAuto.size(0), coarse * (3 * 5 * 7));
    checkClosedLeaf(viaAuto);

    Dune::CpGrid viaLgr;
    { auto raw = g.raw(); viaLgr.processEclipseFormat(raw, false); }
    useConformingBuilder(viaLgr, makeBuilder(g));
    viaLgr.addLgrsUpdateLeafView({{3,5,7}}, {{0,0,0}}, {{3,2,2}}, {"WHOLE"});
    BOOST_CHECK_EQUAL(viaAuto.size(0), viaLgr.size(0));
    BOOST_CHECK_CLOSE(totalVolume(viaAuto), totalVolume(viaLgr), 1e-9);
}

// autoRefine validates its factors: non-positive and even both throw.
BOOST_AUTO_TEST_CASE(autoRefineRejectsNonPositiveAndEvenFactors)
{
    auto g = makeCartesian(2, 2, 1);
    Dune::CpGrid grid;
    { auto raw = g.raw(); grid.processEclipseFormat(raw, false); }
    useConformingBuilder(grid, makeBuilder(g));

    BOOST_CHECK_THROW(grid.autoRefine({0, 3, 5}), std::invalid_argument);  // non-positive
    BOOST_CHECK_THROW(grid.autoRefine({3, -1, 5}), std::invalid_argument); // negative
    BOOST_CHECK_THROW(grid.autoRefine({4, 3, 5}), std::invalid_argument);  // even
    BOOST_CHECK_EQUAL(grid.maxLevel(), 0);                                 // unchanged
}

// mark()/adapt() belong to the Trilinear backend; the Conforming one says so instead of ignoring them.
BOOST_AUTO_TEST_CASE(markAndAdaptRefusedByConforming)
{
    auto g = makeCartesian(2, 2, 1);
    Dune::CpGrid grid;
    { auto raw = g.raw(); grid.processEclipseFormat(raw, false); }
    useConformingBuilder(grid, makeBuilder(g));
    const auto element = *Dune::elements(grid.leafGridView()).begin();
    BOOST_CHECK_THROW(grid.mark(1, element), std::logic_error);
    BOOST_CHECK_EQUAL(grid.getMark(element), 0);
    BOOST_CHECK_THROW(grid.preAdapt(), std::logic_error);
    BOOST_CHECK_THROW(grid.adapt(), std::logic_error);
    BOOST_CHECK_NO_THROW(grid.postAdapt());
    BOOST_CHECK_NO_THROW(grid.syncDistributedGlobalCellIds());
    BOOST_CHECK_EQUAL(grid.maxLevel(), 0);
}
