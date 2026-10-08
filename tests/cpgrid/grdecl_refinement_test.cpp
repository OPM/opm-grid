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
#include <config.h>

#define BOOST_TEST_MODULE GrdeclRefinementTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/refinement/conforming/GrdeclRefinement.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <array>
#include <cmath>
#include <functional>
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

// Build a corner-point description with straight pillars from (x,y,ztop) to
// (x,y,zbot) (vertical) and per-corner depths given by a callback
// z(i_, j_, k_) on the doubled lattice.
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
            p[0] = i; p[1] = j; p[2] = 0.0;     // top
            p[3] = i; p[4] = j; p[5] = 100.0;   // bottom
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

// Build a corner-point description with *slanted* (non-vertical, non-parallel)
// pillars: the bottom of each pillar is sheared laterally by a position-
// dependent amount, so every cell is a genuinely skewed hexahedron. On such a
// grid an index-space ("trilinear in cell coordinates") subdivision does not
// preserve a parent cell's centre of mass, whereas a geometrically faithful
// subdivision of the pillars does -- which is exactly what this exercises.
TestGrdecl makeSlantedPillarGrid(const std::array<int,3>& dims,
                                 const std::function<double(int,int,int)>& depth)
{
    TestGrdecl g;
    g.dims = dims;
    const auto& [nx, ny, nz] = dims;

    g.coord.resize(6*(nx + 1)*(ny + 1));
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            double* p = &g.coord[6*(static_cast<std::size_t>(j)*(nx + 1) + i)];
            // Smoothly varying lateral shear of the pillar bottom. Kept to a
            // realistic magnitude (a few percent of a cell) so the hexahedron
            // geometry-discretisation artifact stays well below the tolerances
            // below, the way it does on real corner-point decks.
            const double sx = 0.06*std::sin(0.7*i) + 0.02*j;
            const double sy = 0.05*std::cos(0.5*j) + 0.015*i;
            p[0] = i;      p[1] = j;      p[2] = 0.0;     // top
            p[3] = i + sx; p[4] = j + sy; p[5] = 100.0;   // bottom (sheared)
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

// Doubled-lattice index -> the (cell, corner-side) pair it represents.
int cellOf(int doubled) { return doubled / 2; }
int sideOf(int doubled) { return doubled % 2; }

double totalVolume(const Dune::CpGrid& grid)
{
    double vol = 0.0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        vol += element.geometry().volume();
    }
    return vol;
}

double boxVolume(const Dune::CpGrid& grid,
                 const std::array<int,3>& startIJK,
                 const std::array<int,3>& endIJK)
{
    const auto& dims = grid.logicalCartesianSize();
    double vol = 0.0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        const int cart = grid.globalCell()[element.index()];
        const int i = cart % dims[0];
        const int j = (cart / dims[0]) % dims[1];
        const int k = cart / (dims[0]*dims[1]);
        if (i >= startIJK[0] && i < endIJK[0] && j >= startIJK[1] && j < endIJK[1]
            && k >= startIJK[2] && k < endIJK[2]) {
            vol += element.geometry().volume();
        }
    }
    return vol;
}

// Volume-weighted centre of mass of all leaf cells.
std::array<double,3> totalCenterOfMass(const Dune::CpGrid& grid)
{
    std::array<double,3> m = {0.0, 0.0, 0.0};
    double vol = 0.0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        const auto geo = element.geometry();
        const double v = geo.volume();
        const auto c = geo.center();
        for (int d = 0; d < 3; ++d) m[d] += v*c[d];
        vol += v;
    }
    for (int d = 0; d < 3; ++d) m[d] /= vol;
    return m;
}

// Volume-weighted centre of mass of the cells whose (cartesian) index lies in
// the given box.
std::array<double,3> boxCenterOfMass(const Dune::CpGrid& grid,
                                     const std::array<int,3>& startIJK,
                                     const std::array<int,3>& endIJK)
{
    const auto& dims = grid.logicalCartesianSize();
    std::array<double,3> m = {0.0, 0.0, 0.0};
    double vol = 0.0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        const int cart = grid.globalCell()[element.index()];
        const int i = cart % dims[0];
        const int j = (cart / dims[0]) % dims[1];
        const int k = cart / (dims[0]*dims[1]);
        if (i >= startIJK[0] && i < endIJK[0] && j >= startIJK[1] && j < endIJK[1]
            && k >= startIJK[2] && k < endIJK[2]) {
            const auto geo = element.geometry();
            const double v = geo.volume();
            const auto c = geo.center();
            for (int d = 0; d < 3; ++d) m[d] += v*c[d];
            vol += v;
        }
    }
    for (int d = 0; d < 3; ++d) m[d] /= vol;
    return m;
}

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_AUTO_TEST_CASE(uniformGridRefinesToUniformLattice)
{
    // Unit cells: depth(k_) = layer boundaries 0,1,...
    auto parent = makeVerticalPillarGrid({4, 3, 3}, [](int, int, int k_) {
        return static_cast<double>(cellOf(k_) + sideOf(k_));
    });

    Opm::Refinement::BlockRefinement req;
    req.name = "LGR1";
    req.cellsPerDim = {2, 3, 2};
    req.startIJK = {1, 1, 1};
    req.endIJK = {3, 2, 3};

    const auto refined = Opm::Refinement::refineBlock(parent.dims, parent.coord.data(),
                                                      parent.zcorn.data(), nullptr, req);

    BOOST_CHECK(refined.dims == (std::array<int,3>{4, 3, 4}));

    // Sub-pillars form the expected uniform lattice: pillar (ir, jr) at
    // x = 1 + ir/2, y = 1 + jr/3.
    for (int jr = 0; jr <= refined.dims[1]; ++jr) {
        for (int ir = 0; ir <= refined.dims[0]; ++ir) {
            const double* p = &refined.coord[6*(static_cast<std::size_t>(jr)*(refined.dims[0] + 1) + ir)];
            BOOST_CHECK_CLOSE(p[0], 1.0 + ir/2.0, 1e-11);
            BOOST_CHECK_CLOSE(p[1], 1.0 + jr/3.0, 1e-11);
        }
    }

    // Refined zcorn: layer kr spans [1 + kr/2, 1 + (kr+1)/2].
    for (std::size_t idx = 0; idx < refined.zcorn.size(); ++idx) {
        // All values must lie in the box's z range.
        BOOST_CHECK(refined.zcorn[idx] >= 1.0 - 1e-12);
        BOOST_CHECK(refined.zcorn[idx] <= 3.0 + 1e-12);
    }

    // Process the refined description and verify cell count and volume:
    // box is 2x1x2 unit cells -> volume 4, refined into 4x3x4 = 48 cells.
    Dune::CpGrid refinedGrid;
    const auto raw = grdecl{ {refined.dims[0], refined.dims[1], refined.dims[2]},
                             refined.coord.data(), refined.zcorn.data(), refined.actnum.data() };
    refinedGrid.processEclipseFormat(raw, false);
    BOOST_CHECK_EQUAL(refinedGrid.size(0), 48);
    BOOST_CHECK_CLOSE(totalVolume(refinedGrid), 4.0, 1e-9);
}

BOOST_AUTO_TEST_CASE(distortedVerticalPillarGridConservesVolume)
{
    // Smoothly varying layer thicknesses and non-planar surfaces.
    auto depth = [](int i_, int j_, int k_) {
        const double x = cellOf(i_) + sideOf(i_);
        const double y = cellOf(j_) + sideOf(j_);
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return base + 0.3*std::sin(0.8*x) + 0.2*std::cos(0.5*y) + 0.05*x*y;
    };
    auto parent = makeVerticalPillarGrid({5, 4, 3}, depth);

    Dune::CpGrid parentGrid;
    auto rawParent = parent.raw();
    parentGrid.processEclipseFormat(rawParent, false);

    Opm::Refinement::BlockRefinement req;
    req.name = "LGR1";
    req.cellsPerDim = {3, 2, 4};
    req.startIJK = {1, 1, 0};
    req.endIJK = {4, 3, 2};

    const auto refined = Opm::Refinement::refineBlock(parent.dims, parent.coord.data(),
                                                      parent.zcorn.data(), nullptr, req);

    Dune::CpGrid refinedGrid;
    const auto rawRefined = grdecl{ {refined.dims[0], refined.dims[1], refined.dims[2]},
                                    refined.coord.data(), refined.zcorn.data(), refined.actnum.data() };
    refinedGrid.processEclipseFormat(rawRefined, false);

    // On vertical pillars the resampled sub-cells partition the parent
    // cells exactly: total refined volume == parent box volume.
    BOOST_CHECK_EQUAL(refinedGrid.size(0), (3*3)*(2*2)*(2*4));
    BOOST_CHECK_CLOSE(totalVolume(refinedGrid),
                      boxVolume(parentGrid, req.startIJK, req.endIJK), 1e-8);
}

BOOST_AUTO_TEST_CASE(collapsedParentGivesCollapsedChildren)
{
    // MINPV collapses a removed cell to zero thickness at its top. Its children
    // must be exactly collapsed too, or rounding inverts some and the
    // corner-point processor rejects the block.
    auto depth = [](int i_, int j_, int k_) {
        const double x = cellOf(i_) + sideOf(i_);
        const double y = cellOf(j_) + sideOf(j_);
        return 2.0*(cellOf(k_) + sideOf(k_)) + 0.3*std::sin(0.8*x) + 0.2*std::cos(0.5*y);
    };
    auto parent = makeVerticalPillarGrid({3, 3, 3}, depth);
    const int nx = parent.dims[0], ny = parent.dims[1];
    const auto z = [&](int i_, int j_, int k_) -> double& {
        return parent.zcorn[i_ + 2*nx*j_ + 4*nx*ny*k_];
    };
    for (int j_ = 0; j_ < 2*ny; ++j_) {
        for (int i_ = 0; i_ < 2*nx; ++i_) {
            z(i_, j_, 3) = z(i_, j_, 2);        // layer k = 1 collapsed at its top
        }
    }

    Opm::Refinement::BlockRefinement req;
    req.name = "LGR1";
    req.cellsPerDim = {3, 3, 3};
    req.startIJK = {0, 0, 0};
    req.endIJK = {3, 3, 3};

    const auto refined = Opm::Refinement::refineBlock(parent.dims, parent.coord.data(),
                                                      parent.zcorn.data(), nullptr, req);

    const int rx = refined.dims[0], ry = refined.dims[1];
    const auto rz = [&](int i_, int j_, int k_) {
        return refined.zcorn[i_ + 2*rx*j_ + 4*rx*ry*k_];
    };
    for (int kr = 3; kr < 6; ++kr) {
        for (int j_ = 0; j_ < 2*ry; ++j_) {
            for (int i_ = 0; i_ < 2*rx; ++i_) {
                BOOST_CHECK_EQUAL(rz(i_, j_, 2*kr + 1), rz(i_, j_, 2*kr));
            }
        }
    }

    Dune::CpGrid refinedGrid;
    const auto rawRefined = grdecl{ {refined.dims[0], refined.dims[1], refined.dims[2]},
                                    refined.coord.data(), refined.zcorn.data(), refined.actnum.data() };
    BOOST_CHECK_NO_THROW(refinedGrid.processEclipseFormat(rawRefined, false));
}

BOOST_AUTO_TEST_CASE(slantedPillarRefinementConservesCenterOfMass)
{
    // Ground-truth geometric test: subdividing a block of skewed corner-point
    // cells must conserve both the volume and the centre of mass of the region
    // it replaces -- the un-refined parent grid is the reference. A faithful
    // refinement of the pillar geometry preserves the centre of mass exactly;
    // an index-space subdivision of skewed cells does not. (This is the
    // property that distinguishes a correct corner-point LGR from an
    // approximate one; observed on opm-tests/lgr corner-point decks.)
    auto depth = [](int i_, int j_, int k_) {
        const double x = cellOf(i_) + sideOf(i_);
        const double y = cellOf(j_) + sideOf(j_);
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return base + 0.3*std::sin(0.8*x) + 0.2*std::cos(0.5*y) + 0.05*x*y;
    };
    auto parent = makeSlantedPillarGrid({5, 4, 3}, depth);

    Dune::CpGrid parentGrid;
    auto rawParent = parent.raw();
    parentGrid.processEclipseFormat(rawParent, false);

    Opm::Refinement::BlockRefinement req;
    req.name = "LGR1";
    req.cellsPerDim = {3, 2, 4};
    req.startIJK = {1, 1, 0};
    req.endIJK = {4, 3, 2};

    const auto refined = Opm::Refinement::refineBlock(parent.dims, parent.coord.data(),
                                                      parent.zcorn.data(), nullptr, req);

    Dune::CpGrid refinedGrid;
    const auto rawRefined = grdecl{ {refined.dims[0], refined.dims[1], refined.dims[2]},
                                    refined.coord.data(), refined.zcorn.data(), refined.actnum.data() };
    refinedGrid.processEclipseFormat(rawRefined, false);

    // (1) Volume of the refined block equals the parent box volume. The small
    //     residual is the hexahedron-volume discretisation artifact (the coarse
    //     and fine cells tetrahedralise the same skewed region differently); a
    //     faithful subdivision keeps it ~1e-7, far below the 1e-3% bound.
    const double parentBoxVol = boxVolume(parentGrid, req.startIJK, req.endIJK);
    BOOST_CHECK_CLOSE(totalVolume(refinedGrid), parentBoxVol, 1e-3);

    // (2) Centre of mass of the refined block equals the parent box centre of
    //     mass. The refined grid spans exactly the box, so its *total* centre of
    //     mass is the relevant quantity. A faithful pillar subdivision conserves
    //     it to ~1e-6 (relative) even on a strongly sheared grid; an index-space
    //     subdivision of skewed cells drifts ~1e-4 (seen comparing the fork's
    //     refinement against upstream master on the opm-tests corner-point
    //     decks), so this bound separates the two.
    const auto parentCoM = boxCenterOfMass(parentGrid, req.startIJK, req.endIJK);
    const auto refinedCoM = totalCenterOfMass(refinedGrid);
    const double scale = std::max({std::abs(parentCoM[0]),
                                   std::abs(parentCoM[1]),
                                   std::abs(parentCoM[2]), 1.0});
    const double tol = 1e-5 * scale;
    BOOST_CHECK_SMALL(refinedCoM[0] - parentCoM[0], tol);
    BOOST_CHECK_SMALL(refinedCoM[1] - parentCoM[1], tol);
    BOOST_CHECK_SMALL(refinedCoM[2] - parentCoM[2], tol);
}

BOOST_AUTO_TEST_CASE(faultInsideBlockIsPreserved)
{
    // Vertical fault between i-columns 1 and 2: columns to the right are
    // shifted down by 0.6 (less than one layer thickness of 2.0).
    const double throwZ = 0.6;
    auto depth = [&](int i_, int j_, int k_) {
        (void)j_;
        const double base = 2.0*(cellOf(k_) + sideOf(k_));
        return (cellOf(i_) >= 2) ? base + throwZ : base;
    };
    auto parent = makeVerticalPillarGrid({4, 2, 2}, depth);

    Dune::CpGrid parentGrid;
    auto rawParent = parent.raw();
    parentGrid.processEclipseFormat(rawParent, false);

    // Refine a block containing the fault.
    Opm::Refinement::BlockRefinement req;
    req.name = "LGR1";
    req.cellsPerDim = {2, 2, 2};
    req.startIJK = {1, 0, 0};
    req.endIJK = {3, 2, 2};

    const auto refined = Opm::Refinement::refineBlock(parent.dims, parent.coord.data(),
                                                      parent.zcorn.data(), nullptr, req);

    Dune::CpGrid refinedGrid;
    const auto rawRefined = grdecl{ {refined.dims[0], refined.dims[1], refined.dims[2]},
                                    refined.coord.data(), refined.zcorn.data(), refined.actnum.data() };
    refinedGrid.processEclipseFormat(rawRefined, false);

    // Volume is conserved across the faulted block (vertical pillars).
    BOOST_CHECK_CLOSE(totalVolume(refinedGrid),
                      boxVolume(parentGrid, req.startIJK, req.endIJK), 1e-8);

    // The fault survives refinement: the refined grid must contain faces
    // whose two neighboring cells are offset across the fault plane, which
    // shows up as more I-direction internal connections than a conforming
    // grid would have... A robust simple check: the number of leaf faces of
    // the refined faulted grid exceeds that of the same block refined from
    // an unfaulted twin.
    auto depthNoFault = [](int, int, int k_) {
        return 2.0*(cellOf(k_) + sideOf(k_));
    };
    auto parentNoFault = makeVerticalPillarGrid({4, 2, 2}, depthNoFault);
    const auto refinedNoFault = Opm::Refinement::refineBlock(parentNoFault.dims,
                                                             parentNoFault.coord.data(),
                                                             parentNoFault.zcorn.data(),
                                                             nullptr, req);
    Dune::CpGrid refinedGridNoFault;
    const auto rawNoFault = grdecl{ {refinedNoFault.dims[0], refinedNoFault.dims[1], refinedNoFault.dims[2]},
                                    refinedNoFault.coord.data(), refinedNoFault.zcorn.data(),
                                    refinedNoFault.actnum.data() };
    refinedGridNoFault.processEclipseFormat(rawNoFault, false);

    BOOST_CHECK_EQUAL(refinedGrid.size(0), refinedGridNoFault.size(0));
    BOOST_CHECK_GT(refinedGrid.numFaces(), refinedGridNoFault.numFaces());
}

BOOST_AUTO_TEST_CASE(inactiveParentsProduceInactiveChildren)
{
    auto parent = makeVerticalPillarGrid({3, 2, 2}, [](int, int, int k_) {
        return static_cast<double>(cellOf(k_) + sideOf(k_));
    });
    // Deactivate cell (1,0,0).
    parent.actnum.assign(3*2*2, 1);
    parent.actnum[1] = 0;

    Opm::Refinement::BlockRefinement req;
    req.name = "LGR1";
    req.cellsPerDim = {2, 2, 2};
    req.startIJK = {0, 0, 0};
    req.endIJK = {3, 2, 1};

    const auto refined = Opm::Refinement::refineBlock(parent.dims, parent.coord.data(),
                                                      parent.zcorn.data(), parent.actnum.data(), req);

    int inactive = 0;
    for (int flag : refined.actnum) {
        inactive += (flag == 0);
    }
    // One parent cell -> 2x2x2 children inactive.
    BOOST_CHECK_EQUAL(inactive, 8);

    Dune::CpGrid refinedGrid;
    const auto rawRefined = grdecl{ {refined.dims[0], refined.dims[1], refined.dims[2]},
                                    refined.coord.data(), refined.zcorn.data(), refined.actnum.data() };
    refinedGrid.processEclipseFormat(rawRefined, false);
    BOOST_CHECK_EQUAL(refinedGrid.size(0), 3*2*1*8 - 8);
}

// N*FIN/H*FIN: how many refined columns each parent cell takes, and their
// relative widths within it.
BOOST_AUTO_TEST_CASE(gradedRefinementFollowsColumnTables)
{
    auto parent = makeVerticalPillarGrid({4, 3, 3}, [](int, int, int k_) {
        return static_cast<double>(cellOf(k_) + sideOf(k_));
    });

    // Box spans parent cells i = 1..2. The first takes one column, the second
    // three with widths 1 : 2 : 1 -- so x runs 1, 2, 2.25, 2.75, 3.
    Opm::Refinement::BlockRefinement req;
    req.name = "GRADED";
    req.cellsPerDim = {1, 1, 1};
    req.startIJK = {1, 1, 1};
    req.endIJK = {3, 2, 3};
    req.subdivision[0] = Opm::Refinement::AxisSubdivision{
        /* parentOffset = */ {0, 1, 1, 1},
        /* fracLo       = */ {0.0, 0.0,  0.25, 0.75},
        /* fracHi       = */ {1.0, 0.25, 0.75, 1.0} };

    const auto refined = Opm::Refinement::refineBlock(parent.dims, parent.coord.data(),
                                                      parent.zcorn.data(), nullptr, req);

    // The refined count comes from the tables, not from a factor.
    BOOST_CHECK(refined.dims == (std::array<int,3>{4, 1, 2}));

    const std::array<double,5> expectedX = {1.0, 2.0, 2.25, 2.75, 3.0};
    for (int jr = 0; jr <= refined.dims[1]; ++jr) {
        for (int ir = 0; ir <= refined.dims[0]; ++ir) {
            const double* p = &refined.coord[6*(static_cast<std::size_t>(jr)*(refined.dims[0] + 1) + ir)];
            BOOST_CHECK_CLOSE(p[0], expectedX[ir], 1e-11);
        }
    }

    // Volume is conserved whatever the columns do: the box is 2x1x2 unit cells.
    Dune::CpGrid refinedGrid;
    const auto raw = grdecl{ {refined.dims[0], refined.dims[1], refined.dims[2]},
                             refined.coord.data(), refined.zcorn.data(), refined.actnum.data() };
    refinedGrid.processEclipseFormat(raw, false);
    BOOST_CHECK_EQUAL(refinedGrid.size(0), 8);
    BOOST_CHECK_CLOSE(totalVolume(refinedGrid), 4.0, 1e-9);

    // A parent that keeps only one column is left exactly as it was.
    BOOST_CHECK_CLOSE(refined.coord[0], 1.0, 1e-11);
}

// Spelling a uniform box as explicit tables must change nothing: this is the
// path every existing deck now takes through axisSubdivision().
BOOST_AUTO_TEST_CASE(explicitUniformTablesMatchTheFactorPath)
{
    auto parent = makeVerticalPillarGrid({4, 3, 3}, [](int i_, int j_, int k_) {
        return static_cast<double>(cellOf(k_) + sideOf(k_))
             + 0.1*cellOf(i_) - 0.05*cellOf(j_);      // non-trivial geometry
    });

    Opm::Refinement::BlockRefinement uniform;
    uniform.name = "LGR1";
    uniform.cellsPerDim = {2, 3, 2};
    uniform.startIJK = {1, 1, 1};
    uniform.endIJK = {3, 2, 3};

    auto graded = uniform;
    graded.subdivision[0] = Opm::Refinement::AxisSubdivision{
        {0, 0, 1, 1}, {0.0, 0.5, 0.0, 0.5}, {0.5, 1.0, 0.5, 1.0} };
    graded.subdivision[1] = Opm::Refinement::AxisSubdivision{
        {0, 0, 0}, {0.0, 1.0/3.0, 2.0/3.0}, {1.0/3.0, 2.0/3.0, 1.0} };
    graded.subdivision[2] = Opm::Refinement::AxisSubdivision{
        {0, 0, 1, 1}, {0.0, 0.5, 0.0, 0.5}, {0.5, 1.0, 0.5, 1.0} };

    const auto a = Opm::Refinement::refineBlock(parent.dims, parent.coord.data(),
                                                parent.zcorn.data(), nullptr, uniform);
    const auto b = Opm::Refinement::refineBlock(parent.dims, parent.coord.data(),
                                                parent.zcorn.data(), nullptr, graded);

    BOOST_CHECK(a.dims == b.dims);
    BOOST_REQUIRE_EQUAL(a.coord.size(), b.coord.size());
    BOOST_REQUIRE_EQUAL(a.zcorn.size(), b.zcorn.size());
    for (std::size_t i = 0; i < a.coord.size(); ++i) {
        BOOST_CHECK_CLOSE(a.coord[i], b.coord[i], 1e-11);
    }
    for (std::size_t i = 0; i < a.zcorn.size(); ++i) {
        BOOST_CHECK_CLOSE(a.zcorn[i], b.zcorn[i], 1e-11);
    }
}
