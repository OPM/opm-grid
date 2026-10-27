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

#define BOOST_TEST_MODULE LevelGridAssemblerTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/cpgrid/refinement/conforming/LevelGridAssembler.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <dune/common/parallel/mpihelper.hh>

#include <array>
#include <cmath>
#include <functional>
#include <map>
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
    auto parent = makeVerticalPillarGrid({4, 3, 2}, depth);

    Dune::CpGrid parentGrid;
    auto rawParent = parent.raw();
    parentGrid.processEclipseFormat(rawParent, false);

    Opm::Refinement::BlockRefinement req;
    req.name = "LGR1";
    req.cellsPerDim = {2, 3, 2};
    req.startIJK = {1, 0, 0};
    req.endIJK = {3, 2, 2};

    auto& storage = parentGrid.currentData();
    auto level = Opm::Refinement::assembleBlockLevelGrid(*storage[0], parent.dims,
                                                         parent.coord.data(), parent.zcorn.data(),
                                                         nullptr, req, /* levelIndex = */ 1,
                                                         /* parentLevel = */ 0,
                                                         storage, Dune::MPIHelper::getCommunicator());
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
    auto parent = makeVerticalPillarGrid({3, 2, 2}, [](int, int, int k_) {
        return static_cast<double>(cellOf(k_) + sideOf(k_));
    });
    parent.actnum.assign(3*2*2, 1);
    parent.actnum[4] = 0; // cell (1,1,0)

    Dune::CpGrid parentGrid;
    auto rawParent = parent.raw();
    parentGrid.processEclipseFormat(rawParent, false);

    Opm::Refinement::BlockRefinement req;
    req.name = "LGR1";
    req.cellsPerDim = {2, 2, 2};
    req.startIJK = {0, 0, 0};
    req.endIJK = {3, 2, 1};

    auto& storage = parentGrid.currentData();
    auto level = Opm::Refinement::assembleBlockLevelGrid(*storage[0], parent.dims,
                                                         parent.coord.data(), parent.zcorn.data(),
                                                         parent.actnum.data(), req, 1,
                                                         /* parentLevel = */ 0,
                                                         storage, Dune::MPIHelper::getCommunicator());
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
    auto parent = makeVerticalPillarGrid({4, 2, 2}, depth);

    Dune::CpGrid parentGrid;
    auto rawParent = parent.raw();
    parentGrid.processEclipseFormat(rawParent, false);

    Opm::Refinement::BlockRefinement req;
    req.name = "LGR1";
    req.cellsPerDim = {2, 2, 2};
    req.startIJK = {1, 0, 0};
    req.endIJK = {3, 2, 2};

    auto& storage = parentGrid.currentData();
    auto level = Opm::Refinement::assembleBlockLevelGrid(*storage[0], parent.dims,
                                                         parent.coord.data(), parent.zcorn.data(),
                                                         nullptr, req, 1,
                                                         /* parentLevel = */ 0,
                                                         storage, Dune::MPIHelper::getCommunicator());
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
