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

  Distributed (MPI) refinement of a rank-interior LGR box: the box lives
  entirely inside one rank's interior (PLAN Track 1 step 6). Only that rank
  refines it; other ranks carry an empty placeholder level grid and a
  purely coarse leaf. Verifies per-rank cell structure and that the
  distributed leaf conserves the total (refined) volume and cell count.
*/
#include <config.h>

#define BOOST_TEST_MODULE DistributedBuilderTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/refinement/conforming/ConformingBlockBuilder.hpp>
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <dune/common/parallel/mpihelper.hh>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <array>
#include <iostream>
#include <cmath>
#include <map>
#include <set>
#include <algorithm>
#include <functional>
#include <memory>
#include <sstream>
#include <string>
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
        g.coord = coord.data(); g.zcorn = zcorn.data();
        g.actnum = actnum.empty() ? nullptr : actnum.data();
        return g;
    }
};

TestGrdecl makeUnitGrid(const std::array<int,3>& dims)
{
    TestGrdecl g;
    g.dims = dims;
    const auto& [nx, ny, nz] = dims;
    g.coord.resize(6*(nx + 1)*(ny + 1));
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            double* p = &g.coord[6*(static_cast<std::size_t>(j)*(nx + 1) + i)];
            p[0] = i; p[1] = j; p[2] = 0.0;
            p[3] = i; p[4] = j; p[5] = static_cast<double>(nz);
        }
    }
    g.zcorn.resize(8*static_cast<std::size_t>(nx)*ny*nz);
    for (int k = 0; k < 2*nz; ++k) {
        const double z = (k + 1) / 2;
        for (std::size_t idx = 0; idx < 4*static_cast<std::size_t>(nx)*ny; ++idx) {
            g.zcorn[k*4*static_cast<std::size_t>(nx)*ny + idx] = z;
        }
    }
    return g;
}

// Unit grid with a vertical fault in front of column iFault; that column and those after it drop by dz.
TestGrdecl makeFaultedGrid(const std::array<int,3>& dims, int iFault, double dz)
{
    auto g = makeUnitGrid(dims);
    const std::size_t nx = dims[0];
    for (std::size_t z = 0; z < g.zcorn.size(); ++z) {
        if (static_cast<int>((z % (2*nx)) / 2) >= iFault) {
            g.zcorn[z] += dz;
        }
    }
    return g;
}

void useConformingBuilder(Dune::CpGrid& grid, std::unique_ptr<Opm::Refinement::Builder> builder)
{
    grid.setLgrBackend(Opm::Refinement::Backend::Conforming);
    grid.setRefinementBuilder(std::move(builder));
}

// Interior leaf cell ids keyed by cell centre (distinct on a unit grid).
std::map<std::array<long,3>, long> leafIdsByCentre(const Dune::CpGrid& grid)
{
    std::map<std::array<long,3>, long> ids;
    for (const auto& e : Dune::elements(grid.leafGridView(), Dune::Partitions::interior)) {
        const auto x = e.geometry().center();
        ids[{std::lround(x[0]*1e6), std::lround(x[1]*1e6), std::lround(x[2]*1e6)}] = grid.globalIdSet().id(e);
    }
    return ids;
}

// Leaf vertex ids keyed by position.
std::map<std::array<long,3>, long> leafVertexIdsByPosition(const Dune::CpGrid& grid)
{
    std::map<std::array<long,3>, long> ids;
    for (const auto& v : Dune::vertices(grid.leafGridView())) {
        const auto x = v.geometry().center();
        ids[{std::lround(x[0]*1e6), std::lround(x[1]*1e6), std::lround(x[2]*1e6)}] = grid.globalIdSet().id(v);
    }
    return ids;
}

// Number of entries of mine that are missing from, or differ in, reference; summed over ranks.
int mismatches(const Dune::CpGrid& grid, const std::map<std::array<long,3>, long>& mine,
               const std::map<std::array<long,3>, long>& reference)
{
    int bad = 0;
    for (const auto& [x, id] : mine) {
        const auto it = reference.find(x);
        bad += it == reference.end() || it->second != id;
    }
    return grid.comm().sum(bad);
}

} // anonymous namespace

BOOST_GLOBAL_FIXTURE(MPIFixture);

#if HAVE_MPI
// Distributed refinement of a rank-interior box. The builder infrastructure
// (rank-interior classification, empty placeholder level grids on non-owning
// ranks, self-communicator level grids so only the owning rank refines
// without deadlock) is in place. The full distributed *leaf* assembly is the
// remaining blocker (it currently faults reading the distributed level-zero
// overlap structure), so this case is disabled pending that fix; the
// rank-interior enforcement is covered by boxTouchingOverlapThrows below.
BOOST_AUTO_TEST_CASE(rankInteriorBoxRefinedInParallel)
{
    const std::array<int,3> dims = {{12, 12, 4}};
    auto g = makeUnitGrid(dims);

    Dune::CpGrid grid;
    grid.createCartesian(dims, {{1.0, 1.0, 1.0}}); // unit cells, matching makeUnitGrid
    if (grid.comm().size() < 2) {
        return;
    }

    std::vector<int> parts(static_cast<std::size_t>(dims[0])*dims[1]*dims[2]);
    const int np = grid.comm().size();
    for (int k = 0; k < dims[2]; ++k) {
        for (int j = 0; j < dims[1]; ++j) {
            for (int i = 0; i < dims[0]; ++i) {
                parts[i + dims[0]*j + dims[0]*dims[1]*k] = (i < 8) ? 0 : (1 + ((i - 8) % (np - 1)));
            }
        }
    }
    grid.loadBalance(parts, false, true, 2);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        g.dims, g.coord, g.zcorn, g.actnum));
    grid.addLgrsUpdateLeafView({{2,2,2}}, {{2,2,1}}, {{5,5,3}}, {"LGR1"});
    BOOST_CHECK_EQUAL(grid.maxLevel(), 1);

    int refined = 0;
    double interiorVolume = 0.0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        if (element.partitionType() != Dune::InteriorEntity) {
            continue;
        }
        interiorVolume += element.geometry().volume();
        if (element.hasFather()) {
            ++refined;
        }
    }
    BOOST_CHECK_EQUAL(grid.comm().sum(refined > 0 ? 1 : 0), 1);
    BOOST_CHECK_EQUAL(grid.comm().sum(refined), 18*8);
    BOOST_CHECK_CLOSE(grid.comm().sum(interiorVolume), 12.0*12.0*4.0, 1e-8);
}

// Mechanics assembles on vertices, so leaf vertices need what flow never asks
// for: one global id per vertex on every rank, and exactly one owner.
BOOST_AUTO_TEST_CASE(rankInteriorLeafVerticesConsistent)
{
    const std::array<int,3> dims = {{12, 12, 4}};
    auto g = makeUnitGrid(dims);

    Dune::CpGrid grid;
    grid.createCartesian(dims, {{1.0, 1.0, 1.0}});
    if (grid.comm().size() < 2) {
        return;
    }
    std::vector<int> parts(static_cast<std::size_t>(dims[0])*dims[1]*dims[2]);
    const int np = grid.comm().size();
    for (int k = 0; k < dims[2]; ++k) {
        for (int j = 0; j < dims[1]; ++j) {
            for (int i = 0; i < dims[0]; ++i) {
                parts[i + dims[0]*j + dims[0]*dims[1]*k] = (i < 8) ? 0 : (1 + ((i - 8) % (np - 1)));
            }
        }
    }
    grid.loadBalance(parts, false, true, 2);
    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        g.dims, g.coord, g.zcorn, g.actnum));
    grid.addLgrsUpdateLeafView({{2,2,2}}, {{2,2,1}}, {{5,5,3}}, {"LGR1"});

    // id, x, y, z, partition type, rank
    std::vector<double> mine;
    const auto& gv = grid.leafGridView();
    const auto& ids = grid.globalIdSet();
    for (const auto& v : Dune::vertices(gv)) {
        const auto x = v.geometry().center();
        mine.insert(mine.end(), {static_cast<double>(ids.id(v)), x[0], x[1], x[2],
                                 static_cast<double>(v.partitionType()),
                                 static_cast<double>(grid.comm().rank())});
    }
    const int n = static_cast<int>(mine.size());
    std::vector<int> counts(np);
    grid.comm().allgather(&n, 1, counts.data());
    std::vector<int> displ(np + 1, 0);
    for (int r = 0; r < np; ++r) {
        displ[r + 1] = displ[r] + counts[r];
    }
    std::vector<double> all(displ[np]);
    grid.comm().allgatherv(mine.data(), n, all.data(), counts.data(), displ.data());

    std::map<long, std::array<double,3>> posOfId;
    std::map<std::array<long,3>, long> idOfPos;
    std::map<long, std::array<int,2>> owners; // interior count, border count
    int idClash = 0, posClash = 0;
    for (std::size_t e = 0; e < all.size(); e += 6) {
        const long id = static_cast<long>(all[e]);
        const std::array<double,3> x {all[e+1], all[e+2], all[e+3]};
        const std::array<long,3> key {std::lround(x[0]*1e6), std::lround(x[1]*1e6), std::lround(x[2]*1e6)};
        const auto [it, fresh] = posOfId.emplace(id, x);
        if (!fresh && (std::abs(it->second[0]-x[0]) + std::abs(it->second[1]-x[1])
                       + std::abs(it->second[2]-x[2])) > 1e-9) {
            ++idClash;
        }
        const auto [pit, pfresh] = idOfPos.emplace(key, id);
        if (!pfresh && pit->second != id) {
            ++posClash;
        }
        const auto type = static_cast<Dune::PartitionType>(static_cast<int>(all[e+4]));
        auto& o = owners[id];
        o[0] += (type == Dune::InteriorEntity);
        o[1] += (type == Dune::BorderEntity);
    }
    int badOwner = 0;
    for (const auto& [id, o] : owners) {
        const bool ok = (o[0] == 1 && o[1] == 0) || (o[0] == 0 && o[1] >= 1);
        badOwner += !ok;
    }
    if (grid.comm().rank() == 0) {
        std::cout << "leaf vertices: " << owners.size() << " ids, " << idOfPos.size()
                  << " positions; one id at two positions " << idClash
                  << ", one position under two ids " << posClash
                  << ", ids without a single owner " << badOwner << std::endl;
    }
    BOOST_CHECK_EQUAL(idClash, 0);
    BOOST_CHECK_EQUAL(posClash, 0);
    BOOST_CHECK_EQUAL(badOwner, 0);
}

// Refined cell ids depend on the LGR and the cell's place in it, not on the partition.
BOOST_AUTO_TEST_CASE(refinedCellIdsIndependentOfPartition)
{
    const std::array<int,3> dims = {{16, 16, 2}};
    auto g = makeUnitGrid(dims);
    const auto refine = [&g](Dune::CpGrid& grid) {
        useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
            g.dims, g.coord, g.zcorn, g.actnum));
        grid.addLgrsUpdateLeafView({{2,2,2}, {3,3,1}, {2,2,1}},
                                   {{2,2,0}, {11,11,0}, {1,1,1}},
                                   {{5,5,2}, {14,14,1}, {3,3,3}},
                                   {"LGR1", "LGR2", "NEST1"},
                                   {"GLOBAL", "GLOBAL", "LGR1"});
    };

    Dune::CpGrid serial(MPI_COMM_SELF);
    serial.createCartesian(dims, {{1.0, 1.0, 1.0}});
    refine(serial);
    BOOST_REQUIRE_EQUAL(serial.maxLevel(), 3);
    const auto& ids = serial.globalIdSet();

    long offset = 0; // one above every level-zero id, then past each LGR's Cartesian box
    for (const auto& e : Dune::elements(serial.levelGridView(0))) {
        offset = std::max(offset, static_cast<long>(ids.id(e)) + 1);
    }
    for (const auto& v : Dune::vertices(serial.levelGridView(0))) {
        offset = std::max(offset, static_cast<long>(ids.id(v)) + 1);
    }
    for (int level = 1; level <= serial.maxLevel(); ++level) {
        int wrong = 0;
        for (const auto& e : Dune::elements(serial.levelGridView(level))) {
            wrong += ids.id(e) != offset + e.getLevelCartesianIdx();
        }
        BOOST_CHECK_EQUAL(wrong, 0);
        const auto& d = serial.currentData()[level]->logicalCartesianSize();
        offset += static_cast<long>(d[0])*d[1]*d[2];
    }
    std::set<long> cellIds, pointIds;
    int leafDiffers = 0;
    for (const auto& e : Dune::elements(serial.leafGridView())) {
        cellIds.insert(ids.id(e));
        leafDiffers += ids.id(e) != ids.id(e.getLevelElem());
    }
    for (const auto& v : Dune::vertices(serial.leafGridView())) {
        pointIds.insert(ids.id(v));
    }
    BOOST_CHECK_EQUAL(leafDiffers, 0);
    BOOST_CHECK_EQUAL(cellIds.size(), static_cast<std::size_t>(serial.size(0)));
    BOOST_CHECK_EQUAL(pointIds.size(), static_cast<std::size_t>(serial.size(3)));
    BOOST_CHECK(std::ranges::none_of(pointIds, [&cellIds](long id) { return cellIds.count(id) > 0; }));

    Dune::CpGrid grid;
    grid.createCartesian(dims, {{1.0, 1.0, 1.0}});
    const int np = grid.comm().size();
    if (np < 2) {
        return;
    }
    // 8x8 quadrants: LGR1 (and NEST1) on rank 0, LGR2 on rank 3 % np.
    std::vector<int> parts(static_cast<std::size_t>(dims[0])*dims[1]*dims[2]);
    for (std::size_t c = 0; c < parts.size(); ++c) {
        const int i = static_cast<int>(c % dims[0]);
        const int j = static_cast<int>((c / dims[0]) % dims[1]);
        parts[c] = (i/8 + 2*(j/8)) % np;
    }
    grid.loadBalance(parts, false, true, 2);
    refine(grid);

    const auto reference = leafIdsByCentre(serial);
    const auto mine = leafIdsByCentre(grid);
    BOOST_CHECK_EQUAL(mismatches(grid, mine, reference), 0);
    BOOST_CHECK_EQUAL(grid.comm().sum(static_cast<int>(mine.size())), static_cast<int>(reference.size()));
    BOOST_CHECK_EQUAL(mismatches(grid, leafVertexIdsByPosition(grid), leafVertexIdsByPosition(serial)), 0);
}

// Corners of the split faces on a faulted box boundary exist only on the leaf; their ids
// must not depend on the partition either.
BOOST_AUTO_TEST_CASE(splitCornerIdsIndependentOfPartition)
{
    const std::array<int,3> dims = {{8, 8, 2}};
    const auto g = makeFaultedGrid(dims, 4, 0.6);
    std::ostringstream deckString;
    deckString << "RUNSPEC\nDIMENS\n 8 8 2 /\nGRID\nCOORD\n";
    for (const double v : g.coord) {
        deckString << ' ' << v;
    }
    deckString << " /\nZCORN\n";
    for (const double v : g.zcorn) {
        deckString << ' ' << v;
    }
    deckString << " /\nPORO\n 128*0.2 /\nCARFIN\n'LGR1' 2 4 1 2 1 2 6 4 4 /\nENDFIN\n";
    const auto deck = Opm::Parser{}.parseString(deckString.str());
    Opm::EclipseState state(deck);
    auto eclGrid = state.getInputGrid();
    const auto build = [&](Dune::CpGrid& grid) {
        grid.setLgrBackend(Opm::Refinement::Backend::Conforming);
        grid.processEclipseFormat(&eclGrid, &state, false, false, false);
    };
    const auto refine = [](Dune::CpGrid& grid) {
        grid.addLgrsUpdateLeafView({{2,2,2}}, {{1,0,0}}, {{4,2,2}}, {"LGR1"}); // high side on the fault
    };

    Dune::CpGrid serial(MPI_COMM_SELF);
    build(serial);
    refine(serial);
    BOOST_REQUIRE_GT(static_cast<std::size_t>(serial.size(3)),
                     serial.currentData().back()->cornerHistorySize()); // split corners exist

    Dune::CpGrid grid;
    const int np = grid.comm().size();
    if (np < 2) {
        return;
    }
    build(grid);
    std::vector<int> parts(static_cast<std::size_t>(dims[0])*dims[1]*dims[2]);
    for (std::size_t c = 0; c < parts.size(); ++c) {
        const int j = static_cast<int>((c / dims[0]) % dims[1]);
        parts[c] = (j < 4) ? 0 : 1 + (j - 4) % (np - 1);
    }
    grid.loadBalance(parts, false, true, 2);
    refine(grid);

    BOOST_CHECK_EQUAL(mismatches(grid, leafIdsByCentre(grid), leafIdsByCentre(serial)), 0);
    BOOST_CHECK_EQUAL(mismatches(grid, leafVertexIdsByPosition(grid), leafVertexIdsByPosition(serial)), 0);
}

BOOST_AUTO_TEST_CASE(boxTouchingOverlapThrows)
{
    const std::array<int,3> dims = {{12, 4, 2}};
    auto g = makeUnitGrid(dims);

    Dune::CpGrid grid;
    grid.createCartesian(dims, {{1.0, 1.0, 1.0}}); // unit cells, matching makeUnitGrid

    const int np = grid.comm().size();
    if (np != 2) {
        return;
    }

    // Split exactly at i=6; a box straddling i=6 crosses the rank boundary.
    std::vector<int> parts(static_cast<std::size_t>(dims[0])*dims[1]*dims[2]);
    for (int k = 0; k < dims[2]; ++k) {
        for (int j = 0; j < dims[1]; ++j) {
            for (int i = 0; i < dims[0]; ++i) {
                parts[i + dims[0]*j + dims[0]*dims[1]*k] = (i < 6) ? 0 : 1;
            }
        }
    }
    grid.loadBalance(parts, false, true, 2);

    useConformingBuilder(grid, std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
        g.dims, g.coord, g.zcorn, g.actnum));

    // Box i in [5,8) straddles the i=6 partition boundary -> must throw.
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView({{2,2,2}}, {{5,1,0}}, {{8,3,2}}, {"LGR1"}),
                      std::runtime_error);  // rethrown collectively on every rank
}
#endif // HAVE_MPI

// A deck-built Conforming grid refines after loadBalance() from the input rank 0 retained.
BOOST_AUTO_TEST_CASE(retainedInputReachesEveryRank)
{
    const std::string deckString = R"(RUNSPEC
DIMENS
 10 6 3 /
GRID
CARFIN
'LGR1' 2 3 3 3 1 1 4 1 1 /
ENDFIN
DX
 180*1 /
DY
 180*1 /
DZ
 180*1 /
TOPS
 60*0 /
PORO
 180*0.2 /
)";
    const auto deck = Opm::Parser{}.parseString(deckString);
    Opm::EclipseState state(deck);
    auto eclGrid = state.getInputGrid();

    Dune::CpGrid grid;
    if (grid.comm().size() < 2) {
        return;
    }
    grid.setLgrBackend(Opm::Refinement::Backend::Conforming);
    grid.processEclipseFormat(&eclGrid, &state, false, false, false);
    // The box (i = 1..2) stays inside rank 0's block i < 6, clear of its overlap.
    const int np = grid.comm().size();
    std::vector<int> parts(10*6*3);
    for (std::size_t c = 0; c < parts.size(); ++c) {
        const int i = static_cast<int>(c % 10);
        parts[c] = (i < 6) ? 0 : 1 + (i - 6) % (np - 1);
    }
    grid.loadBalance(parts, false, true, 2);
    grid.addLgrsUpdateLeafView({{2, 1, 1}}, {{1, 2, 0}}, {{3, 3, 1}}, {"LGR1"});

    int interior = 0;
    int refined = 0;
    for (const auto& element : Dune::elements(grid.leafGridView(), Dune::Partitions::interior)) {
        ++interior;
        refined += element.hasFather() ? 1 : 0;
    }
    BOOST_CHECK_EQUAL(grid.comm().sum(interior), 180 - 2 + 4);
    BOOST_CHECK_EQUAL(grid.comm().sum(refined), 4);
}
