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

#define BOOST_TEST_MODULE RefinedDistributionTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <algorithm>
#include <array>
#include <map>
#include <cmath>
#include <cstdint>
#include <string>

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

// Sends each cell's volume and centroid from the global leaf; the receiver compares with its own.
class GeometryCheckHandle
{
public:
    using DataType = double;
    bool contains(int, int codim) const { return codim == 0; }
    bool fixedSize(int, int) const { return true; }
    template <class T> std::size_t size(const T&) const { return 4; }
    template <class B, class T> void gather(B& buffer, const T& e) const
    {
        const auto geometry = e.geometry();
        buffer.write(geometry.volume());
        for (const auto x : geometry.center()) {
            buffer.write(x);
        }
    }
    template <class B, class T> void scatter(B& buffer, const T& e, std::size_t)
    {
        const auto geometry = e.geometry();
        double value = 0.0;
        buffer.read(value);
        double mismatch = std::abs(value - geometry.volume());
        for (const auto x : geometry.center()) {
            buffer.read(value);
            mismatch = std::max(mismatch, std::abs(value - x));
        }
        mismatches += (mismatch > 1e-10) ? 1 : 0;
        ++received;
    }
    int mismatches = 0;
    int received = 0;
};

} // anonymous namespace

// Refined before load balancing, the leaf is distributed whole and balanced by refined cells.
BOOST_AUTO_TEST_CASE(refinedLeafIsDistributedAndBalanced)
{
    // Half the grid refined 2x2x2: 288 coarse + 288*8 refined cells, unit volume each coarse cell.
    const std::string deckString = R"(RUNSPEC
DIMENS
 12 12 4 /
GRID
CARFIN
'LGR1' 1 6 1 12 1 4 12 24 8 /
ENDFIN
DX
 576*1 /
DY
 576*1 /
DZ
 576*1 /
TOPS
 144*0 /
PORO
 576*0.2 /
)";
    const auto deck = Opm::Parser{}.parseString(deckString);
    Opm::EclipseState state(deck);
    auto eclGrid = state.getInputGrid();

    Dune::CpGrid grid;
    grid.setLgrBackend(Opm::Refinement::Backend::Conforming);
    grid.processEclipseFormat(&eclGrid, &state, false, false, false);
    grid.addLgrsUpdateLeafView({{2, 2, 2}}, {{0, 0, 0}}, {{6, 12, 4}}, {"LGR1"});

    GeometryCheckHandle handle;
    if (grid.comm().size() > 1) {
        grid.loadBalance(handle, Dune::EdgeWeightMethod::uniformEdgeWgt, nullptr, {}, /*serialPartitioning=*/false);
    }

    int interior = 0;
    double volume = 0.0;
    for (const auto& element : Dune::elements(grid.leafGridView(), Dune::Partitions::interior)) {
        ++interior;
        volume += element.geometry().volume();
    }
    const auto& comm = grid.comm();
    BOOST_CHECK_EQUAL(comm.sum(interior), 288 + 288*8);
    BOOST_CHECK_CLOSE(comm.sum(volume), 576.0, 1e-10);
    if (comm.size() > 1) {
        BOOST_CHECK_EQUAL(comm.sum(handle.mismatches), 0);
        BOOST_CHECK_EQUAL(handle.received, grid.leafGridView().size(0));
        // Unweighted, the ranks holding the refined half would own about 8 times the cells.
        const double average = (288.0 + 288*8) / comm.size();
        BOOST_CHECK_LT(comm.max(interior), 1.3 * average);
        BOOST_CHECK_GT(comm.min(interior), 0.7 * average);
    }

    // Every cell keeps the id it has in the serially refined grid.
    Dune::CpGrid serial(MPI_COMM_SELF);
    serial.setLgrBackend(Opm::Refinement::Backend::Conforming);
    serial.processEclipseFormat(&eclGrid, &state, false, false, false);
    serial.addLgrsUpdateLeafView({{2, 2, 2}}, {{0, 0, 0}}, {{6, 12, 4}}, {"LGR1"});
    const auto byCentre = [](const Dune::CpGrid& g) {
        std::map<std::array<long,3>, long> ids;
        for (const auto& e : Dune::elements(g.leafGridView(), Dune::Partitions::interior)) {
            const auto x = e.geometry().center();
            ids[{std::lround(x[0]*1e6), std::lround(x[1]*1e6), std::lround(x[2]*1e6)}] = g.globalIdSet().id(e);
        }
        return ids;
    };
    const auto reference = byCentre(serial);
    int differ = 0;
    for (const auto& [x, id] : byCentre(grid)) {
        const auto it = reference.find(x);
        differ += it == reference.end() || it->second != id;
    }
    BOOST_CHECK_EQUAL(comm.sum(differ), 0);

    // And every vertex.
    const auto vertexIds = [](const Dune::CpGrid& g) {
        std::map<std::array<long,3>, long> ids;
        for (const auto& v : Dune::vertices(g.leafGridView())) {
            const auto x = v.geometry().center();
            ids[{std::lround(x[0]*1e6), std::lround(x[1]*1e6), std::lround(x[2]*1e6)}] = g.globalIdSet().id(v);
        }
        return ids;
    };
    const auto vertexReference = vertexIds(serial);
    int vertexDiffer = 0;
    for (const auto& [x, id] : vertexIds(grid)) {
        const auto it = vertexReference.find(x);
        vertexDiffer += it == vertexReference.end() || it->second != id;
    }
    BOOST_CHECK_EQUAL(comm.sum(vertexDiffer), 0);

    // So does stableCellId(), although the distributed leaf is flat.
    const auto stableByCentre = [](const Dune::CpGrid& g) {
        const auto sid = g.stableCellId();
        std::map<std::array<long,3>, std::int64_t> ids;
        for (const auto& e : Dune::elements(g.leafGridView(), Dune::Partitions::interior)) {
            const auto x = e.geometry().center();
            ids[{std::lround(x[0]*1e6), std::lround(x[1]*1e6), std::lround(x[2]*1e6)}] = sid[e.index()];
        }
        return ids;
    };
    const auto stableReference = stableByCentre(serial);
    int stableDiffer = 0;
    for (const auto& [x, id] : stableByCentre(grid)) {
        const auto it = stableReference.find(x);
        stableDiffer += it == stableReference.end() || it->second != id;
    }
    BOOST_CHECK_EQUAL(comm.sum(stableDiffer), 0);
}
