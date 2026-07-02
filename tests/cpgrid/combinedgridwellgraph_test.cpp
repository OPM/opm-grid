// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>

#include <dune/common/version.hh>

#define BOOST_TEST_MODULE CombinedGridWellGraph
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/GridEnums.hpp>
#include <opm/grid/common/ZoltanGraphFunctions.hpp>
#include <opm/grid/utility/OpmWellType.hpp>

#if HAVE_OPM_COMMON
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#endif

#include <algorithm>
#include <numeric>

#if HAVE_OPM_COMMON
namespace {
    // create Wells, we only use well name and cell locations
    auto createConnection (int i, int j, int k)
    {
        return Opm::Connection(i,j,k,0, 0,Opm::Connection::State::OPEN,
                                   Opm::Connection::Direction::Z,
                                   Opm::Connection::CTFKind::DeckValue, 0,
                                   5.,Opm::Connection::CTFProperties(),0,false);
    }
    auto createWell (const std::string& name)
    {
        using namespace Opm;
        return Dune::cpgrid::OpmWellType(name,name,0,0,0,0,0.,WellType(),
                   Well::ProducerCMode(),Connection::Order(),UnitSystem(),
                   0.,0.,false,false,0,Well::GasInflowEquation());
    };
} // end anonymous namespace
#endif

BOOST_AUTO_TEST_CASE(CombinedGridWellGraph)
{
    ///! beware, test functionality, not implementation

    /// construct CombinedGridWellGraph
    Dune::CpGrid grid;
    std::array<int, 3> dims { 3, 3, 1 };
    std::array<double, 3> size { 1., 1., 1. };
    grid.createCartesian(dims, size);

    auto wellCon = std::make_shared<Opm::WellConnections>(); // do not confuse with Dune::cpgrid::WellConnections
    std::vector<Dune::cpgrid::OpmWellType> wells;
    wellCon->add(createConnection(0,0,0));
    wellCon->add(createConnection(1,0,0));
    wellCon->add(createConnection(2,0,0));
    wells.push_back(createWell("first_row"));
    wells[0].updateConnections(wellCon,true);

    wellCon = std::make_shared<Opm::WellConnections>(); // reset
    wellCon->add(createConnection(0,0,0));
    wellCon->add(createConnection(1,1,0));
    wellCon->add(createConnection(2,2,0));
    wells.push_back(createWell("diag"));
    wells[1].updateConnections(wellCon,true);

    // std::vector<OpmWellType> wells;
    std::unordered_map<std::string, std::set<int>> possibleFutureConnections;
    std::vector<double> transmissibilities{0., 1., 2., 0., 0., 3., 4., 0., 0., 5., 6., 0., // X
                                           0., 0., 0.,11.,12.,13.,14.,15.,16., 0., 0., 0., // Y
                                           0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}; // Z
    Dune::cpgrid::CombinedGridWellGraph gridWellGraph(grid, &wells, possibleFutureConnections,
        transmissibilities.data(), false, Dune::EdgeWeightMethod::defaultTransEdgeWgt);
    const auto& wellEdges = gridWellGraph.getWellsGraph();

    /// log transmissibilities do not affect well connections
    {
        Dune::cpgrid::CombinedGridWellGraph logGridWellGraph(grid, &wells, possibleFutureConnections,
            transmissibilities.data(), false, Dune::EdgeWeightMethod::logTransEdgeWgt);
        const auto& logWellEdges = logGridWellGraph.getWellsGraph();
        BOOST_REQUIRE(wellEdges == logWellEdges);
    }

    /// connections of every well are interconnected
    BOOST_REQUIRE(wellEdges.size() == 9); // all vertices
    std::set<int> edges{1, 2, 4, 8}; // both wells have the cell 0
    BOOST_REQUIRE(wellEdges[0].size() == 4);
    BOOST_REQUIRE(wellEdges[0] == edges);
    edges = std::set<int>{0,2};
    BOOST_REQUIRE(wellEdges[1].size() == 2);
    BOOST_REQUIRE(wellEdges[1] == edges);
    edges = std::set<int>{0,1};
    BOOST_REQUIRE(wellEdges[2].size() == 2);
    BOOST_REQUIRE(wellEdges[2] == edges);
    BOOST_REQUIRE(wellEdges[3].size() == 0);
    edges = std::set<int>{0,8};
    BOOST_REQUIRE(wellEdges[4].size() == 2);
    BOOST_REQUIRE(wellEdges[4] == edges);
    BOOST_REQUIRE(wellEdges[5].size() == 0);
    BOOST_REQUIRE(wellEdges[6].size() == 0);
    BOOST_REQUIRE(wellEdges[7].size() == 0);
    edges = std::set<int>{0,4};
    BOOST_REQUIRE(wellEdges[8].size() == 2);
    BOOST_REQUIRE(wellEdges[8] == edges);

    /// well edge weight is by defaul equal to the sum of all grid edges
    auto wellEdgeWeight = calculateWellEdgeWeight<float>(grid, gridWellGraph);
    // the sum (102) gets multipied by 1e18 to accomodate partitioners that use integral weights (Metis)
    BOOST_REQUIRE(std::abs(wellEdgeWeight / 1e18 - 102.) < 1e-5);

    // a positive multiplier changes the weight from sum_of_faces to multiplier_\times_average_face
    gridWellGraph.setMultiplyWellConnectivities(21.); // there are 42 grid edges in total
    wellEdgeWeight = calculateWellEdgeWeight<float>(grid, gridWellGraph);
    BOOST_REQUIRE(std::abs(wellEdgeWeight / 1e18 - 51.) < 1e-5);

    /// edges from well connections overwrite grid edges
    int neighborCounter{0};
    std::vector<uint> gID(9), nborGID(32, 0); // 32=2x(12 inner faces + 4 new connections due to wells)
    ZOLTAN_ID_PTR nborGIDData = nborGID.data(), gIDData = gID.data();
    std::vector<float> ewgts(32, 0.);
    for (int i=0; i<9; ++i) {
        gID[i] = i;
        fillNBORGIDAndWeightsForSpecificCellAndIncrementNeighborCounterForGridWithWells(gridWellGraph, i, gIDData, neighborCounter, nborGIDData, ewgts.data(), wellEdgeWeight);
    }

    BOOST_REQUIRE(neighborCounter == 32);
    // check if the well edge weights are at correct places
    for (const auto i : std::vector<int>{0,1,2,3, 5,6, 8,9, 14,15, 28,29}) {
        BOOST_REQUIRE(ewgts[i] == wellEdgeWeight);
    }
    const auto weightSum = std::accumulate(ewgts.begin(), ewgts.end(), 0.);
    // sum of edge weights is 102, there are 6 well connections in total and grid edges {1,2} were overwritten
    BOOST_REQUIRE(std::abs(weightSum/1e18 - 2*(102 + 6*wellEdgeWeight/1e18 - 1 - 2)) < 1e-3);
}

bool init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    boost::unit_test::unit_test_main(&init_unit_test_func,
        argc, argv);
}
