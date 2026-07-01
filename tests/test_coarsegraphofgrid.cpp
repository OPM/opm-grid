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

#define BOOST_TEST_MODULE GraphRepresentationOfGrid
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/WellConnections.hpp>
#include <opm/grid/GraphOfGrid.hpp>
#include <opm/grid/CoarseGraphOfGrid.hpp>
#include <opm/grid/GraphOfGridWrappers.hpp>
#include <opm/grid/utility/OpmWellType.hpp>

#if HAVE_OPM_COMMON
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#endif

#include <algorithm>

// Create "Transmissibility graph" from grid and set all transmissibilities to 1.0 
std::unique_ptr<Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>> createAdjacencyMatrix(const Dune::CpGrid& grid) 
{
    size_t numCells = grid.numCells();
    auto gridView = grid.leafGridView();

    using ElementMapper =
      Dune::MultipleCodimMultipleGeomTypeMapper<typename Dune::CpGrid::LeafGridView>;
    const auto elemMapper = ElementMapper { gridView, Dune::mcmgElementLayout() };

    Dune::MatrixIndexSet op;
    op.resize(numCells, numCells);

    // Create adjacency for transmissibility graph
    for (const auto& elem : elements(gridView, Dune::Partitions::interiorBorder)) {
        auto d = elemMapper.index(elem);
        op.add(d, d);
        for (const auto& is : intersections(gridView, elem)) {
            if (!is.neighbor()) {
                continue;
            }

            const auto J = static_cast<unsigned int>(elemMapper.index(is.outside()));
            op.add(d, J);
        }
    }

    // Allocate the BCRSMatrix wrapped in a unique_ptr
    auto graph = std::make_unique<Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>>();
    op.exportIdx(*graph);
    *graph = 1.0;

    return graph;
}

// basic test to check if the graph was constructed correctly. Copied from test_graphofgrid.cpp
BOOST_AUTO_TEST_CASE(SimpleGraph)
{
    Dune::CpGrid grid;
    std::array<int,3> dims{2,2,2};
    std::array<double,3> size{2.,2.,2.};
    grid.createCartesian(dims,size);
    if (grid.size(0)==0)
        return;

    Dune::cpgrid::WellConnections wells;

    auto graph = createAdjacencyMatrix(grid);
    double strongConnectionThreshold = 1.1; // no nodes will be merged
    int maxNodeSize = 2;
    bool allowDistWells = false;
    Opm::CoarseGraphOfGrid<Dune::CpGrid> cgog(grid, Dune::EdgeWeightMethod::uniformEdgeWgt,
                                              graph.get(), strongConnectionThreshold,
                                              maxNodeSize, allowDistWells, wells);

    BOOST_REQUIRE(cgog.cSize()==8); // number of graph vertices
    BOOST_REQUIRE(cgog.getCoarseNodes()[0].size()==1); // weight of all nodes is one 

    const std::vector< std::map<int, double> > edges = cgog.getCoarseEdges();

    BOOST_REQUIRE(edges[0].size()==3);
    auto edgeL = edges[2];

    BOOST_REQUIRE(edgeL.size()==3); // neighbors of vertex 2 are: 0, 3, 6
    BOOST_REQUIRE(edgeL[0]==1.0); 
    BOOST_REQUIRE(edgeL[3]==1.0);
    BOOST_REQUIRE(edgeL[6]==1.0);
    BOOST_REQUIRE_THROW(edgeL.at(4),std::out_of_range); // not a neighbor

    BOOST_REQUIRE_THROW(edges.at(10),std::logic_error);
}

// basic test to check if nodes get merged.
//Similar to SimpleGraphWithVertexContraction in test_graphofgrid.cpp
BOOST_AUTO_TEST_CASE(MergeTwoNodes)
{
    Dune::CpGrid grid;
    std::array<int,3> dims{2,2,2};
    std::array<double,3> size{2.,2.,2.};
    grid.createCartesian(dims,size);
    if (grid.size(0)==0)
        return;

    Dune::cpgrid::WellConnections wells;

    auto graph = createAdjacencyMatrix(grid);
    (*graph)[0][1] = 2.0; // Set 0->1 connection to 2
    double strongConnectionThreshold = 1.1; // Node 0 and 1 will be merged
    int maxNodeSize = 2;
    bool allowDistWells = false;
    Opm::CoarseGraphOfGrid<Dune::CpGrid> cgog(grid, Dune::EdgeWeightMethod::uniformEdgeWgt,
                                              graph.get(), strongConnectionThreshold,
                                              maxNodeSize, allowDistWells, wells);

    BOOST_REQUIRE(cgog.cSize()==7); // number of graph vertices
    BOOST_REQUIRE(cgog.getCoarseNodes()[0].size()==2); // weight of node 0 is two 

    const auto edges = cgog.getCoarseEdges();

    auto edge0 = edges[0];
    BOOST_REQUIRE(edge0.size()==4);

    auto map_to_coarse = cgog.getMapToCoarse();

    BOOST_REQUIRE(map_to_coarse[0]==0);
    BOOST_REQUIRE(map_to_coarse[1]==0);
    BOOST_REQUIRE(map_to_coarse[2]==1);

    BOOST_REQUIRE(edge0[map_to_coarse[2]]==1.0);
    BOOST_REQUIRE(edge0[map_to_coarse[3]]==1.0);
    BOOST_REQUIRE(edge0[map_to_coarse[4]]==1.0);
    BOOST_REQUIRE(edge0[map_to_coarse[5]]==1.0);
}

// basic test to check if the coarse graph gets correct edge weights
BOOST_AUTO_TEST_CASE(CoarseEdgeWeights)
{
    Dune::CpGrid grid;
    std::array<int,3> dims{3,3,1};
    std::array<double,3> size{2.,2.,2.};
    grid.createCartesian(dims,size);
    if (grid.size(0)==0)
        return;

    Dune::cpgrid::WellConnections wells;

    auto graph = createAdjacencyMatrix(grid);
    (*graph)[0][1] = 2.0;
    (*graph)[0][3] = 2.0;
    (*graph)[7][8] = 0.0; // Fault between cell 7 and 8
    (*graph)[8][7] = 0.0;

    double strongConnectionThreshold = 1.1; 
    int maxNodeSize = 3;
    bool allowDistWells = false;
    Opm::CoarseGraphOfGrid<Dune::CpGrid> cgog(grid, Dune::EdgeWeightMethod::uniformEdgeWgt,
                                              graph.get(), strongConnectionThreshold,
                                              maxNodeSize, allowDistWells, wells);

    BOOST_REQUIRE(cgog.cSize()==7); // number of graph vertices
    BOOST_REQUIRE(cgog.getCoarseNodes()[0].size()==3); // Node weight of node 0 is 3 

    const std::vector< std::map<int, double> > edges = cgog.getCoarseEdges();
    auto edge0 = edges[0];

    BOOST_REQUIRE(edge0.size()==3); // The merged node has three edges

    auto map_to_coarse = cgog.getMapToCoarse();
    BOOST_REQUIRE(edge0[map_to_coarse[2]]==1.0);
    BOOST_REQUIRE(edge0[map_to_coarse[4]]==2.0); //Node 0 has two connections to the centre node
    BOOST_REQUIRE(edge0[map_to_coarse[6]]==1.0);

    auto edge7 = edges[map_to_coarse[7]];
    auto edge8 = edges[map_to_coarse[8]];

    BOOST_REQUIRE(edge7.size()==2); // Only two edges because of fault
    BOOST_REQUIRE(edge8.size()==1); // Only one edges because of fault

    BOOST_REQUIRE_THROW(edge7.at(map_to_coarse[8]),std::out_of_range); // not a neighbor
    BOOST_REQUIRE_THROW(edge8.at(map_to_coarse[7]),std::out_of_range); // not a neighbor
}


// basic test to check if the graph uses maxNodeSize correctly
BOOST_AUTO_TEST_CASE(MaxNodeSize)
{
    Dune::CpGrid grid;
    std::array<int,3> dims{3,3,1};
    std::array<double,3> size{2.,2.,2.};
    grid.createCartesian(dims,size);
    if (grid.size(0)==0)
        return;

    Dune::cpgrid::WellConnections wells;

    auto graph = createAdjacencyMatrix(grid);
    (*graph)[0][1] = 2.0;
    (*graph)[0][3] = 3.0; // This is lager, so should be merged
    (*graph)[6][7] = 2.0;
    (*graph)[7][8] = 2.0;

    double strongConnectionThreshold = 1.1;
    int maxNodeSize = 2;
    bool allowDistWells = false;
    Opm::CoarseGraphOfGrid<Dune::CpGrid> cgog(grid, Dune::EdgeWeightMethod::uniformEdgeWgt,
                                              graph.get(), strongConnectionThreshold,
                                              maxNodeSize, allowDistWells, wells);

    BOOST_REQUIRE(cgog.cSize()==7); // number of graph vertices
    BOOST_REQUIRE(cgog.getCoarseNodes()[0].size()==2); // weight of node 0 is two

    auto map_to_coarse = cgog.getMapToCoarse();
    BOOST_REQUIRE(map_to_coarse[0]==map_to_coarse[3]); // node 0 merged with node 3
    BOOST_REQUIRE(map_to_coarse[0]!=map_to_coarse[1]); // node 1 not merged with node 0 and 3

    BOOST_REQUIRE(cgog.getCoarseNodes()[map_to_coarse[6]].size()==2); // weight of node 6 is two
    BOOST_REQUIRE(map_to_coarse[6]==map_to_coarse[7]); // node 6 merged with node 7
    BOOST_REQUIRE(map_to_coarse[7]!=map_to_coarse[8]); // node 8 not merged with node 6 and 7
}

#if HAVE_OPM_COMMON
namespace {
    // create Wells, we only use well name and cell locations
    // copied from test_graphofgrid.cpp
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

#if HAVE_MPI && HAVE_OPM_COMMON
// Check if merging wells work
// Similar to addWellConnections in test_graphofgrid.cpp
BOOST_AUTO_TEST_CASE(MergeWells)
{
    Dune::CpGrid grid;
    std::array<int,3> dims{2,2,2};
    std::array<double,3> size{1.,1.,1.};
    grid.createCartesian(dims,size);
    if (grid.size(0)==0)
        return;

    auto wellCon = std::make_shared<Opm::WellConnections>(); // do not confuse with Dune::cpgrid::WellConnections
    wellCon->add(createConnection(0,0,0)); // 0
    wellCon->add(createConnection(0,1,0)); // 2
    wellCon->add(createConnection(0,1,1)); // 6
    std::vector<Dune::cpgrid::OpmWellType> wells;
    wells.push_back(createWell("first"));
    wells[0].updateConnections(wellCon,true);

    wellCon = std::make_shared<Opm::WellConnections>(); //reset
    wellCon->add(createConnection(0,0,1)); // 4
    wellCon->add(createConnection(1,1,0)); // 3
    wells.push_back(createWell("second"));
    wells[1].updateConnections(wellCon,true);

    wellCon = std::make_shared<Opm::WellConnections>(); //reset
    wellCon->add(createConnection(0,0,1)); // 4
    wellCon->add(createConnection(1,0,1)); // 5
    wells.push_back(createWell("third")); // intersects with second
    wells[2].updateConnections(wellCon,true);

    Dune::cpgrid::WellConnections wellConnections(wells,std::unordered_map<std::string, std::set<int>>(),grid);

    auto graph = createAdjacencyMatrix(grid);
    (*graph)[0][1] = 2.0; // should not be merged, due to maxNodeSize = 2

    double strongConnectionThreshold = 1.1;
    int maxNodeSize = 2;
    bool allowDistWells = false;
    Opm::CoarseGraphOfGrid<Dune::CpGrid> cgog(grid, Dune::EdgeWeightMethod::uniformEdgeWgt,graph.get(),
                                              strongConnectionThreshold,
                                              maxNodeSize, allowDistWells, wellConnections);

    BOOST_REQUIRE(cgog.cSize()==4);

    int err;
    int nVer = getCoarseGraphNumVertices(&cgog,&err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(nVer == 4);
    std::vector<uint> gIDs(nVer);
    std::vector<float> objWeights(nVer);
    getCoarseGraphVerticesList(&cgog, 1, 1, gIDs.data(), nullptr, 1, objWeights.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    std::ranges::sort(gIDs);
    BOOST_REQUIRE(objWeights[0]==3.0 && objWeights[1]==1.0 && objWeights[2]==3.0 && objWeights[3]==1.0);
    std::vector<int> numEdges(nVer);
    getCoarseGraphNumEdges(&cgog, 1, 1, nVer, gIDs.data(), nullptr, numEdges.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(numEdges[0]==3 && numEdges[1]==2 && numEdges[2]==3 && numEdges[3]==2);
    int nEdges = 10; // sum of numEdges[i]
    std::vector<uint> nborGIDs(nEdges);
    std::vector<int> nborProc(nEdges);
    std::vector<float> edgeWeights(nEdges);
    getCoarseGraphEdgeList(&cgog, 1, 1, nVer, gIDs.data(), nullptr, numEdges.data(), nborGIDs.data(), nborProc.data(), 1, edgeWeights.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);

    // check all edgeWeights. Note that nborGIDs are not sorted
    for (int i=0; i<3; ++i)
    {
        switch (nborGIDs[i])
        {
            case 1: BOOST_REQUIRE(edgeWeights[i]==1); break;
            case 2: BOOST_REQUIRE(edgeWeights[i]==3); break;
            case 3: BOOST_REQUIRE(edgeWeights[i]==1); break;
            default: throw("CoarseGraph was constructed badly.");
        }
    }
    for (int i=3; i<5; ++i)
    {
        switch (nborGIDs[i])
        {
            case 0: BOOST_REQUIRE(edgeWeights[i]==1); break;
            case 2: BOOST_REQUIRE(edgeWeights[i]==2); break;
            default: throw("CoarseGraph was constructed badly.");
        }
    }
    for (int i=5; i<8; ++i)
    {
        switch (nborGIDs[i])
        {
            case 0: BOOST_REQUIRE(edgeWeights[i]==3); break;
            case 1: BOOST_REQUIRE(edgeWeights[i]==2); break;
            case 3: BOOST_REQUIRE(edgeWeights[i]==2); break;
            default: throw("CoarseGraph was constructed badly.");
        }
    }
    for (int i=8; i<10; ++i)
    {
        switch (nborGIDs[i])
        {
            case 0: BOOST_REQUIRE(edgeWeights[i]==1); break;
            case 2: BOOST_REQUIRE(edgeWeights[i]==2); break;
            default: throw("CoarseGraph was constructed badly.");
        }
    }
}

// Check if merging multiple intersecting wells work
// Similar to IntersectingWells in test_graphofgrid.cpp
BOOST_AUTO_TEST_CASE(MergeWellsMoreIntersecting)
{
    Dune::CpGrid grid;
    std::array<int,3> dims{5,4,3};
    std::array<double,3> size{1.,1.,1.};
    grid.createCartesian(dims,size);
    if (grid.size(0)==0)
        return;

    auto wellCon = std::make_shared<Opm::WellConnections>(); // do not confuse with Dune::cpgrid::WellConnections
    wellCon->add(createConnection(0,0,0)); // 0
    wellCon->add(createConnection(1,0,0)); // 1
    wellCon->add(createConnection(2,0,0)); // 2
    wellCon->add(createConnection(3,0,0)); // 3
    wellCon->add(createConnection(4,0,0)); // 4
    std::vector<Dune::cpgrid::OpmWellType> wells;
    wells.push_back(createWell("first"));
    wells[0].updateConnections(wellCon,true);

    wellCon = std::make_shared<Opm::WellConnections>(); //reset
    wellCon->add(createConnection(2,2,2)); // 52
    wellCon->add(createConnection(2,2,1)); // 32
    wellCon->add(createConnection(2,2,0)); // 12
    wells.push_back(createWell("second"));
    wells[1].updateConnections(wellCon,true);

    wellCon = std::make_shared<Opm::WellConnections>(); //reset
    wellCon->add(createConnection(4,3,2)); // 59
    wellCon->add(createConnection(3,1,2)); // 48
    wellCon->add(createConnection(2,3,1)); // 37
    wells.push_back(createWell("third")); // 
    wells[2].updateConnections(wellCon,true);

    wellCon = std::make_shared<Opm::WellConnections>(); //reset
    wellCon->add(createConnection(2,3,1)); // 37
    wellCon->add(createConnection(3,3,1)); // 38
    wellCon->add(createConnection(4,3,1)); // 39
    wellCon->add(createConnection(4,2,1)); // 34
    wells.push_back(createWell("forth"));
    wells[3].updateConnections(wellCon,true);

    wellCon = std::make_shared<Opm::WellConnections>(); //reset
    wellCon->add(createConnection(2,0,0)); // 2
    wellCon->add(createConnection(3,1,0)); // 8
    wells.push_back(createWell("fifth"));
    wells[4].updateConnections(wellCon,true);

    wellCon = std::make_shared<Opm::WellConnections>(); //reset
    wellCon->add(createConnection(2,0,0)); // 2
    wellCon->add(createConnection(3,3,1)); // 38
    wells.push_back(createWell("sixth"));
    wells[5].updateConnections(wellCon,true);
    Dune::cpgrid::WellConnections wellConnections(wells,std::unordered_map<std::string, std::set<int>>(),grid);

    auto graph = createAdjacencyMatrix(grid);

    double strongConnectionThreshold = 1.1;
    int maxNodeSize = 2;
    bool allowDistWells = false;
    Opm::CoarseGraphOfGrid<Dune::CpGrid> cgog(grid, Dune::EdgeWeightMethod::uniformEdgeWgt,
                                              graph.get(), strongConnectionThreshold,
                                              maxNodeSize, allowDistWells, wellConnections);

    auto map_to_coarse = cgog.getMapToCoarse();
    BOOST_REQUIRE(cgog.cSize() == 47);

    int err;
    int nVer = getCoarseGraphNumVertices(&cgog,&err);
    BOOST_REQUIRE(nVer == 47);
    std::vector<uint> gIDs(nVer);
    std::vector<float> objWeights(nVer);
    getCoarseGraphVerticesList(&cgog, 1, 1, gIDs.data(), nullptr, 1, objWeights.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);

    for (int i=0; i<nVer; ++i) {
        switch (gIDs[i])
        {
            case 0:  BOOST_REQUIRE(objWeights[i]==12.); break;
            case 7:  BOOST_REQUIRE(objWeights[i]==3.); break;
            default: BOOST_REQUIRE(objWeights[i]==1.); // ordinary vertex
        }
    }

    std::vector<int> numEdges(nVer);
    getCoarseGraphNumEdges(&cgog, 1, 1, nVer, gIDs.data(), nullptr, numEdges.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(numEdges[7]==12);
    BOOST_REQUIRE(numEdges[0]==26);
    BOOST_REQUIRE(numEdges[42]==3);

    int nEdges = 232;
    std::vector<uint> nborGIDs(nEdges);
    std::vector<int> nborProc(nEdges);
    std::vector<float> edgeWeights(nEdges);
    getCoarseGraphEdgeList(&cgog, 1, 1, nVer, gIDs.data(), nullptr, numEdges.data(), nborGIDs.data(), nborProc.data(), 1, edgeWeights.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);

    int checked=0;
    for (int i=0; i<26; ++i)
    {
        switch (nborGIDs[i])
        {
            // neighboring two well cells adds up the edge weight
            case 3: case 4: case 23: case 27: case 42: case 46:
                BOOST_REQUIRE(edgeWeights[i]==2.);
                ++checked;
                break;
            default: BOOST_REQUIRE(edgeWeights[i]==1.);
        }
    }
    BOOST_REQUIRE(checked==6);
}

// Check if merging wells and transmissibility works
BOOST_AUTO_TEST_CASE(MergeWithTransAndWells)
{
    Dune::CpGrid grid;
    std::array<int,3> dims{4,4,1};
    std::array<double,3> size{1.,1.,1.};
    grid.createCartesian(dims,size);
    if (grid.size(0)==0)
        return;

    auto wellCon = std::make_shared<Opm::WellConnections>(); // do not confuse with Dune::cpgrid::WellConnections
    wellCon->add(createConnection(0,0,0)); // 0
    wellCon->add(createConnection(1,0,0)); // 1
    wellCon->add(createConnection(2,0,0)); // 2
    std::vector<Dune::cpgrid::OpmWellType> wells;
    wells.push_back(createWell("first"));
    wells[0].updateConnections(wellCon,true);

    wellCon = std::make_shared<Opm::WellConnections>(); //reset
    wellCon->add(createConnection(0,2,0)); // 8 
    wellCon->add(createConnection(1,2,0)); // 9
    wellCon->add(createConnection(2,2,0)); // 10
    wells.push_back(createWell("second"));
    wells[1].updateConnections(wellCon,true);

    wellCon = std::make_shared<Opm::WellConnections>(); //reset
    wellCon->add(createConnection(3,2,0)); // 11 
    wellCon->add(createConnection(3,3,0)); // 15
    wells.push_back(createWell("third"));
    wells[2].updateConnections(wellCon,true);

    Dune::cpgrid::WellConnections wellConnections(wells,std::unordered_map<std::string, std::set<int>>(),grid);

    auto graph = createAdjacencyMatrix(grid);
    (*graph)[0][4] = 2.0; // should be merged
    (*graph)[4][8] = 2.0; // should be merged and add second well too
    (*graph)[2][3] = 2.0; // should be merged
    (*graph)[3][7] = 2.0; // should be merged
    (*graph)[11][10] = 2.0; // not merged because 10 aready merged. If graph[10][11]=2, would be merged

    // How coarsening happens in this example 
    //
    // +---+---+---+---+       +---+---+---+---+
    // | w1| w1| w1|   |       |               |
    // +---+---+---+---+       +   +---+---+   +
    // |   |   |   |   |       |   |   |   |   |
    // +---+---+---+---+   ->  +   +---+---+---+
    // | w2| w2| w2| w3|       |           |   |
    // +---+---+---+---+       +---+---+---+   +
    // |   |   |   | w3|       |   |   |   |   |
    // +---+---+---+---+       +---+---+---+---+

    int maxNodeSize = 14;
    bool allowDistWells = false;
    Opm::CoarseGraphOfGrid<Dune::CpGrid> cgog(grid, Dune::EdgeWeightMethod::uniformEdgeWgt,graph.get(),
                                              1.1, maxNodeSize, allowDistWells, wellConnections);

    BOOST_REQUIRE(cgog.getCoarseNodes()[0].size()==9);
    auto map_to_coarse = cgog.getMapToCoarse();

    int err;
    int nVer = getCoarseGraphNumVertices(&cgog,&err);
    BOOST_REQUIRE(nVer == 7);
    std::vector<uint> gIDs(nVer);
    std::vector<float> objWeights(nVer);
    getCoarseGraphVerticesList(&cgog, 1, 1, gIDs.data(), nullptr, 1, objWeights.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);

    for (int i=0; i<nVer; ++i) {
        switch (gIDs[i])
        {
            case 0:  BOOST_REQUIRE(objWeights[i]==9.); break;
            case 3:  BOOST_REQUIRE(objWeights[i]==2.); break;
            default: BOOST_REQUIRE(objWeights[i]==1.); // ordinary vertex
        }
    }

    std::vector<int> numEdges(nVer);
    getCoarseGraphNumEdges(&cgog, 1, 1, nVer, gIDs.data(), nullptr, numEdges.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(numEdges[0]==6);
    BOOST_REQUIRE(numEdges[1]==2);
    BOOST_REQUIRE(numEdges[2]==2);
    BOOST_REQUIRE(numEdges[3]==2);
    BOOST_REQUIRE(numEdges[4]==2);
    BOOST_REQUIRE(numEdges[5]==3);
    BOOST_REQUIRE(numEdges[6]==3);

    int nEdges = 20;
    std::vector<uint> nborGIDs(nEdges);
    std::vector<int> nborProc(nEdges);
    std::vector<float> edgeWeights(nEdges);
    getCoarseGraphEdgeList(&cgog, 1, 1, nVer, gIDs.data(), nullptr, numEdges.data(), nborGIDs.data(), nborProc.data(), 1, edgeWeights.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);

    // Check edgeWeights for the first node
    for (int i=0; i<6; ++i)
    {
        switch (nborGIDs[i])
        {
            case 1:
                BOOST_REQUIRE(edgeWeights[i]==3.); break;
            case 2:
                BOOST_REQUIRE(edgeWeights[i]==3.); break;
            case 3:
                BOOST_REQUIRE(edgeWeights[i]==2.); break;
            default: BOOST_REQUIRE(edgeWeights[i]==1.);
        }
    }
}
#endif

bool
init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    boost::unit_test::unit_test_main(&init_unit_test_func,
                                     argc, argv);
}
