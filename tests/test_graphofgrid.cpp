// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2024 Equinor ASA.

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
#include <opm/grid/CpGrid.hpp>

#include <opm/grid/GraphOfGrid.hpp>
#include <opm/grid/GraphOfGridWrappers.hpp>

#include <opm/grid/utility/OpmWellType.hpp>
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

// basic test to check if the graph was constructed correctly
BOOST_AUTO_TEST_CASE(SimpleGraph)
{
    Dune::CpGrid grid;
    std::array<int,3> dims{2,2,2};
    std::array<double,3> size{2.,2.,2.};
    grid.createCartesian(dims,size);
    Opm::GraphOfGrid gog(grid);
    if (grid.size(0)==0) // in prallel run, non-root ranks are empty
        return;

    BOOST_REQUIRE(gog.size()==8); // number of graph vertices
    BOOST_REQUIRE(gog.numEdges(0)==3); // each vertex has 3 neighbors

    auto edgeL = gog.edgeList(2);
    BOOST_REQUIRE(edgeL.size()==3); // neighbors of vertex 2 are: 0, 3, 6
    BOOST_REQUIRE(edgeL[0]==1.);
    BOOST_REQUIRE(edgeL[3]==1.);
    BOOST_REQUIRE(edgeL[6]==1.);
    BOOST_REQUIRE_THROW(edgeL.at(4),std::out_of_range); // not a neighbor

    BOOST_REQUIRE_THROW(gog.edgeList(10),std::logic_error); // vertex 10 is not in the graph
}

// test vertex contraction on a simple graph
BOOST_AUTO_TEST_CASE(SimpleGraphWithVertexContraction)
{
    Dune::CpGrid grid;
    std::array<int,3> dims{2,2,2};
    std::array<double,3> size{2.,2.,2.};
    grid.createCartesian(dims,size);
    Opm::GraphOfGrid gog(grid);
    if (grid.size(0)==0)
        return;

    auto edgeL = gog.edgeList(3); // std::map<int,float>(gID,edgeWeight)
    BOOST_REQUIRE(edgeL[1]==1);
    BOOST_REQUIRE_THROW(edgeL.at(0),std::out_of_range);
    gog.contractVertices(0,1);
    BOOST_REQUIRE(gog.size()==7);
    edgeL = gog.edgeList(3);
    BOOST_REQUIRE_THROW(edgeL.at(1),std::out_of_range);
    BOOST_REQUIRE(edgeL[0]==1);
    edgeL = gog.edgeList(0);
    BOOST_REQUIRE(edgeL.size()==4);
    BOOST_REQUIRE(edgeL[2]==1); // neighbor of 0
    BOOST_REQUIRE(edgeL[3]==1); // neighbor of 1
    BOOST_REQUIRE_THROW(edgeL.at(1),std::out_of_range); // removed vertex, former neighbor of 0

    gog.contractVertices(0,2);
    BOOST_REQUIRE(gog.size()==6);
    BOOST_REQUIRE(gog.getVertex(0).weight==3.);
    edgeL = gog.edgeList(0);
    BOOST_REQUIRE(edgeL.size()==4);
    BOOST_REQUIRE(edgeL[3]==2);
    BOOST_REQUIRE(gog.edgeList(3).size()==2);
    BOOST_REQUIRE(gog.edgeList(3).at(0)==2);
    // contracting vertices removes higher ID from the graph
    // (when well is added, IDs removed from the graph are stored in the well)
    BOOST_REQUIRE_THROW(gog.getVertex(1),std::logic_error);

    auto v5e = gog.getVertex(5).edges;
    BOOST_REQUIRE(v5e==gog.edgeList(5));
    BOOST_REQUIRE(v5e==gog.edgeList(6)); // 5 and 6 have the same neighbors (1, 2 got merged)
    BOOST_REQUIRE(v5e!=gog.edgeList(7));

}

BOOST_AUTO_TEST_CASE(SimpleGraphWithTransmissibilities)
{
    Dune::CpGrid grid;
    std::array<int,3> dims{3,3,1};
    std::array<double,3> size{1.,1.,1.};
    grid.createCartesian(dims,size);
    // boundary faces should not appear in the graph, give them -1
    // other faces get value 10*ID1+ID2, where ID1<ID2 are cell global IDs
    std::vector<double> transmissiblities(24,-1);
    transmissiblities[grid.cellFace(0,1)] =  1;
    transmissiblities[grid.cellFace(1,1)] = 12;
    transmissiblities[grid.cellFace(3,1)] = 34;
    transmissiblities[grid.cellFace(4,1)] = 45;
    transmissiblities[grid.cellFace(6,1)] = 67;
    transmissiblities[grid.cellFace(7,1)] = 78;
    transmissiblities[grid.cellFace(0,3)] =  3;
    transmissiblities[grid.cellFace(1,3)] = 14;
    transmissiblities[grid.cellFace(2,3)] = 25;
    transmissiblities[grid.cellFace(3,3)] = 36;
    transmissiblities[grid.cellFace(4,3)] = 47;
    transmissiblities[grid.cellFace(5,3)] = 58;
    Opm::GraphOfGrid gog(grid,transmissiblities.data());
    if (grid.size(0)==0)
        return;

    int checked=0;
    double sum=0;
    BOOST_REQUIRE(gog.size()==9);
    for (int i=0; i<9; ++i)
    {
        const auto& edges = gog.edgeList(i);
        for (const auto& v : edges)
        {
            const double transm = i<v.first ? 10*i+v.first : i+10*v.first;
            BOOST_CHECK(transm==v.second);
            ++checked;
            sum += transm;
        }
    }
    BOOST_REQUIRE(checked==24); // each face gets checked twice
    BOOST_CHECK(sum==840); // 2*sum(transmissibilities) excluding boundaries

    gog.addWell(std::set<int>{0,1,3});
    gog.addWell(std::set<int>{2,6,7});
    BOOST_REQUIRE(gog.size()==5);
    {
        const auto& edges = gog.edgeList(0);
        checked=0;
        BOOST_REQUIRE(edges.size()==2);
        for (const auto& v : edges)
        {
            switch(v.first)
            {
                case 2: BOOST_CHECK(v.second==12+36); break;
                case 4: BOOST_CHECK(v.second==14+34); break;
                default: throw("Well 0 has a wrong edge.");
            }
        }
    }
    {
        const auto& edges = gog.edgeList(2);
        checked=0;
        BOOST_REQUIRE(edges.size()==4);
        for (const auto& v : edges)
        {
            switch(v.first)
            {
                case 0: BOOST_CHECK(v.second==12+36); break;
                case 4: BOOST_CHECK(v.second==47); break;
                case 5: BOOST_CHECK(v.second==25); break;
                case 8: BOOST_CHECK(v.second==78); break;
                default: throw("Well 2 has a wrong edge.");
            }
        }
    }
}

#if HAVE_MPI
BOOST_AUTO_TEST_CASE(WrapperForZoltan)
{
    Dune::CpGrid grid;
    std::array<int,3> dims{5,4,3};
    std::array<double,3> size{1.,1.,1.};
    grid.createCartesian(dims,size);
    Opm::GraphOfGrid gog(grid);
    if (grid.size(0)==0)
        return;

    int err;
    int nVer = getGraphOfGridNumVertices(&gog,&err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(nVer == 60);

    std::vector<uint> gIDs(nVer);
    std::vector<float> objWeights(nVer);
    getGraphOfGridVerticesList(&gog, 1, 1, gIDs.data(), nullptr, 1, objWeights.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(objWeights[18]==1); // all weights are 1 at this point

    std::vector<int> numEdges(nVer);
    getGraphOfGridNumEdges(&gog, 1, 1, nVer, gIDs.data(), nullptr, numEdges.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    int nEdges=0;
    for (int i=0; i<nVer; ++i)
    {
        switch (gIDs[i])
        {
            case 0:  BOOST_REQUIRE(numEdges[i]==3); break;
            case 9:  BOOST_REQUIRE(numEdges[i]==4); break;
            case 37: BOOST_REQUIRE(numEdges[i]==5); break;
            case 26: BOOST_REQUIRE(numEdges[i]==6); break;
        }
        nEdges += numEdges[i];
    }
    BOOST_REQUIRE(nEdges==266);

    std::vector<uint> nborGIDs(nEdges);
    std::vector<int> nborProc(nEdges);
    std::vector<float> edgeWeights(nEdges);
    getGraphOfGridEdgeList(&gog, 1, 1, nVer, gIDs.data(), nullptr, numEdges.data(), nborGIDs.data(), nborProc.data(), 1, edgeWeights.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE_MESSAGE(nborProc[145]==0, "Implementation detail: default process in GraphofGrid is 0");
    BOOST_REQUIRE(edgeWeights[203]==1.); // all are 1., no vertices were contracted

    numEdges[16] = 8;
    getGraphOfGridEdgeList(&gog, 1, 1, nVer, gIDs.data(), nullptr, numEdges.data(), nborGIDs.data(), nborProc.data(), 1, edgeWeights.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_FATAL);
}

BOOST_AUTO_TEST_CASE(GraphWithWell)
{
    Dune::CpGrid grid;
    std::array<int,3> dims{5,4,3};
    std::array<double,3> size{1.,1.,1.};
    grid.createCartesian(dims,size);
    Opm::GraphOfGrid gog(grid);
    if (grid.size(0)==0)
        return;

    std::unordered_map<std::string, std::set<int>> wells{
        {"shape L on the front face", {5,10,15,35,55} },
        {"lying 8 on the right face", {20,1,41,22,3,43,24} },
        {"disconnected vertices", {58,12} } };
    addFutureConnectionWells(gog,wells);
    BOOST_REQUIRE(gog.getWells().size()==3);
    int err;
    int nVer = getGraphOfGridNumVertices(&gog,&err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(nVer == 49);

    std::vector<uint> gIDs(nVer);
    std::vector<float> objWeights(nVer);
    getGraphOfGridVerticesList(&gog, 1, 1, gIDs.data(), nullptr, 1, objWeights.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    for (int i=0; i<nVer; ++i)
    {
        switch (gIDs[i])
        {
            case 1:  BOOST_REQUIRE(objWeights[i]==7.); break;
            case 5:  BOOST_REQUIRE(objWeights[i]==5.); break;
            case 12: BOOST_REQUIRE(objWeights[i]==2.); break;
            default: BOOST_REQUIRE(objWeights[i]==1.); // ordinary vertex
        }
    }
}

BOOST_AUTO_TEST_CASE(IntersectingWells)
{
    Dune::CpGrid grid;
    std::array<int,3> dims{5,4,3};
    std::array<double,3> size{1.,1.,1.};
    grid.createCartesian(dims,size);
    Opm::GraphOfGrid gog(grid);
    if (grid.size(0)==0)
        return;

    std::array<std::set<int>,3> wells{std::set<int>{0,1,2,3,4},
                                      std::set<int>{52,32,12},
                                      std::set<int>{59,48,37}};
                        // later add  std::set<int>{37,38,39,34},
                        //                         {2,8} and {2,38}
    for (const auto& w : wells)
    {
        gog.addWell(w,false);
    }
    BOOST_REQUIRE(gog.getWells().size()==3);

    int err;
    int nVer = getGraphOfGridNumVertices(&gog,&err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(nVer == 52);

    gog.addWell(std::set<int>{37,38,39,34}); // intersects with previous
    BOOST_REQUIRE(gog.getWells().size()==3);
    nVer = getGraphOfGridNumVertices(&gog,&err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(nVer == 49);

    gog.addWell(std::set<int>{2,8});
    nVer = getGraphOfGridNumVertices(&gog,&err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(nVer == 48);

    gog.addWell(std::set<int>{2,38}); // joins two wells
    BOOST_REQUIRE(gog.getWells().size()==2);
    nVer = getGraphOfGridNumVertices(&gog,&err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(nVer == 47);

    gog.addWell(std::set<int>{8,38});
    nVer = getGraphOfGridNumVertices(&gog,&err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(nVer == 47);

    std::vector<uint> gIDs(nVer);
    std::vector<float> objWeights(nVer);
    getGraphOfGridVerticesList(&gog, 1, 1, gIDs.data(), nullptr, 1, objWeights.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);

    for (int i=0; i<nVer; ++i)
    {
        switch (gIDs[i])
        {
            case 0:  BOOST_REQUIRE(objWeights[i]==12.); break;
            case 12: BOOST_REQUIRE(objWeights[i]==3.); break;
            default: BOOST_REQUIRE(objWeights[i]==1.); // ordinary vertex
        }
    }

    int nOut = 3;
    std::vector<int> numEdges(nOut);
    std::vector<uint> gID{12,0,54};
    getGraphOfGridNumEdges(&gog, 1, 1, nOut, gID.data(), nullptr, numEdges.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(numEdges[0]==12);
    BOOST_REQUIRE(numEdges[1]==26);
    BOOST_REQUIRE(numEdges[2]==3);

    int nEdges = 41;
    std::vector<uint> nborGIDs(nEdges);
    std::vector<int> nborProc(nEdges);
    std::vector<float> edgeWeights(nEdges);
    getGraphOfGridEdgeList(&gog, 1, 1, nOut, gID.data(), nullptr, numEdges.data(), nborGIDs.data(), nborProc.data(), 1, edgeWeights.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);

    // neighbors of the well with cells 12, 32, 52
    int checked = 0;
    for (int i=0; i<12; ++i)
    {
        BOOST_REQUIRE(edgeWeights[i]==1);
        switch (nborGIDs[i])
        {
            case  7: case 11: case 13: case 17:
            case 27: case 31: case 33: case  0: // 37 is a well with ID 0
            case 47: case 51: case 53: case 57:
                ++checked;
        }
    }
    BOOST_REQUIRE(checked==12);

    // neighbors of the well with cells 0,1,2,3,4,8,34,37,38,39,48,59
    checked=0;
    for (int i=12; i<38; ++i)
    {
        switch (nborGIDs[i])
        {
            // neighboring two well cells adds up the edge weight
            case 7: case 9: case 28: case 33: case 54: case 58:
                BOOST_REQUIRE(edgeWeights[i]==2.);
                ++checked;
                break;
            default: BOOST_REQUIRE(edgeWeights[i]==1.);
        }
    }
    BOOST_REQUIRE(checked==6);
    checked=0;

    // neighbors of the cell with global ID 54
    for (int i=38; i<41; ++i)
    {
        switch (nborGIDs[i])
        {
            case 0: // contains cells 34 and 59
                ++checked;
                BOOST_REQUIRE(edgeWeights[i]==2.);
                break;
            case 49: case 53:
                ++checked;
                BOOST_REQUIRE(edgeWeights[i]==1.);
        }
    }
    BOOST_REQUIRE(checked==3);

    const auto& wellList = gog.getWells();
    BOOST_REQUIRE(wellList.size()==2);
    std::set<int> well1{12,32,52};
    std::set<int> well2{0,1,2,3,4,8,34,37,38,39,48,59};
    if (wellList.begin()->size()==3)
    {
        BOOST_REQUIRE( *wellList.begin()==well1 );
        BOOST_REQUIRE( *wellList.rbegin()==well2 );
    }
    else
    {
        BOOST_REQUIRE( *wellList.begin()==well2 );
        BOOST_REQUIRE( *wellList.rbegin()==well1 );
    }
}
#endif // HAVE_MPI

BOOST_AUTO_TEST_CASE(WellWithBuffers)
{
    Dune::CpGrid grid;
    std::array<int,3> dims{5,1,1};
    std::array<double,3> size{1.,1.,1.};
    grid.createCartesian(dims,size);
    Opm::GraphOfGrid gog(grid);
    if (grid.size(0)==0)
        return;

    gog.addNeighboringCellsToWells(); // adding buffer to zero wells does nothing
    BOOST_REQUIRE(gog.size()==5);
    std::set<int> well{0,1};
    gog.addWell(well);

    // buffers of negative or zero size are ignored
    gog.addNeighboringCellsToWells(0);
    BOOST_REQUIRE(gog.size()==4);
    gog.addNeighboringCellsToWells(-4);
    BOOST_REQUIRE(gog.size()==4);

    gog.addNeighboringCellsToWells(); // no arg is 1 layer
    BOOST_REQUIRE(gog.size()==3);
    gog.addNeighboringCellsToWells(2);
    BOOST_REQUIRE(gog.size()==1);
}

BOOST_AUTO_TEST_CASE(NeighboringWellsWithBuffers)
{
    Dune::CpGrid grid;
    std::array<int,3> dims{6,1,1};
    std::array<double,3> size{1.,1.,1.};
    grid.createCartesian(dims,size);
    Opm::GraphOfGrid gog(grid);
    if (grid.size(0)==0)
        return;

    std::set<int> well0{0,1};
    std::set<int> well1{2,3};
    gog.addWell(well0);
    gog.addWell(well1);
    BOOST_REQUIRE(gog.size()==4);
    gog.addNeighboringCellsToWells();
    BOOST_REQUIRE(gog.size()==2);
    BOOST_REQUIRE(gog.getWells().size()==1);
    BOOST_REQUIRE(*gog.getWells().begin()==(std::set<int>{0,1,2,3,4}));
}

BOOST_AUTO_TEST_CASE(WellsWithIntersectingBuffers)
{
    Dune::CpGrid grid;
    std::array<int,3> dims{6,1,1};
    std::array<double,3> size{1.,1.,1.};
    grid.createCartesian(dims,size);
    Opm::GraphOfGrid gog(grid);
    if (grid.size(0)==0)
        return;

    std::set<int> well0{0,1};
    std::set<int> well1{3,4};
    gog.addWell(well0);
    gog.addWell(well1);
    BOOST_REQUIRE(gog.size()==4);
    gog.addNeighboringCellsToWells();
    BOOST_REQUIRE(gog.size()==1);
    BOOST_REQUIRE(gog.getWells().size()==1);
    BOOST_REQUIRE(*gog.getWells().begin()==(std::set<int>{0,1,2,3,4,5}));
}

BOOST_AUTO_TEST_CASE(WellsWithIntersectingBuffers2)
{
    Dune::CpGrid grid;
    std::array<int,3> dims{6,4,1};
    std::array<double,3> size{1.,1.,1.};
    grid.createCartesian(dims,size);
    Opm::GraphOfGrid gog(grid);
    if (grid.size(0)==0)
        return;

    std::set<int> well0{0,1};
    std::set<int> well1{5,11};
    std::set<int> well2{3,8,9,14};
    std::set<int> well3{18,22};
    gog.addWell(well0);
    gog.addWell(well1);
    gog.addWell(well2);
    gog.addWell(well3);
    BOOST_REQUIRE(gog.size()==18);
    gog.addNeighboringCellsToWells();
    BOOST_REQUIRE(gog.size()==2);
    const auto& wells = gog.getWells();
    BOOST_REQUIRE(wells.size()==2);
    for (const auto& w : wells)
    {
        if (*w.begin()==0)
        {
            BOOST_REQUIRE(w==(std::set<int>{0,1,2,3,4,5,6,7,8,9,10,11,13,14,15,17,20}));
        }
        else
        {
            BOOST_REQUIRE(w==(std::set<int>{12,16,18,19,21,22,23}));
        }
    }
    // adding one layer contracts everything into one vertex, another layer does nothing
    gog.addNeighboringCellsToWells(2);
    BOOST_REQUIRE(gog.size()==1);
}

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

#if HAVE_MPI
// Create yet another small grid with wells and test graph properties.
// This time wells are supplied via OpmWellType interface
BOOST_AUTO_TEST_CASE(addWellConnections)
{
    // create a grid
    Dune::CpGrid grid;
    std::array<int,3> dims{2,2,2};
    std::array<double,3> size{1.,1.,1.};
    grid.createCartesian(dims,size);
    Opm::GraphOfGrid gog(grid);
    if (grid.size(0)==0)
        return;
    BOOST_REQUIRE(gog.size()==8);

    auto wellCon = std::make_shared<Opm::WellConnections>(); // do not confuse with Dune::cpgrid::WellConnections
    wellCon->add(createConnection(0,0,0));
    wellCon->add(createConnection(0,1,0));
    wellCon->add(createConnection(0,1,1));
    std::vector<Dune::cpgrid::OpmWellType> wells;
    wells.push_back(createWell("first"));
    wells[0].updateConnections(wellCon,true);

    wellCon = std::make_shared<Opm::WellConnections>(); //reset
    wellCon->add(createConnection(0,0,1));
    wellCon->add(createConnection(1,1,0));
    wells.push_back(createWell("second"));
    wells[1].updateConnections(wellCon,true);

    wellCon = std::make_shared<Opm::WellConnections>(); //reset
    wellCon->add(createConnection(0,0,1));
    wellCon->add(createConnection(1,0,1));
    wells.push_back(createWell("third")); // intersects with second
    wells[2].updateConnections(wellCon,true);

    Dune::cpgrid::WellConnections wellConnections(wells,std::unordered_map<std::string, std::set<int>>(),gog.getGrid());
    BOOST_REQUIRE(wellConnections.size()==3);
    BOOST_REQUIRE(wellConnections[0]==(std::set<int>{0,2,6}));
    BOOST_REQUIRE(wellConnections[1].size()==2);
    BOOST_REQUIRE(wellConnections[1]==(std::set<int>{3,4}));
    BOOST_REQUIRE(wellConnections[2].size()==2);
    BOOST_REQUIRE(wellConnections[2]==(std::set<int>{4,5}));

    Opm::addWellConnections(gog,wellConnections,true);
    BOOST_REQUIRE(gog.size()==4);
    BOOST_REQUIRE(gog.getWells().size()==2); // second and third got merged (in gog)

    int err;
    int nVer = getGraphOfGridNumVertices(&gog,&err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(nVer == 4);
    std::vector<uint> gIDs(nVer);
    std::vector<float> objWeights(nVer);
    getGraphOfGridVerticesList(&gog, 1, 1, gIDs.data(), nullptr, 1, objWeights.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    std::sort(gIDs.begin(),gIDs.end());
    BOOST_REQUIRE(gIDs[0]==0 && gIDs[1]==1 && gIDs[2]==3 && gIDs[3]==7);
    std::vector<int> numEdges(nVer);
    getGraphOfGridNumEdges(&gog, 1, 1, nVer, gIDs.data(), nullptr, numEdges.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(numEdges[0]==3 && numEdges[1]==2 && numEdges[2]==3 && numEdges[3]==2);
    int nEdges = 10; // sum of numEdges[i]
    std::vector<uint> nborGIDs(nEdges);
    std::vector<int> nborProc(nEdges);
    std::vector<float> edgeWeights(nEdges);
    getGraphOfGridEdgeList(&gog, 1, 1, nVer, gIDs.data(), nullptr, numEdges.data(), nborGIDs.data(), nborProc.data(), 1, edgeWeights.data(), &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);

    // check all edgeWeights. Note that nborGIDs are not sorted
    for (int i=0; i<3; ++i)
    {
        switch (nborGIDs[i])
        {
            case 1: BOOST_REQUIRE(edgeWeights[i]==1); break;
            case 3: BOOST_REQUIRE(edgeWeights[i]==3); break;
            case 7: BOOST_REQUIRE(edgeWeights[i]==1); break;
            default: throw("GraphOfGrid was constructed badly.");
        }
    }
    for (int i=3; i<5; ++i)
    {
        switch (nborGIDs[i])
        {
            case 0: BOOST_REQUIRE(edgeWeights[i]==1); break;
            case 3: BOOST_REQUIRE(edgeWeights[i]==2); break;
            default: throw("GraphOfGrid was constructed badly.");
        }
    }
    for (int i=5; i<8; ++i)
    {
        switch (nborGIDs[i])
        {
            case 0: BOOST_REQUIRE(edgeWeights[i]==3); break;
            case 1: BOOST_REQUIRE(edgeWeights[i]==2); break;
            case 7: BOOST_REQUIRE(edgeWeights[i]==2); break;
            default: throw("GraphOfGrid was constructed badly.");
        }
    }
    for (int i=8; i<10; ++i)
    {
        switch (nborGIDs[i])
        {
            case 0: BOOST_REQUIRE(edgeWeights[i]==1); break;
            case 3: BOOST_REQUIRE(edgeWeights[i]==2); break;
            default: throw("GraphOfGrid was constructed badly.");
        }
    }

}
#endif // HAVE_MPI

BOOST_AUTO_TEST_CASE(gIDtoRankCorrection)
{
    // create a grid with wells
    Dune::CpGrid grid;
    std::array<int,3> dims{2,3,2};
    std::array<double,3> size{1.,1.,1.};
    grid.createCartesian(dims,size);
    Opm::GraphOfGrid gog(grid);
    if (grid.size(0)==0)
        return;

    // well needs at least 2 cells for vertex contraction
    gog.addWell(std::set<int>{});
    gog.addWell(std::set<int>{1});
    BOOST_REQUIRE(gog.getWells().size()==0);

    gog.addWell(std::set<int>{0,1,2});
    gog.addWell(std::set<int>{5,8,11});
    const auto& wells = gog.getWells();
    BOOST_REQUIRE(wells.size()==2);

    std::vector<int>gIDtoRank(12,1);
    gIDtoRank[0]=0; // well {0,1,2}
    gIDtoRank[8]=2; // inside well {5,8,11}, to be rewritten unless skipped
    extendGIDtoRank(gog,gIDtoRank,1); // skip wells on rank 1
    BOOST_CHECK(gIDtoRank[8]==2);
    for (int i=0; i<12; ++i)
    {
        if (i<3)
            BOOST_CHECK(gIDtoRank[i]==0);
        else if (i!=8)
            BOOST_CHECK(gIDtoRank[i]==1);
    }
    extendGIDtoRank(gog,gIDtoRank);
    BOOST_CHECK(gIDtoRank[8]==1);
}

// getWellRanks takes wellConnections and vector gIDtoRank mapping cells to their ranks
// and returns a vector of well ranks
BOOST_AUTO_TEST_CASE(test_getWellRanks)
{
    // create a grid with wells
    Dune::CpGrid grid;
    std::array<int,3> dims{1,2,4};
    std::array<double,3> size{1.,1.,1.};
    grid.createCartesian(dims,size);
    if (grid.size(0)==0)
        return;

    auto wellCon = std::make_shared<Opm::WellConnections>();
    wellCon->add(createConnection(0,0,0));
    wellCon->add(createConnection(0,1,0));
    wellCon->add(createConnection(0,1,1));
    std::vector<Dune::cpgrid::OpmWellType> wells;
    wells.push_back(createWell("first"));
    wells[0].updateConnections(wellCon,true);

    wellCon = std::make_shared<Opm::WellConnections>(); // reset
    wellCon->add(createConnection(0,0,2));
    wellCon->add(createConnection(0,1,2));
    wells.push_back(createWell("second"));
    wells[1].updateConnections(wellCon,true);

    wells.push_back(createWell("third"));

    std::vector<int> gIDtoRank{4,4,1,4,3,3,2,2};
    std::unordered_map<std::string, std::set<int>> futureConnections;
    futureConnections.emplace("third",std::set<int>{6,7});
    Dune::cpgrid::WellConnections wellConnections(wells,futureConnections,grid);
    auto wellRanks = Opm::getWellRanks(gIDtoRank,wellConnections);
    BOOST_REQUIRE(wellRanks.size()==3);
    BOOST_CHECK(wellRanks[0]==4);
    BOOST_CHECK(wellRanks[1]==3);
    BOOST_CHECK(wellRanks[2]==2);
}

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
