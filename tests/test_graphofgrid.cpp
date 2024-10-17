#include <config.h>

#include <dune/common/version.hh>

#define BOOST_TEST_MODULE GraphRepresentationOfGrid
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <opm/grid/CpGrid.hpp>

#include <opm/grid/GraphOfGrid.hpp>

// basic test to check if the graph was constructed correctly
BOOST_AUTO_TEST_CASE(SimpleGraph)
{
	Dune::CpGrid grid;
	std::array<int,3> dims{2,2,2};
	std::array<double,3> size{2.,2.,2.};
	grid.createCartesian(dims,size);
	Opm::GraphOfGrid gog(grid);

	BOOST_REQUIRE(gog.size()==8); // number of graph vertices
	BOOST_REQUIRE(gog.numEdges(0)==3); // each vertex has 3 neighbors

	auto edgeL = gog.edgeList(2);
	BOOST_REQUIRE(edgeL.size()==3); // neighbors of vertex 2 are: 0, 3, 6
	BOOST_REQUIRE(edgeL[0]==1.);
	BOOST_REQUIRE(edgeL[3]==1.);
	BOOST_REQUIRE(edgeL[6]==1.);
	BOOST_REQUIRE(edgeL[4]==0.); // not a neighbor (edgeL's size increased)

	edgeL = gog.edgeList(10); // vertex 10 is not in the graph
	BOOST_REQUIRE(edgeL.size()==0);
}

// test vertex contraction on a simple graph
BOOST_AUTO_TEST_CASE(SimpleGraphWithVertexContraction)
{
	Dune::CpGrid grid;
	std::array<int,3> dims{2,2,2};
	std::array<double,3> size{2.,2.,2.};
	grid.createCartesian(dims,size);
	Opm::GraphOfGrid gog(grid);

	auto edgeL = gog.edgeList(2);
	gog.contractVertices(0,1);
	BOOST_REQUIRE(gog.size()==7);
	edgeL = gog.edgeList(0);
	BOOST_REQUIRE(edgeL.size()==4);
	BOOST_REQUIRE(edgeL[2]==1); // neighbor of 0
	BOOST_REQUIRE(edgeL[3]==1); // neighbor of 1

	gog.contractVertices(0,2);
	BOOST_REQUIRE(gog.size()==6);
	BOOST_REQUIRE(gog.getVertex(0).weight==3.);
	edgeL = gog.edgeList(0);
	BOOST_REQUIRE(edgeL.size()==4);
	BOOST_REQUIRE(edgeL[3]==2);
	BOOST_REQUIRE(gog.edgeList(3).size()==2);
	BOOST_REQUIRE(gog.edgeList(3)[0]==2);
	BOOST_REQUIRE_MESSAGE(gog.getVertex(1).weight==0.,"Implementation detail: \
		Looking up vertex not present in GraphofGrid gives dummy with weight 0.");

    auto v5e = gog.getVertex(5).edges;
	BOOST_REQUIRE(v5e==gog.edgeList(5));
	BOOST_REQUIRE(v5e==gog.edgeList(6)); // 5 and 6 have the same neighbors (1, 2 got merged)
	BOOST_REQUIRE(v5e!=gog.edgeList(7));

}

#include <opm/grid/GraphOfGridWrappers.hpp>

BOOST_AUTO_TEST_CASE(WrapperForZoltan)
{
	Dune::CpGrid grid;
	std::array<int,3> dims{5,4,3};
	std::array<double,3> size{1.,1.,1.};
	grid.createCartesian(dims,size);
	Opm::GraphOfGrid gog(grid);

    int err;
    int nVer = getGraphOfGridNumVertices(&gog,&err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(nVer == 60);

    int gIDs[nVer], lIDs[0];
    float objWeights[nVer];
    getGraphOfGridVerticesList(&gog, 1, 1, gIDs, lIDs, 1, objWeights, &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(objWeights[18]==1); // all weights are 1 at this point

    int numEdges[nVer];
    getGraphOfGridNumEdges(&gog, 1, 1, nVer, gIDs, lIDs, numEdges, &err);
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

    int nborGIDs[nEdges];
    int nborProc[nEdges];
    float edgeWeights[nEdges];
    getGraphOfGridEdgeList(&gog, 1, 1, nVer, gIDs, lIDs, numEdges, nborGIDs, nborProc, 1, edgeWeights, &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE_MESSAGE(nborProc[145]==0, "Implementation detail: default process in GraphofGrid is 0");
    BOOST_REQUIRE(edgeWeights[203]==1.); // all are 1., no vertices were contracted

    numEdges[16] = 8;
    std::string message("Expecting an error message from getGraphOfGridEdgeList, the vertex "
                        + gIDs[16] + std::string(" has a wrong number of edges."));
    Opm::OpmLog::info(message);
    getGraphOfGridEdgeList(&gog, 1, 1, nVer, gIDs, lIDs, numEdges, nborGIDs, nborProc, 1, edgeWeights, &err);
    BOOST_REQUIRE(err==ZOLTAN_FATAL);
}

BOOST_AUTO_TEST_CASE(GraphWithWell)
{
	Dune::CpGrid grid;
	std::array<int,3> dims{5,4,3};
	std::array<double,3> size{1.,1.,1.};
	grid.createCartesian(dims,size);
	Opm::GraphOfGrid gog(grid);

    std::unordered_map<std::string, std::set<int>> wells{
    	{"shape L on the front face", {5,10,15,35,55} },
    	{"lying 8 on the right face", {20,1,41,22,3,43,24} },
    	{"disconnected vertices", {58,12} } };
    addFutureConnectionWells(gog,wells);
    int err;
    int nVer = getGraphOfGridNumVertices(&gog,&err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(nVer == 49);

	int gIDs[nVer], lIDs[0];
    float objWeights[nVer];
    getGraphOfGridVerticesList(&gog, 1, 1, gIDs, lIDs, 1, objWeights, &err);
    BOOST_REQUIRE(err=ZOLTAN_OK);
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

	std::array<std::set<int>,3> wells{std::set<int>{0,1,2,3,4},
	                                  std::set<int>{52,32,12},
	                                  std::set<int>{59,48,37}};
			            // later add  std::set<int>{37,38,39,34},
			            //                         {2,8} and {2,38}
	for (const auto& w : wells)
		gog.addWell(w,false);

    int err;
    int nVer = getGraphOfGridNumVertices(&gog,&err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(nVer == 52);

	gog.addWell(std::set<int>{37,38,39,34}); // intersects with previous
    nVer = getGraphOfGridNumVertices(&gog,&err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(nVer == 49);

    gog.addWell(std::set<int>{2,8});
    nVer = getGraphOfGridNumVertices(&gog,&err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(nVer == 48);

    gog.addWell(std::set<int>{2,38});
    nVer = getGraphOfGridNumVertices(&gog,&err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(nVer == 47);

    gog.addWell(std::set<int>{8,38});
    nVer = getGraphOfGridNumVertices(&gog,&err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(nVer == 47);

	int gIDs[nVer], lIDs[0];
    float objWeights[nVer];
    getGraphOfGridVerticesList(&gog, 1, 1, gIDs, lIDs, 1, objWeights, &err);
    BOOST_REQUIRE(err=ZOLTAN_OK);

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
    int numEdges[nOut];
    int gID[nOut]; gID[0]=12; gID[1]=0; gID[2]=54;
    getGraphOfGridNumEdges(&gog, 1, 1, nOut, gID, lIDs, numEdges, &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(numEdges[0]==12);
    BOOST_REQUIRE(numEdges[1]==26);
    BOOST_REQUIRE(numEdges[2]==3);

    int nEdges = 41;
    int nborGIDs[nEdges];
    int nborProc[nEdges];
    float edgeWeights[nEdges];
    getGraphOfGridEdgeList(&gog, 1, 1, nOut, gID, lIDs, numEdges, nborGIDs, nborProc, 1, edgeWeights, &err);
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
