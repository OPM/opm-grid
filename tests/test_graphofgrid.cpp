#include <config.h>

#include <dune/common/version.hh>

#define BOOST_TEST_MODULE GraphRepresentationOfGrid
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <opm/grid/CpGrid.hpp>

#include <opm/grid/GraphOfGrid.cpp>

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
	BOOST_REQUIRE(edgeL[4]==0.); // not a neighbor (beware! size got increased)
}

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
	BOOST_REQUIRE(edgeL[2]==1);
	BOOST_REQUIRE(edgeL[3]==1);

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

	auto it = gog.begin();
	BOOST_REQUIRE(it->first==0); // gID are ordered, 0 is first
	++it;
	BOOST_REQUIRE(it->first==3); // 1, 2 got merged, 3 is next
}

#include <opm/grid/GraphOfGridWrappers.cpp>

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
    BOOST_REQUIRE(gIDs[6]==6);
    BOOST_REQUIRE(gIDs[59]==59);
    BOOST_REQUIRE(objWeights[18]==1);

    int numEdges[nVer];
    getGraphOfGridNumEdges(&gog, 1, 1, nVer, gIDs, lIDs, numEdges, &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(numEdges[0]==3);
    BOOST_REQUIRE(numEdges[9]==4);
    BOOST_REQUIRE(numEdges[37]==5);
    BOOST_REQUIRE(numEdges[26]==6);

    int nEdges=0;
    for (int i=0; i<nVer; ++i)
    	nEdges += numEdges[i];
    BOOST_REQUIRE(nEdges==266);

    int nborGIDs[nEdges];
    int nborProc[nEdges];
    float edgeWeights[nEdges];
    getGraphOfGridEdgeList(&gog, 1, 1, nVer, gIDs, lIDs, numEdges, nborGIDs, nborProc, 1, edgeWeights, &err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE_MESSAGE(nborProc[145]==0, "Implementation detail: default process in GraphofGrid is 0");
    BOOST_REQUIRE(edgeWeights[203]==1.); // all are 1., no vertices were contracted

    numEdges[16] = 8;
    std::cerr << "Expecting an error message from getGraphOfGridEdgeList, the vertex 16 has a wrong number of edges." <<std::endl;
    getGraphOfGridEdgeList(&gog, 1, 1, nVer, gIDs, lIDs, numEdges, nborGIDs, nborProc, 1, edgeWeights, &err);
    BOOST_REQUIRE(err==ZOLTAN_FATAL);
}

BOOST_AUTO_TEST_CASE(GrapWithWell)
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
    addFutureConnectionWell(gog,wells);
    int err;
    int nVer = getGraphOfGridNumVertices(&gog,&err);
    BOOST_REQUIRE(err==ZOLTAN_OK);
    BOOST_REQUIRE(nVer == 49);

	int gIDs[nVer], lIDs[0];
    float objWeights[nVer];
    getGraphOfGridVerticesList(&gog, 1, 1, gIDs, lIDs, 1, objWeights, &err);
    BOOST_REQUIRE(err=ZOLTAN_OK);
    BOOST_REQUIRE(gIDs[48]==59); // last ID
    BOOST_REQUIRE(gIDs[47]==57); // pre-last, 58 was contracted to 12
    BOOST_REQUIRE(gIDs[4]==5); // 3 was contracted to 1
    BOOST_REQUIRE(objWeights[1]==7.); // well with 7 vertices
    BOOST_REQUIRE(objWeights[4]==5.); // well with 5 vertices
    BOOST_REQUIRE(gIDs[19]==25);
    BOOST_REQUIRE(objWeights[19]==1.); // ordinary vertex

	// Well well1("SomeName1",std::array<int>{0,1,2,3,4});
	// Well well2("SomeName2",std::array<int>{52,32,12});
	// Well well3("SomeName3",std::array<int>{59,48,37});
	// Well well4("SomeName4",std::array<int>{37,38,39,34}); // merges with well3
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
