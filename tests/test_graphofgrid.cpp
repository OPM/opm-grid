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
	BOOST_REQUIRE(gog.edgeList(3)[0]==2);
	BOOST_REQUIRE_MESSAGE(gog.getVertex(1).weight==0.,"Implementation detail: \
		Looking up vertex not present in GraphofGrid gives dummy with weight 0.");

}

BOOST_AUTO_TEST_CASE(GrapWithWell)
{
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
