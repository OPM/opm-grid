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
	// vertex 1 and vertex 2 are neighbors
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
