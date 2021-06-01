#include <config.h>

#define BOOST_TEST_MODULE CartGridTest
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif

#include <opm/grid/common/SubGridView.hpp>

// Warning suppression for Dune includes.
#include <opm/grid/utility/platform_dependent/disable_warnings.h>

#include <dune/common/unused.hh>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/GridHelpers.hpp>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/yaspgrid.hh>

#include <opm/grid/cpgrid/dgfparser.hh>


#define DISABLE_DEPRECATED_METHOD_CHECK 1
using Dune::referenceElement; //grid check assume usage of Dune::Geometry
#include <dune/grid/test/gridcheck.hh>


// Re-enable warnings.
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>

#include <cmath>
#include <iostream>

#if HAVE_MPI
struct MPIError
{
    MPIError(std::string s, int e) : errorstring(std::move(s)), errorcode(e){}
    std::string errorstring;
    int errorcode;
};

void MPI_err_handler(MPI_Comm*, int* err_code, ...)
{
    std::vector<char> err_string(MPI_MAX_ERROR_STRING);
    int err_length;
    MPI_Error_string(*err_code, err_string.data(), &err_length);
    std::string s(err_string.data(), err_length);
    std::cerr << "An MPI Error ocurred:" << std::endl << s << std::endl;
    throw MPIError(s, *err_code);
}
#endif

template <class GridView>
void testGridInteriorIteration( const GridView& gridView, const int nElem )
{
    typedef typename GridView::template Codim<0>::Iterator ElemIterator;
    typedef typename GridView::IntersectionIterator IsIt;
    typedef typename GridView::template Codim<0>::Geometry Geometry;

    int numElem = 0;
    ElemIterator elemIt = gridView.template begin<0>();
    ElemIterator elemEndIt = gridView.template end<0>();
    for (; elemIt != elemEndIt; ++elemIt) {
        if (elemIt->partitionType() != Dune::InteriorEntity) {
            continue;
        }
        const Geometry& elemGeom = elemIt->geometry();
        BOOST_CHECK_CLOSE(elemGeom.volume(), 1.0, 1e-8);

        typename Geometry::LocalCoordinate local( 0.5 );
        typename Geometry::GlobalCoordinate global = elemGeom.global( local );
        typename Geometry::GlobalCoordinate center = elemGeom.center();
        BOOST_CHECK_SMALL((center - global).two_norm(), 1e-12);

        int numIs = 0;
        IsIt isIt = gridView.ibegin(*elemIt);
        IsIt isEndIt = gridView.iend(*elemIt);
        for (; isIt != isEndIt; ++isIt, ++ numIs)
        {
            const auto& intersection = *isIt;
            const auto& isGeom = intersection.geometry();
            BOOST_CHECK_CLOSE(isGeom.volume(), 1.0, 1e-8);

            if (intersection.neighbor())
            {
                BOOST_CHECK_EQUAL(numIs, intersection.indexInInside());
                BOOST_CHECK_CLOSE(intersection.outside().geometry().volume(), 1.0, 1e-8);
                BOOST_CHECK_CLOSE(intersection.inside().geometry().volume(), 1.0, 1e-8);
            }
        }

        BOOST_CHECK_EQUAL(numIs, 2 * GridView::dimension);

        ++ numElem;
    }

    BOOST_CHECK_EQUAL(numElem, nElem);
}


template <class Grid>
auto getSeeds(const Grid& grid, const std::vector<int>& indices)
{
    assert(std::is_sorted(indices.begin(), indices.end()));
    using EntitySeed = typename Grid::template Codim<0>::Entity::EntitySeed;
    std::vector<EntitySeed> seeds(indices.size());
    auto it = grid.template leafbegin<0>();
    int previous = 0;
    for (std::size_t c = 0; c < indices.size(); ++c) {
        std::advance(it, indices[c] - previous);
        seeds[c] = it->seed();
        previous = indices[c];
    }
    return seeds;
}


template <class Grid>
void testGrid(Grid& grid, const std::string& name, const std::size_t nElem, const std::size_t nVertices)
{
    typedef typename Grid::LeafGridView GridView;

    std::cout << name << std::endl;

    testGridInteriorIteration( grid.leafGridView(), nElem );

    std::cout << "create vertex mapper\n";
    Dune::MultipleCodimMultipleGeomTypeMapper<GridView> mapper(grid.leafGridView(), Dune::mcmgVertexLayout());

    BOOST_CHECK_EQUAL(mapper.size(), nVertices);

    Dune::SubGridView<Grid> sgv(grid, getSeeds(grid, {0, 1, 2}));
    testGridInteriorIteration(sgv, 3);

}


#if HAVE_ECL_INPUT
BOOST_AUTO_TEST_CASE(deck)
{
    // ------------ Test grid from deck. ------------
    const char* deckString =
R"(
RUNSPEC
METRIC
DIMENS
2 2 2 /
GRID
DXV
2*1 /
DYV
2*1 /
DZ
8*1 /
TOPS
8*100.0 /
)";

    const auto deck = Opm::Parser{}.parseString(deckString);

    Dune::CpGrid grid;
    const int* actnum = deck.hasKeyword("ACTNUM") ? deck.getKeyword("ACTNUM").getIntData().data() : nullptr;
    Opm::EclipseGrid ecl_grid(deck , actnum);

    grid.processEclipseFormat(&ecl_grid, nullptr, false, false, false);
    testGrid( grid, "CpGrid_ecl", 8, 27 );
}
#endif // HAVE_ECL_INPUT


BOOST_AUTO_TEST_CASE(dgf)
{
    // ------------ Test grid from dgf. ------------
    std::stringstream dgfFile;
    // create grid with 4 cells in each direction
    dgfFile << "DGF" << std::endl;
    dgfFile << "Interval" << std::endl;
    dgfFile << "0 0 0" << std::endl;
    dgfFile << "4 4 4" << std::endl;
    dgfFile << "4 4 4" << std::endl;
    dgfFile << "#" << std::endl;

    Dune::GridPtr< Dune::CpGrid > gridPtr( dgfFile );
    testGrid( *gridPtr, "CpGrid_dgf", 64, 125 );
}


BOOST_AUTO_TEST_CASE(yasp)
{
    // ------------ Test YaspGrid. ------------

    Dune::YaspGrid<3, Dune::EquidistantCoordinates<double, 3>> yaspGrid({4.0, 4.0, 4.0}, {4, 4, 4});
    testGrid(yaspGrid, "YaspGrid", 64, 125);
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
#if HAVE_MPI
    // register a throwing error handler to allow for
    // debugging with "catch throw" in gdb
    MPI_Errhandler handler;
    MPI_Comm_create_errhandler(MPI_err_handler, &handler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, handler);
#endif
    return boost::unit_test::unit_test_main(&init_unit_test, argc, argv);
}
