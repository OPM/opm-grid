#include <config.h>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE DistributedCpGridTests
#include <boost/test/unit_test.hpp>

#include <dune/grid/CpGrid.hpp>
#include <dune/geometry/referenceelements.hh>
#include <dune/common/fvector.hh>

BOOST_AUTO_TEST_CASE(distribute)
{

    int m_argc = boost::unit_test::framework::master_test_suite().argc;
    char** m_argv = boost::unit_test::framework::master_test_suite().argv;
    Dune::MPIHelper::instance(m_argc, m_argv);
    Dune::CpGrid grid;
    std::array<int, 3> dims={{10, 10, 10}};
    std::array<double, 3> size={{ 1.0, 1.0, 1.0}};

    grid.createCartesian(dims, size);
    std::vector<int> cell_indices, face_indices, point_indices;
    std::vector<Dune::CpGrid::Traits::Codim<0>::Geometry::GlobalCoordinate > cell_centers, face_centers, point_centers;
    int cell_size = grid.leafView().size(0);
    int face_size = grid.leafView().size(1);
    int point_size = grid.leafView().size(3);
    
    const Dune::CpGrid::LeafIndexSet& ix = grid.leafIndexSet();
    for (Dune::CpGrid::Codim<0>::LeafIterator it = grid.leafbegin<0>();
             it != grid.leafend<0>(); ++it) {
        Dune::GeometryType gt = it->type () ;
        const Dune::GenericReferenceElement<Dune::CpGrid::ctype, 3>& ref=
            Dune::GenericReferenceElements<Dune::CpGrid::ctype, 3>::general(gt);

        cell_indices.push_back(ix.index(*it));
        cell_centers.push_back(it->geometry().center());
        for(Dune::CpGrid::LeafIntersectionIterator iit=grid.leafView().ibegin(*it), 
                endiit = grid.leafView().iend(*it); iit!=endiit; ++iit)
        {
            //            face_indices.push_back(ix.index(*it->subEntity<1>(iit->indexInInside())));
            face_centers.push_back(iit->geometry().center());
            for(int i=0; i<4; ++i){
                point_indices.push_back(ix.subIndex(*it, ref.subEntity(iit->indexInInside(),1,i,3), 3));
                //ref.subEntity(iit->indexInInside(),1,i,dim).geometry().center();
            }
        }
    }
    
    grid.scatterGrid();

    BOOST_REQUIRE(cell_size  == grid.leafView().size(0));
    BOOST_REQUIRE(face_size  == grid.leafView().size(1));
    BOOST_REQUIRE(point_size == grid.leafView().size(3));
    
    int cell_index=0, face_index=0, point_index=0;

    const Dune::CpGrid::LeafIndexSet& ix1 = grid.leafIndexSet();
    BOOST_REQUIRE(&ix!=&ix1);
    
    for (Dune::CpGrid::Codim<0>::LeafIterator it = grid.leafbegin<0>();
             it != grid.leafend<0>(); ++it) {
        Dune::GeometryType gt = it->type () ;
        const Dune::GenericReferenceElement<Dune::CpGrid::ctype, 3>& ref=
            Dune::GenericReferenceElements<Dune::CpGrid::ctype, 3>::general(gt);

        BOOST_REQUIRE(cell_indices[cell_index]==ix1.index(*it));
        BOOST_REQUIRE(cell_centers[cell_index++]==it->geometry().center());
        for(Dune::CpGrid::LeafIntersectionIterator iit=grid.leafView().ibegin(*it), 
                endiit = grid.leafView().iend(*it); iit!=endiit; ++iit)
        {
            //BOOST_REQUIRE(face_indices[face_index]==ix1.index(*it->subEntity<1>(iit->indexInInside())));
            BOOST_REQUIRE(face_centers[face_index++]==iit->geometry().center());
            for(int i=0; i<4; ++i){
                BOOST_REQUIRE(point_indices[point_index++]==ix1.subIndex(*it, ref.subEntity(iit->indexInInside(),1,i,3), 3));
                //ref.subEntity(iit->indexInInside(),1,i,dim).geometry().center();
            }
        }
    }
    
}