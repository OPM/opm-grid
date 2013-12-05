#include <config.h>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE DistributedCpGridTests
#include <boost/test/unit_test.hpp>

#include <dune/grid/CpGrid.hpp>
#include <dune/geometry/referenceelements.hh>
#include <dune/common/fvector.hh>

class MPIError {
public:
  /** @brief Constructor. */
  MPIError(std::string s, int e) : errorstring(s), errorcode(e){}
  /** @brief The error string. */
  std::string errorstring;
  /** @brief The mpi error code. */
  int errorcode;
};

void MPI_err_handler(MPI_Comm *comm, int *err_code, ...){
  char *err_string=new char[MPI_MAX_ERROR_STRING];
  int err_length;
  MPI_Error_string(*err_code, err_string, &err_length);
  std::string s(err_string, err_length);
  std::cerr << "An MPI Error ocurred:"<<std::endl<<s<<std::endl;
  delete[] err_string;
  throw MPIError(s, *err_code);
}

class DummyDataHandle
{
public:
    typedef double DataType;
    bool fixedsize()
    {
        return true;
    }
    
    template<class T>
    std::size_t size(const T& t)
    {
        return 1;
    }
    template<class B, class T>
    void gather(B& buffer, const T& t)
    {
        buffer.write(100.0);
    }
    template<class B, class T>
    void scatter(B& buffer, const T& t, std::size_t s)
    {
        double val;
        std::cout<<"Scattering ";
        for(std::size_t i=0; i<s; ++i)
        {            
            buffer.read(val);
            std::cout<<val<<" to "<<t.index()<<" "<<i;
        }
        std::cout<<std::endl;
    }
    bool contains(int dim, int codim)
    {
        return dim==3 && (codim<=1 || codim==3);
    }
};

BOOST_AUTO_TEST_CASE(distribute)
{

    int m_argc = boost::unit_test::framework::master_test_suite().argc;
    char** m_argv = boost::unit_test::framework::master_test_suite().argv;
    Dune::MPIHelper::instance(m_argc, m_argv);
    MPI_Errhandler handler;
    MPI_Errhandler_create(MPI_err_handler, &handler);
    MPI_Errhandler_set(MPI_COMM_WORLD, handler);
    int procs;
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
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

    if(procs==1)
    {
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
    }
    DummyDataHandle data;

    grid.loadBalance(data);

    
    if(procs==1)
    {
        // Check whether the scattered grid is identical to the orinal one.
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
    }else
    {
        grid.communicate(data, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
    }
}
