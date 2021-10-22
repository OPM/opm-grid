#include <config.h>

#include <dune/common/version.hh>

#define BOOST_TEST_MODULE CPGridAquiferTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/CommunicationUtils.hpp>
#include <opm/parser/eclipse/Parser/ErrorGuard.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#if HAVE_MPI
class MPIError {
public:
  /** @brief Constructor. */
  MPIError(std::string s, int e) : errorstring(s), errorcode(e){}
  /** @brief The error string. */
  std::string errorstring;
  /** @brief The mpi error code. */
  int errorcode;
};

void MPI_err_handler(MPI_Comm *, int *err_code, ...){
  char *err_string=new char[MPI_MAX_ERROR_STRING];
  int err_length;
  MPI_Error_string(*err_code, err_string, &err_length);
  std::string s(err_string, err_length);
  std::cerr << "An MPI Error ocurred:"<<std::endl<<s<<std::endl;
  delete[] err_string;
  throw MPIError(s, *err_code);
}
#endif

struct CommandLineDataFileMpiInit
{
    void setup()
    {
        auto& argv = boost::unit_test::framework::master_test_suite().argv;
        auto& argc = boost::unit_test::framework::master_test_suite().argc;
        Dune::MPIHelper::instance(argc, argv);
        if(argc>1)
            caseFileName = std::string(argv[argc-1]);
    }
    void teardown()
    {
    }
    static const Dune::MPIHelper& helper(){
        auto& argv = boost::unit_test::framework::master_test_suite().argv;
        auto& argc = boost::unit_test::framework::master_test_suite().argc;
        return Dune::MPIHelper::instance(argc, argv);
    }
    static std::string caseFileName;
};

std::string CommandLineDataFileMpiInit::caseFileName;

BOOST_TEST_GLOBAL_FIXTURE( CommandLineDataFileMpiInit);


BOOST_AUTO_TEST_CASE(CpGridAquiferTest)
{    
    if(CommandLineDataFileMpiInit::helper().size()==1)
    {
         std::cerr<<"Test should be run parallel"<<std::endl;
        return;
    }

    Opm::Parser parser;
    Opm::ParseContext context;
    Opm::ErrorGuard guard;
    BOOST_REQUIRE(CommandLineDataFileMpiInit::caseFileName.size());
    const auto deck = parser.parseFile(CommandLineDataFileMpiInit::caseFileName, context,
                                       guard);
    Dune::CpGrid grid;
    Opm::EclipseGrid ecl_grid(deck);
    grid.processEclipseFormat(&ecl_grid, nullptr, false, false, false);

    auto aquifer_cells = grid.sortedNumAquiferCells();
    for( auto& cell: aquifer_cells)
    {
        cell = grid.globalCell()[cell];
    }
    grid.loadBalance();
    auto load_balanced_aquifer_cells = grid.sortedNumAquiferCells();

    for( auto&& cell: load_balanced_aquifer_cells)
    {
        cell = grid.globalCell()[cell];
    }

    auto [gathered_aquifer_cells, offset] = Opm::allGatherv(load_balanced_aquifer_cells, CommandLineDataFileMpiInit::helper().getCollectiveCommunication());

    offset.size();
    std::sort(aquifer_cells.begin(), aquifer_cells.end());
    std::sort(gathered_aquifer_cells.begin(), gathered_aquifer_cells.end());

    BOOST_TEST(aquifer_cells == gathered_aquifer_cells);
}
