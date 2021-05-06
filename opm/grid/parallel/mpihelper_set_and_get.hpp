// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 * File: opm/grid/parallel/mpihelper_set_and_get.hpp
 * 
 * This MPI helper class has been derived from the Dune library MPIHelper classes
 * defined in the file dune-common/opm/grid/parallel/mpihelper_set_and_get.hpp
 * 
 * See LICENCE.md in this directory for licence terms
 * 
 * 
 */


#ifndef DUNE_MPIHELPER
#define DUNE_MPIHELPER

#if HAVE_MPI
#include <cassert>
#endif

#if HAVE_MPI
#include <mpi.h>
#endif

#include <dune/common/parallel/collectivecommunication.hh>
#if HAVE_MPI
#include <dune/common/parallel/mpicollectivecommunication.hh>
#include <dune/common/stdstreams.hh>
#endif
#include <dune/common/visibility.hh>

// typedef MPI_Comm MPICommunicator;

namespace Opm
{
  
  
  /**
   * @brief A real mpi helper that can set as well as get the MPI communicator.
   *
   * This helper should be used for parallel programs.
   */
  class MPISetAndGetHelper
  {
  public:
    
    /**
     * @brief The type of the mpi communicator.
     */
     typedef MPI_Comm MPICommunicator;

    /** \brief get the default communicator
     *
     *  Return a communicator to exchange data with all processes
     *
     *  \returns MPI_COMM_WORLD
     */
    static MPICommunicator getCommunicator ()
    {
      return DUNE_MPI_COMM_;
    }
    
    /** \brief set the default communicator
     * 
     *  Return a communicator to exchange data with all processes.
     *  This function also sets @rank_ and size_ to values from the input
     *  communicator
     *
     *  \returns MPI_COMM_WORLD
     */
    static void setCommunicator (MPICommunicator DUNE_MPI_COMM)
    {
            
       DUNE_MPI_COMM_ = DUNE_MPI_COMM;
       MPI_Comm_rank(DUNE_MPI_COMM_,&rank_);
       MPI_Comm_size(DUNE_MPI_COMM_,&size_);
      
    }

    /** \brief get a local communicator
     *
     *  Returns a communicator to exchange data with the local process only
     *
     *  \returns MPI_COMM_SELF
     */
    static MPICommunicator getLocalCommunicator ()
    {
      return MPI_COMM_SELF;
    }

    static Dune::CollectiveCommunication<MPICommunicator>
    getCollectiveCommunication()
    {
      return Dune::CollectiveCommunication<MPICommunicator>(getCommunicator());
    }
    /**
     * @brief Get the singleton instance of the helper.
     *
     * This method has to be called with the same arguments
     * that the main method of the program was called:
     * \code
     * int main(int argc, char** argv){
     *   MPISetAndGetHelper::instance(argc, argv);
     *   // program code comes here
     *   ...
     * }
     * \endcode
     * @param argc The number of arguments provided to main.
     * @param argv The arguments provided to main.
     */
    DUNE_EXPORT static MPISetAndGetHelper& instance(int& argc, char**& argv)
    {
      // create singleton instance
      static MPISetAndGetHelper singleton (argc, argv);
      return singleton;
    }

    /**
     * @brief return rank of process
     */
    static int rank ()  { return rank_; }
    /**
     * @brief return number of processes
     */
    static int size ()  { return size_; }

  private:
    static int rank_;
    static int size_;
    static MPICommunicator DUNE_MPI_COMM_;
    bool initializedHere_;
    void prevent_warning(int){}

    //! \brief calls MPI_Init with argc and argv as parameters
    MPISetAndGetHelper(int& argc, char**& argv)
    : initializedHere_(false)
    {
      int wasInitialized = -1;
      MPI_Initialized( &wasInitialized );
      if(!wasInitialized)
      {
        rank_ = -1;
        size_ = -1;
        static int is_initialized = MPI_Init(&argc, &argv);
        prevent_warning(is_initialized);
        initializedHere_ = true;
        
        MPISetAndGetHelper::setCommunicator (MPI_COMM_WORLD) ; 
      }

      // MPISetAndGetHelper::setCommunicator (MPI_COMM_WORLD) ; 
      // MPI_Comm_rank(MPI_COMM_WORLD,&rank_);
      // MPI_Comm_size(MPI_COMM_WORLD,&size_);

      //assert( rank_ >= 0 );
      //assert( size_ >= 1 );

      std::cout << "Called  MPI_Init on p=" << rank_ << "!" << std::endl;
    }
    //! \brief calls MPI_Finalize
    ~MPISetAndGetHelper()
    {
      int wasFinalized = -1;
      MPI_Finalized( &wasFinalized );
      if(!wasFinalized && initializedHere_)
      {
        MPI_Finalize();
        std::cout << "Called MPI_Finalize on p=" << rank_ << "!" <<std::endl;
      }

    }
    MPISetAndGetHelper(const MPISetAndGetHelper&);
    MPISetAndGetHelper& operator=(const MPISetAndGetHelper);
  };
  


} // end namespace Dune
#endif
