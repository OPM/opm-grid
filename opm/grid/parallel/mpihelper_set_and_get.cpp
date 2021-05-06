// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:


#include <dune/common/parallel/collectivecommunication.hh>

#include <dune/common/parallel/mpicollectivecommunication.hh>
#include "mpihelper_set_and_get.hpp"

namespace Opm
{
 
  // define our static member variables of class MPISetAndGetHelper
  MPISetAndGetHelper::MPICommunicator MPISetAndGetHelper::DUNE_MPI_COMM_ {MPI_COMM_WORLD};
  int MPISetAndGetHelper::rank_ = 0;
  int MPISetAndGetHelper::size_ = 0 ;



}
