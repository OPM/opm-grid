/*
  Copyright 2024 OPM-OP AS

  This file is part of The Open Porous Media project  (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef DUNE_CPGRID_METISPARTITION_HEADER
#define DUNE_CPGRID_METISPARTITION_HEADER

#include <unordered_set>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/ZoltanGraphFunctions.hpp>
#include <opm/grid/common/WellConnections.hpp>
#include <opm/grid/common/ZoltanPartition.hpp>
#include <opm/grid/common/GridPartitioning.hpp>

// We want to use METIS, but if METIS is installed together with Scotch, then METIS uses some artifacts from Scotch.
// For this type, we need to include scotch.h.
#if IS_SCOTCH_METIS_HEADER
extern "C" {
  #include <scotch.h>
}
// And also, we need to set the version number here manually to 5
#ifndef SCOTCH_METIS_VERSION
#define SCOTCH_METIS_VERSION 5
#endif /* SCOTCH_METIS_VERSION */
#endif

#if defined(HAVE_METIS) && HAVE_MPI
extern "C" {
  #include <metis.h>
}
#endif // defined(HAVE_METIS) && HAVE_MPI

#if defined(HAVE_METIS) && HAVE_MPI
namespace Dune
{
namespace cpgrid
{
#if defined(REALTYPEWIDTH)
  using real_t = ::real_t;
#else
  using real_t = double;
#endif

#if defined(IDXTYPEWIDTH)
  using idx_t = ::idx_t;
#elif IS_SCOTCH_METIS_HEADER
  using idx_t = SCOTCH_Num;
#else
  using idx_t = int;
#endif

#if IS_SCOTCH_METIS_HEADER
  // NOTE: scotchmetis does not define a return type for METIS functions
  #define METIS_OK 1
#endif

/// \brief Partition a CpGrid using METIS
///
/// This function will extract graph information
/// from the grid, and the wells and use it to partition the grid.

/// @param grid The grid to partition
/// @param wells The wells of the eclipse If null wells will be neglected.
/// @param possibleFutureConnections Possible future connections of wells that might get added through an ACTIONX.
///                                  The grid will then be partitioned such that these connections are on the same
///                                  partition. If NULL, they will be neglected.
/// @param transmissibilities The transmissibilities associated with the
///             faces
/// @paramm cc  The MPI communicator to use for the partitioning.
///             The will be partitioned among the partiticipating processes.
/// @param edgeWeightMethod The method used to calculate the weights associated
///             with the edges of the graph (uniform, transmissibilities, log thereof)
/// @param root The process number that holds the global grid.
/// @param imbalanceTol Set the imbalance tolerance used by METIS, i.e. the entries of the parameter ubvec
/// @param allowDistributedWells Allow the perforation of a well to be distributed to the
///        interior region of multiple processes.
/// @return A tuple consisting of a vector that contains for each local cell of the original grid the
///         the number of the process that owns it after repartitioning,
///         a vector containing a pair of name  and a boolean indicating whether this well has
///         perforated cells local to the process of all wells,
///         vector containing information for each exported cell (global id
///         of cell, process id to send to, attribute there), a vector containing
///         information for each imported cell (global index, process id that sends, attribute here, local index
///         here), and a WellConnections object containing information about the well connections
///         (if argument wells was not null and this is the root rank this will contain connections in
///          form of global indices)

std::tuple<std::vector<int>,
           std::vector<std::pair<std::string, bool>>,
           std::vector<std::tuple<int, int, char>>,
           std::vector<std::tuple<int, int, char, int>>,
           WellConnections>
metisSerialGraphPartitionGridOnRoot(const CpGrid& grid,
                                    const std::vector<OpmWellType> * wells,
                                    const std::unordered_map<std::string, std::set<std::array<int,3>>>* possibleFutureConnections,
                                    const double* transmissibilities,
                                    const Communication<MPI_Comm>& cc,
                                    EdgeWeightMethod edgeWeightsMethod,
                                    int root,
                                    real_t imbalanceTol,
                                    bool allowDistributedWells,
                                    const std::map<std::string,std::string>& params);
}
}


#endif // defined(HAVE_METIS) && HAVE_MPI
#endif // header guard
