/*
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services.
  Copyright 2015 Statoil AS

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
#ifndef DUNE_CPGRID_ZOLTANPARTITION_HEADER
#define DUNE_CPGRID_ZOLTANPARTITION_HEADER

#include <unordered_set>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/ZoltanGraphFunctions.hpp>
#include <opm/grid/common/WellConnections.hpp>

#if HAVE_MPI
namespace Dune
{
namespace cpgrid
{
/// \brief Creates all the information needed from the partitioning
///
/// Note that if wells is not null and allowDitributedWells is false,
/// then the partitioning is postprocessed such that all cells that a
/// perforated by a well are part of the interior of exactly one process.
/// Otherwise they might spreaded between multiple processes.
///
/// \param grid The grid
/// \param cc  The MPI communicator used for the partitioning.
/// \param wells Pointer to vector with all possible wells (all perforations) of the problem.
///              nullptr is possible
/// \param gridAndWells Graph representing grid and wells (must match the other parameters!)
/// \param root The rank that has the whole grid before loadbalancing.
/// \param numExport The number of exported cells (created e.g. by Zoltan)
/// \param numImport The number of imported cells (created e.g. by Zoltan)
/// \param exportGlobalIds C-array of the global ids of exported cells (created e.g. by Zoltan)
/// \param exportLocalIds C-array of the local ids of exported cells (created e.g. by Zoltan)
/// \param exportToPart C-array with target partition of exported cells (created e.g. by Zoltan)
/// \param importGlobalIds C-array of the global ids of exported cells (created e.g. by Zoltan)
/// \param allowDistributedWells Whether wells might be distibuted to the interior of multiple processes.
/// \return  A tuple consisting of a vector that contains for each local cell of the original grid the
///         the number of the process that owns it after repartitioning,
///         a vector containing a pair of name  and a boolean indicating whether this well has
///         perforated cells local to the process of all wells,
///         vector containing information for each exported cell (global id
///         of cell, process id to send to, attribute there), and a vector containing
///         information for each imported cell (global index, process id that sends, attribute here, local index
///         here), and a WellConnections object containing information about the well connections
///         (if argument wells was not null and this is the root rank this will contain connections in
///          form of global indices)
template<class Id>
std::tuple<std::vector<int>, std::vector<std::pair<std::string,bool>>,
           std::vector<std::tuple<int,int,char> >,
           std::vector<std::tuple<int,int,char,int> >,
           WellConnections>
makeImportAndExportLists(const Dune::CpGrid& cpgrid,
                         const Dune::CollectiveCommunication<MPI_Comm>& cc,
                         const std::vector<Dune::cpgrid::OpmWellType> * wells,
                         const Dune::cpgrid::CombinedGridWellGraph* gridAndWells,
                         int root,
                         int numExport,
                         int numImport,
                         const Id* exportLocalGids,
                         const Id* exportGlobalGids,
                         const int* exportToPart,
                         const Id* importGlobalGids,
                         bool allowDistributedWells = false);

template<class Id>
std::tuple<int, std::vector<Id> >
scatterExportInformation(int numExport, const Id* exportGlobalGids,
                         const int* exportToPart, int root,
                         const Dune::CollectiveCommunication<MPI_Comm>& cc);
} // end namespace cpgrid
} // end namespace Dune
#endif //HAVE_MPI
#if defined(HAVE_ZOLTAN) && defined(HAVE_MPI)
namespace Dune
{
namespace cpgrid
{
/// \brief Partition a CpGrid using Zoltan
///
/// This function will extract Zoltan's graph information
/// from the grid, and the wells and use it to partition the grid.
/// In case the global grid is available on all processes, it
/// will nevertheless only use the information on the root process
/// to partition it as Zoltan cannot identify this situation.
/// @param grid The grid to partition
/// @param wells The wells of the eclipse If null wells will be neglected.
/// @param transmissibilities The transmissibilities associated with the
///             faces
/// @paramm cc  The MPI communicator to use for the partitioning.
///             The will be partitioned among the partiticipating processes.
/// @param edgeWeightMethod The method used to calculate the weights associated
///             with the edges of the graph (uniform, transmissibilities, log thereof)
/// @param root The process number that holds the global grid.
/// @param zoltanImbalanceTol Set the imbalance tolerance used by Zoltan
/// \param allowDistributedWells Allow the perforation of a well to be distributed to the
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

std::tuple<std::vector<int>,std::vector<std::pair<std::string,bool>>,
           std::vector<std::tuple<int,int,char> >,
           std::vector<std::tuple<int,int,char,int> >,
           WellConnections>
zoltanGraphPartitionGridOnRoot(const CpGrid& grid,
                               const std::vector<OpmWellType> * wells,
                               const double* transmissibilities,
                               const CollectiveCommunication<MPI_Comm>& cc,
                               EdgeWeightMethod edgeWeightsMethod, int root,
                               const double zoltanImbalanceTol,
                               bool allowDistributedWells);

/// \brief Partition a CpGrid using Zoltan serially only on rank 0
///
/// This function will extract Zoltan's graph information
/// from the grid, and the wells and use it to partition the grid.
/// In case the global grid is available on all processes, it
/// will nevertheless only use the information on the root process
/// to partition it as Zoltan cannot identify this situation.
/// @param grid The grid to partition
/// @param wells The wells of the eclipse If null wells will be neglected.
/// @param transmissibilities The transmissibilities associated with the
///             faces
/// @paramm cc  The MPI communicator to use for the partitioning.
///             The will be partitioned among the partiticipating processes.
/// @param edgeWeightMethod The method used to calculate the weights associated
///             with the edges of the graph (uniform, transmissibilities, log thereof)
/// @param root The process number that holds the global grid.
/// @param zoltanImbalanceTol Set the imbalance tolerance used by Zoltan
/// \param allowDistributedWells Allow the perforation of a well to be distributed to the
///        interior region of multiple processes.
/// @return A tuple consisting of a vector that contains for each local cell of the original grid the
///         the number of the process that owns it after repartitioning,
///         a set of names of wells that should be defunct in a parallel
///         simulation, vector containing information for each exported cell (global id
///         of cell, process id to send to, attribute there), a vector containing
///         information for each imported cell (global index, process id that sends, attribute here, local index
///         here), and a WellConnections object containing information about the well connections
///         (if argument wells was not null and this is the root rank this will contain connections in
///          form of global indices)
///
/// @note This function will only do *serial* partioning.
std::tuple<std::vector<int>, std::vector<std::pair<std::string,bool>>,
           std::vector<std::tuple<int,int,char> >,
           std::vector<std::tuple<int,int,char,int> >,
           WellConnections>
zoltanSerialGraphPartitionGridOnRoot(const CpGrid& grid,
                               const std::vector<OpmWellType> * wells,
                               const double* transmissibilities,
                               const CollectiveCommunication<MPI_Comm>& cc,
                               EdgeWeightMethod edgeWeightsMethod, int root,
                               const double zoltanImbalanceTol,
                               bool allowDistributedWells);
}
}
#endif // HAVE_ZOLTAN
#endif // header guard
