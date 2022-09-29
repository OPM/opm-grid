/*
  Copyright 2016 Dr. Blatt - HPC-Simulation-Software & Services.
  Copyright 2016 Statoil AS

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
#ifndef DUNE_CPGRID_WELL_CONNECTIONS_HEADER_INCLUDED
#define DUNE_CPGRID_WELL_CONNECTIONS_HEADER_INCLUDED

#include <set>
#include <unordered_set>
#include <vector>
#include <functional>

#ifdef HAVE_MPI
#include <opm/grid/utility/platform_dependent/disable_warnings.h>
#include "mpi.h"
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>
#endif

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 7)
#include <dune/common/parallel/communication.hh>
#else
#include <dune/common/parallel/mpicollectivecommunication.hh>
#endif


#include <opm/grid/utility/OpmParserIncludes.hpp>
#include <opm/grid/CpGrid.hpp>

namespace Dune
{
namespace cpgrid
{

/// \brief A class calculating and representing all connections of wells.
///
/// Wells are identified by their position as exported by the wells method
/// of the eclipse parser.
/// For each well the container stores at the well index all indices of cells
/// that the well perforates.
class WellConnections
{

public:
    /// \brief The const iterator type.
    typedef std::vector<std::set<int> >::const_iterator const_iterator;

    /// \brief The iterator type (always const).
    typedef const_iterator iterator;

    WellConnections() = default;

    /// \brief Constructor
    /// \param wells The eclipse information about the wells
    /// \param cartesianSize The logical cartesian size of the grid.
    /// \param cartesian_to_compressed Mapping of cartesian index
    ///        compressed cell index. The compressed index is used
    ///        to represent the well conditions.
    WellConnections(const std::vector<OpmWellType>& wells,
                    const std::array<int, 3>& cartesianSize,
                    const std::vector<int>& cartesian_to_compressed);

    /// \brief Constructor
    /// \param wells The eclipse information about the wells
    /// \param cpGrid The corner point grid
    WellConnections(const std::vector<OpmWellType>& wells,
                    const Dune::CpGrid& cpGrid);

    /// \brief Initialze the data of the container
    /// \param schedule The eclipse information
    /// \param cartesianSize The logical cartesian size of the grid.
    /// \param cartesian_to_compressed Mapping of cartesian index
    ///        compressed cell index. The compressed index is used
    ///        to represent the well conditions.
    void init(const std::vector<OpmWellType>& wells,
              const std::array<int, 3>& cartesianSize,
              const std::vector<int>& cartesian_to_compressed);

    /// \brief Access all connections of a well
    /// \param i The index of the well (position of the well in the
    ///          eclipse schedule.
    /// \return The set of compressed indices of cells perforated by the well.
    const std::set<int>& operator[](std::size_t i) const
    {
        return well_indices_[i];
    }

    /// \brief Get a begin iterator
    const_iterator begin() const
    {
        return well_indices_.begin();
    }

    /// \brief Get the end iterator
    const_iterator end() const
    {
        return well_indices_.end();
    }

    /// \breif Get the number of wells
    std::size_t size() const
    {
        return well_indices_.size();
    }
private:
    /// Stores at index i all cells that are perforated by
    /// the well at position i of the eclipse schedule.
    std::vector<std::set<int> > well_indices_;
};


#ifdef HAVE_MPI
/// \brief Determines the wells that have perforate cells for each process.
///
/// On the root process omputes for all processes all indices of wells that
/// will perforate local cells.
/// Note that a well might perforate local cells of multiple processes
///
/// \param parts The partition number for each cell
/// \param well The eclipse information about the wells.
/// \param cpGrid The unbalanced grid we compute on.
/// \return On the rank that has the global grid a vector with the well
///         indices for process i at index i.
std::vector<std::vector<int> >
perforatingWellIndicesOnProc(const std::vector<int>& parts,
                  const std::vector<Dune::cpgrid::OpmWellType>& wells,
                  const CpGrid& cpgrid);

/// \brief Computes wells assigned to processes.
///
/// Computes for all processes all indices of wells that
/// will be assigned to this process.
/// \param parts The partition number for each cell
/// \param gid Functor that turns cell index to global id.
/// \param eclipseState The eclipse information
/// \param well_connecton The information about the perforations of each well.
/// \param exportList List of cells to be exported. Each entry is a tuple of the
///                   global cell index, the process rank to export to, and the
///                   attribute on that rank (assumed to be owner)
/// \param exportList List of cells to be imported. Each entry is a tuple of the
///                   global cell index, the process rank that exports it, and the
///                   attribute on this rank (assumed to be owner)
/// \param cc Information about the parallelism together with the decomposition.
/// \return On rank 0 a vector with the well indices for process i
///         at index i.
std::vector<std::vector<int> >
postProcessPartitioningForWells(std::vector<int>& parts,
                                std::function<int(int)> gid,
                                const std::vector<OpmWellType>&  wells,
                                const WellConnections& well_connections,
                                std::vector<std::tuple<int,int,char>>& exportList,
                                std::vector<std::tuple<int,int,char,int>>& importList,
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 7)
                                const Communication<MPI_Comm>& cc);
#else
                                const CollectiveCommunication<MPI_Comm>& cc);
#endif

/// \brief Computes whether wells are perforating cells on this process.
/// \param wells_on_proc well indices assigned to each process
/// \param eclipseState The eclipse information
/// \param cc The communicator
/// \param root The rank of the process that has the complete partitioning
///             information.
/// \return Vector of pairs of well name and a boolean indicating whether the
///         well with this name perforates cells here. Sorted by well name!
std::vector<std::pair<std::string,bool>>
computeParallelWells(const std::vector<std::vector<int> >& wells_on_proc,
                     const std::vector<OpmWellType>&  wells,
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 7)
                     const Communication<MPI_Comm>& cc,
#else
                     const CollectiveCommunication<MPI_Comm>& cc,
#endif
                     int root);
#endif
} // end namespace cpgrid
} // end namespace Dune
#endif
