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

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include <dune/common/parallel/mpicollectivecommunication.hh>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
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
    /// \param eclipseState The eclipse information
    /// \param cartesianSize The logical cartesian size of the grid.
    /// \param cartesian_to_compressed Mapping of cartesian index
    ///        compressed cell index. The compressed index is used
    ///        to represent the well conditions.
    WellConnections(const Opm::EclipseState& eclipseState,
                    const std::array<int, 3>& cartesianSize,
                    const std::vector<int>& cartesian_to_compressed);

    /// \brief Initialze the data of the container
    /// \param eclipseState The eclipse information
    /// \param cartesianSize The logical cartesian size of the grid.
    /// \param cartesian_to_compressed Mapping of cartesian index
    ///        compressed cell index. The compressed index is used
    ///        to represent the well conditions.
    void init(const Opm::EclipseState& eclipseState,
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

/// \brief Computes wells assigned to processes.
///
/// Computes for all processes all indices of wells that
/// will be assigned to this process.
/// \param parts The partition number for each cell
/// \param eclipseState The eclipse information
/// \param well_connecton The informatio about the perforations of each well.
/// \param no_procs The number of processes.
std::vector<std::vector<int> >
postProcessPartitioningForWells(std::vector<int>& parts,
                                const Opm::EclipseState& eclipseState,
                                const WellConnections& well_connections,
                                std::size_t no_procs);

#ifdef HAVE_MPI
/// \brief Computes that names that of all wells not handled by this process
/// \param wells_on_proc well indices assigned to each process
/// \param eclipseState The eclipse information
/// \param cc The communicator
/// \param root The rank of the process that has the complete partitioning
///             information.
std::unordered_set<std::string>
computeDefunctWellNames(const std::vector<std::vector<int> >& wells_on_proc,
                        const Opm::EclipseState& eclipseState,
                        const CollectiveCommunication<MPI_Comm>& cc,
                        int root);
#endif
} // end namespace cpgrid
} // end namespace Dune
#endif
