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

class WellConnections
{
public:
    typedef std::vector<std::set<int> >::const_iterator const_iterator;
    typedef const_iterator iterator;

    WellConnections(const Opm::EclipseStateConstPtr eclipseState,
                    const std::array<int, 3>& cartesianSize,
                    const std::vector<int>& cartesian_to_compressed);

    const std::set<int>& operator[](std::size_t i) const
    {
        return well_indices_[i];
    }

    const_iterator begin() const
    {
        return well_indices_.begin();
    }

    const_iterator end() const
    {
        return well_indices_.end();
    }

    std::size_t size() const
    {
        return well_indices_.size();
    }
private:
    std::vector<std::set<int> > well_indices_;
};

std::vector<std::vector<int> >
postProcessPartitioningForWells(std::vector<int>& parts,
                                const Opm::EclipseStateConstPtr eclipseState,
                                const WellConnections& well_connections,
                                std::size_t no_procs);

#ifdef HAVE_MPI
std::unordered_set<std::string>
computeDefunctWellNames(const std::vector<std::vector<int> >& wells_on_proc,
                        const Opm::EclipseStateConstPtr eclipseState,
                        const CollectiveCommunication<MPI_Comm>& cc,
                        int root);
#endif
} // end namespace cpgrid
} // end namespace Dune
#endif
