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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <map>

#include <dune/grid/common/WellConnections.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/CompletionSet.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>

namespace Dune
{
namespace cpgrid
{
WellConnections::WellConnections(const Opm::EclipseStateConstPtr eclipseState,
                                 const std::array<int, 3>& cartesianSize,
                                 const std::vector<int>& cartesian_to_compressed)
{
    init(eclipseState, cartesianSize, cartesian_to_compressed);
}

void WellConnections::init(const Opm::EclipseStateConstPtr eclipseState,
                           const std::array<int, 3>& cartesianSize,
                           const std::vector<int>& cartesian_to_compressed)
{
    std::vector<const Opm::Well*>  wells  = eclipseState->getSchedule()->getWells();
    int last_time_step = eclipseState->getSchedule()->getTimeMap()->size()-1;
    well_indices_.resize(wells.size());

    // We assume that we know all the wells.
    int index=0;
    for (const auto well : wells) {
        std::set<int>& well_indices = well_indices_[index];
        Opm::CompletionSetConstPtr completionSet = well->getCompletions(last_time_step);
        for (size_t c=0; c<completionSet->size(); c++) {
            Opm::CompletionConstPtr completion = completionSet->get(c);
            int i = completion->getI();
            int j = completion->getJ();
            int k = completion->getK();
            int cart_grid_idx = i + cartesianSize[0]*(j + cartesianSize[1]*k);
            int compressed_idx = cartesian_to_compressed[cart_grid_idx];
            if ( compressed_idx >= 0 ) // Ignore completions in inactive cells.
            {
                well_indices.insert(compressed_idx);
            }
        }
        ++index;
    }
}

std::vector<std::vector<int> >
postProcessPartitioningForWells(std::vector<int>& parts,
                                const Opm::EclipseStateConstPtr eclipseState,
                                const WellConnections& well_connections,
                                std::size_t no_procs)
{
    std::vector<const Opm::Well*>  wells  = eclipseState->getSchedule()->getWells();
    // Contains for each process the indices of the wells assigned to it.
    std::vector<std::vector<int> > well_indices_on_proc(no_procs);

    if( ! well_connections.size() )
    {
        // No wells to be processed
        return well_indices_on_proc;
    }

    // prevent memory allocation
    for(auto& well_indices : well_indices_on_proc)
    {
        well_indices.reserve(wells.size());
    }

    // Check that all completions of a well have ended up on one process.
    // If that is not the case for well then move them manually to the
    // process that already has the most completions on it.
    int well_index = 0;

    for (const auto* well: wells) {
        const auto& well_completions = well_connections[well_index];
        std::map<int,std::size_t> no_completions_on_proc;
        for ( auto completion_index: well_completions )
        {
            ++no_completions_on_proc[parts[completion_index]];
        }

        int owner = no_completions_on_proc.begin()->first;

        if ( no_completions_on_proc.size() > 1 )
        {
            // partition with the most completions on it becomes new owner
            int new_owner = std::max_element(no_completions_on_proc.begin(),
                                             no_completions_on_proc.end(),
                                             [](const std::pair<int,std::size_t>& p1,
                                                const std::pair<int,std::size_t>& p2){
                                                 return ( p1.second > p2.second );
                                             })->first;
            std::cout << "Manually moving well " << well->name() << " to partition "
                      << new_owner << std::endl;

            for ( auto completion_cell : well_completions )
            {
                parts[completion_cell] = new_owner;
            }

            owner = new_owner;
        }

        well_indices_on_proc[owner].push_back(well_index);
        ++well_index;
    }
    return well_indices_on_proc;

}

#ifdef HAVE_MPI
std::unordered_set<std::string>
computeDefunctWellNames(const std::vector<std::vector<int> >& wells_on_proc,
                        const Opm::EclipseStateConstPtr eclipseState,
                        const CollectiveCommunication<MPI_Comm>& cc,
                        int root)
{
    std::vector<const Opm::Well*>  wells  = eclipseState->getSchedule()->getWells();
    std::vector<int> my_well_indices;
    const int well_information_tag = 267553;

    if( root == cc.rank() )
    {
        std::vector<MPI_Request> reqs(cc.size(), MPI_REQUEST_NULL);
        my_well_indices = wells_on_proc[root];
        for ( int i=0; i < cc.size(); ++i )
        {
            if(i==root)
            {
                continue;
            }
            MPI_Isend(const_cast<int*>(wells_on_proc[i].data()),
                      wells_on_proc[i].size(),
                      MPI_INT, i, well_information_tag, cc, &reqs[i]);
        }
        std::vector<MPI_Status> stats(reqs.size());
        MPI_Waitall(reqs.size(), reqs.data(), stats.data());
    }
    else
    {
        MPI_Status stat;
        MPI_Probe(root, well_information_tag, cc, &stat);
        int msg_size;
        MPI_Get_count(&stat, MPI_INT, &msg_size);
        my_well_indices.resize(msg_size);
        MPI_Recv(my_well_indices.data(), msg_size, MPI_INT, root,
                 well_information_tag, cc, &stat);
    }

    // Compute defunct wells in parallel run.
    std::vector<int> defunct_wells(wells.size(), true);

    for(auto well_index : my_well_indices)
    {
        defunct_wells[well_index] = false;
    }

    // We need to use well names as only they are consistent.
    std::unordered_set<std::string> defunct_well_names;

    for(auto defunct = defunct_wells.begin(); defunct != defunct_wells.end(); ++defunct)
    {
        if ( *defunct )
        {
            defunct_well_names.insert(wells[defunct-defunct_wells.begin()]->name());
        }
    }

    return defunct_well_names;
}
#endif
} // end namespace cpgrid
} // end namespace Dune
