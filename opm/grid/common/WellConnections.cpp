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

#include <opm/grid/common/WellConnections.hpp>

#include <opm/grid/utility/OpmParserIncludes.hpp>

namespace Dune
{
namespace cpgrid
{
WellConnections::WellConnections(const std::vector<const OpmWellType*>& wells,
                                 const std::array<int, 3>& cartesianSize,
                                 const std::vector<int>& cartesian_to_compressed)
{
    init(wells, cartesianSize, cartesian_to_compressed);
}

void WellConnections::init(const std::vector<const OpmWellType*>& wells,
                           const std::array<int, 3>& cartesianSize,
                           const std::vector<int>& cartesian_to_compressed)
{
#if HAVE_ECL_INPUT
    well_indices_.resize(wells.size());

    // We assume that we know all the wells.
    int index=0;
    for (const auto well : wells) {
        std::set<int>& well_indices = well_indices_[index];
        const auto& connectionSet = well->getConnections( );
        for (size_t c=0; c<connectionSet.size(); c++) {
            const auto& connection = connectionSet.get(c);
            int i = connection.getI();
            int j = connection.getJ();
            int k = connection.getK();
            int cart_grid_idx = i + cartesianSize[0]*(j + cartesianSize[1]*k);
            int compressed_idx = cartesian_to_compressed[cart_grid_idx];
            if ( compressed_idx >= 0 ) // Ignore connections in inactive cells.
            {
                well_indices.insert(compressed_idx);
            }
        }
        ++index;
    }
#endif
}

std::vector<std::vector<int> >
postProcessPartitioningForWells(std::vector<int>& parts,
                                const std::vector<const OpmWellType*>& wells,
                                const WellConnections& well_connections,
                                std::size_t no_procs)
{
    // Contains for each process the indices of the wells assigned to it.
    std::vector<std::vector<int> > well_indices_on_proc(no_procs);

#if HAVE_ECL_INPUT
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

    // Check that all connections of a well have ended up on one process.
    // If that is not the case for well then move them manually to the
    // process that already has the most connections on it.
    int well_index = 0;

    for (const auto* well: wells) {
        const auto& connections = well_connections[well_index];
        std::map<int,std::size_t> no_connections_on_proc;
        for ( auto connection_index: connections )
        {
            ++no_connections_on_proc[parts[connection_index]];
        }

        int owner = no_connections_on_proc.begin()->first;

        if ( no_connections_on_proc.size() > 1 )
        {
            // partition with the most connections on it becomes new owner
            int new_owner = std::max_element(no_connections_on_proc.begin(),
                                             no_connections_on_proc.end(),
                                             [](const std::pair<int,std::size_t>& p1,
                                                const std::pair<int,std::size_t>& p2){
                                                 return ( p1.second > p2.second );
                                             })->first;
            std::cout << "Manually moving well " << well->name() << " to partition "
                      << new_owner << std::endl;

            for ( auto connection_cell : connections )
            {
                parts[connection_cell] = new_owner;
            }

            owner = new_owner;
        }

        well_indices_on_proc[owner].push_back(well_index);
        ++well_index;
    }
#endif

    return well_indices_on_proc;

}

#ifdef HAVE_MPI
std::unordered_set<std::string>
computeDefunctWellNames(const std::vector<std::vector<int> >& wells_on_proc,
                        const std::vector<const OpmWellType*>& wells,
                        const CollectiveCommunication<MPI_Comm>& cc,
                        int root)
{
    // We need to use well names as only they are consistent.
    std::unordered_set<std::string> defunct_well_names;

#if HAVE_ECL_INPUT
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

    for(auto defunct = defunct_wells.begin(); defunct != defunct_wells.end(); ++defunct)
    {
        if ( *defunct )
        {
            defunct_well_names.insert(wells[defunct-defunct_wells.begin()]->name());
        }
    }
#endif

    return defunct_well_names;
}
#endif
} // end namespace cpgrid
} // end namespace Dune
