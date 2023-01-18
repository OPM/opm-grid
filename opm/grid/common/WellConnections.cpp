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
#include <opm/grid/common/WellConnections.hpp>

#include <opm/grid/cpgrid/CpGridData.hpp>

#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <dune/common/parallel/mpitraits.hh>
#include <dune/istl/owneroverlapcopy.hh>

#include <map>

namespace
{

struct Less
{
    template<typename T>
    bool operator()( const T& t1, const T& t2 )
    {
        return std::get<0>(t1) < std::get<0>(t2);
    }
    template<typename T>
    bool operator()( const T& t, int i )
    {
        return std::get<0>(t) < i;
    }
    template<typename T>
    bool operator()( int i, const T& t )
    {
        return i < std::get<0>(t);
    }
};
}

namespace Dune
{
namespace cpgrid
{
WellConnections::WellConnections(const std::vector<OpmWellType>& wells,
                                 const std::array<int, 3>& cartesianSize,
                                 const std::vector<int>& cartesian_to_compressed)
{
    init(wells, cartesianSize, cartesian_to_compressed);
}

WellConnections::WellConnections(const std::vector<OpmWellType>& wells,
                                 const Dune::CpGrid& cpGrid)
{
    const auto& cpgdim = cpGrid.logicalCartesianSize();
    // create compressed lookup from cartesian.
    std::vector<int> cartesian_to_compressed(cpgdim[0]*cpgdim[1]*cpgdim[2], -1);

    for( int i=0; i < cpGrid.numCells(); ++i )
    {
        cartesian_to_compressed[cpGrid.globalCell()[i]] = i;
    }
    init(wells, cpgdim, cartesian_to_compressed);
}

void WellConnections::init([[maybe_unused]] const std::vector<OpmWellType>& wells,
                           [[maybe_unused]] const std::array<int, 3>& cartesianSize,
                           [[maybe_unused]] const std::vector<int>& cartesian_to_compressed)
{
#if HAVE_ECL_INPUT
    well_indices_.resize(wells.size());

    // We assume that we know all the wells.
    int index=0;
    for (const auto& well : wells) {
        std::set<int>& well_indices = well_indices_[index];
        const auto& connectionSet = well.getConnections( );
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

#ifdef HAVE_MPI
std::vector<std::vector<int> >
perforatingWellIndicesOnProc(const std::vector<int>& parts,
                             const std::vector<Dune::cpgrid::OpmWellType>& wells,
                             const CpGrid& cpGrid)
{
    auto numProcs = cpGrid.comm().size();
    std::vector<std::vector<int> > wellIndices(numProcs);

    if (cpGrid.numCells())
    {
        // root process that has global cells
        WellConnections wellConnections(wells, cpGrid);
        if (!wellConnections.size())
        {
            return wellIndices;
        }
        // prevent memory allocation
        for (auto &well_indices : wellIndices) {
            well_indices.reserve(wells.size());
        }

        for (std::size_t wellIndex = 0; wellIndex < wells.size(); ++wellIndex) {
            const auto &connections = wellConnections[wellIndex];
            std::map<int, std::size_t> connectionsOnProc;
            for (const auto& connection_index : connections) {
                ++connectionsOnProc[parts[connection_index]];
            }

            for (const auto& entry: connectionsOnProc)
            {
                if (entry.second > 0) // Should be unnecessary
                {
                    wellIndices[entry.first].push_back(wellIndex);
                }
            }
        }
    }
    return wellIndices;
}
std::vector<std::vector<int> >
postProcessPartitioningForWells(std::vector<int>& parts,
                                [[maybe_unused]] std::function<int(int)> gid,
                                [[maybe_unused]] const std::vector<OpmWellType>& wells,
                                [[maybe_unused]] const WellConnections& well_connections,
                                [[maybe_unused]] std::vector<std::tuple<int,int,char>>& exportList,
                                [[maybe_unused]] std::vector<std::tuple<int,int,char,int>>& importList,
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 7)
                                const Communication<MPI_Comm>& cc)
#else
                                const CollectiveCommunication<MPI_Comm>& cc)
#endif
{
    auto no_procs = cc.size();
    auto noCells = parts.size();
    std::vector<std::size_t> cellsPerProc(no_procs);
    cc.allgather(&noCells, 1, cellsPerProc.data());

    // Contains for each process the indices of the wells assigned to it.
    std::vector<std::vector<int> > well_indices_on_proc(no_procs);

#if HAVE_ECL_INPUT
    const auto& mpiType =  MPITraits<std::size_t>::getType();
    std::map<int, std::vector<int>> addCells, removeCells;
    using AttributeSet = CpGridData::AttributeSet;

    if (noCells && well_connections.size()) {

        // prevent memory allocation
        for (auto &well_indices : well_indices_on_proc) {
            well_indices.reserve(wells.size());
        }

        // Check that all connections of a well have ended up on one process.
        // If that is not the case for well then move them manually to the
        // process that already has the most connections on it.
        int well_index = 0;

        for (const auto &well : wells) {
            const auto &connections = well_connections[well_index];
            if (connections.empty()) {
                // No connections, nothing to move or worry about.
                ++well_index;
                continue;
            }
            std::map<int, std::size_t> no_connections_on_proc;
            for (auto connection_index : connections) {
                ++no_connections_on_proc[parts[connection_index]];
            }
            assert(!no_connections_on_proc.empty());
            int owner = no_connections_on_proc.begin()->first;

            if (no_connections_on_proc.size() > 1) {
                // partition with the most connections on it becomes new owner
                int new_owner = std::max_element(no_connections_on_proc.begin(),
                                                 no_connections_on_proc.end(),
                                                 [](const auto& p1, const auto& p2)
                                                 { return (p1.second < p2.second); })
                                    ->first;
                std::cout << "Manually moving well " << well.name()
                          << " to partition " << new_owner << std::endl;

                // all cells moving to new_owner. Might already contain cells from
                // previous wells.
                auto &add = addCells[new_owner];
                auto addOldSize = add.size(); // remember beginning of this well

                for (auto connection_cell : connections) {
                    const auto &global = gid(connection_cell);
                    auto old_owner = parts[connection_cell];
                    if (old_owner != new_owner) // only parts might be moved
                    {
                        removeCells[old_owner].push_back(global);
                        add.push_back(global);
                        parts[connection_cell] = new_owner;
                    }
                }
                auto oldEnd = add.begin()+addOldSize;
                std::sort(oldEnd, add.end()); // we need ascending order
                auto exportCandidate =  exportList.begin();

                for (auto movedCell = oldEnd; movedCell != add.end(); ++movedCell) {
                    exportCandidate = std::lower_bound(exportCandidate, exportList.end(), *movedCell,
                                                       Less());
                    assert(exportCandidate != exportList.end() && std::get<0>(*exportCandidate) == *movedCell);
                    std::get<1>(*exportCandidate) = new_owner;
                }
                owner = new_owner;
            }

            well_indices_on_proc[owner].push_back(well_index);
            ++well_index;
        }
        auto sorter = [](std::pair<const int, std::vector<int>> &pair) {
                          auto &vec = pair.second;
                          std::sort(vec.begin(), vec.end());
                      };
        std::for_each(addCells.begin(), addCells.end(), sorter);
        std::for_each(removeCells.begin(), removeCells.end(), sorter);
    }

    // setup receives for each process that owns cells of the original grid
    auto noSource = std::count_if(cellsPerProc.begin(), cellsPerProc.end(),
                                  [](const std::size_t &i) { return i > 0; });
    std::vector<MPI_Request> requests(noSource, MPI_REQUEST_NULL);
    std::vector<std::vector<std::size_t>> sizeBuffers(cc.size());
    auto begin = cellsPerProc.begin();
    auto req = requests.begin();
    int tag = 7823;

    for (auto it = begin, end = cellsPerProc.end(); it != end; ++it) {
        auto otherRank = it - begin;
        if ( *it > 0 ) {
            sizeBuffers[otherRank].resize(2);
            MPI_Irecv(sizeBuffers[otherRank].data(), 2, mpiType, otherRank, tag, cc, &(*req));
            ++req;
        }
    }

    // Send the sizes
    if (!parts.empty()) {
        for (int otherRank = 0; otherRank < cc.size(); ++otherRank)
        {
            std::size_t sizes[2] = {0, 0};
            auto candidate = addCells.find(otherRank);
            if (candidate != addCells.end())
                sizes[0] = candidate->second.size();

            candidate = removeCells.find(otherRank);
            if (candidate != removeCells.end())
                sizes[1] = candidate->second.size();
            MPI_Send(sizes, 2, mpiType, otherRank, tag, cc);
        }
    }
    std::vector<MPI_Status> statuses(requests.size());
    MPI_Waitall(requests.size(), requests.data(), statuses.data());

    auto messages = std::count_if(
                                  sizeBuffers.begin(), sizeBuffers.end(),
                                  [](const std::vector<std::size_t> &v) { return v.size() ? v[0] + v[1] : 0; });
    requests.resize(messages);
    ++tag;

    req = requests.begin();
    std::vector<std::vector<std::size_t>> cellIndexBuffers(messages); // receive buffers for indices of each rank.
    auto cellIndexBufferIt = cellIndexBuffers.begin();

    for (auto it = begin, end = cellsPerProc.end(); it != end; ++it) {
        auto otherRank = it - begin;
        const auto& sizeBuffer = sizeBuffers[otherRank];
        if ( sizeBuffer.size() >= 2 && (sizeBuffer[0] + sizeBuffer[1])) {
            auto &cellIndexBuffer = *cellIndexBufferIt;
            cellIndexBuffer.resize(sizeBuffer[0] + sizeBuffer[1]);
            MPI_Irecv(cellIndexBuffer.data(), cellIndexBuffer.size(), mpiType, otherRank, tag, cc,
                      &(*req));
            ++req;
            ++cellIndexBufferIt;
        }
    }

    // Send data if we have cells.
    if (!parts.empty()) {
        for (int otherRank = 0; otherRank < cc.size(); ++otherRank)
        {
            std::vector<std::size_t> buffer;
            auto candidate = addCells.find(otherRank);
            if (candidate != addCells.end())
                buffer.insert(buffer.end(), candidate->second.begin(),
                              candidate->second.end());

            candidate = removeCells.find(otherRank);
            if (candidate != removeCells.end())
                buffer.insert(buffer.end(), candidate->second.begin(),
                              candidate->second.end());

            if (!buffer.empty()) {
                MPI_Send(buffer.data(), buffer.size(), mpiType, otherRank, tag, cc);
            }
        }
    }
    statuses.resize(messages);
    // Wait for the messages
    MPI_Waitall(requests.size(), requests.data(), statuses.data());

    // unpack data
    auto status = statuses.begin();
    for (const auto &cellIndexBuffer : cellIndexBuffers) {
        int otherRank = status->MPI_SOURCE;
        if (!cellIndexBuffer.empty()) {
            // add cells that moved here
            auto noAdded = sizeBuffers[otherRank][0];
            importList.reserve(importList.size() + noAdded);
            auto middle = importList.end();
            std::vector<std::tuple<int, int, char>> addToImport;
            std::size_t offset = 0;
            for (; offset != noAdded; ++offset)
                importList.emplace_back(cellIndexBuffer[offset], otherRank,
                                        AttributeSet::owner, -1);
            std::inplace_merge(importList.begin(), middle, importList.end(),
                               Less());

            // remove cells that moved to another process
            auto noRemoved = sizeBuffers[otherRank][1];
            if (noRemoved) {
                std::vector<std::tuple<int, int, char, int>> tmp(importList.size());
                auto newEnd =
                    std::set_difference(importList.begin(), importList.end(),
                                        cellIndexBuffer.begin() + noAdded,
                                        cellIndexBuffer.end(), tmp.begin(), Less());
                tmp.resize(newEnd - tmp.begin());
                importList.swap(tmp);
            }
        }
        ++status;
    }
#endif

    return well_indices_on_proc;

}

std::vector<std::pair<std::string,bool>>
computeParallelWells([[maybe_unused]] const std::vector<std::vector<int> >& wells_on_proc,
                     [[maybe_unused]] const std::vector<OpmWellType>& wells,
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 7)
                     [[maybe_unused]] const Communication<MPI_Comm>& cc,
#else
                     [[maybe_unused]] const CollectiveCommunication<MPI_Comm>& cc,
#endif
                     [[maybe_unused]] int root)
{
    // We need to use well names as only they are consistent.
    std::vector<std::pair<std::string,bool>> parallel_wells;

#if HAVE_ECL_INPUT
    std::vector<int> my_well_indices;
    std::vector<std::string> globalWellNames;
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
        // Broadcast well names
        // 1. Compute packed size and broadcast
        std::size_t sizes[2] = {wells.size(),0};
        int wellMessageSize = 0;
        MPI_Pack_size(2, MPITraits<std::size_t>::getType(), cc, &wellMessageSize);
        for(const auto& well: wells)
        {
            int size;
            MPI_Pack_size(well.name().size() + 1, MPI_CHAR, cc, &size); // +1 for '\0' delimiter
            sizes[1] += well.name().size() + 1;
            wellMessageSize += size;
        }
        MPI_Bcast(&wellMessageSize, 1, MPI_INT, root, cc);
        // 2. Send number of wells and their names in one message
        globalWellNames.reserve(wells.size());
        std::vector<char> buffer(wellMessageSize);
        int pos = 0;
        MPI_Pack(&sizes, 2, MPITraits<std::size_t>::getType(), buffer.data(), wellMessageSize, &pos, cc);
        for(const auto& well: wells)
        {
            MPI_Pack(well.name().c_str(), well.name().size()+1, MPI_CHAR, buffer.data(),
                     wellMessageSize, &pos, cc); // +1 for '\0' delimiter
            globalWellNames.push_back(well.name());
        }
        MPI_Bcast(buffer.data(), wellMessageSize, MPI_PACKED, root, cc);
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

        // 1. receive broadcasted message Size
        int wellMessageSize;
        MPI_Bcast(&wellMessageSize, 1, MPI_INT, root, cc);

        // 2. Receive number of wells and their names in one message
        globalWellNames.reserve(wells.size());
        std::vector<char> buffer(wellMessageSize);
        MPI_Bcast(buffer.data(), wellMessageSize, MPI_PACKED, root, cc);
        std::size_t sizes[2];
        int pos = 0;
        MPI_Unpack(buffer.data(), wellMessageSize, &pos, &sizes, 2, MPITraits<std::size_t>::getType(), cc);
        // unpack all string at once
        std::vector<char> cstr(sizes[1]);
        MPI_Unpack(buffer.data(), wellMessageSize, &pos, cstr.data(), sizes[1], MPI_CHAR, cc);
        pos = 0;
        for(std::size_t i = 0; i < sizes[0]; ++i)
        {
            globalWellNames.emplace_back(cstr.data()+pos);
            pos += globalWellNames.back().size() + 1; // +1 because of '\0' delimiter
        }
    }

    // Compute wells active/inactive in parallel run.
    // boolean indicates whether the well perforates local cells
    parallel_wells.reserve(wells.size());

    for(const auto& well_name: globalWellNames)
    {
        parallel_wells.emplace_back(well_name, false);
    }

    for(auto well_index : my_well_indices)
    {
        parallel_wells[well_index].second = true;
    }

    std::sort(parallel_wells.begin(), parallel_wells.end());
#ifndef NDEBUG
    std::string last;
    for(const auto& wpair: parallel_wells)
    {
        assert(last != wpair.first);
        last = wpair.first;
    }
#endif
#endif

    return parallel_wells;
}
#endif
} // end namespace cpgrid
} // end namespace Dune
