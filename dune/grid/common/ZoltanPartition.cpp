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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <dune/grid/common/ZoltanPartition.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#if defined(HAVE_ZOLTAN) && defined(HAVE_MPI)
namespace Dune
{
namespace cpgrid
{
std::pair<std::vector<int>, std::vector<int> >
zoltanGraphPartitionGridOnRoot(const CpGrid& cpgrid,
                               const Opm::EclipseStateConstPtr eclipseState,
                               const double* transmissibilities,
                               const CollectiveCommunication<MPI_Comm>& cc,
                               int root)
{
    int rc = ZOLTAN_OK - 1;
    float ver = 0;
    struct Zoltan_Struct *zz;
    int changes, numGidEntries, numLidEntries, numImport, numExport;
    ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
    int *importProcs, *importToPart, *exportProcs, *exportToPart;
    int argc=0;
    char** argv = 0 ;
    rc = Zoltan_Initialize(argc, argv, &ver);
    zz = Zoltan_Create(cc);
    if ( rc != ZOLTAN_OK )
    {
        OPM_THROW(std::runtime_error, "Could not initialize Zoltan!");
    }

    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
    Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
    Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "3");
    Zoltan_Set_Param(zz, "CHECK_GRAPH", "2");
    Zoltan_Set_Param(zz,"EDGE_WEIGHT_DIM","0");
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
    Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", ".35");  /* 0-remove all, 1-remove none */

    bool pretendEmptyGrid = cc.rank()!=root;
    std::shared_ptr<CombinedGridWellGraph> grid_and_wells;

    if( eclipseState )
    {
        Zoltan_Set_Param(zz,"EDGE_WEIGHT_DIM","1");
        grid_and_wells.reset(new CombinedGridWellGraph(cpgrid, eclipseState,
                                                       transmissibilities, pretendEmptyGrid));
        Dune::cpgrid::setCpGridZoltanGraphFunctions(zz, *grid_and_wells,
                                                    pretendEmptyGrid);
    }
    else
    {
        Dune::cpgrid::setCpGridZoltanGraphFunctions(zz, cpgrid, pretendEmptyGrid);
    }

    rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
                             &changes,        /* 1 if partitioning was changed, 0 otherwise */
                             &numGidEntries,  /* Number of integers used for a global ID */
                             &numLidEntries,  /* Number of integers used for a local ID */
                             &numImport,      /* Number of vertices to be sent to me */
                             &importGlobalGids,  /* Global IDs of vertices to be sent to me */
                             &importLocalGids,   /* Local IDs of vertices to be sent to me */
                             &importProcs,    /* Process rank for source of each incoming vertex */
                             &importToPart,   /* New partition for each incoming vertex */
                             &numExport,      /* Number of vertices I must send to other processes*/
                             &exportGlobalGids,  /* Global IDs of the vertices I must send */
                             &exportLocalGids,   /* Local IDs of the vertices I must send */
                             &exportProcs,    /* Process to which I send each of the vertices */
                             &exportToPart);  /* Partition to which each vertex will belong */
    int                         size = cpgrid.numCells();
    int                         rank  = cc.rank();
    std::vector<int>            parts(size, rank);
    std::vector<std::vector<int> > wells_on_proc;

    for ( int i=0; i < numExport; ++i )
    {
        parts[exportLocalGids[i]] = exportProcs[i];
    }

    if( eclipseState && ! pretendEmptyGrid )
    {
        wells_on_proc = grid_and_wells->postProcessPartitioningForWells(parts);
#ifndef NDEBUG
        int index = 0;
        for( auto well : grid_and_wells->getWellsGraph() )
        {
            int part=parts[index];
            std::set<std::pair<int,int> > cells_on_other;
            for( auto vertex : well )
            {
                if( part != parts[vertex] )
                {
                    cells_on_other.insert(std::make_pair(vertex, parts[vertex]));
                }
            }
            if ( cells_on_other.size() )
            {
                OPM_THROW(std::domain_error, "Well is distributed between processes, which should not be the case!");
            }
            ++index;
        }
#endif
    }
    // free space allocated for zoltan.
    Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);
    Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, &importProcs, &importToPart);
    Zoltan_Destroy(&zz);

    cc.broadcast(&parts[0], parts.size(), root);
    std::vector<int> my_well_indices;

    if( pretendEmptyGrid)
    {
        std::vector<MPI_Request> reqs(cc.size(), MPI_REQUEST_NULL);
        my_well_indices = wells_on_proc[root];
        for ( int i=0; i < cc.size(); ++i )
        {
            if(i==root)
            {
                continue;
            }
            MPI_Isend(wells_on_proc[i].data(), wells_on_proc[i].size(),
                      MPI_INT, i, 267553, cc, &reqs[i]);
        }
        std::vector<MPI_Status> stats(reqs.size());
        MPI_Waitall(reqs.size(), reqs.data(), stats.data());
    }
    else
    {
        MPI_Status stat;
        MPI_Probe(root, 267553, cc, &stat);
        int msg_size;
        MPI_Get_count(&stat, MPI_INT, &msg_size);
        my_well_indices.resize(msg_size);
        MPI_Recv(my_well_indices.data(), msg_size, MPI_INT, root, 267553,
                 cc, &stat);
    }

    return std::make_pair(parts, my_well_indices);
}
}
}
#endif // HAVE_ZOLTAN
