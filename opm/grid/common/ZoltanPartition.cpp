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
#include <opm/grid/common/ZoltanPartition.hpp>

#include <opm/grid/utility/OpmParserIncludes.hpp>

#if defined(HAVE_ZOLTAN) && defined(HAVE_MPI)
namespace Dune
{
namespace cpgrid
{
std::tuple<std::vector<int>, std::unordered_set<std::string>,
           std::vector<std::tuple<int,int,char> >,
           std::vector<std::tuple<int,int,char,int> > >
zoltanGraphPartitionGridOnRoot(const CpGrid& cpgrid,
                               const std::vector<OpmWellType> * wells,
                               const double* transmissibilities,
                               const CollectiveCommunication<MPI_Comm>& cc,
                               EdgeWeightMethod edgeWeightsMethod, int root)
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
    Zoltan_Set_Param(zz, "IMBALANCE_TOL", "1.1");
    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
    Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
    Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
    Zoltan_Set_Param(zz, "CHECK_GRAPH", "2");
    Zoltan_Set_Param(zz,"EDGE_WEIGHT_DIM","0");
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
    Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", ".35");  /* 0-remove all, 1-remove none */

    // For the load balancer one process has the whole grid and
    // all others an empty partition before loadbalancing.
    bool partitionIsEmpty     = cc.rank()!=root;

    std::shared_ptr<CombinedGridWellGraph> grid_and_wells;

    if( wells )
    {
        Zoltan_Set_Param(zz,"EDGE_WEIGHT_DIM","1");
        grid_and_wells.reset(new CombinedGridWellGraph(cpgrid,
                                                       wells,
                                                       transmissibilities,
                                                       partitionIsEmpty,
                                                       edgeWeightsMethod));
        Dune::cpgrid::setCpGridZoltanGraphFunctions(zz, *grid_and_wells,
                                                    partitionIsEmpty);
    }
    else
    {
        Dune::cpgrid::setCpGridZoltanGraphFunctions(zz, cpgrid, partitionIsEmpty);
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
    // List entry: (global) index, process to export to, attribute there
    std::vector<std::tuple<int,int,char>> myExportList(numExport);
    // List entry: process to import from, global index, attribute here, local index
    // (determined later)
    std::vector<std::tuple<int,int,char,int>>myImportList(numImport);
    myExportList.reserve(1.2*myExportList.size());
    myImportList.reserve(1.2*myImportList.size());
    using AttributeSet = Dune::OwnerOverlapCopyAttributeSet::AttributeSet;

    for ( int i=0; i < numExport; ++i )
    {
        parts[exportLocalGids[i]] = exportProcs[i];
        myExportList[i] = std::make_tuple(exportGlobalGids[i], exportProcs[i], static_cast<char>(AttributeSet::owner));
    }

    for ( int i=0; i < numImport; ++i )
    {
        myImportList[i] = std::make_tuple(importGlobalGids[i], importProcs[i], static_cast<char>(AttributeSet::owner),-1);
    }
    // Add cells that stay here to the lists. Somehow I could not persuade Zoltan to do this.
    for ( std::size_t i = 0; i < parts.size(); ++i)
    {
        if ( parts[i] == rank )
        {
            myExportList.emplace_back(i, rank, static_cast<char>(AttributeSet::owner) );
            myImportList.emplace_back(i, rank, static_cast<char>(AttributeSet::owner), -1 );
        }
    }

    std::inplace_merge(myImportList.begin(), myImportList.begin() + numImport, myImportList.end());
    std::inplace_merge(myExportList.begin(), myExportList.begin() + numExport, myExportList.end());
    // free space allocated for zoltan.
    Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);
    Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, &importProcs, &importToPart);
    Zoltan_Destroy(&zz);

    if( wells )
    {
        wells_on_proc =
            postProcessPartitioningForWells(parts,
                                            cpgrid.globalCell(),
                                            *wells,
                                            grid_and_wells->getWellConnections(),
                                            myExportList, myImportList,
                                            cc);

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

    std::unordered_set<std::string> defunct_well_names;

    if( wells )
    {
        defunct_well_names = computeDefunctWellNames(wells_on_proc,
                                                     *wells,
                                                     cc,
                                                     root);
    }

    return std::make_tuple(parts, defunct_well_names, myExportList, myImportList);
}
}
}
#endif // HAVE_ZOLTAN
