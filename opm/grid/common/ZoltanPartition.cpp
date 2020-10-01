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
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/Entity.hpp>
#include <algorithm>

#if defined(HAVE_ZOLTAN) && defined(HAVE_MPI)
namespace Dune
{
namespace cpgrid
{

namespace {
void setDefaultZoltanParameters(Zoltan_Struct* zz) {
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
}


std::tuple<std::vector<int>, std::vector<std::pair<std::string,bool>>,
           std::vector<std::tuple<int,int,char> >,
           std::vector<std::tuple<int,int,char,int> > >
makeImportAndExportLists(const CpGrid& cpgrid,
                         const CollectiveCommunication<MPI_Comm>& cc,
                         const std::vector<OpmWellType> * wells,
                         const CombinedGridWellGraph& gridAndWells,
                         int root,
                         int numExport,
                         int numImport,
                         const ZOLTAN_ID_PTR exportLocalGids,
                         const ZOLTAN_ID_PTR exportGlobalGids,
                         const int* exportToPart,
                         const ZOLTAN_ID_PTR importGlobalGids) {
    int                         size = cpgrid.numCells();
    int                         rank  = cc.rank();
    std::vector<int>            parts(size, rank);
    std::vector<std::vector<int> > wellsOnProc;

    // List entry: process to export to, (global) index, process rank, attribute there (not needed?)
    std::vector<std::tuple<int,int,char>> myExportList(numExport);
    // List entry: process to import from, global index, process rank, attribute here, local index
    // (determined later)
    std::vector<std::tuple<int,int,char,int>> myImportList(numImport);
    myExportList.reserve(1.2*myExportList.size());
    myImportList.reserve(1.2*myImportList.size());
    using AttributeSet = CpGridData::AttributeSet;

    for ( int i=0; i < numExport; ++i )
    {
        parts[exportLocalGids[i]] = exportToPart[i];
        myExportList[i] = std::make_tuple(exportGlobalGids[i], exportToPart[i], static_cast<char>(AttributeSet::owner));
    }

    for ( int i=0; i < numImport; ++i )
    {
        myImportList[i] = std::make_tuple(importGlobalGids[i], root, static_cast<char>(AttributeSet::owner),-1);
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




    if( wells )
    {
        auto gidGetter = [&cpgrid](int i) { return cpgrid.globalIdSet().id(createEntity<0>(cpgrid, i, true));};
        wellsOnProc =
            postProcessPartitioningForWells(parts,
                                            gidGetter,
                                            *wells,
                                            gridAndWells.getWellConnections(),
                                            myExportList, myImportList,
                                            cc);


#ifndef NDEBUG
        int index = 0;
        for( auto well : gridAndWells.getWellsGraph() )
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

    std::vector<std::pair<std::string,bool>> parallel_wells;
    if( wells )
    {
        parallel_wells = computeParallelWells(wellsOnProc,
                                              *wells,
                                              cc,
                                              root);
    }
    return std::make_tuple(parts, parallel_wells, myExportList, myImportList);
}
} // anon namespace

std::tuple<std::vector<int>, std::vector<std::pair<std::string,bool>>,
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
    setDefaultZoltanParameters(zz);

    // For the load balancer one process has the whole grid and
    // all others an empty partition before loadbalancing.
    bool partitionIsEmpty     = cc.rank()!=root;

    std::shared_ptr<CombinedGridWellGraph> gridAndWells;

    if( wells )
    {
        Zoltan_Set_Param(zz,"EDGE_WEIGHT_DIM","1");
        gridAndWells.reset(new CombinedGridWellGraph(cpgrid,
                                                       wells,
                                                       transmissibilities,
                                                       partitionIsEmpty,
                                                       edgeWeightsMethod));
        Dune::cpgrid::setCpGridZoltanGraphFunctions(zz, *gridAndWells,
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

    auto importExportLists = makeImportAndExportLists(cpgrid,
                                     cc,
                                     wells,
                                     *gridAndWells,
                                     root,
                                     numExport,
                                     numImport,
                                     exportLocalGids,
                                     exportGlobalGids,
                                     exportProcs,
                                     importGlobalGids);
    Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);
    Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, &importProcs, &importToPart);
    Zoltan_Destroy(&zz);

    return importExportLists;
}



class ZoltanSerialPartitioner
{
public:
    ZoltanSerialPartitioner(const CpGrid& _cpgrid,
                            const std::vector<OpmWellType>* _wells,
                            const double* _transmissibilities,
                            const CollectiveCommunication<MPI_Comm>& _cc,
                            EdgeWeightMethod _edgeWeightsMethod,
                            int _root)
        : cpgrid(_cpgrid)
        , wells(_wells)
        , transmissibilities(_transmissibilities)
        , cc(_cc)
        , edgeWeightsMethod(_edgeWeightsMethod)
        , root(_root)
    {
        if (wells) {
            const bool partitionIsEmpty = cc.rank() != root;
            gridAndWells.reset(
                new CombinedGridWellGraph(cpgrid, wells, transmissibilities, partitionIsEmpty, edgeWeightsMethod));
        }
    }

    std::tuple<std::vector<int>,
               std::vector<std::pair<std::string, bool>>,
               std::vector<std::tuple<int, int, char>>,
               std::vector<std::tuple<int, int, char, int>>>
    partition()
    {
        MPI_Barrier(cc);

        // Initialize Zoltan and perform partitioning.
        int rc = ZOLTAN_OK;
        if (cc.rank() == root) {
            rc = callZoltan();
            cc.broadcast<int>(&rc, 1, root);
        } else {
            cc.broadcast<int>(&rc, 1, root);
        }
        if (rc != ZOLTAN_OK) {
            OPM_THROW(std::runtime_error, "Could not initialize Zoltan, or Zoltan partitioning failed.");
        }

        // Build and communicate import/export data.
        // 1. Send number of exports/imports.
        if (cc.rank() == root) {
            numberOfExportedVerticesPerProcess.resize(cc.size(), 0);
            for (int i = 0; i < numExport; ++i) {
                ++numberOfExportedVerticesPerProcess[exportToPart[i]];
            }
            int dummyForRoot = 0;
            cc.scatter<int>(numberOfExportedVerticesPerProcess.data(), &dummyForRoot, 1, root);
        } else {
            cc.scatter<int>(nullptr, &numImport, 1, root);
        }

        // 2. Build the imports/exports themselves.
        std::string error;
        if (cc.rank() == root) {
            offsets.resize(cc.size() + 1, 0);
            std::partial_sum(numberOfExportedVerticesPerProcess.begin(),
                             numberOfExportedVerticesPerProcess.end(),
                             offsets.begin() + 1);
            globalIndicesToSend.resize(numExport, 0);
            const int commSize = cc.size();
            std::vector<int> currentIndex(commSize, 0);
            for (int i = 0; i < numExport; ++i) {
                if (exportToPart[i] >= commSize) {
                    std::ostringstream oss;
                    oss << "Something wrong with Zoltan decomposition. "
                        << "Debug information: exportToPart[i] = " << exportToPart[i] << ", "
                        << "currentIndex.size() = " << currentIndex.size();
                    error = oss.str();
                    break;
                }
                const auto index = currentIndex[exportToPart[i]]++ + offsets[exportToPart[i]];
                if (index >= numExport) {
                    std::ostringstream oss;
                    oss << "Something wrong with Zoltan decomposition. "
                        << "index " << index << ", "
                        << "globalIndicesToSend.size() = " << globalIndicesToSend.size()
                        << "\n\noffsets[exportToPart[i]] = " << offsets[exportToPart[i]]
                        << "\n\ncurrentIndex[exportToPart[i]] = " << currentIndex[exportToPart[i]]
                        << "\n\nexportToPart[i] = " << exportToPart[i];
                    error = oss.str();
                    break;
                }
                globalIndicesToSend[index] = exportGlobalGids[i];
            }
        } else {
            importGlobalGidsVector.resize(numImport, 0);
        }
        // Check for errors
        int ok = error.empty();
        if (cc.rank() == root) {
            cc.broadcast<int>(&ok, 1, root);
        } else {
            cc.broadcast<int>(&ok, 1, root);
        }
        if (!ok) {
            OPM_THROW(std::runtime_error, error);
        }

        // 3. Communicate the imports/exports.
        if (cc.rank() == root) {
            std::vector<unsigned int> dummyIndicesForRoot(1, 0);
            cc.scatterv<unsigned int>(globalIndicesToSend.data(),
                                      numberOfExportedVerticesPerProcess.data(),
                                      offsets.data(),
                                      dummyIndicesForRoot.data(),
                                      0,
                                      root);
        } else {
            cc.scatterv<unsigned int>(nullptr, nullptr, nullptr, importGlobalGidsVector.data(), numImport, root);
            importGlobalGids = importGlobalGidsVector.data();
        }

        auto importExportLists = makeImportAndExportLists(cpgrid,
                                                          cc,
                                                          wells,
                                                          *gridAndWells,
                                                          root,
                                                          numExport,
                                                          numImport,
                                                          exportLocalGids,
                                                          exportGlobalGids,
                                                          exportToPart,
                                                          importGlobalGids);

        return importExportLists;
    }

    ~ZoltanSerialPartitioner()
    {
        // free space allocated for zoltan.
        if (cc.rank() == root) {
            Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);
            Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, &importProcs, &importToPart);
            Zoltan_Destroy(&zz);
        }
    }

private:
    // Methods

    int callZoltan()
    {
        int argc = 0;
        char** argv = 0;
        float ver = 0;

        int rc = Zoltan_Initialize(argc, argv, &ver);
        zz = Zoltan_Create(MPI_COMM_SELF);
        if (rc != ZOLTAN_OK) {
            return rc;
        }

        setDefaultZoltanParameters(zz);
        Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", std::to_string(cc.size()).c_str());

        // For the load balancer one process has the whole grid and
        // all others an empty partition before loadbalancing.
        bool partitionIsEmpty = cc.rank() != root;

        if (wells) {
            Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "1");
            Dune::cpgrid::setCpGridZoltanGraphFunctions(zz, *gridAndWells, partitionIsEmpty);
        } else {
            Dune::cpgrid::setCpGridZoltanGraphFunctions(zz, cpgrid, partitionIsEmpty);
        }

        rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
                                 &changes, /* 1 if partitioning was changed, 0 otherwise */
                                 &numGidEntries, /* Number of integers used for a global ID */
                                 &numLidEntries, /* Number of integers used for a local ID */
                                 &numImport, /* Number of vertices to be sent to me */
                                 &importGlobalGids, /* Global IDs of vertices to be sent to me */
                                 &importLocalGids, /* Local IDs of vertices to be sent to me */
                                 &importProcs, /* Process rank for source of each incoming vertex */
                                 &importToPart, /* New partition for each incoming vertex */
                                 &numExport, /* Number of vertices I must send to other processes*/
                                 &exportGlobalGids, /* Global IDs of the vertices I must send */
                                 &exportLocalGids, /* Local IDs of the vertices I must send */
                                 &exportProcs, /* Process to which I send each of the vertices */
                                 &exportToPart); /* Partition to which each vertex will belong */

        // This is very important: by default, zoltan sets numImport to -1,
        // as we are only running Zoltan in serial.
        // In order to make sense of the rest of  the code, this must be set
        // to 0.
        numImport = 0;
        return rc;
    }

    // Data members
    const CpGrid& cpgrid;
    const std::vector<OpmWellType>* wells;
    const double* transmissibilities;
    const CollectiveCommunication<MPI_Comm>& cc;
    EdgeWeightMethod edgeWeightsMethod;
    int root;
    std::string errorOnRoot;

    struct Zoltan_Struct* zz = nullptr;
    int changes = 0;
    int numGidEntries = 0;
    int numLidEntries = 0;
    int numImport = 0;
    int numExport = 0;
    ZOLTAN_ID_PTR importGlobalGids = nullptr;
    ZOLTAN_ID_PTR importLocalGids = nullptr;
    ZOLTAN_ID_PTR exportGlobalGids = nullptr;
    ZOLTAN_ID_PTR exportLocalGids = nullptr;
    int *importProcs, *importToPart, *exportProcs, *exportToPart;
    std::unique_ptr<CombinedGridWellGraph> gridAndWells;
    std::vector<unsigned int> importGlobalGidsVector;
    std::vector<int> numberOfExportedVerticesPerProcess;
    std::vector<unsigned int> globalIndicesToSend;
    std::vector<int> offsets;
};


std::tuple<std::vector<int>,
           std::vector<std::pair<std::string, bool>>,
           std::vector<std::tuple<int, int, char>>,
           std::vector<std::tuple<int, int, char, int>>>
zoltanSerialGraphPartitionGridOnRoot(const CpGrid& cpgrid,
                                     const std::vector<OpmWellType>* wells,
                                     const double* transmissibilities,
                                     const CollectiveCommunication<MPI_Comm>& cc,
                                     EdgeWeightMethod edgeWeightsMethod,
                                     int root)
{
    ZoltanSerialPartitioner partitioner(cpgrid, wells, transmissibilities, cc, edgeWeightsMethod, root);
    return partitioner.partition();
}




} // namespace cpgrid
} // namespace Dune
#endif // HAVE_ZOLTAN
