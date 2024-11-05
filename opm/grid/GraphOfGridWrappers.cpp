// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2024 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>
#include "GraphOfGridWrappers.hpp"

namespace Opm {

#if HAVE_MPI
int getGraphOfGridNumVertices(void* pGraph, int *err)
{
    const GraphOfGrid<Dune::CpGrid>&  gog = *static_cast<const GraphOfGrid<Dune::CpGrid>*>(pGraph);
    int size = gog.size();
    *err = ZOLTAN_OK;
    return size;
}

void getGraphOfGridVerticesList(void* pGraph,
               [[maybe_unused]] int dimGlobalID,
               [[maybe_unused]] int dimLocalID,
                                ZOLTAN_ID_PTR gIDs,
               [[maybe_unused]] ZOLTAN_ID_PTR lIDs,
                                int weightDim,
                                float *objWeights,
                                int *err)
{
    assert(dimGlobalID==1); // ID is a single int
    assert(weightDim==1); // vertex weight is a single float
    const GraphOfGrid<Dune::CpGrid>& gog = *static_cast<const GraphOfGrid<Dune::CpGrid>*>(pGraph);
    int i=0;
    for (const auto& v : gog)
    {
        gIDs[i] = v.first;
        // lIDs are left unused
        objWeights[i] = v.second.weight;
        ++i;
    }
    *err = ZOLTAN_OK;
}

void getGraphOfGridNumEdges(void *pGraph,
           [[maybe_unused]] int dimGlobalID,
           [[maybe_unused]] int dimLocalID,
                            int numCells,
                            ZOLTAN_ID_PTR gIDs,
           [[maybe_unused]] ZOLTAN_ID_PTR lIDs,
                            int *numEdges,
                            int *err)
{
    assert(dimGlobalID==1); // ID is a single int
    const GraphOfGrid<Dune::CpGrid>& gog = *static_cast<const GraphOfGrid<Dune::CpGrid>*>(pGraph);
    for (int i=0; i<numCells; ++i)
    {
        int nE = gog.numEdges(gIDs[i]);
        if (nE== -1)
        {
            std::ostringstream ostr;
            ostr << "getGraphOfGridNumEdges error: Vertex with ID " << gIDs[i] << " is not in graph.";
            OpmLog::error(ostr.str());
            *err = ZOLTAN_FATAL;
            return;
        }
        numEdges[i] = nE;
    }
    *err = ZOLTAN_OK;
}

void getGraphOfGridEdgeList(void *pGraph,
           [[maybe_unused]] int dimGlobalID,
           [[maybe_unused]] int dimLocalID,
                            int numCells,
                            ZOLTAN_ID_PTR gIDs,
           [[maybe_unused]] ZOLTAN_ID_PTR lIDs,
                            int *numEdges,
                            ZOLTAN_ID_PTR nborGIDs,
                            int *nborProc,
                            int weightDim,
                            float *edgeWeights,
                            int *err)
{
    assert(dimGlobalID==1); // ID is a single int
    assert(weightDim==1); // edge weight is a single float
    const GraphOfGrid<Dune::CpGrid>&  gog = *static_cast<const GraphOfGrid<Dune::CpGrid>*>(pGraph);
    int id=0;
    for (int i=0; i<numCells; ++i)
    {
        const auto& eList = gog.edgeList(gIDs[i]);
        if ((int)eList.size()!=numEdges[i])
        {
            std::ostringstream ostr;
            ostr << "getGraphOfGridEdgeList error: Edge number disagreement"
                 << " between Zoltan (" << numEdges[i] << ") and Graph ("
                 << eList.size() << ") for vertex with ID " << gIDs[i] << std::endl;
            OpmLog::error(ostr.str());
            *err = ZOLTAN_FATAL;
            return;
        }
        for (const auto& e : eList)
        {
            nborGIDs[id]= e.first;
            nborProc[id]= gog.getVertex(e.first).nproc;
            edgeWeights[id]= e.second;
            ++id;
        }
    }
    *err = ZOLTAN_OK;
}

template<typename Zoltan_Struct>
void setGraphOfGridZoltanGraphFunctions(Zoltan_Struct *zz,
                      const GraphOfGrid<Dune::CpGrid>& gog,
                                                  bool pretendNull)
{
    GraphOfGrid<Dune::CpGrid>* pGraph = const_cast<GraphOfGrid<Dune::CpGrid>*>(&gog);
    if (pretendNull)
    {
        Zoltan_Set_Num_Obj_Fn(zz, Dune::cpgrid::getNullNumCells, pGraph);
        Zoltan_Set_Obj_List_Fn(zz, Dune::cpgrid::getNullVertexList, pGraph);
        Zoltan_Set_Num_Edges_Multi_Fn(zz, Dune::cpgrid::getNullNumEdgesList, pGraph);
        Zoltan_Set_Edge_List_Multi_Fn(zz, Dune::cpgrid::getNullEdgeList, pGraph);
    }
    else
    {
        Zoltan_Set_Num_Obj_Fn(zz, getGraphOfGridNumVertices, pGraph);
        Zoltan_Set_Obj_List_Fn(zz, getGraphOfGridVerticesList, pGraph);
        Zoltan_Set_Num_Edges_Multi_Fn(zz, getGraphOfGridNumEdges, pGraph);
        Zoltan_Set_Edge_List_Multi_Fn(zz, getGraphOfGridEdgeList, pGraph);
    }
}
#endif // HAVE_MPI

void addFutureConnectionWells (GraphOfGrid<Dune::CpGrid>& gog,
    const std::unordered_map<std::string, std::set<int>>& wells,
                                                     bool checkWellIntersections)
{
    // create compressed lookup from cartesian.
    const auto& grid = gog.getGrid();
    const auto& cpgdim = grid.logicalCartesianSize();
    std::vector<int> cartesian_to_compressed(cpgdim[0]*cpgdim[1]*cpgdim[2], -1);
    for( int i=0; i < grid.numCells(); ++i )
    {
        cartesian_to_compressed[grid.globalCell()[i]] = i;
    }

    for (const auto& w: wells)
    {
        std::set<int> wellsgID;
        for (const int& cell : w.second)
        {
            int gID = cartesian_to_compressed.at(cell);
            assert(gID!=-1); // well should be an active cell
            wellsgID.insert(gID);
        }
        gog.addWell(wellsgID,checkWellIntersections);
    }
}

void addWellConnections (GraphOfGrid<Dune::CpGrid>& gog,
               const Dune::cpgrid::WellConnections& wells,
                                               bool checkWellIntersections)
{
    for (const auto& w : wells)
    {
        gog.addWell(w,checkWellIntersections);
    }
}

void extendGIDtoRank (const GraphOfGrid<Dune::CpGrid>& gog,
                                     std::vector<int>& gIDtoRank,
                                            const int& thisRank)
{
    for (const auto& w : gog.getWells())
    {
        auto wellID = *w.begin();
        if (gIDtoRank[wellID]!=thisRank)
        {
            for (const auto& gID : w)
            {
                gIDtoRank.at(gID) = gIDtoRank.at(wellID);
            }
        }
    }
}

template<typename TheTuple>
void extendImportExportList (const GraphOfGrid<Dune::CpGrid>& gog,
                                       std::vector<TheTuple>& cellList,
                                                          int skippedRank,
                                      const std::vector<int>& gIDtoRank)
{
    // using TheTuple = std::tuple<int,int,char>; or std::tuple<int,int,char,int>
    using CellList = std::vector<TheTuple>;
    // make a list of wells for easy identification. Contains ID, begin, end
    using iter = std::set<int>::const_iterator;
    std::unordered_map<int,std::tuple<iter,iter>> wellMap;
    for (const auto& well : gog.getWells())
    {
        if (gIDtoRank.size()>0)
        {
            auto wellID = *well.begin();
            if (gIDtoRank.at(wellID)!=skippedRank)
            {
                wellMap[wellID] = std::make_tuple(well.begin(),well.end());
            }
        }
        else
        {
            wellMap[*well.begin()] = std::make_tuple(well.begin(),well.end());
        }
    }
    if (wellMap.size()==0)
    {
        return;
    }

    CellList addToList;
    // iterate once through the original cellList
    for (const auto& cellProperties : cellList)
    {
        // if a cell is in any well, add cells of the well to cellList
        auto pWell = wellMap.find(std::get<0>(cellProperties));
        if (pWell!=wellMap.end())
        {
            if (std::get<1>(cellProperties)!=skippedRank)
            {
                const auto& [begin,end] = std::pair(std::get<0>(pWell->second),std::get<1>(pWell->second));
                for (auto pgID = begin; pgID!=end; ++pgID)
                {
                    // cells in one well have the same attributes (except ID)
                    if (*pgID!=std::get<0>(cellProperties)) // avoid adding cell that is already in the list
                    {
                        TheTuple wellCell = cellProperties;
                        std::get<0>(wellCell) = *pgID;
                        addToList.push_back(wellCell);
                    }
                }
            }
            wellMap.erase(pWell);
            if (wellMap.empty())
            {
                break;
            }
        }
    }

    // add new cells to the cellList and sort it. It is assumed that cellList starts sorted.
    auto compareTuple = [](const auto& a, const auto& b){return std::get<0>(a)<std::get<0>(b);};
    std::sort(addToList.begin(),addToList.end(),compareTuple);
    auto origSize = cellList.size();
    auto totsize = origSize+addToList.size();
    cellList.reserve(totsize);
    cellList.insert(cellList.end(),addToList.begin(),addToList.end());
    std::inplace_merge(cellList.begin(),cellList.begin()+origSize,cellList.end(),compareTuple);
}

std::vector<int> getWellRanks(const std::vector<int>& gIDtoRank,
                 const Dune::cpgrid::WellConnections& wellConnections)
{
    std::vector<int> wellIndices(wellConnections.size());
    for (std::size_t wellIndex = 0; wellIndex < wellConnections.size(); ++wellIndex)
    {
        int wellID = *(wellConnections[wellIndex].begin());
        wellIndices[wellIndex] = gIDtoRank[wellID];
    }
    return wellIndices;
}

#if HAVE_MPI
std::vector<std::pair<std::string,bool>>
wellsOnThisRank(const std::vector<Dune::cpgrid::OpmWellType>& wells,
                const std::vector<int>& wellRanks,
                const Dune::cpgrid::CpGridDataTraits::Communication& cc,
                int root)
{
    auto numProcs = cc.size();
    std::vector<std::vector<int>> wells_on_proc(numProcs);
    for (std::size_t i=0; i<wellRanks.size(); ++i)
    {
        wells_on_proc.at(wellRanks[i]).push_back(i);
    }
    return Dune::cpgrid::computeParallelWells(wells_on_proc, wells, cc, root);
}

template<class Id>
std::tuple<std::vector<int>,
           std::vector<std::pair<std::string,bool>>,
           std::vector<std::tuple<int,int,char> >,
           std::vector<std::tuple<int,int,char,int> > >
makeImportAndExportLists(const GraphOfGrid<Dune::CpGrid>& gog,
                         const Dune::Communication<MPI_Comm>& cc,
                         const std::vector<Dune::cpgrid::OpmWellType> * wells,
                         const Dune::cpgrid::WellConnections& wellConnections,
                         int root,
                         int numExport,
                         int numImport,
        [[maybe_unused]] const Id* exportLocalGids,
                         const Id* exportGlobalGids,
                         const int* exportToPart,
                         const Id* importGlobalGids)
{
    const auto& cpgrid = gog.getGrid();
    int size = cpgrid.numCells();
    int rank  = cc.rank();
    std::vector<int> gIDtoRank(size, rank);
    std::vector<std::vector<int> > wellsOnProc;

    // List entry: process to export to, (global) index, process rank, attribute there (not needed?)
    std::vector<std::tuple<int,int,char>> myExportList(numExport);
    // List entry: process to import from, global index, process rank, attribute here, local index (determined later)
    std::vector<std::tuple<int,int,char,int>> myImportList(numImport);
    myExportList.reserve(1.2*myExportList.size());
    myImportList.reserve(1.2*myImportList.size());
    using AttributeSet = Dune::cpgrid::CpGridData::AttributeSet;

    for ( int i=0; i < numExport; ++i )
    {
        gIDtoRank[exportGlobalGids[i]] = exportToPart[i];
        myExportList[i] = std::make_tuple(exportGlobalGids[i], exportToPart[i], static_cast<char>(AttributeSet::owner));
    }
    // partitioner sees only one cell per well, modify remaining
    extendGIDtoRank(gog,gIDtoRank,rank);

    for ( int i=0; i < numImport; ++i )
    {
        myImportList[i] = std::make_tuple(importGlobalGids[i], root, static_cast<char>(AttributeSet::owner),-1);
    }


    // Add cells that stay here to the lists. Somehow I could not persuade Zoltan to do this.
    for ( std::size_t i = 0; i < gIDtoRank.size(); ++i)
    {
        if ( gIDtoRank[i] == rank )
        {
            myExportList.emplace_back(i, rank, static_cast<char>(AttributeSet::owner) );
            myImportList.emplace_back(i, rank, static_cast<char>(AttributeSet::owner), -1 );
        }
    }
    std::inplace_merge(myImportList.begin(), myImportList.begin() + numImport, myImportList.end());
    std::inplace_merge(myExportList.begin(), myExportList.begin() + numExport, myExportList.end());

    // Complete the lists by adding cells that were contracted in the graph.
    // All cells on the rank "rank" were added before.
    extendImportExportList(gog,myImportList,rank,gIDtoRank);
    extendImportExportList(gog,myExportList,rank,gIDtoRank);

    std::vector<std::pair<std::string,bool>> parallel_wells;
    if( wells )
    {
        auto wellRanks = getWellRanks(gIDtoRank, wellConnections);
        parallel_wells = wellsOnThisRank(*wells, wellRanks, cc, root);
    }
    return std::make_tuple(gIDtoRank, parallel_wells, myExportList, myImportList);
}

namespace {
void setDefaultZoltanParameters(Zoltan_Struct* zz)
{
    Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
    Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "0");
    Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
    Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM","1");
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1");
    Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", ".35");  /* 0-remove all, 1-remove none */
    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
#ifndef NDEBUG
    Zoltan_Set_Param(zz, "CHECK_GRAPH", "2");
#else
    Zoltan_Set_Param(zz, "CHECK_GRAPH", "0");
#endif
}

} // anon namespace

std::tuple<std::vector<int>, std::vector<std::pair<std::string,bool>>,
           std::vector<std::tuple<int,int,char> >,
           std::vector<std::tuple<int,int,char,int> >,
           Dune::cpgrid::WellConnections>
zoltanPartitioningWithGraphOfGrid(const Dune::CpGrid& grid,
                                  const std::vector<Dune::cpgrid::OpmWellType> * wells,
                                  const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections,
                 [[maybe_unused]] const double* transmissibilities,
                                  const Dune::cpgrid::CpGridDataTraits::Communication& cc,
                 [[maybe_unused]] Dune::EdgeWeightMethod edgeWeightsMethod,
                                  int root,
                                  const double zoltanImbalanceTol,
                                  const std::map<std::string,std::string>& params)
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
    Zoltan_Set_Param(zz, "IMBALANCE_TOL", std::to_string(zoltanImbalanceTol).c_str());
    for (const auto& [key, value] : params)
        Zoltan_Set_Param(zz, key.c_str(), value.c_str());

    // root process has the whole grid, other ranks nothing
    bool partitionIsEmpty = cc.rank()!=root;

    // prepare graph and contract well cells
    GraphOfGrid gog(grid);
    Dune::cpgrid::WellConnections wellConnections(*wells,possibleFutureConnections,grid);
    addWellConnections(gog,wellConnections);

    // call partitioner
    setGraphOfGridZoltanGraphFunctions(zz,gog,partitionIsEmpty);
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

    // arrange output into tuples and add well cells
    auto importExportLists = makeImportAndExportLists(gog,
                                                      cc,
                                                      wells,
                                                      wellConnections,
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

    // add wellConnections to the importExportLists and return it
    auto result = std::tuple(std::move(std::get<0>(importExportLists)),
                             std::move(std::get<1>(importExportLists)),
                             std::move(std::get<2>(importExportLists)),
                             std::move(std::get<3>(importExportLists)),
                             wellConnections);
    return result;
}
#endif // HAVE_MPI

// explicit template instantiations
template void extendImportExportList(const GraphOfGrid<Dune::CpGrid>&,
                                     std::vector<std::tuple<int,int,char> >& ,
                                     int,
                                     const std::vector<int>&);
template void extendImportExportList(const GraphOfGrid<Dune::CpGrid>&,
                                     std::vector<std::tuple<int,int,char,int> >& ,
                                     int,
                                     const std::vector<int>&);

} // end namespace Opm
