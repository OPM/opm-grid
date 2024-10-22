// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

#ifndef GRAPH_OF_GRID_WRAPPERS_HEADER
#define GRAPH_OF_GRID_WRAPPERS_HEADER

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/grid/GraphOfGrid.hpp>

namespace Opm {
/*
  This file contains wrappers for GraphOfGrid that satisfy interface
  requirements of graph partitioners like Zoltan and (TODO!) Metis.

  Additionally, parsing wells is done here.
*/
  #define ZOLTAN_OK 1
  #define ZOLTAN_FATAL 0
  namespace {
    using ZOLTAN_ID_PTR = int*;
  }

  /// \brief callback function for ZOLTAN_NUM_OBJ_FN
  ///
  /// returns the number of vertices in the graph
  int getGraphOfGridNumVertices(void* pGraph, int *err)
  {
    const GraphOfGrid<Dune::CpGrid>&  gog = *static_cast<const GraphOfGrid<Dune::CpGrid>*>(pGraph);
    int size = gog.size();
    *err = ZOLTAN_OK;
    return size;
  }

  /// \brief callback ftion for ZOLTAN_OBJ_LIST_FN
  ///
  /// fills the vector gIDs with vertex global IDs
  ///  and the vector objWeights with their weights
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

  /// \brief callback ftion for ZOLTAN_NUM_EDGES_MULTI_FN
  ///
  /// takes the list of global IDs (gIDs) and fills (consecutively)
  /// vector numEdges with the number of their edges
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

  /// \brief callback ftion for ZOLTAN_EDGE_LIST_MULTI_FN
  ///
  /// takes the list of global IDs (gIDs) and fills (consecutively):
  /// vector nborGIDs with the list of neighbors (all into 1 vector),
  /// vector nborProc with neighbors' process numbers,
  /// vector edgeWeights with edge weights.
  /// The vector numEdges provides the number of edges for each gID
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
      auto eList = gog.edgeList(gIDs[i]);
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

  /// \brief Adds well to the GraphOfGrid
  ///
  /// Translates wells' cartesian ID to global ID used in the graph.
  /// Adding the well contracts vertices of the well into one vertex.
  ///
  /// checkWellIntersections==true makes the algorithm check if wells
  /// intersect and if their cell IDs are present in the graph.
  /// Setting it to false makes the algorithm faster but leaves user
  /// responsible for keeping wells disjoint.
  void addFutureConnectionWells (GraphOfGrid<Dune::CpGrid>& gog,
     const std::unordered_map<std::string, std::set<int>>& wells,
     bool checkWellIntersections=true)
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

  /// \brief Add well cells' global IDs to the list
  ///
  /// Output of the partitioning is missing vertices that were contracted.
  /// This function fills in omitted gIDs and gives them the properties
  /// (like process number and ownership) of their representative cell (well ID).
  template<typename TheTuple>
  void extendImportExportList (const GraphOfGrid<Dune::CpGrid>& gog,
                 std::vector<TheTuple>& cellList)
  {
    // using TheTuple = std::tuple<int,int,char>; or std::tuple<int,int,char,int>
    using CellList = std::vector<TheTuple>;
    // make a list of wells for easy identification. Contains ID, begin, end
    using iter = std::set<int>::const_iterator;
    std::list<std::tuple<int,iter,iter>> wellList;
    for (const auto& well : gog.getWells())
    {
      wellList.push_front(std::make_tuple(*well.begin(),well.begin(),well.end()));
    }

    CellList addToList;
    // iterate once through the original cellList
    for (const auto& cell : cellList)
    {
      // if a cell is in any well, add cells of the well to cellList
      // and remove the well from the wellList to quicken the search
      for (auto wID=wellList.begin(); wID!=wellList.end(); )
      {
        if (std::get<0>(*wID) == std::get<0>(cell))
        {
          for (auto pgID = std::get<1>(*wID); pgID!=std::get<2>(*wID); ++pgID)
          {
            // cells in one well have the same attributes (except ID)
            if (*pgID!=std::get<0>(cell))
            {
              TheTuple wellCell = cell;
              std::get<0>(wellCell) = *pgID;
              addToList.push_back(wellCell);
            }
          }
          wID = wellList.erase(wID);
        }
        else
        {
          ++wID;
        }
      }
      if (wellList.empty())
      {
        break;
      }
    }
    int totsize = cellList.size()+addToList.size();
    cellList.reserve(totsize);
    cellList.insert(cellList.end(),addToList.begin(),addToList.end());
  }

} // end namespace Opm

#endif // GRAPH_OF_GRID_WRAPPERS_HEADER
