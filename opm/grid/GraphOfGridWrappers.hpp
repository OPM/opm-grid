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
  /// \brief callback ftion for ZOLTAN_NUM_OBJ_FN
  /// returns the number of vertices in the graph
  int getGraphOfGridNumVertices(void* pGraph, int *err)
  {
    const GraphOfGrid<Dune::CpGrid>&  gog = *static_cast<const GraphOfGrid<Dune::CpGrid>*>(pGraph);
    int size = gog.size();
    *err = ZOLTAN_OK;
    return size;
  }

  /// \brief callback ftion for ZOLTAN_OBJ_LIST_FN
  /// fills the vector gIDs with vertex global IDs
  ///  and the vector objWeights with their weights
  void getGraphOfGridVerticesList(void* pGraph,
                                  int dimGlobalID,
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
  /// takes the list of global IDs (gIDs) and fills (consecutively)
  /// vector numEdges with the number of their edges
  void getGraphOfGridNumEdges(void *pGraph,
                              int dimGlobalID,
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
  /// takes the list of global IDs (gIDs) and fills (consecutively):
  /// vector nborGIDs with the list of neighbors (all into 1 vector),
  /// vector nborProc with neighbors' process numbers,
  /// vector edgeWeights with edge weights.
  /// The vector numEdges provides the number of edges for each gID
  void getGraphOfGridEdgeList(void *pGraph, int dimGlobalID,
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

  // Wells:
  /// \brief Adds well to the GraphOfGrid
  /// Adding the well contracts vertices of the well into one vertex.
  void addFutureConnectionWells (GraphOfGrid<Dune::CpGrid>& gog,
     const std::unordered_map<std::string, std::set<int>>& wells)
  {
    for (const auto& w : wells)
      gog.addWell(w.second);
  }

} // end namespace Opm

#endif // GRAPH_OF_GRID_WRAPPERS_HEADER