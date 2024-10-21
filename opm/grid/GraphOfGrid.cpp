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

#include "GraphOfGrid.hpp"

namespace Opm{

template<typename Grid>
void GraphOfGrid<Grid>::createGraph ()
{
  // load vertices (grid cells) into graph
  for (auto it=grid.template leafbegin<0>(); it!=grid.template leafend<0>(); ++it)
  {
    VertexProperties vertex;
    // get vertex's global ID
    int gID = grid.globalIdSet().id(*it);

    // iterate over vertex's faces and store neighbors' IDs
    for (int face_lID=0; face_lID<grid.numCellFaces(gID); ++face_lID)
    {
      const int face  = grid.cellFace(gID, face_lID);
      int otherCell   = grid.faceCell(face, 0);
      if (otherCell == -1) // -1 means no cell, face is at boundary
        continue;
      if (otherCell == gID)
        otherCell = grid.faceCell(face, 1);
      if (otherCell == -1)
        continue;
      WeightType weight = 1; // default edge weight
      vertex.edges.try_emplace(otherCell,weight);
    }

    graph.try_emplace(gID,vertex);
  }

}

template<typename Grid>
int GraphOfGrid<Grid>::contractVertices (int gID1, int gID2)
{
  // ensure gID1<gID2
  if (gID1==gID2)
      return gID1;
  if (gID2<gID1)
    std::swap(gID1,gID2);

  // check if the gIDs are in the graph or a well
  // do nothing if the vertex is not there
  auto pgID1 = find(gID1);
  auto pgID2 = find(gID2);
  if (pgID1==graph.end() || pgID2==graph.end())
    return -1;

  gID1 = pgID1->first;
  gID2 = pgID2->first;
  // ensure that gID1<gID2
  if (gID1==gID2)
    return gID1;
  if (gID2<gID1)
    std::swap(gID1,gID2);

  // add up vertex weights
  graph[gID1].weight += graph[gID2].weight;

  // Merge the list of neighbors,
  // for common neighbors add up edge weights.
  // Remove the edge between gID1, gID2.
  auto& v1e = graph[gID1].edges;
  v1e.erase(gID2);
  for (const auto& edge : graph[gID2].edges)
  {
    if (v1e.find(edge.first)==v1e.end())
    { // new edge
      if (edge.first != gID1)
      {
        v1e.insert(edge);
        // remap neighbor's edge
        graph[edge.first].edges.erase(gID2);
        graph[edge.first].edges.emplace(gID1,edge.second);
      }
    }
    else
    { // common neighbor, add edge weight
      v1e[edge.first] += edge.second;
      graph[edge.first].edges.erase(gID2);
      graph[edge.first].edges[gID1] += edge.second;
    }
  }

  // erase the second vertex to conclude contraction
  graph.erase(gID2);
  return gID1;
}

template<typename Grid>
int GraphOfGrid<Grid>::wellID (int gID) const
{
  for (const auto& w : wells)
  {
    auto pgID = w.find(gID);
    if (pgID!=w.end())
      return *(w.begin()); // well entries are ordered
  }
  return -1;
}

template<typename Grid>
void GraphOfGrid<Grid>::addWell (const std::set<int>& well, bool checkIntersection)
{
  if (well.size()<2)
    return;
  int wID = *(well.begin());

  if (checkIntersection)
  {
    std::set<int> newWell;
    for (int gID : well)
    {
      // check if the cell is already in some well
      if (newWell.find(gID)!=newWell.end())
        continue;
      for (auto w=wells.begin(); w!=wells.end(); ++w)
      {
        if (w->find(gID)!=w->end())
        {
          // gID is in another well => remap it and join wells
          if (wID==gID)
            wID = *(w->begin());
          gID = *(w->begin());
          newWell.insert(w->begin(),w->end());
          wells.erase(w);
          break; // wells are assumed disjoint, each gID has max 1 match
        }
      }
      wID = contractVertices(wID,gID);
      assert(wID!=-1 && "Added well vertex was not found in the grid (or its wells).");
    }
    newWell.insert(well.begin(),well.end());
    wells.push_front(newWell);
  }
  else
  {
    for (int gID : well)
      wID = contractVertices(wID,gID);
    wells.emplace_front(well);
  }
}

template class GraphOfGrid<Dune::CpGrid>;
} // namespace Opm
