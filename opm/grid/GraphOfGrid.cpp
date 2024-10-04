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

#include <config.h>

#include <opm/grid/CpGrid.hpp>

namespace Opm {


/// \brief A class storing a graph representation of the grid
///
/// Stores the list of all cell global IDs and for each cell
/// a list of global IDs of its neighbors.
/// In addition, weights of graph vertices and edges are stored.
/// 
/// Features edge contractions, which adds weights of merged vertices
/// and of edges to every shared neighbor. Intended use is for loadbalancing
/// to ensure that no well is split between processes.
template<typename Grid>
class GraphOfGrid{
  using WeightType = float;

  struct VertexProperties
  {
    WeightType weight = 1;
    std::map<int,WeightType> edges; // neighbor's global ID and edge's weight
  };

public:
  GraphOfGrid (const Grid& grid_)
    : grid(grid_)
  {
    createGraph();
  }

  /// \brief Number of graph vertices
  int size () const
  {
    return graph.size();
  }

  /// \brief Number of vertices for given vertex
  // returns -1 if vertex with such global ID is not in the graph
  int numEdges (int gID) const
  {
    auto pgID = graph.find(gID);
    if (pgID == graph.end())
      return -1;
    else
      return pgID->second.edges.size();
  }

private:
  /// \brief Create a graph representation of the grid
  void createGraph ();

  const Grid& grid;
  std::map<int, VertexProperties> graph;
};

  template<typename Grid>
  void GraphOfGrid<Grid>::createGraph ()
  {
    // load vertices (grid cell IDs) into graph
    for (auto it=grid.template leafbegin<0>(); it!=grid.template leafend<0>(); ++it)
    {
      VertexProperties vertex;
      // get vertex's global ID
      int gID = grid.globalIdSet().id(*it); //it->something;

      // iterate over vertex's faces and store neighbors' IDs
      for (int face_lID=0; face_lID<grid.numCellFaces(gID); ++face_lID)
      {
        const int face  = grid.cellFace(gID, face_lID);
        int otherCell   = grid.faceCell(face, 0);
        if ( otherCell == gID || otherCell == -1 ) // -1 means no cell, face is at boundary
        {
            otherCell = grid.faceCell(face, 1);
            if ( otherCell == gID || otherCell == -1 )
                continue;
        }
        WeightType weight = 1;
        vertex.edges.try_emplace(otherCell,weight);
      }

      graph.try_emplace(gID,vertex);
    }

  }
} // namespace Opm