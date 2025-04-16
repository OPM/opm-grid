// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2024 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#ifndef OPM_GRAPH_OF_GRID_HEADER
#define OPM_GRAPH_OF_GRID_HEADER

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
    using EdgeList = std::map<int,WeightType>;

    struct VertexProperties
    {
        int nproc = 0; // number of processor
        WeightType weight = 1; // vertex weight
        EdgeList edges;
    };

public:
    explicit GraphOfGrid (const Grid& grid_,
                          const double* transmissibilities=nullptr,
                          const Dune::EdgeWeightMethod edgeWeightMethod=Dune::EdgeWeightMethod::defaultTransEdgeWgt,
                          int level = -1)
        : grid(grid_)
    {
        createGraph(transmissibilities,edgeWeightMethod, level);
    }

    const Grid& getGrid() const
    {
        return grid;
    }

    /// \brief Number of graph vertices
    int size () const
    {
        return graph.size();
    }

    auto begin() const
    {
        return graph.begin();
    }
    auto end() const
    {
        return graph.end();
    }

    /// \brief Get iterator to the vertex with this global ID
    /// or ID of the well containing it
    auto find(int gID) const
    {
        // search the graph first, and then wells
        auto pgID = graph.find(gID);
        if (pgID == graph.end())
        {
            gID = wellID(gID);
            if (gID == -1)
            {
                return graph.end();
            }
            pgID = graph.find(gID);
        }
        return pgID;
    }

    /// \brief Return properties of vertex of given ID.
    ///
    /// If no such vertex exists, returns vertex with
    /// process -1, weight 0, and empty edgeList.
    /// If the vertex is in a well, return the well's vertex.
    const VertexProperties& getVertex (int gID) const
    {
        auto pgID = find(gID);
        if (pgID == graph.end())
        {
            OPM_THROW(std::logic_error, "GraphOfGrid::getVertex: gID is not in the graph!");
        }
        return pgID->second;
    }

    /// \brief Number of vertices for given vertex
    ///
    // returns -1 if vertex with such global ID is not in the graph (or wells)
    int numEdges (int gID) const
    {
        auto pgID = find(gID);
        if (pgID == graph.end())
        {
            return -1;
        }
        else
        {
            return pgID->second.edges.size();
        }
    }

    /// \brief List of neighbors for given vertex
    const EdgeList& edgeList(int gID) const
    {
        // get iterator to the vertex or the well containing it
        auto pgID = find(gID);
        if (pgID==graph.end())
        {
            OPM_THROW(std::logic_error, "GraphOfGrid::edgeList: gID is not in the graph!");
        }
        return pgID->second.edges;
    }

    /// \brief Contract two vertices
    ///
    /// Vertex weights are added, and edges are merged. Edge weights
    /// for their common neighbors are added up.
    /// Returns global ID of the resulting vertex, which is smaller ID.
    /// If either gID is in a well, well's ID can be returned if it is smaller.
    int contractVertices (int gID1, int gID2);

    /// \brief Register the well to the list of wells
    ///
    /// If checkIntersection==true, it checks if any of well's cells is
    /// in another well(s) and merges them together.
    /// checkIntersection==false skips those (possibly expensive) checks
    /// but leaves it to user to guarantee that wells are disjoint and
    /// that all cell global IDs are in the graph
    void addWell (const std::set<int>& well, bool checkIntersection=true);

    /// \brief Return the list of wells
    const auto& getWells () const
    {
        return wells;
    }

    /// \brief Contract a layer of verices around each well into it
    ///
    /// Representing a well by one node guarantees that the well won't
    /// be split over several processes. Giving the well an extra layer
    /// of cells distances that well from the subdomain boundary.
    void addNeighboringCellsToWells ();
    void addNeighboringCellsToWells (int layers)
    {
        for (int i=0; i<layers; ++i)
        {
            addNeighboringCellsToWells();
        }
    }

private:
    /// \brief Create a graph representation of the grid
    ///
    /// If transmissibilities are not supplied, edge weight=1
    /// \param level Integer representing the level grid to be considered.
    ///        Default leaf grid view (current_view_data_) set to -1.
    void createGraph (const double* transmissibilities=nullptr,
                      const Dune::EdgeWeightMethod edgeWeightMethod=Dune::EdgeWeightMethod::defaultTransEdgeWgt,
                      int level = -1);

    /// \brief Identify the well containing the cell with this global ID
    ///
    /// returns the smallest cell-ID in the well or
    /// returns -1 if no well contains given gID
    int wellID (int gID) const;

    /// \brief Merge the given well with any overlapping existing wells and contract all vertices.
    ///
    /// This method takes a set of cell indices representing a well and checks if any of them already belong to
    /// a previously registered well in 'wells'. If overlaps are found, the overlapping wells are merged,
    /// and the corresponding vertices are contracted.
    /// The final merged set of vertices is added to the front of 'wells'.
    ///
    /// \param well A set of cell indices representing a well to be contracted and added into 'wells'.
    void mergeWellIndices(const std::set<int>& well);

    /// \brief Contract all vertices in the given well into a single vertex and register the well.
    ///
    /// \param well A set of cell indices representing a well to be contracted and added into 'wells'.
    void contractWellAndAdd(const std::set<int>& well);

    const Grid& grid;
    std::unordered_map<int, VertexProperties> graph; // <gID, VertexProperties>
    std::list<std::set<int>> wells;
};

} // namespace Opm

#endif // OPM_GRAPH_OF_GRID_HEADER
