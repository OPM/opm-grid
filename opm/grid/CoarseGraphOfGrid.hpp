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

#ifndef OPM_COARSE_GRAPH_OF_GRID_HEADER
#define OPM_COARSE_GRAPH_OF_GRID_HEADER

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/WellConnections.hpp>
#include <queue>

namespace Opm {

struct WgtIdx {

    double wgt;
    int idx;

    bool operator<(const WgtIdx& other) const {
    return wgt < other.wgt;
    }
};

/// \brief A class storing a Coarse graph representation of the grid
///
/// Similar to GraphOfGrid, but here nodes are merged if
/// the transmissibility on the edge connecting them is below
/// a certain threshold.
template<typename Grid>
class CoarseGraphOfGrid{

    using TransGraph = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    using Row = typename TransGraph::row_type;

public:

    explicit CoarseGraphOfGrid (const Grid& grid_,
                                const Dune::EdgeWeightMethod edgeWeightMethod,
                                TransGraph* tg,
                                double coarseThreshold,
                                int coarsePartitionMaxNodeSize,
                                bool allowDistributedWells,
                                const Dune::cpgrid::WellConnections& wellConn)
        : grid(grid_), transGraph(tg)
    {
        createCoarseGraph(edgeWeightMethod, coarseThreshold, coarsePartitionMaxNodeSize,
                          allowDistributedWells, wellConn);
    }

    const Grid& getGrid() const
    {
        return grid;
    }

    /// \brief Number of graph vertices
    int cSize() const
    {
        return coarseNodes.size();
    }

    std::vector<std::vector<int>> getCoarseNodes() const
    {
        return coarseNodes;
    }

    std::vector<std::map<int, double> > getCoarseEdges() const
    {
        return cedges;
    }

    /// \brief returns map of global cell id to coarse vertex id
    std::vector<int> getMapToCoarse() const
    {
        return map_to_coarse_;
    }

    /// \brief Return the list of wells
    const auto& getWells () const
    {
        return wells;
    }
private:


    /// \brief Coarsen the graph by doing a (d)epth (f)irst (s)earch with priority (q)ueue.
    ///
    /// Recurseive dfs that merges strongly connected nodes (transmissibility) in the
    /// partitioning graph.
    /// A priority queue is used to make sure the largest connections are prioritised in the
    /// merging of vertices.
    void dfsq(Row row, std::priority_queue<WgtIdx> &q, int v, int master,
              double w, int maxNode, std::vector<bool>& visited,
              std::vector<int>& cnode, std::vector<std::tuple<int,int,double> >& edges);


    /// \brief Merge vertices that share a common well 
    void mergeWellCellsForCoarseGraph(std::vector<int>& hasWell,
                                      std::vector<std::vector<int>>& wellPerf,
                                      const Dune::cpgrid::WellConnections& wells);

    /// \brief Coarsen the graph by doing a (d)epth (f)irst (s)earch with priority (q)ueue
    /// and merges (w)ells.
    ///
    /// Similar to dfsq, but merges wells before strongly connected nodes.
    void dfsqw(Row row, std::priority_queue<WgtIdx> &q, int v, int master,
               double w, int maxNode, std::vector<bool>& visited,
               std::vector<int>& cnode, std::vector<std::tuple<int,int,double> >& edges,
               std::vector<int>& hasWell, std::vector<std::vector<int>>& wellPerf);

    /// \brief Create the coarse graph merging all vertices connected with a large transmissibility.
    ///
    /// Creates a coarsened graph by merging vertices connected with a
    /// transmissibility larger than coarseThreshold.
    /// Coarse vertices will not be larger than coarsePartitionMaxNodeSize.
    /// Also meges wells if allowDistributedWells == true.
    void createCoarseGraph(const Dune::EdgeWeightMethod edgeWeightMethod,
                           double coarseThreshold,
                           int coarsePartitionMaxNodeSize,
                           bool allowDistributedWells,
                           const Dune::cpgrid::WellConnections& wells);

    const Grid& grid;
    std::list<std::set<int>> wells;

    Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>* transGraph;
    std::vector<int> map_to_coarse_;
    std::vector<std::map<int, double> > cedges;
    std::vector<std::vector<int>> coarseNodes;
    
};

} // namespace Opm

#endif // OPM_COARSE_GRAPH_OF_GRID_HEADER
