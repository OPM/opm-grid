// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2026 Equinor ASA.

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

#include <config.h>
#include "CoarseGraphOfGrid.hpp"

#include <numeric>

namespace Opm {

template<typename Grid>
void CoarseGraphOfGrid<Grid>::mergeWellCellsForCoarseGraph(std::vector<int>& hasWell,
                                                           std::vector<std::vector<int>>& wellPerf,
                                                           const Dune::cpgrid::WellConnections& wellConn)
{
    int wellId = 0;

    for (const auto& well : wellConn) {

        bool cellInMultWells = false;
        std::set<int> otherWell;
        std::vector<int> perfs;

        // Loop over all connections in the well. Create list perf of connections
        for (int idx : well) {

            if (!cellInMultWells) {
                //Check if connection is intersected by other well
                if (hasWell[idx]!=-1) {
                    cellInMultWells = true;
                    otherWell.insert( hasWell[idx] );
                } else {
                    hasWell[idx] = wellId;
                    perfs.push_back(idx);
                }
            }
            else {
                if (hasWell[idx]!=-1) {
                    otherWell.insert( hasWell[idx] );
                }
            }
        }
        // If current well intersects with other wells. Add the well connections together
        if (cellInMultWells) {
            int minOtherWell = *otherWell.begin();
            for (int idx : well) {
                if (hasWell[idx]!=minOtherWell) {
                    hasWell[idx] = minOtherWell;
                    wellPerf[minOtherWell].push_back(idx);
                }
            }
            // If more than one intersection add all wells to smallest wellId and delete others.
            for (int otherWellId : otherWell) {
                if (otherWellId != minOtherWell) {
                    for (int idx : wellPerf[otherWellId]) {
                        if (hasWell[idx]!=minOtherWell) {
                            hasWell[idx] = minOtherWell;
                            wellPerf[minOtherWell].push_back(idx);
                        }
                    }
                    wellPerf[otherWellId].clear();
                }
            }
        }
        else {
            wellPerf.push_back(perfs);
            wellId++;
        }
    }
}

template<typename Grid>
void CoarseGraphOfGrid<Grid>::dfsqw(const Row& row, std::priority_queue<WgtIdx> &q, int v, int master,
                                    double w, int maxNode, std::vector<bool>& visited,
                                    std::vector<std::vector<std::tuple<int,int,double>>>& gEdges,
                                    const std::vector<int>& hasWell, const std::vector<std::vector<int>>& wellPerf)
{

    auto& current_cnode = coarseNodes.back();
    auto& current_edges = gEdges.back();

    std::vector<int> wellIdxs;
    if (hasWell[v] == -1) {
        visited[v] = true;
        map_to_coarse_[v] = master;
        current_cnode.push_back(v);

        // Add all neighboring vertices of v with transmissibility larger than w to the queue.
        auto col = row.begin();
        for (; col != row.end(); ++col) {
            int nab = col.index();
            double wgt = (*transGraph)[v][nab];
            if ( wgt > w) {
                if (!visited[nab]) {
                    q.push({wgt, nab});
                }
            }
        }
    } else {
        // If node has a well, merge all connections to  master node.
        int wellId = hasWell[v];
        const std::vector<int>& perfs = wellPerf[wellId];

        for (const auto& idx : perfs) {
            visited[idx] = true;
            map_to_coarse_[idx] = master;
            current_cnode.push_back(idx);
            wellIdxs.push_back(idx);
        }
        for (const auto& idx : perfs) {
            auto wrow = (*transGraph)[idx];
            auto col = wrow.begin();
            for (; col != wrow.end(); ++col) {
                int nab = col.index();
                double wgt = (*transGraph)[idx][nab];
                if ( wgt > w) {
                    if (!visited[nab]) {
                        q.push({wgt, nab});
                    }
                }
            }
        }
    }

    // Only merge more vertices if the current coarse node is smaller than maxNode.
    if ( (int)current_cnode.size() < maxNode ) {
        if (!q.empty()) {

            // Find strongest connection in queue q not already merged to current_cnode.
            auto strongCon = q.top();
            int nab = strongCon.idx;
            q.pop();
            while (visited[nab] && !q.empty()) {
                strongCon = q.top();
                nab = strongCon.idx;
                q.pop();
            }
            // Call dfsq reflexively on strongest connection in queue
            if (!visited[nab])
                dfsqw((*transGraph)[nab],q,nab,master,w,maxNode,visited,gEdges,hasWell,wellPerf);
        }
    } else {
        q = std::priority_queue<WgtIdx>();
    }

    // Add connection between v and nab if v and nab are not merged
    if (wellIdxs.size() > 0) {
        for (const auto& idx : wellIdxs) {
            auto wrow = (*transGraph)[idx];
            auto col = wrow.begin();
            for (; col != wrow.end(); ++col) {
                int nab = col.index();
                if (map_to_coarse_[v]!=map_to_coarse_[nab]) {
                    current_edges.push_back({v,nab,(*transGraph)[idx][nab]});
                }
            }
        }
    }
    else {
        auto col = row.begin();
        for (; col != row.end(); ++col) {
            int nab = col.index();
            if (map_to_coarse_[v]!=map_to_coarse_[nab]) {
                current_edges.push_back({v,nab,(*transGraph)[v][nab]});
            }
        }
    }
}

template<typename Grid>
void CoarseGraphOfGrid<Grid>::createCoarseGraph(const Dune::EdgeWeightMethod edgeWeightMethod,
                                                double coarseThreshold,
                                                int coarsePartitionMaxNodeSize,
                                                bool allowDistributedWells,
                                                int root,
                                                const Dune::cpgrid::WellConnections& wellConn)
{
    int N = grid.size(0);
    const auto& rank = grid.comm().rank();

    // List to keep track of visited vertices in original fine graph
    std::vector<bool> visited(N, false);
    // List to know if vertex i is perferated by a well. hasWell[i]==-1 means no, hasWell[i]!=-1 yes
    std::vector<int> hasWell(N,-1);
    // For each well i, wellPerf[i] lists vertices perferated by well i.
    // If wells intersect they are merged.
    std::vector<std::vector<int>> wellPerf;

    // Only add data to hasWell and wellPerf if allowDistributedWells==false
    if (!allowDistributedWells) {
        mergeWellCellsForCoarseGraph(hasWell, wellPerf, wellConn);
    }
    // Map from index of fine graph to index of coarse graph.
    map_to_coarse_.resize(N, -1);

    // Vector that describes the coarse graph.
    // Each coarse vertex v has vector<(int,int,double)> = gEdges[v],
    // where each tuple (idx,nab,wgt) is index idx of fine graph vertex idx contained in v,
    // (idx,nab) is edge in fine graph, and wgt is transmissibility between idx and nab.
    std::vector<std::vector<std::tuple<int,int,double> >> gEdges;

    // Counter for coarse graph index
    int newV = 0;
    // Keep track of largest coarse vertex.
    int biggest = 0;

    if (rank == root) {
        // Loop over all idecies in fine graph
        for (int v = 0; v < N; ++v) {

            // if idx v is not visited, add it to coarse graph
            if (!visited[v]) {

                std::priority_queue<WgtIdx> q;

                // Allocate coarse node v in gEdges and coarseNodes
                gEdges.emplace_back();
                coarseNodes.emplace_back();

                // Call depth first search from row transGraph[v]
                dfsqw((*transGraph)[v],q,v,newV,coarseThreshold,
                      coarsePartitionMaxNodeSize,visited,
                      gEdges,hasWell,wellPerf);

                newV++;

                if ((int)coarseNodes.back().size() > biggest)
                    biggest = coarseNodes.back().size();
            }
        }
        std::cout << "Coarse partitioning graph size: " << coarseNodes.size() <<" Largest node: "<< biggest << std::endl;

        // Construct the coarse graph from gEdges
        //
        // Loop over coarse nodes with cIdx=0,...,newV -1
        for (const std::vector<std::tuple<int,int,double> >& es : gEdges ) {

            // ce represents the edges of cIdx. If cIdx has connection with cNab,
            // wgt=ce[cNab] exists and is greater than 0.
            std::map<int, double> ce;

            // Loop over edges of all vertices from fine graph contained in coarse node cIdx
            for (const std::tuple<int,int,double>& fe : es) {

                // Coarse index of the fine graph edge
                int coarseNab = map_to_coarse_[std::get<1>(fe)];
                // Transmissibility of fine edge
                double transVal = std::get<2>(fe);
                // Besed on edgeWeightMethod, choose weight of coarse edgeWeight
                double weight = edgeWeightMethod == 0 ? 1.0 : transVal;

                // Only add fine edge to coarse graph if transmissibility is non-zero. Overlap cells
                // between partitions are not added if connection has non-zero transmissibility.
                if (transVal > 0) {
                    if ( ce.count(coarseNab) == 1 ) {
                        ce[coarseNab] += weight;
                    } else {
                        ce.insert({coarseNab,weight});
                    }
                }
            }
            cedges.push_back(ce);
        }
    }
}

template class CoarseGraphOfGrid<Dune::CpGrid>;

} // namespace Opm
