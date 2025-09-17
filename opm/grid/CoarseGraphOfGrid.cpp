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

#include <config.h>
#include "CoarseGraphOfGrid.hpp"

#include <numeric>


namespace Opm {

template<typename Grid>
void CoarseGraphOfGrid<Grid>::dfsq(Row row, std::priority_queue<WgtIdx> &q, int v, int master,
                                   double w, int maxNode, std::vector<bool>& visited,
                                   std::vector<int>& cnode, std::vector<std::tuple<int,int,double> >& edges)
{
    visited[v] = true;
    map_to_coarse_[v] = master;
    cnode.push_back(v);

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

    // Only merge more vertices if the current coarse node is smaller than maxNode.
    if ( (int)cnode.size() < maxNode ) {
        if (!q.empty()) {

            // Find strongest connection in queue q not already merged to cnode.
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
                dfsq((*transGraph)[nab],q,nab,master,w,maxNode,visited,cnode,edges);
        }
    } else {
        q = std::priority_queue<WgtIdx>();
    }

    // Add connection between v and nab if v and nab are not merged
    col = row.begin();
    for (; col != row.end(); ++col) {
        int nab = col.index();
        if (map_to_coarse_[v]!=map_to_coarse_[nab]) {
            edges.push_back({v,nab,(*transGraph)[v][nab]});
        }
    }
}

template<typename Grid>
void CoarseGraphOfGrid<Grid>::mergeWellCellsForCoarseGraph(std::vector<int>& hasWell,
                                                           std::vector<std::vector<int>>& wellPerf,
                                                           const Dune::cpgrid::WellConnections& wellConn)
{
    int wellId = 0;

    for (const auto& well : wellConn) {

        bool cellInMultWells = false;
        std::unordered_set<int> otherWell;
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
            int minOtherWell = * std::min_element(otherWell.begin(),otherWell.end());
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
void CoarseGraphOfGrid<Grid>::dfsqw(Row row, std::priority_queue<WgtIdx> &q, int v, int master,
                              double w, int maxNode, std::vector<bool>& visited,
                              std::vector<int>& cnode, std::vector<std::tuple<int,int,double> >& edges,
                              std::vector<int>& hasWell, std::vector<std::vector<int>>& wellPerf)
{
    std::vector<int> wellIdxs;
    if (hasWell[v] == -1) {
        visited[v] = true;
        map_to_coarse_[v] = master;
        cnode.push_back(v);

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
        // If node has a well, murge all connections to  master node.
        int wellId = hasWell[v];
        std::vector<int> perfs = wellPerf[wellId];

        for (const auto& idx : perfs) {
            visited[idx] = true;
            map_to_coarse_[idx] = master;
            cnode.push_back(idx);
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
    if ( (int)cnode.size() < maxNode ) {
        if (!q.empty()) {

            // Find strongest connection in queue q not already merged to cnode.
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
                dfsqw((*transGraph)[nab],q,nab,master,w,maxNode,visited,cnode,edges,hasWell,wellPerf);
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
                    edges.push_back({v,nab,(*transGraph)[idx][nab]});
                }
            }
        }
    }
    else {
        auto col = row.begin();
        for (; col != row.end(); ++col) {
            int nab = col.index();
            if (map_to_coarse_[v]!=map_to_coarse_[nab]) {
                edges.push_back({v,nab,(*transGraph)[v][nab]});
            }
        }
    }
}

template<typename Grid>
void CoarseGraphOfGrid<Grid>::createCoarseGraph(const Dune::EdgeWeightMethod edgeWeightMethod,
                                                double coarseThreshold,
                                                int coarsePartitionMaxNodeSize,
                                                bool allowDistributedWells,
                                                const Dune::cpgrid::WellConnections& wellConn)
{
    int N = grid.size(0);
    const auto& rank = grid.comm().rank();

    std::vector<bool> visited(N, false);
    std::vector<int> hasWell(N,-1);
    std::vector<std::vector<int>> wellPerf;
    if (!allowDistributedWells) {
        mergeWellCellsForCoarseGraph(hasWell, wellPerf, wellConn);
    }
    map_to_coarse_.resize(N, -1);
    std::vector<int> c2f;

    std::vector<std::vector<std::tuple<int,int,double> >> gEdges;

    int newV = 0;
    int biggest = 0;

    if (rank == 0) {
        for (int v = 0; v < N; ++v) {

            if (!visited[v]) {

                std::priority_queue<WgtIdx> q;
                c2f.push_back(v);
                std::vector<int> cnode;
                std::vector<std::tuple<int,int,double> > edges;

                if (allowDistributedWells) {
                    dfsq((*transGraph)[v],q,v,newV,coarseThreshold,
                         coarsePartitionMaxNodeSize,visited,cnode,edges);
                } else {
                    dfsqw((*transGraph)[v],q,v,newV,coarseThreshold,
                          coarsePartitionMaxNodeSize,visited,cnode,edges,
                          hasWell, wellPerf);
                }

                newV++;
                gEdges.push_back(edges);
                coarseNodes.push_back(cnode);
                if ((int)cnode.size() > biggest)
                    biggest = cnode.size();
            }
        }
        std::cout << "Coarse partitioning graph size: " << coarseNodes.size() <<" Largest node: "<< biggest << std::endl;

        for (std::vector<std::tuple<int,int,double> > es : gEdges ) {
            std::map<int, double> ce;
            for (std::tuple<int,int,double> fe : es) {
                int coarseNab = map_to_coarse_[std::get<1>(fe)];
                double transVal = std::get<2>(fe);
                double weight = edgeWeightMethod == 0 ? 1.0 : transVal;
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
