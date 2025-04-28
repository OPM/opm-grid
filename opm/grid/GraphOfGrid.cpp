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
#include "GraphOfGrid.hpp"

#include <numeric>

namespace Opm {

template<typename Grid>
void GraphOfGrid<Grid>::createGraph (const double* transmissibilities,
                                     const Dune::EdgeWeightMethod edgeWeightMethod,
                                     [[maybe_unused]] int level)
{
    // Find the lowest positive transmissibility in the grid.
    // This includes boundary faces, even though they will not appear in the graph.
    WeightType logMinTransm = std::numeric_limits<WeightType>::max();
    if (transmissibilities && edgeWeightMethod==Dune::EdgeWeightMethod::logTransEdgeWgt)
    {
        for (int face = 0; face < grid.numFaces(); ++face)
        {
            WeightType transm = transmissibilities[face];
            if (transm > 0 && transm < logMinTransm)
            {
                logMinTransm = transm;
            }
        }
        if (logMinTransm == std::numeric_limits<WeightType>::max()) {
            OPM_THROW(std::domain_error, "All transmissibilities are negative, zero, or bigger than the limit of the WeightType.");
        }
        logMinTransm = std::log(logMinTransm);
    }

    const auto& rank = grid.comm().rank();
    // load vertices (grid cells) into graph
    graph.reserve(grid.size(0));
    for (auto it=grid.template leafbegin<0>(); it!=grid.template leafend<0>(); ++it)
    {
        VertexProperties vertex;
        vertex.nproc = rank;
        // get vertex's global ID
        int gID = grid.globalIdSet().id(*it);

        // iterate over vertex's faces and store neighbors' IDs
        for (int face_lID=0; face_lID<grid.numCellFaces(gID); ++face_lID)
        {
            const int face  = grid.cellFace(gID, face_lID);
            int otherCell   = grid.faceCell(face, 0);
            if (otherCell == -1) // -1 means no cell, face is at boundary
            {
                continue;
            }
            if (otherCell == gID)
            {
                otherCell = grid.faceCell(face, 1);
            }
            if (otherCell == -1)
            {
                continue;
            }
            WeightType weight;
            if (transmissibilities) {
                switch (edgeWeightMethod) {
                case 0:
                    weight = 1.;
                    break;
                case 1:
                    weight = transmissibilities[face];
                    break;
                case 2:
                    weight = 1 + std::log(transmissibilities[face]) - logMinTransm;
                    break;
                default:
                    OPM_THROW(std::invalid_argument, "GraphOfGrid recognizes only EdgeWeightMethod of value 0, 1, or 2.");
                }
            } else {
                weight = 1.;
            }
            vertex.edges.try_emplace(otherCell, weight);
        }

        graph.try_emplace(gID, vertex);
    }

}

// CpGrid Specialization
template<>
void GraphOfGrid<Dune::CpGrid>::createGraph (const double* transmissibilities,
                                             const Dune::EdgeWeightMethod edgeWeightMethod,
                                             int level)
{
    // Find the lowest positive transmissibility in the grid.
    // This includes boundary faces, even though they will not appear in the graph.
    WeightType logMinTransm = std::numeric_limits<WeightType>::max();
    if (transmissibilities && edgeWeightMethod==Dune::EdgeWeightMethod::logTransEdgeWgt)
    {
        for (int face = 0; face < grid.numFaces(level); ++face)
        {
            WeightType transm = transmissibilities[face];
            if (transm > 0 && transm < logMinTransm)
            {
                logMinTransm = transm;
            }
        }
        if (logMinTransm == std::numeric_limits<WeightType>::max()) {
            OPM_THROW(std::domain_error, "All transmissibilities are negative, zero, or bigger than the limit of the WeightType.");
        }
        logMinTransm = std::log(logMinTransm);
    }

    const auto& rank = grid.comm().rank();
    // load vertices (grid cells) into graph
    graph.reserve(grid.numCells(level));

    // Select data according to level/leaf grid to be distributed
    bool validLevel = (level>-1) && (level <= grid.maxLevel());
    auto it = validLevel?  grid.template lbegin<0>(level) :  grid.template leafbegin<0>();
    auto itEnd = validLevel? grid.template lend<0>(level) : grid.template leafend<0>();
            
    for (; it!=itEnd; ++it)
    {
        VertexProperties vertex;
        vertex.nproc = rank;
        // get vertex's global ID
        int gID = validLevel? grid.currentData()[level]->globalIdSet().id(*it) : grid.globalIdSet().id(*it);

        // iterate over vertex's faces and store neighbors' IDs
       
        // To access cell_to_face_, face_to_cell_ local ID is needed
        for (int face_lID=0; face_lID<grid.numCellFaces(it->index(), level); ++face_lID)
        {
            const int face  = grid.cellFace(it->index(), face_lID, level);
            int otherCell   = grid.faceCell(face, 0, level);
            if (otherCell == -1) // -1 means no cell, face is at boundary
            {
                continue;
            }
            if (otherCell == it->index())
            {
                otherCell = grid.faceCell(face, 1, level);
            }
            if (otherCell == -1)
            {
                continue;
            }
            WeightType weight;
            if (transmissibilities) {
                switch (edgeWeightMethod) {
                case 0:
                    weight = 1.;
                    break;
                case 1:
                    weight = transmissibilities[face];
                    break;
                case 2:
                    weight = 1 + std::log(transmissibilities[face]) - logMinTransm;
                    break;
                default:
                    OPM_THROW(std::invalid_argument, "GraphOfGrid recognizes only EdgeWeightMethod of value 0, 1, or 2.");
                }
            } else {
                weight = 1.;
            }
            // lookup otherCell gID
            const auto& otherCellElem = validLevel? Dune::cpgrid::Entity<0>(*(grid.currentData()[level]), otherCell, true) :
                Dune::cpgrid::Entity<0>(*(grid.currentData().back()), otherCell, true);
            int otherCellgID = validLevel? grid.currentData()[level]->globalIdSet().id(otherCellElem) : grid.globalIdSet().id(otherCellElem);
            vertex.edges.try_emplace(otherCellgID /*otherCell*/, weight);
        }

        graph.try_emplace(gID /*it->index()*/, vertex);
    }

}

template<typename Grid>
int GraphOfGrid<Grid>::contractVertices (int gID1, int gID2)
{
    // check if the gIDs are in the graph or a well
    // do nothing if the vertex is not there
    auto pgID1 = find(gID1);
    auto pgID2 = find(gID2);
    if (pgID1==graph.end() || pgID2==graph.end())
    {
        return -1;
    }

    gID1 = pgID1->first;
    gID2 = pgID2->first;
    // ensure that gID1<gID2
    if (gID1==gID2)
    {
        return gID1;
    }
    if (gID2<gID1)
    {
        std::swap(gID1, gID2);
    }

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
                graph[edge.first].edges.emplace(gID1, edge.second);
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
        {
            return *(w.begin()); // well entries are ordered
        }
    }
    return -1;
}

template<typename Grid>
void GraphOfGrid<Grid>::mergeWellIndices(const std::set<int>& well)
{
    if (well.empty())
        return;
    
    int wellIdx = *(well.begin());
    std::set<int> newWell;
    for (int idx : well)
    {
        // check if the cell is already in some well
        if (newWell.find(idx) != newWell.end()) {
            continue;
        }
        for (auto w = wells.begin(); w != wells.end();) {
            if (w->find(idx) != w->end()) {
                // idx is in another well => remap it and join wells
                if ( wellIdx == idx ) {
                    wellIdx = *(w->begin());
                }
                idx = *(w->begin());
                newWell.insert(w->begin(), w->end());
                w = wells.erase(w);
                break; // GraphOfGrid::wells are constructed to be disjoint, each idx has max 1 match
            }
            else {
                ++w;
            }
        }
        wellIdx = contractVertices(wellIdx, idx);
        assert( wellIdx!=-1 && "Added well vertex was not found in the grid (or its wells).");
    }
    newWell.insert(well.begin(), well.end());
    wells.push_front(newWell);
}

template<typename Grid>
void GraphOfGrid<Grid>::contractWellAndAdd(const std::set<int>& well)
{
    if (well.empty())
        return;
    
    int wID = *(well.begin());
    std::accumulate(well.begin(), well.end(), wID,
                    [this](const auto wId, const auto gID)
                    { return contractVertices(wId, gID); });
    wells.emplace_front(well);
}


template<typename Grid>
void GraphOfGrid<Grid>::addWell (const std::set<int>& well, bool checkIntersection)
{
    if (well.size()<2)
        return;

    if (checkIntersection) {
        mergeWellIndices(well);
    }
    else {
        contractWellAndAdd(well);
    }
}

template<typename Grid>
void GraphOfGrid<Grid>::addNeighboringCellsToWells ()
{
    // mark all cells that will be added to wells (addding them one
    // by one would require recursive checks for neighboring wells)
    std::vector<std::set<int>> buffer(wells.size());
    int i=0;
    for (auto& w : wells)
    {
        buffer[i].insert(*w.begin()); // intersects with its well
        for (const auto& v : edgeList(*w.begin()))
        {
            buffer[i].insert(v.first);
        }
        ++i;
    }

    // intersecting wells and buffers will be merged
    for (const auto& b : buffer)
    {
        addWell(b);
    }
}

template class GraphOfGrid<Dune::CpGrid>;

} // namespace Opm
