/*
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services.
  Copyright 2015 NTNU
  Copyright 2015 Statoil AS

  This file is part of The Open Porous Media project  (OPM).

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
*/
#ifndef DUNE_CPGRID_ZOLTAN_GRAPH_FUNCTIONS_HEADER
#define DUNE_CPGRID_ZOLTAN_GRAPH_FUNCTIONS_HEADER

#include <opm/grid/utility/OpmParserIncludes.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/WellConnections.hpp>

#if defined(HAVE_ZOLTAN) && defined(HAVE_MPI)

#include <mpi.h>

// Zoltan redefines HAVE_MPI. Therefore we need to back it up, undef, and
// redifine it after the header is included
#undef HAVE_MPI
#include <zoltan.h>
#undef HAVE_MPI
#define HAVE_MPI 1

namespace Dune
{
namespace cpgrid
{
/// \brief Get the number of cells of the grid.
///
/// The cells are the vertices of the graph.
/// \return The number of vertices of the graph representing the grid.
inline int getCpGridNumCells(void* cpGridPointer, int* err)
{
    const Dune::CpGrid&  grid = *static_cast<const Dune::CpGrid*>(cpGridPointer);
    *err = ZOLTAN_OK;
    return grid.numCells();
}

/// \brief Get the list of vertices of the graph of the grid.
void getCpGridVertexList(void* cpGridPointer, int numGlobalIds,
                         int numLocalIds, ZOLTAN_ID_PTR gids,
                         ZOLTAN_ID_PTR lids, int wgtDim,
                         float *objWgts, int *err);

/// \brief Get the number of edges the graph of the grid.
void getCpGridNumEdgesList(void *cpGridPointer, int sizeGID, int sizeLID,
                           int numCells,
                           ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                           int *numEdges, int *err);

/// \brief Get the list of edges of the graph of the grid.
void getCpGridEdgeList(void *cpGridPointer, int sizeGID, int sizeLID,
                       int numCells, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                       int *num_edges,
                       ZOLTAN_ID_PTR nborGID, int *nborProc,
                       int wgt_dim, float *ewgts, int *err);

/// \brief Get a list of vertices with zero enties
void getNullVertexList(void* cpGridPointer, int numGlobalIds,
                         int numLocalIds, ZOLTAN_ID_PTR gids,
                         ZOLTAN_ID_PTR lids, int wgtDim,
                         float *objWgts, int *err);

/// \brief Get zero as the number of edges the graph of the grid.
void getNullNumEdgesList(void *cpGridPointer, int sizeGID, int sizeLID,
                           int numCells,
                           ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                           int *numEdges, int *err);

/// \brief Get a list of edges of size zero.
void getNullEdgeList(void *cpGridPointer, int sizeGID, int sizeLID,
                       int numCells, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                       int *num_edges,
                       ZOLTAN_ID_PTR nborGID, int *nborProc,
                       int wgt_dim, float *ewgts, int *err);

/// \brief Get always zero as the number of cells of the grid.
///
/// The cells are the vertices of the graph.
inline int getNullNumCells(void* cpGridPointer, int* err)
{
    (void) cpGridPointer;
    *err = ZOLTAN_OK;
    return 0;
}

/// \brief Get the number of edges the graph of the grid and the wells.
void getCpGridWellsNumEdgesList(void *cpGridWellsPointer, int sizeGID, int sizeLID,
                           int numCells,
                           ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                           int *numEdges, int *err);

/// \brief Get the list of edges of the graph of the grid and the wells
void getCpGridWellsEdgeList(void *cpGridWellsPointer, int sizeGID, int sizeLID,
                       int numCells, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                       int *num_edges,
                       ZOLTAN_ID_PTR nborGID, int *nborProc,
                       int wgt_dim, float *ewgts, int *err);

/// \brief A graph repesenting a grid together with the well completions.
///
/// The edges of the graph are formed by the superset of the edges representing
/// the faces of the grid and the ones that represent the connections within the
/// wells. If a well has completions
/// on cell i and cell j, then there is an edge from i to j and j to i in the graph.
/// Even for shut wells the connections will exist.
class CombinedGridWellGraph
{
public:
    typedef std::vector<std::set<int> > GraphType;

    /// \brief Create a graph representing a grid together with the wells.
    /// \param grid The grid.
    /// \param eclipseState The eclipse state to extract the well information from.
    /// \param pretendEmptyGrid True if we should pretend the grid and wells are empty.
    CombinedGridWellGraph(const Dune::CpGrid& grid,
                          const std::vector<const OpmWellType*> * wells,
                          const double* transmissibilities,
                          bool pretendEmptyGrid);

    /// \brief Access the grid.
    const Dune::CpGrid& getGrid() const
    {
        return grid_;
    }

    const GraphType& getWellsGraph() const
    {
        return wellsGraph_;
    }

    double transmissibility(int face_index) const
    {
        return transmissibilities_ ? (1.0e18*transmissibilities_[face_index]) : 1;
    }

    const WellConnections& getWellConnections() const
    {
        return well_indices_;
    }
private:

    void addCompletionSetToGraph()
    {
        for(const auto& well_indices: well_indices_)
        {
            for( auto well_idx = well_indices.begin(); well_idx != well_indices.end();
                 ++well_idx)
            {
                auto well_idx2 = well_idx;
                for( ++well_idx2; well_idx2 != well_indices.end();
                     ++well_idx2)
                {
                    wellsGraph_[*well_idx].insert(*well_idx2);
                    wellsGraph_[*well_idx2].insert(*well_idx);
                }
            }
        }
    }


    const Dune::CpGrid& grid_;
    GraphType wellsGraph_;
    const double* transmissibilities_;
    WellConnections well_indices_;
};


/// \brief Sets up the call-back functions for ZOLTAN's graph partitioning.
/// \param zz The struct with the information for ZOLTAN.
/// \param grid The grid to partition.
/// \param pretendNull If true, we will pretend that the grid has zero cells.
void setCpGridZoltanGraphFunctions(Zoltan_Struct *zz, const Dune::CpGrid& grid,
                                   bool pretendNull=false);

void setCpGridZoltanGraphFunctions(Zoltan_Struct *zz,
                                   const CombinedGridWellGraph& graph,
                                   bool pretendNull);
} // end namespace cpgrid
} // end namespace Dune

#endif // HAVE_ZOLTAN
#endif // header guard
