/*
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services.
  Copyright 2015 NTNU

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
#include <config.h>
#include <dune/grid/common/ZoltanGraphFunctions.hpp>
#ifdef HAVE_ZOLTAN
namespace Dune
{
namespace cpgrid
{
void getCpGridVertexList(void* cpGridPointer, int numGlobalIdEntries,
                         int numLocalIdEntries, ZOLTAN_ID_PTR gids,
                         ZOLTAN_ID_PTR lids, int wgtDim,
                         float *objWgts, int *err)
{
    (void) wgtDim; (void) objWgts;
    const Dune::CpGrid&  grid = *static_cast<const Dune::CpGrid*>(cpGridPointer);
    auto& globalIdSet         =  grid.globalIdSet();
    auto& localIdSet          =  grid.localIdSet();

    if ( numGlobalIdEntries != numLocalIdEntries || numGlobalIdEntries != 1 )
    {
        std::cerr<<"numGlobalIdEntries="<<numGlobalIdEntries<<" numLocalIdEntries="<<numLocalIdEntries<<" grid cells="
                 <<grid.numCells()<<std::endl;
        *err = ZOLTAN_FATAL;
        return;
    }
    int idx = 0;
    for (auto cell = grid.leafbegin<0>(), cellEnd = grid.leafend<0>();
         cell != cellEnd; ++cell)
    {
        gids[idx]   = globalIdSet.id(*cell);
        lids[idx++] = localIdSet.id(*cell);
    }
    *err = ZOLTAN_OK;
}

void getCpGridNumEdgesList(void *cpGridPointer, int sizeGID, int sizeLID,
                           int numCells,
                           ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                           int *numEdges, int *err)
{
    (void) globalID;
    const Dune::CpGrid&  grid = *static_cast<const Dune::CpGrid*>(cpGridPointer);
    if ( sizeGID != 1 || sizeLID != 1 || numCells != grid.numCells() )
    {
        *err = ZOLTAN_FATAL;
        return;
    }
    for( int i = 0; i < numCells;  i++ )
    {
        // For the graph there is an edge only if the face has two neighbors.
        // Therefore we need to check each face
        int edges = 0;
        int lid   = localID[i];
        for ( int local_face = 0; local_face < grid.numCellFaces(static_cast<int>(localID[i])); ++local_face )
        {
            const int face = grid.cellFace(lid, local_face);
            if ( grid.faceCell(face, 0) != -1 && grid.faceCell(face, 1) != -1 )
            {
                ++edges;
            }
        }
        numEdges[i] = edges;
    }
    std::cout<<std::endl;
    *err = ZOLTAN_OK;
}

void getCpGridEdgeList(void *cpGridPointer, int sizeGID, int sizeLID,
                       int numCells, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                       int *numEdges,
                       ZOLTAN_ID_PTR nborGID, int *nborProc,
                       int wgtDim, float *ewgts, int *err)
{
    (void) wgtDim; (void) globalID; (void) numEdges; (void) ewgts;
    const Dune::CpGrid&  grid = *static_cast<const Dune::CpGrid*>(cpGridPointer);
    if ( sizeGID != 1 || sizeLID != 1 || numCells != grid.numCells() )
    {
        *err = ZOLTAN_FATAL;
        return;
    }
    int oldidx = 0;
    int idx = 0;

    for( int cell = 0; cell < numCells;  cell++ )
    {
        const int currentCell = localID[cell];
        for ( int local_face = 0 ; local_face < grid.numCellFaces(static_cast<int>(localID[cell])); ++local_face )
        {
            const int face  = grid.cellFace(currentCell, local_face);
            int otherCell   = grid.faceCell(face, 0);
            if ( otherCell == currentCell || otherCell == -1 )
            {
                otherCell = grid.faceCell(face, 1);
                if ( otherCell == currentCell || otherCell == -1 )
                {
                    continue;
                }
                else
                {
                    nborGID[idx++] = globalID[otherCell];
                    continue;
                }
            }
            nborGID[idx++] = globalID[otherCell];
        }
        assert(numEdges[cell] == idx - oldidx);
        oldidx = idx;
    }

    const int myrank = grid.comm().rank();

    for ( int i = 0; i < idx; ++i )
    {
        nborProc[i] = myrank;
    }
#ifdef DEBUG
    // The above relies heavily on the grid not being distributed already.
    // Therefore we check here that all cells are owned by us.
    GlobalLookupIndexSet<Dune::CpGrid::ParallelIndexSet> globalIdxSet(getCellIndexSet(),
                                                                      grid.numCells());
    for ( int cell = 0; cell < numCells;  cell++ )
    {
        if ( globaIdxSet.pair(cell).second.attribute() !=
             Dune::CpGrid::ParallelIndexSet::AttributeSet::owner )
        {
            *err = ZOLTAN_FATAL;
        }
    }
#endif
}

void setCpGridZoltanGraphFunctions(Zoltan_Struct *zz, Dune::CpGrid& grid)
{
      Zoltan_Set_Num_Obj_Fn(zz, getCpGridNumCells, &grid);
      Zoltan_Set_Obj_List_Fn(zz, getCpGridVertexList, &grid);
      Zoltan_Set_Num_Edges_Multi_Fn(zz, getCpGridNumEdgesList, &grid);
      Zoltan_Set_Edge_List_Multi_Fn(zz, getCpGridEdgeList, &grid);
}
} // end namespace cpgrid
} // end namespace Dune
#endif // HAVE_ZOLTAN
