/*
  Copyright 2014, 2015 Dr. Markus Blatt - HPC-Simulation-Software & Services.
  Copyright 2014 Statoil AS
  Copyright 2015 NTNU

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
*/

#include <config.h>

#include <dune/grid/cpgrid/GridHelpers.hpp>

namespace Opm
{
// Interface functions using CpGrid

namespace UgGridHelpers
{

int numCells(const Dune::CpGrid& grid)
{
    return grid.numCells();
}

int numFaces(const  Dune::CpGrid& grid)
{
    return grid.numFaces();
}

int dimensions(const Dune::CpGrid&)
{
    return Dune::CpGrid::dimension;
}

int numCellFaces(const Dune::CpGrid& grid)
{
    return grid.numCellFaces();    
}

const int* cartDims(const Dune::CpGrid& grid)
{
    return &(grid.logicalCartesianSize()[0]);
}

const int*  globalCell(const Dune::CpGrid& grid)
{
    return &(grid.globalCell()[0]);
}

CellCentroidTraits<Dune::CpGrid>::IteratorType
beginCellCentroids(const Dune::CpGrid& grid)
{
    return CellCentroidTraits<Dune::CpGrid>::IteratorType(grid, 0);
}

double cellCentroidCoordinate(const Dune::CpGrid& grid, int cell_index,
                              int coordinate)
{
    return grid.cellCentroid(cell_index)[coordinate];
}

FaceCentroidTraits<Dune::CpGrid>::IteratorType
beginFaceCentroids(const Dune::CpGrid& grid)
{
    return FaceCentroidTraits<Dune::CpGrid>::IteratorType(grid, 0);
}

const double* cellCentroid(const Dune::CpGrid& grid, int cell_index)
{
    return &(grid.cellCentroid(cell_index)[0]);
}

double cellVolume(const  Dune::CpGrid& grid, int cell_index)
{
    return grid.cellVolume(cell_index);
}

CellVolumeIterator beginCellVolumes(const Dune::CpGrid& grid)
{
    return CellVolumeIterator(grid, 0);
}

CellVolumeIterator endCellVolumes(const Dune::CpGrid& grid)
{
    return CellVolumeIterator(grid, numCells(grid));
}

const FaceCentroidTraits<Dune::CpGrid>::ValueType&
faceCentroid(const Dune::CpGrid& grid, int face_index)
{
    return grid.faceCentroid(face_index);
}

Dune::cpgrid::Cell2FacesContainer cell2Faces(const Dune::CpGrid& grid)
{
    return Dune::cpgrid::Cell2FacesContainer(&grid);
}

FaceCellTraits<Dune::CpGrid>::Type
faceCells(const Dune::CpGrid& grid)
{
    return Dune::cpgrid::FaceCellsContainerProxy(&grid);
}

Face2VerticesTraits<Dune::CpGrid>::Type
face2Vertices(const Dune::CpGrid& grid)
{
    return Dune::cpgrid::FaceVerticesContainerProxy(&grid);
}

const double* vertexCoordinates(const Dune::CpGrid& grid, int index)
{
    return &(grid.vertexPosition(index)[0]);
}

const double* faceNormal(const Dune::CpGrid& grid, int face_index)
{
    return &(grid.faceNormal(face_index)[0]);
}

double faceArea(const Dune::CpGrid& grid, int face_index)
{
    return grid.faceArea(face_index);
}

int faceTag(const Dune::CpGrid& grid,
            const Dune::cpgrid::Cell2FacesRow::iterator& cell_face)
{
    return grid.faceTag(cell_face);
}
} // end namespace UgGridHelpers

} // end namespace Opm
