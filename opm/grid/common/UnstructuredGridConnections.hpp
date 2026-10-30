/*
  Copyright 2026 Equinor ASA

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
#ifndef OPM_UNSTRUCTUREDGRIDCONNECTIONS_HEADER_INCLUDED
#define OPM_UNSTRUCTUREDGRIDCONNECTIONS_HEADER_INCLUDED

#include <opm/grid/UnstructuredGrid.h>

#include <utility>
#include <vector>

namespace Dune::unstructured
{

/// Connect a cell to the cell below it when there is space between them that is not a cell, as where MINPV removed a cell
/// and left empty space, with PINCH.
///
/// The bottom side of the first cell of a pair is a face at the boundary of the grid, and the top side of the second cell
/// as well. This is what the CpGrid does with the pairs that MinpvProcessor gives: the face of the first cell is given the
/// second cell as its other cell, with the geometry of the face of the first cell. The second cell keeps its own top face
/// and gets the face as one more, with the tag of its top side.
///
/// A pair is left out if the first cell has no bottom face at the boundary, which is the case if the cells are
/// connected already, or if a cell is not in the grid.
///
/// \param grid   A three dimensional grid with the tags of the sides of the cells, which is changed.
/// \param cells  Pairs of the indices of the cells in the grid: the cell above and the cell below.
/// \return The number of connections that were made.
int connectCellsAcrossGaps(UnstructuredGrid& grid, const std::vector<std::pair<int, int>>& cells);

} // namespace Dune::unstructured

#endif // OPM_UNSTRUCTUREDGRIDCONNECTIONS_HEADER_INCLUDED
