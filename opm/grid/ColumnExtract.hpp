/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#include <opm/grid/UnstructuredGrid.h>
#include <vector>

struct UnstructuredGrid;

namespace Opm {

/// Extract each column of the grid.
///  \note Assumes the pillars of the grid are all vertically aligned.
///  \param grid The grid from which to extract the columns.
///  \param columns will for each (i, j) where (i, j) represents a non-empty column,
////        contain the cell indices contained in the column
///         centered at (i, j) in the second variable, and i+jN in the first variable.
void extractColumn(const UnstructuredGrid& grid, std::vector<std::vector<int> >& columns);

} // namespace Opm
