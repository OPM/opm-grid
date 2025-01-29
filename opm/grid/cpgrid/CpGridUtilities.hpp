/*
  Copyright 2025 ....

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

#ifndef OPM_CPGRIDUTILITIES_HEADER_INCLUDED
#define OPM_CPGRIDUTILITIES_HEADER_INCLUDED

#include <opm/grid/CpGrid.hpp>

namespace Opm
{

/// @brief Retrieves Cartesian indices for a specified Local Grid Refinement (LGR) level in a Dune::CpGrid.
///
/// This function extracts the mapping between active cell indices and Cartesian indices for a given LGR name.
/// If the specified name is not found, an exception is thrown.
///
/// @param grid The Dune::CpGrid
/// @param lgr_name The name of the LGR whose indices are to be retrieved.
/// @return A tuple containing:
///   - A std::vector<int> mapping cell indices to their corresponding Cartesian indices in the LGR.
///   - A std::unordered_map<int, int> mapping Cartesian indices back to cell indices (handles inactive parent cells).
///   - A std::vector<std::array<int, 3>> storing the (i, j, k) Cartesian coordinates for active cells.
std::tuple<std::vector<int>, std::unordered_map<int, int>, std::vector<std::array<int, 3>>>
lgrIJK(const Dune::CpGrid& grid, const std::string& lgr_name);

} // namespace Opm

#endif // OPM_CPGRIDUTILITIES_HEADER_INCLUDED
