/*
  Copyright 2025 Equinor ASA.

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

#ifndef OPM_NESTEDREFINEMENTUTILITIES_HEADER_INCLUDED
#define OPM_NESTEDREFINEMENTUTILITIES_HEADER_INCLUDED

#include <array>
#include <string>
#include <vector>

namespace Opm
{

bool isNameInTheList(const std::vector<std::string>& lgr_parent_grid_names,
                     const std::string& parent_grid_name);

std::vector<int>
collectIndicesWithSameParentGridName(const std::vector<std::string>& lgr_parent_grid_names,
                                     const std::string& parent_grid_name);

std::tuple<std::vector<std::array<int,3>>,
           std::vector<std::array<int,3>>,
           std::vector<std::array<int,3>>,
           std::vector<std::string>>
filterLgrDataPerParentGridName(const std::vector<std::array<int,3>>& cells_per_dim_vec,
                               const std::vector<std::array<int,3>>& startIJK_vec,
                               const std::vector<std::array<int,3>>& endIJK_vec,
                               const std::vector<std::string>& lgr_name_vec,
                               const std::vector<std::string>& lgr_parent_grid_name_vec,
                               const std::string& parent_grid_name);
    
} // namespace Opm

#endif // OPM_NESTEDREFINEMENTUTILITIES_HEADER_INCLUDED
