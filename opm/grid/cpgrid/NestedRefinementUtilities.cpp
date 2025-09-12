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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <algorithm>
#include <array>
#include <cassert>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <tuple>

namespace Opm
{

bool isNameInTheList(const std::vector<std::string>& lgr_parent_grid_names,
                     const std::string& grid_name)
{
    return (std::find(lgr_parent_grid_names.begin(), lgr_parent_grid_names.end(),
                      grid_name) != lgr_parent_grid_names.end());
}

std::vector<int> getLgrDataIndicesByParentGrid(const std::vector<std::string>& lgr_parent_grid_names,
                                               const std::string& parent_grid_name)
{
    if (!isNameInTheList(lgr_parent_grid_names, parent_grid_name)) {
        throw std::invalid_argument("Parent grid name does not exist.\n");
    }

    std::vector<int> indices;

    for (std::size_t i = 0; i < lgr_parent_grid_names.size(); ++i) {
        if (lgr_parent_grid_names[i] == parent_grid_name) {
            indices.push_back(i);
        }
    }
    return indices;
}

std::tuple<std::vector<std::array<int,3>>,
           std::vector<std::array<int,3>>,
           std::vector<std::array<int,3>>,
           std::vector<std::string>>
filterLgrDataPerParentGridName(const std::vector<std::array<int,3>>& cells_per_dim_vec,
                               const std::vector<std::array<int,3>>& startIJK_vec,
                               const std::vector<std::array<int,3>>& endIJK_vec,
                               const std::vector<std::string>& lgr_name_vec,
                               const std::vector<std::string>& lgr_parent_grid_name_vec,
                               const std::string& parent_grid_name)
{
    const auto lgrs_with_given_parent_grid = getLgrDataIndicesByParentGrid(lgr_parent_grid_name_vec,
                                                                           parent_grid_name);
    // Assume all vector sizes are equal
    const auto& size =  lgrs_with_given_parent_grid.size();
    std::vector<std::array<int,3>> filtered_cells_per_dim_vec{};
    std::vector<std::array<int,3>> filtered_startIJK_vec{};
    std::vector<std::array<int,3>> filtered_endIJK_vec{};
    std::vector<std::string> filtered_lgr_name_vec{};
    filtered_cells_per_dim_vec.reserve(size);
    filtered_startIJK_vec.reserve(size);
    filtered_endIJK_vec.reserve(size);
    filtered_lgr_name_vec.reserve(size);

    for (const auto& lgr_idx : lgrs_with_given_parent_grid) {
        filtered_cells_per_dim_vec.push_back(cells_per_dim_vec[lgr_idx]);
        filtered_startIJK_vec.push_back(startIJK_vec[lgr_idx]);
        filtered_endIJK_vec.push_back(endIJK_vec[lgr_idx]);
        filtered_lgr_name_vec.push_back(lgr_name_vec[lgr_idx]);
    }
    return std::make_tuple(filtered_cells_per_dim_vec, filtered_startIJK_vec, filtered_endIJK_vec, filtered_lgr_name_vec);
}

bool areParentGridsAvailableBeforeTheirLgrs(const std::map<std::string,int>& existing_grid_names,
                                            const std::vector<std::string>& new_lgr_names,
                                            const std::vector<std::string>& new_lgrs_parent_grid_names)
{
    assert(new_lgr_names.size() == new_lgrs_parent_grid_names.size());

    for (std::size_t i = 0; i < new_lgrs_parent_grid_names.size(); ++i) {
        const std::string& parent = new_lgrs_parent_grid_names[i];

        // Case 1: Parent already exists
        if (existing_grid_names.find(parent) != existing_grid_names.end()) {
            continue;
        }

        // Case 2: Parent must appear earlier among the new LGRs
        bool foundEarlier = std::find(new_lgr_names.begin(),
                                      new_lgr_names.begin() + i,
                                      parent) != new_lgr_names.begin() + i;

        if (!foundEarlier) {
            return false; // parent not found before child
        }
    }
    return true; // all parent grids valid (exist before their LGRs)
}

} // namespace Opm
