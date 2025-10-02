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

/// @brief Check whether a grid name exists in a list of parent grid names.
///
/// @param [in] lgr_parent_grid_names  List of parent grid names (std::string)
///                                    used for refinement.
/// @param [in] grid_name              Grid name to search for.
/// @return true if grid_name is found in lgr_parent_grid_names,
///         false otherwise.
bool isNameInTheList(const std::vector<std::string>& lgr_parent_grid_names,
                     const std::string& grid_name);

/// @brief Group LGRs by their parent grid, identified through indices.
///
/// @param [in] lgr_parent_grid_names  List of parent grid names (std::string)
///                                    used for refinement.
/// @param [in] parent_grid_name       Parent grid name to match.
/// @return A vector of indices corresponding to LGRs associated with the
///         given parent grid. These indices can be used to access LGR data
///         and trigger refinement.
std::vector<int>
getLgrDataIndicesByParentGrid(const std::vector<std::string>& lgr_parent_grid_names,
                              const std::string& parent_grid_name);

/// @brief Retrieve LGR data associated with a specific parent grid.
///
/// Filters LGR information (cell refinements, block bounds, and LGR names)
/// based on the given parent grid name.
///
/// @param [in] cells_per_dim_vec        Vector of triplets specifying the number of
///                                      refined cells in each dimension for each patch.
/// @param [in] startIJK_vec             Vector of triplets specifying the starting
///                                      Cartesian indices of each patch.
/// @param [in] endIJK_vec               Vector of triplets specifying the ending
///                                      Cartesian indices of each patch.
/// @param [in] lgr_name_vec             Vector of LGR/level names.
/// @param [in] lgr_parent_grid_name_vec Vector of parent grid names corresponding
///                                      to each LGR.
/// @param [in] parent_grid_name         Parent grid name to filter by.
/// @return A tuple containing:
///         - Vector of cell refinement triplets (cells per dimension),
///         - Vector of starting index triplets,
///         - Vector of ending index triplets,
///         - Vector of LGR names,
///         all restricted to the specified parent grid.
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

/// @brief Check whether each parent grid exists before its child LGRs.
///
/// Ensures that a parent grid is either already present in the set of existing
/// grids or appears earlier in the list of new LGRs, before its children.
///
/// @param [in] existing_grid_names         Map of existing grid names.
/// @param [in] new_lgr_names               Names of new LGRs to be created.
/// @param [in] new_lgrs_parent_grid_names  Corresponding parent grid names for each LGR.
/// @return true if all parent grids exist before their child LGRs, false otherwise.
bool areParentGridsAvailableBeforeTheirLgrs(const std::map<std::string,int>& existing_grid_names,
                                            const std::vector<std::string>& new_lgr_names,
                                            const std::vector<std::string>& new_lgrs_parent_grid_names);

} // namespace Opm

#endif // OPM_NESTEDREFINEMENTUTILITIES_HEADER_INCLUDED
