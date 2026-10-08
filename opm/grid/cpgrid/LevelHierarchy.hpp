//===========================================================================
//
// File: CpGridData.hpp
//
// Created: Sep 17 21:11:41 2013
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
//            Markus Blatt        <markus@dr-blatt.de>
//            Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// Comment: Major parts of this file originated in dune/grid/CpGrid.hpp
//          and got transfered here during refactoring for the parallelization.
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010, 2013, 2022-2023 Equinor ASA.
  Copyright 2013 Dr. Blatt - HPC-Simulation-Software & Services

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

#ifndef OPM_LEVELHIERARCHY_HEADER
#define OPM_LEVELHIERARCHY_HEADER

#include <array>
#include <memory>
#include <tuple>
#include <vector>

namespace Dune
{
namespace cpgrid
{

class CpGridData;

/// How one CpGridData relates to the other levels of a refined CpGrid.
struct LevelHierarchy
{
    /** Mark elements to be refined **/
    std::vector<int> mark;
    /** Level of the current CpGridData (0 when it's "GLOBAL", 1,2,.. for LGRs). */
    int level{0};
    /** Copy of (CpGrid object).data_ associated with the CpGridData object. */
    std::vector<std::shared_ptr<CpGridData>>* level_data_ptr{};
    // SUITABLE FOR ALL LEVELS EXCEPT FOR LEAFVIEW
    /** Map between level and leafview cell indices. Only cells (from that level) that appear in leafview count. -1 when the cell vanished.*/
    std::vector<int> level_to_leaf_cells;  // In entry 'level cell index', we store 'leafview cell index'
    /** Parent cells and their children. Entry is {-1, {}} when cell has no children.*/  // {level LGR, {child0, child1, ...}}
    std::vector<std::tuple<int,std::vector<int>>> parent_to_children_cells;
    /** Amount of children cells per parent cell in each direction. */  // {# children in x-direction, ... y-, ... z-}
    std::array<int,3> cells_per_dim;
    // SUITABLE ONLY FOR LEAFVIEW
    /** Relation between leafview and (possible different) level(s) cell indices. */  // {level, cell index in that level}
    std::vector<std::array<int,2>> leaf_to_level_cells;
    /** Corner history. corner_history[ corner index ] = {level where the corner was born, its index there }, {-1,-1} otherwise. */
    std::vector<std::array<int,2>> corner_history;
    // SUITABLE FOR ALL LEVELS INCLUDING LEAFVIEW
    /** Child cells and their parents. Entry is {-1,-1} when cell has no father. */  // {level parent cell, parent cell index}
    std::vector<std::array<int,2>> child_to_parent_cells;
    /** Level-grid or Leaf-grid cell to parent cell and refined-cell-in-parent-cell index (number between zero and total amount
        of children per parent (cells_per_dim[0]_*cells_per_dim[1]*cells_per_dim[2])). Entry is -1 when cell has no father. */
    std::vector<int> cell_to_idxInParentCell;
    /** To keep track of refinement processes */
    int refinement_max_level{0};
};

} // namespace cpgrid
} // namespace Dune

#endif // OPM_LEVELHIERARCHY_HEADER
