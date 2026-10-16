/*
  Copyright 2026 Equinor ASA.

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
#ifndef OPM_GRID_REFINEMENT_REFINEDIDS_HEADER_INCLUDED
#define OPM_GRID_REFINEMENT_REFINEDIDS_HEADER_INCLUDED

#include <dune/common/parallel/mpihelper.hh>

#include <array>
#include <memory>
#include <vector>

namespace Dune::cpgrid { class CpGridData; }

namespace Opm::Refinement
{

/// Leaf ids returned by assignRefinedIds().
struct RefinedLeafIds
{
    std::vector<int> cells;
    std::vector<int> points;
    int end = 0; ///< One above every cell and point id.
};

/// Global ids of refined entities that do not depend on rank or partition, above every
/// level-zero id: cells by LGR and LGR-local Cartesian index; points by the level they are
/// born on and their index in that level grid, which is built whole from the global input;
/// leaf-only corners of faulted split faces by box and coordinate order (splitKeys, one
/// {box, ordinal} per leaf corner past the corner history).
/// Sets the mappings of the level grids storage[1..numBoxes]. Collective on comm.
RefinedLeafIds assignRefinedIds(const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& storage,
                                int numBoxes,
                                const std::vector<std::array<int,2>>& leafToLevel,
                                const std::vector<std::array<int,2>>& leafCornerHistory,
                                const std::vector<std::array<int,2>>& splitKeys,
                                Dune::MPIHelper::MPICommunicator comm);

} // namespace Opm::Refinement

#endif // OPM_GRID_REFINEMENT_REFINEDIDS_HEADER_INCLUDED
