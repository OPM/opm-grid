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

#include <memory>
#include <vector>

namespace Dune::cpgrid { class CpGridData; }

namespace Opm::Refinement
{

/// Global id of a refined cell: base + the Cartesian sizes of the LGRs before
/// its own + its LGR-local Cartesian index. The same in serial and on any partition.
class RefinedCellIds
{
public:
    /// \param levels Level grids, level zero first; entries from numLevels on are ignored.
    /// \param base   One above every level-zero id (levelZeroIdEnd).
    RefinedCellIds(const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& levels,
                   int numLevels, int base);

    int operator()(int level, int levelCell) const;

    /// One above every refined cell id.
    int end() const { return end_; }

private:
    std::vector<const std::vector<int>*> cartesian_;
    std::vector<int> offset_;
    int end_ = 0;
};

/// One above every level-zero cell and point id held by this process.
int levelZeroIdEnd(const Dune::cpgrid::CpGridData& level0);

/// Id mappings for refined level grids [1, numLevels) and, if withLeaf, the
/// leaf at data[numLevels]: refined cells from RefinedCellIds, points born on a
/// refined level moved above the refined cells.
void setRefinedIdMappings(const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& data,
                          int numLevels, int base, bool withLeaf);

} // namespace Opm::Refinement

#endif // OPM_GRID_REFINEMENT_REFINEDIDS_HEADER_INCLUDED
