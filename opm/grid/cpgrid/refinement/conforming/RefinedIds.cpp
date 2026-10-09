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

#include "config.h"

#include <opm/grid/cpgrid/refinement/conforming/RefinedIds.hpp>

#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/cpgrid/refinement/GridStateWriter.hpp>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace Opm::Refinement
{

RefinedCellIds::RefinedCellIds(const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& levels,
                               int numLevels, int base)
    : cartesian_(numLevels, nullptr)
    , offset_(numLevels, 0)
{
    std::int64_t next = base;
    for (int level = 1; level < numLevels; ++level) {
        cartesian_[level] = &levels[level]->globalCell();
        offset_[level] = static_cast<int>(next);
        const auto& dims = levels[level]->logicalCartesianSize();
        next += std::int64_t(dims[0]) * dims[1] * dims[2];
        if (next > std::numeric_limits<int>::max()) {
            throw std::overflow_error("Refined cell ids exceed the range of int.");
        }
    }
    end_ = static_cast<int>(next);
}

int RefinedCellIds::operator()(int level, int levelCell) const
{
    return offset_[level] + (*cartesian_[level])[levelCell];
}

int levelZeroIdEnd(const Dune::cpgrid::CpGridData& level0)
{
    // Faces are not entities with ids in CpGrid; their ids are only unique among faces.
    return static_cast<int>(level0.globalIdSet().getMaxGlobalId()) + 1;
}

void setRefinedIdMappings(const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& data,
                          int numLevels, int base, bool withLeaf)
{
    const RefinedCellIds cellIds(data, numLevels, base);
    const int shift = cellIds.end() - base;
    const auto& ids0 = data.front()->globalIdSet();
    const int last = withLeaf ? numLevels : numLevels - 1;
    int maxPointId = cellIds.end() - 1;
    for (int g = 1; g <= last; ++g) {
        auto& grid = *data[g];
        // The local id set is never replaced, so this stays idempotent.
        const auto& raw = grid.localIdSet();
        std::vector<int> cells(grid.size(0)), faces(grid.numFaces()), points(grid.size(3));
        for (int c = 0; c < grid.size(0); ++c) {
            const Dune::cpgrid::Entity<0> e(grid, c, true);
            const int level = e.level();
            cells[c] = (level == 0) ? static_cast<int>(ids0.id(e.getLevelElem()))
                                    : cellIds(level, e.getLevelElem().index());
        }
        for (int f = 0; f < grid.numFaces(); ++f) {
            faces[f] = grid.size(0) + f; // what the local id set would give
        }
        // Leaf corners of faulted split faces belong to no level and have no history.
        const int withHistory = (g == numLevels) ? static_cast<int>(grid.cornerHistorySize())
                                                 : grid.size(3);
        for (int p = 0; p < withHistory; ++p) {
            const auto id = raw.id(Dune::cpgrid::Entity<3>(grid, p, true));
            points[p] = static_cast<int>(id < base ? id : id + shift);
            maxPointId = std::max(maxPointId, points[p]);
        }
        for (int p = withHistory; p < grid.size(3); ++p) {
            points[p] = ++maxPointId;
        }
        GridStateWriter::setGlobalIdMapping(grid, std::move(cells), std::move(faces), std::move(points));
    }
}

} // namespace Opm::Refinement
