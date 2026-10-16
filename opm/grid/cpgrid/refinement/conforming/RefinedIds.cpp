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

#include <dune/common/parallel/communication.hh>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace Opm::Refinement
{

namespace
{

using Dune::cpgrid::CpGridData;

class IdLayout
{
public:
    IdLayout(const std::vector<std::shared_ptr<CpGridData>>& storage, int numBoxes, int base,
             const std::vector<int>& pointCounts, const std::vector<int>& splitCounts)
        : storage_(storage)
        , cellOffset_(numBoxes + 1, 0)
        , pointOffset_(numBoxes + 1, 0)
        , splitOffset_(numBoxes, 0)
    {
        std::int64_t next = base;
        const auto advance = [&next](std::int64_t n) {
            const auto start = next;
            next += n;
            if (next > std::numeric_limits<int>::max()) {
                throw std::overflow_error("Refined global ids exceed the range of int.");
            }
            return static_cast<int>(start);
        };
        for (int level = 1; level <= numBoxes; ++level) {
            const auto& d = storage[level]->logicalCartesianSize();
            cellOffset_[level] = advance(std::int64_t(d[0]) * d[1] * d[2]);
        }
        for (int level = 1; level <= numBoxes; ++level) {
            pointOffset_[level] = advance(pointCounts[level]);
        }
        for (int box = 0; box < numBoxes; ++box) {
            splitOffset_[box] = advance(splitCounts[box]);
        }
        end_ = static_cast<int>(next);
    }

    int cell(int level, int levelCell) const
    {
        return level == 0 ? levelZeroId<0>(levelCell)
                          : cellOffset_[level] + storage_[level]->globalCell()[levelCell];
    }

    /// A point by where it was born: {0, level-zero corner} or {level, corner there}.
    int born(const std::array<int,2>& history) const
    {
        return history[0] == 0 ? levelZeroId<3>(history[1]) : pointOffset_[history[0]] + history[1];
    }

    int point(int level, int levelPoint) const
    {
        const auto& grid = *storage_[level];
        if (level == 0) {
            return levelZeroId<3>(levelPoint);
        }
        const auto history = grid.cornerHistorySize() ? grid.getCornerHistory(levelPoint)
                                                      : std::array<int,2>{-1, -1};
        return history[0] < 0 ? pointOffset_[level] + levelPoint : born(history);
    }

    int split(const std::array<int,2>& key) const { return splitOffset_[key[0]] + key[1]; }
    int end() const { return end_; }

private:
    template <int codim>
    int levelZeroId(int index) const
    {
        const auto& level0 = *storage_[0];
        return static_cast<int>(level0.globalIdSet().id(Dune::cpgrid::Entity<codim>(level0, index, true)));
    }

    const std::vector<std::shared_ptr<CpGridData>>& storage_;
    std::vector<int> cellOffset_, pointOffset_, splitOffset_;
    int end_ = 0;
};

} // anonymous namespace

RefinedLeafIds assignRefinedIds(const std::vector<std::shared_ptr<CpGridData>>& storage,
                                int numBoxes,
                                const std::vector<std::array<int,2>>& leafToLevel,
                                const std::vector<std::array<int,2>>& leafCornerHistory,
                                const std::vector<std::array<int,2>>& splitKeys,
                                Dune::MPIHelper::MPICommunicator comm)
{
    // Faces are not entities with ids in CpGrid; cells and points bound every id.
    int base = static_cast<int>(storage[0]->globalIdSet().getMaxGlobalId()) + 1;
    std::vector<int> pointCounts(numBoxes + 1, 0);
    std::vector<int> splitCounts(numBoxes, 0);
    for (int level = 1; level <= numBoxes; ++level) {
        pointCounts[level] = storage[level]->size(3);
    }
    for (const auto& key : splitKeys) {
        if (key[0] < 0 || key[0] >= numBoxes) {
            throw std::logic_error("A split-face corner belongs to no refined box.");
        }
        splitCounts[key[0]] = std::max(splitCounts[key[0]], key[1] + 1);
    }
    // A box is held whole by one rank; the others see it empty.
    const Dune::Communication<Dune::MPIHelper::MPICommunicator> cc(comm);
    if (cc.size() > 1) {
        base = cc.max(base);
        cc.max(pointCounts.data(), numBoxes + 1);
        if (numBoxes > 0) {
            cc.max(splitCounts.data(), numBoxes);
        }
    }
    const IdLayout ids(storage, numBoxes, base, pointCounts, splitCounts);

    for (int level = 1; level <= numBoxes; ++level) {
        auto& grid = *storage[level];
        std::vector<int> cells(grid.size(0)), faces(grid.numFaces()), points(grid.size(3));
        for (int c = 0; c < grid.size(0); ++c) {
            cells[c] = ids.cell(level, c);
        }
        for (int f = 0; f < grid.numFaces(); ++f) {
            faces[f] = grid.size(0) + f; // what the local id set would give
        }
        for (int p = 0; p < grid.size(3); ++p) {
            points[p] = ids.point(level, p);
        }
        GridStateWriter::setGlobalIdMapping(grid, std::move(cells), std::move(faces), std::move(points));
    }

    RefinedLeafIds leaf;
    leaf.end = ids.end();
    leaf.cells.reserve(leafToLevel.size());
    for (const auto& [level, index] : leafToLevel) {
        leaf.cells.push_back(ids.cell(level, index));
    }
    leaf.points.reserve(leafCornerHistory.size() + splitKeys.size());
    for (const auto& history : leafCornerHistory) {
        leaf.points.push_back(ids.born(history));
    }
    for (const auto& key : splitKeys) {
        leaf.points.push_back(ids.split(key));
    }
    return leaf;
}

} // namespace Opm::Refinement
