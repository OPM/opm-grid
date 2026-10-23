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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <opm/grid/cpgrid/refinement/conforming/RefinedDistribution.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/utility/OpmLog.hpp>

#include <cstddef>
#include <string>

namespace Opm::Refinement
{

std::vector<int> leafPartitionFromLevelZero(const Dune::CpGrid& grid,
                                            const std::vector<int>& level0Part)
{
    const auto& level0 = *grid.currentData().front();
    const auto& dims = level0.logicalCartesianSize();
    std::vector<int> level0Compressed(static_cast<std::size_t>(dims[0]) * dims[1] * dims[2], -1);
    const auto& level0Global = level0.globalCell();
    for (std::size_t c = 0; c < level0Global.size(); ++c) {
        level0Compressed[level0Global[c]] = static_cast<int>(c);
    }
    // A refined leaf cell carries its parent's Cartesian index.
    const auto& leafGlobal = grid.currentLeafData().globalCell();
    std::vector<int> leafPart(leafGlobal.size(), 0);
    for (std::size_t c = 0; c < leafGlobal.size(); ++c) {
        const int parent = level0Compressed[leafGlobal[c]];
        if (parent >= 0 && static_cast<std::size_t>(parent) < level0Part.size()) {
            leafPart[c] = level0Part[parent];
        }
    }
    return leafPart;
}

void prepareLeafForScatter(Dune::cpgrid::CpGridData& leaf,
                           const std::vector<int>& leafPart,
                           const Dune::cpgrid::CpGridDataTraits::Communication& comm)
{
    if (comm.rank() == 0) {
        std::vector<int> counts(comm.size(), 0);
        for (const int p : leafPart) {
            if (p >= 0 && p < comm.size()) {
                ++counts[p];
            }
        }
        std::string msg = "Refined leaf of " + std::to_string(leafPart.size()) + " cells, per rank:";
        for (std::size_t r = 0; r < counts.size(); ++r) {
            msg += " " + std::to_string(counts[r]);
        }
        Opm::OpmLog::info(msg);
    }
#if HAVE_MPI
    // Rank 0 owns the whole serial leaf, so local index = global id keeps the export list ordered.
    using AttributeSet = Dune::cpgrid::CpGridDataTraits::AttributeSet;
    using LocalIndex = Dune::cpgrid::CpGridDataTraits::ParallelIndexSet::LocalIndex;
    auto& indexSet = leaf.cellIndexSet();
    indexSet.beginResize();
    for (int i = 0, n = leaf.size(0); i < n; ++i) {
        indexSet.add(i, LocalIndex(i, AttributeSet::owner, true));
    }
    indexSet.endResize();
#else
    static_cast<void>(leaf);
#endif
}

} // namespace Opm::Refinement
