/*
  Copyright 2025 Equinor ASA.

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

#ifndef OPM_GRID_COMMON_CARTESIANINDEXMAPPERCOLLECTION_HPP
#define OPM_GRID_COMMON_CARTESIANINDEXMAPPERCOLLECTION_HPP

#include <opm/grid/common/LeafCartesianIndexMapper.hpp>
#include <opm/grid/common/LevelCartesianIndexMapper.hpp>

namespace Opm
{
// Interface class to access the local Cartesian grid of each level (when refinement) and leaf grid.
template< class Grid >
class CartesianIndexMapperCollection
{
public:
    // Dimension of the grid.
    static constexpr int dimension = Grid :: dimension;

    // Constructor taking a grid.
    explicit CartesianIndexMapperCollection(const Grid&)
    {}

    // Get level Cartesian index mapper
    LevelCartesianIndexMapper<Grid> getLevelMapper(int level) const
    {
        return LevelCartesianIndexMapper<Grid>{};
    }

    // Get leaf Cartesian index mapper
    LeafCartesianIndexMapper<Grid> getLeafMapper() const
    {
        return LeafCartesianIndexMapper<Grid>{};
    }
private:
    const Grid& grid_;
};

} // end namespace Opm
#endif
