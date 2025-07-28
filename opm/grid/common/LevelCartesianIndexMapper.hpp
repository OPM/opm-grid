/*
  Copyright 2024, 2025 Equinor ASA.

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

#ifndef OPM_GRID_COMMON_LEVELCARTESIANINDEXMAPPER_HPP
#define OPM_GRID_COMMON_LEVELCARTESIANINDEXMAPPER_HPP

#include <array>

namespace Opm
{
// Interface class to access the local Cartesian grid of a level grid.
template< class Grid >
class LevelCartesianIndexMapper
{
public:
    // Dimension of the grid.
    static constexpr int dimension = Grid::dimension;

    // Constructor taking a grid.
    explicit LevelCartesianIndexMapper( const Grid& grid,
                                        int level)
    {}

    // Return the number of cells in each direction (Cartesian dimensions) of the level grid.
    const std::array<int, dimension>& cartesianDimensions() const
    {
        static std::array<int, dimension> a;
        return a;
    }

    // Return total number of cells in the level grid.
    int cartesianSize() const
    {
        return 0;
    }

    // Return number of active cells in the level grid.
    int compressedSize() const
    {
        return 0;
    }

    // Return local/level Cartesian index of a cell in the level grid.
    int cartesianIndex( const int /* levelCompressedElementIndex */) const
    {
        return 0;
    }

    // Compute local/level Cartesian coordinate, i.e. IJK, for a given cell, on the level grid.
    void cartesianCoordinate(const int /* levelCompressedElementIndex*/,
                             std::array<int,dimension>& /* levelCoords */) const
    {
    }
private:
    const Grid& grid_;
    int level_;
};

} // end namespace Opm
#endif
