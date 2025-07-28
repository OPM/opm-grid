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
#ifndef OPM_GRID_COMMON_LEAFCARTESIANINDEXMAPPER_HPP
#define OPM_GRID_COMMON_LEAFCARTESIANINDEXMAPPER_HPP

#include <array>

namespace Opm
{
// Interface class to access the leaf Cartesian grid.
// Relevant for globally refined (corner-point) grids.
template<class Grid>
class LeafCartesianIndexMapper
{
public:
    static constexpr int dimension = Grid::dimension;

    explicit LeafCartesianIndexMapper(const Grid&)
    {
        // TODO?: throw if grid has not been globally refined
    }

    LeafCartesianIndexMapper() = delete;

    // Return the number of cells in each direction (Cartesian dimensions) of the leaf grid.
    const std::array<int, dimension>& cartesianDimensions() const
    {
        static std::array<int, dimension> a;
        return a;
    }

    // Return total number of (inactive/active) cells in the leaf grid.
    int cartesianSize() const
    {
        return 0;
    }

    // Return number of active cells in the leaf grid.
    int compressedSize() const
    {
        return 0;
    }

    // Return leaf Cartesian index of a cell in the leaf grid.
    int cartesianIndex( const int /* leafCompressedElementIndex */) const
    {
        return 0;
    }

    // Compute leaf Cartesian coordinate, i.e. IJK, for a given cell, in the leaf grid.
    void cartesianCoordinate(const int /* leafCompressedElementIndex*/,
                             std::array<int,dimension>& /* leeafCoords */) const
    {
    }

private:
    const Grid& grid_;
};

}

#endif
