//===========================================================================
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
#ifndef OPM_GRID_CPGRID_LEVELCARTESIANINDEXMAPPER_HPP
#define OPM_GRID_CPGRID_LEVELCARTESIANINDEXMAPPER_HPP

#include <opm/grid/common/LevelCartesianIndexMapper.hpp>
#include <opm/grid/CpGrid.hpp>

namespace Dune
{
class CpGrid;
}

namespace Opm
{
// Interface class to access the local Cartesian grid of each level grid (when refinement).
// Further documentation in opm/grid/common/LevelCartesianIndexMapper.hpp
//
// Specialization for CpGrid
template<>
class LevelCartesianIndexMapper<Dune::CpGrid>
{
public:
    static constexpr int dimension = 3 ;

    explicit LevelCartesianIndexMapper(const Dune::CpGrid& grid,
                                       int level)
        : grid_{ &grid }
    {
        if ((level < 0) || (level > grid_->maxLevel())) {
            OPM_THROW(std::invalid_argument, "Invalid level.\n");
        }
        level_ = level;
    }

    LevelCartesianIndexMapper() = delete;

    const std::array<int,3>& cartesianDimensions() const
    {
        return grid_->currentData()[level_]->logicalCartesianSize();
    }

    int cartesianSize() const
    {
        int size = cartesianDimensions()[ 0 ];
        for( int d=1; d<dimension; ++d )
            size *= cartesianDimensions()[ d ];
        return size;
    }

    int compressedSize() const
    {
        return grid_->currentData()[level_]->size(0);
    }

    int cartesianIndex( const int levelCompressedElementIndex) const
    {
        assert(  levelCompressedElementIndex >= 0 && levelCompressedElementIndex <  grid_->currentData()[level_]->size(0) );
        return grid_->currentData()[level_]->globalCell()[levelCompressedElementIndex];
    }

    void cartesianCoordinate(const int levelCompressedElementIndex,
                             std::array<int,dimension>& levelCoords) const
    {
        grid_->currentData()[level_]->getIJK( levelCompressedElementIndex, levelCoords);
    }

private:
    const Dune::CpGrid* grid_;
    int level_;
};

}

#endif
