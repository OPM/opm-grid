//===========================================================================
//
// File: LevelCartesianIndexMapper.hpp
//
// Created: Tue October 01  09:44:00 2024
//
// Author(s): Antonella Ritorto <antonella.ritorto@opm-op.com>
//
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2024 Equinor ASA.

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
#ifndef OPM_CPGRIDLEVELCARTESIANINDEXMAPPER_HH
#define OPM_CPGRIDLEVELCARTESIANINDEXMAPPER_HH

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
    static const int dimension = 3 ;

    explicit LevelCartesianIndexMapper(const Dune::CpGrid& grid) : grid_{ &grid }
    {}

    const std::array<int,3>& cartesianDimensions(int level) const
    {
        return grid_->currentData()[level]->logicalCartesianSize();
    }

    int cartesianSize(int level) const
    {
        return computeCartesianSize(level);
    }

    int compressedSize(int level) const
    {
        validLevel(level);
        return grid_->currentData()[level]->size(0);
    }

    int cartesianIndex( const int compressedElementIndex, const int level) const
    {
        validLevel(level);
        assert(  compressedElementIndex >= 0 && compressedElementIndex <  grid_->currentData()[level]->size(0) );
        return grid_->currentData()[level]->globalCell()[compressedElementIndex];
    }

    void cartesianCoordinate(const int compressedElementIndexOnLevel, std::array<int,dimension>& coordsOnLevel, int level) const
    {
        validLevel(level);
        grid_->currentData()[level]->getIJK( compressedElementIndexOnLevel, coordsOnLevel);
    }

private:
    const Dune::CpGrid* grid_;

    int computeCartesianSize(int level) const
    {
        int size = cartesianDimensions(level)[ 0 ];
        for( int d=1; d<dimension; ++d )
            size *= cartesianDimensions(level)[ d ];
        return size;
    }

    void validLevel(int level) const
    {
        if ((level < 0) || (level > grid_->maxLevel())) {
            throw std::invalid_argument("Invalid level.\n");
        }
    }
};

}

#endif
