//===========================================================================
//
// File: levelcartesianindexmapper.hh
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
#ifndef OPM_POLYHEDRALGRIDLEVELCARTESIANINDEXMAPPER_HH
#define OPM_POLYHEDRALGRIDLEVELCARTESIANINDEXMAPPER_HH

#include <opm/grid/common/LevelCartesianIndexMapper.hpp>
#include <opm/grid/polyhedralgrid.hh>


namespace Dune
{
template<int dim, int dimworld, typename coord_t>
class PolyhedralGrid;

}

namespace Opm
{
// Interface class to access the local Cartesian grid of each level grid (when refinement).
// Further documentation in opm/grid/common/LevelCartesianIndexMapper.hpp
//
// Specialization for PolyhedralGrid
template<int dim, int dimworld, typename coord_t>
class LevelCartesianIndexMapper<Dune::PolyhedralGrid< dim, dimworld, coord_t >>
{
    typedef Dune::PolyhedralGrid< dim, dimworld, coord_t >  Grid;
public:
    static const int dimension = 3 ;

    explicit LevelCartesianIndexMapper(const Grid& grid) : grid_{ &grid }
    {}

    const std::array<int,3>& cartesianDimensions(int level) const
    {
        throwIfLevelPositive(level);
        return grid_->logicalCartesianSize();
    }

    int cartesianSize(int level) const
    {
        throwIfLevelPositive(level);
        return computeCartesianSize(0);
    }

    int compressedSize(int level) const
    {
        throwIfLevelPositive(level);
        return grid_->size(0);
    }

    int cartesianIndex( const int compressedElementIndex, const int level) const
    {
        throwIfLevelPositive(level);
        assert( compressedElementIndex >= 0 && compressedElementIndex < compressedSize(0) );
        return grid_->globalCell()[ compressedElementIndex ];
    }

    void cartesianCoordinate(const int compressedElementIndex, std::array<int,dimension>& coords, int level) const
    {
        throwIfLevelPositive(level);

        int gc = cartesianIndex( compressedElementIndex, 0);
        auto cartesianDimensions = grid_->logicalCartesianSize();
        if( dimension >=2 )
        {
            for( int d=0; d<dimension-2; ++d )
            {
                coords[d] = gc % cartesianDimensions[d];  gc /= cartesianDimensions[d];
            }

            coords[dimension-2] = gc % cartesianDimensions[dimension-2];
            coords[dimension-1] = gc / cartesianDimensions[dimension-1];
        }
        else
            coords[ 0 ] = gc ;
    }

private:
    const Grid* grid_;

    int computeCartesianSize(int level) const
    {
        int size = cartesianDimensions(level)[ 0 ];
        for( int d=1; d<dimension; ++d )
            size *= cartesianDimensions(level)[ d ];
        return size;
    }

    void throwIfLevelPositive(int level) const
    {
        if (level) {
            throw std::invalid_argument("Invalid level.\n");
        }
    }
};

}

#endif
