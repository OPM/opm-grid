//===========================================================================
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
#ifndef OPM_GRID_CPGRID_LEAFCARTESIANINDEXMAPPER_HPP
#define OPM_GRID_CPGRID_LEAFCARTESIANINDEXMAPPER_HPP

#include <opm/grid/common/LeafCartesianIndexMapper.hpp>
#include <opm/grid/CpGrid.hpp>

namespace Dune
{
class CpGrid;
}

namespace Opm
{
// Interface class to access the leaf Cartesian grid.
// Relevant for globally refined grids.
//
// Specialization for CpGrid
template<>
class LeafCartesianIndexMapper<Dune::CpGrid>
{
public:
    static constexpr int dimension = 3 ;

    explicit LeafCartesianIndexMapper(const Dune::CpGrid& grid)
        : grid_{ &grid }
    {}

    LeafCartesianIndexMapper() = delete;

    const std::array<int,3>& cartesianDimensions() const
    {
        return grid_->currentData().back()->logicalCartesianSize();
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
        return grid_->currentData().back()->size(0);
    }

    int cartesianIndex( const int leafCompressedElementIndex) const
    {
        assert(  leafCompressedElementIndex >= 0 && leafCompressedElementIndex <  grid_->currentData().back()->size(0) );
        return grid_->currentData().back()->globalCell()[leafCompressedElementIndex];
    }

    void cartesianCoordinate(const int leafCompressedElementIndex,
                             std::array<int,dimension>& leafCoords) const
    {
        grid_->currentData().back()->getIJK( leafCompressedElementIndex, leafCoords);
    }

private:
    const Dune::CpGrid* grid_;
};

}

#endif
