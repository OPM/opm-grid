//===========================================================================
//
// File: LevelsCartesianIndexMapper.hh
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
#ifndef OPM_CPGRIDLEVELSCARTESIANINDEXMAPPER_HH
#define OPM_CPGRIDLEVELSCARTESIANINDEXMAPPER_HH

#include <dune/grid/common/mcmgmapper.hh>

#include <opm/grid/common/LevelsCartesianIndexMapper.hh>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/Entity.hpp>

#include <memory>

namespace Dune
{
class CpGrid;
}

namespace Opm
{

template<>
class LevelsCartesianIndexMapper<Dune::CpGrid>
{
public:
    static const int dimension = 3 ;

    explicit LevelsCartesianIndexMapper(const Dune::CpGrid& grid) : grid_{std::make_unique<Dune::CpGrid>(grid)}
    {}

    const std::array<int,3>& levelCartesianDimensions(int level) const
    {
        return grid_->currentData()[level]->logicalCartesianSize(); // add logicalCartesianSize in CpGridData.hpp/cpp
    }

    int levelCartesianSize(int level) const
    {
        return computeLevelCartesianSize(level);
    }

    int compressedSize() const
    {
        return grid_->globalCell().size();
    }


    int levelCompressedSize(int level) const
    {
        validLevel(level);
        return grid_->currentData()[level]->size(0);
    }

    int levelCartesianIndex( const int compressedElementIndex, const int level) const
    {
        validLevel(level);
        assert(  compressedElementIndex >= 0 && compressedElementIndex <  grid_->currentData()[level]->size(0) );
        return grid_->currentData()[level]->globalCell()[compressedElementIndex];
    }

    void levelCartesianCoordinate(const int compressedElementIndexOnLevel, std::array<int,dimension>& coordsOnLevel, int level) const
    {
        validLevel(level);
        grid_->currentData()[level]->getIJK( compressedElementIndexOnLevel, coordsOnLevel);
    }

    /** To be deleted */
    /// Additional methods related to LGRs
    int compressedLevelZeroSize() const
    {
        return (*grid_->currentData()[0]).size(0);
    }

    int cartesianIndexLevel( const int compressedElementIndex, const int level) const
    {
        if ((level < 0) || (level > grid_->maxLevel())) {
            throw std::invalid_argument("Invalid level.\n");
        }
        assert(  compressedElementIndex >= 0 && compressedElementIndex <  grid_->currentData()[level]->size(0) );
        return grid_->currentData()[level]->globalCell()[compressedElementIndex];
    }

    void cartesianCoordinateLevel(const int compressedElementIndexOnLevel, std::array<int,dimension>& coordsOnLevel, int level) const
    {
        if ((level < 0) || (level > grid_->maxLevel())) {
            throw std::invalid_argument("Invalid level.\n");
        }
        grid_->currentData()[level]->getIJK( compressedElementIndexOnLevel, coordsOnLevel);
    }
    /// Additional methods related to LGRs. END
    /** end */

private:
    std::unique_ptr<Dune::CpGrid> grid_;

    int computeLevelCartesianSize(int level) const
    {
        int size = levelCartesianDimensions(level)[ 0 ];
        for( int d=1; d<dimension; ++d )
            size *= levelCartesianDimensions(level)[ d ];
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
