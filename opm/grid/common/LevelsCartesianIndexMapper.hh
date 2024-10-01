//===========================================================================
//
// File: LevelsCartesianIndexMapper.hh
//
// Created: Tue October 01  11:44:00 2024
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

#ifndef OPM_LEVELSCARTESIANINDEXMAPPER_HEADER
#define OPM_LEVELSCARTESIANINDEXMAPPER_HEADER

#include <array>

namespace Opm
{
template< class Grid >
class LevelsCartesianIndexMapper
{
public:

    static const int dimension = Grid :: dimension ;


    explicit LevelsCartesianIndexMapper( const Grid& )
    {}

    const std::array<int, dimension>& levelCartesianDimensions(int level) const
    {
        static std::array<int, dimension> a;
        return a;
    }

    int levelCartesianSize(int level) const
    {
        return 0;
    }

    int levelCompressedSize(int level) const
    {
        return 0;
    }

    void levelCartesianCoordinate(const int /* compressedElementIndexOnLevel */,
                                  std::array<int,dimension>& /* coordsOnLevel */,
                                  int /*level*/) const
    {
    }
    int levelCartesianIndex( const int /* compressedElementIndex */ , const int level) const
    {
        return 0;
    }

    /** to be deleted */
    int compressedLevelZeroSize() const
    {
        return 0;
    }
    void cartesianCoordinateLevel(const int /* compressedElementIndexOnLevel */,
                                  std::array<int,dimension>& /* coordsOnLevel */, int /*level*/) const
    {
    }
    int cartesianIndexLevel( const int /* compressedElementIndex */ , const int level) const
    {
        return 0;
    }
    /** to be deleted - END */
};

} // end namespace Opm
#endif
