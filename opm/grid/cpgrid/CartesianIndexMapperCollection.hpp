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

#ifndef OPM_GRID_CPGRID_CARTESIANINDEXMAPPERCOLLECTION_HPP
#define OPM_GRID_CPGRID_CARTESIANINDEXMAPPERCOLLECTION_HPP

#include <opm/grid/common/CartesianIndexMapperCollection.hpp>
#include <opm/grid/cpgrid/LeafCartesianIndexMapper.hpp>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>
#include <opm/grid/CpGrid.hpp>

namespace Dune
{
class CpGrid;
}

namespace Opm
{
// Interface class to access the local Cartesian grid of each level (when refinement) and leaf grid.
// Specialization for CpGrid
template<>
class CartesianIndexMapperCollection<Dune::CpGrid>
{
public:
    // Dimension of the grid.
    static constexpr int dimension = 3;

    // Constructor taking a grid.
    explicit CartesianIndexMapperCollection(const Dune::CpGrid& grid) : grid_{ &grid }
    {}

    CartesianIndexMapperCollection() = delete;

    // Get level Cartesian index mapper
    LevelCartesianIndexMapper<Dune::CpGrid> getLevelMapper(int level) const
    {
        return LevelCartesianIndexMapper<Dune::CpGrid>{*grid_, level};
    }

    // Get leaf Cartesian index mapper
    LeafCartesianIndexMapper<Dune::CpGrid> getLeafMapper() const
    {
        return LeafCartesianIndexMapper<Dune::CpGrid>{*grid_};
    }
private:
    const Dune::CpGrid* grid_;
};

} // end namespace Opm
#endif
