//===========================================================================
//
// File: LookUpCellCentroid.hh
//
// Created: Wed July 26 10:48:00 2023
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
  Copyright 2023 Equinor ASA.

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

#include <dune/grid/common/mcmgmapper.hh>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>

#include <type_traits>

namespace Dune
{
class CpGrid;
}

namespace Opm
{

class EclipseGrid;

/// LookUpCellCentroid struct - To search cell centroids via element index
///
/// Instead of using a specialitation for Dune::CpGrid, we implement to member
/// functions: for Dune:CpGrid and for other Grid types.
template <typename Grid, typename GridView>
struct LookUpCellCentroid
{
    /// \brief:     Constructor taking a GridView, CartesianMapper
    /// \param [in] GridView
    /// \param [in] CartesianIndexMapper
    /// \param [in] EclipseGrid
    explicit LookUpCellCentroid(const GridView& gridView,
                                const Dune::CartesianIndexMapper<Grid>& cartMapper,
                                const Opm::EclipseGrid* eclgrid) :
        gridView_(gridView),
        cartMapper_(&cartMapper),
        eclGrid_(eclgrid)
    {
    }

    /// \brief: getCentroidFromEclGrid
    ///
    ///         For grids different from Dune::CpGrid, it takes an element index, and
    ///         returns its cell centroid, from an EclipseGrid.
    ///
    /// \param [in] elemIdx      Element Index.
    /// \return     centroid     Centroid of the element, computed as in Eclipse.
    std::array<double,3> getCentroidFromEclGrid(std::size_t elemIdx) const;

    /// \brief: getCentroidFromCpGrid
    ///
    ///         For  Dune::CpGrid, it takes an element index, and
    ///         returns its cell centroid, computed as in EclipseGrid.
    ///
    /// \param [in] elemIdx      Element Index.
    /// \return     centroid     Centroid of the element, computed as in Eclipse.
    std::array<double,3> getCentroidFromCpGrid(std::size_t elemIdx) const;

    const GridView& gridView_;
    const Dune::CartesianIndexMapper<Grid>* cartMapper_;
    const Opm::EclipseGrid* eclGrid_;

}; // end LookUpCellCentroid struct
}
// end namespace Opm

template<typename Grid, typename GridView>
std::array<double,3> Opm::LookUpCellCentroid<Grid,GridView>::getCentroidFromEclGrid(std::size_t elemIdx) const
{
    return eclGrid_ -> getCellCenter(cartMapper_ -> cartesianIndex(elemIdx));
}

template<typename Grid, typename GridView>
std::array<double,3> Opm::LookUpCellCentroid<Grid,GridView>::getCentroidFromCpGrid(std::size_t elemIdx) const
{
    static_assert(std::is_same_v<Grid,Dune::CpGrid>);
    return gridView_.grid().getEclCentroid(elemIdx);
}
