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
/// Instead of using a specialitation for Dune::CpGrid, we implement std::enable_if
/// to overload methods with different definitions: for Dune:CpGrid and for other
/// Grid types. An auxiliary defualt template parameter (GridType = Grid) is added
/// to deal with the dependent names at template instantiation.
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

    /// \brief: Call operator
    ///
    ///         For grids different from Dune::CpGrid, it takes an element index, and
    ///         returns its cell centroid, from an EclipseGrid.
    ///
    /// \tparam     GridType     Auxiliary type to overload the method, distinguishing
    ///                          general grids from CpGrid, with std::enable_if.
    ///                          Default: GridType = Grid.
    /// \param [in] elemIdx      Element Index.
    /// \return     centroid     Centroid of the element, computed as in Eclipse.
    template<typename GridType = Grid>
    typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>, std::array<double,3>>
    operator()(std::size_t elemIdx) const;

    /// \brief: Call operator
    ///
    ///         For Dune::CpGrid, it returns a function, taking an integer,
    ///         returning cell centroid, computed as in Eclipse.
    ///
    /// \tparam    GridType      Auxiliary type to overload the method, distinguishing
    ///                          general grids from CpGrid, with std::enable_if.
    ///                          Default: GridType = Grid.
    // \param [in] elemIdx       Element Index.
    /// \return    centroid      Element centroid, computed as in Eclipse.
    template<typename GridType = Grid>
    typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>, std::array<double,3>>
    operator()(std::size_t elemIdx) const;


    const GridView& gridView_;
    const Dune::CartesianIndexMapper<Grid>* cartMapper_;
    const Opm::EclipseGrid* eclGrid_;

}; // end LookUpCellCentroid struct
}
// end namespace Opm


template<typename Grid, typename GridView>
template<typename GridType>
typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>, std::array<double,3>>
Opm::LookUpCellCentroid<Grid,GridView>::operator()(std::size_t elemIdx) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    const auto centroid = this -> eclGrid_ -> getCellCenter(this -> cartMapper_->cartesianIndex(elemIdx));
    return centroid;
}

template<typename Grid, typename GridView>
template<typename GridType>
typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,std::array<double,3>>
Opm::LookUpCellCentroid<Grid,GridView>::operator()(std::size_t elemIdx) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    const auto centroid =  this -> gridView_.grid().getEclCentroid(elemIdx); // Warning! might need to be changed, due to issue in simulators
    return centroid;
}
