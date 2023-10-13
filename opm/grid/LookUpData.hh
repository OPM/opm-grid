//===========================================================================
//
// File: LookUpData.hh
//
// Created: Tue May 23 14:44:00 2023
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

#include <opm/grid/cpgrid/Entity.hpp>

#include <type_traits>
#include <vector>

namespace Dune
{
class CpGrid;

namespace cpgrid
{
template<int codim> class Entity;
}

}

namespace Opm
{

/// LookUpData class - To search data via element index
///
/// Instead of using a specialitation for Dune::CpGrid, we implement std::enable_if
/// to overload methods with different definitions: for Dune:CpGrid and for other
/// Grid types. An auxiliary defualt template parameter (GridType = Grid) is added
/// to deal with the dependent names at template instantiation.
template <typename Grid, typename GridView>
class LookUpData
{
public:
    /// \brief:     Constructor taking a GridView
    /// \param [in] GridView
    explicit LookUpData(const  GridView& gridView) :
        gridView_(gridView),
        elemMapper_(gridView, Dune::mcmgElementLayout())
    {
    }

    /// \brief: Call operator taking an EntityObject and a FeatureVector.
    ///
    ///         Return feature of the entity, via (ACTIVE) INDEX
    ///         For general grids, the feature vector is given for the gridView_.
    ///         For CpGrid, the feature vector is given for level 0.
    ///
    /// \tparam     EntityType  Element type.
    /// \tparam     FeatureType Type of the property of the element, e.g. int, double, float, etc.
    /// \param [in] element     EntityType object.
    /// \param [in] feature_vec Vector with each entry, the feature of an element of
    ///                         the gridView_ [for general grids], or level 0 for CpGrid.
    /// \return feature of the given element.
    template<typename EntityType, typename FeatureType>
    FeatureType operator()(const EntityType& elem, const std::vector<FeatureType>& feature_vec) const;

    /// \brief: Call operator taking an Index and a FeatureVector.
    ///
    ///         Return feature of the entity, via (ACTIVE) INDEX
    ///         For general grids, the feature vector is given for the gridView_.
    ///         For CpGrid, the feature vector is given for level 0.
    /// 
    /// \tparam     FeatureType       Type of the property of the element, e.g. int, double, float, etc.
    /// \param [in] element index
    /// \param [in] feature_vec       Vector with each entry, the feature of an element of
    ///                               the gridView_ [for general grids], or level 0 for CpGrid.
    /// \return feature of the element of the gridView_ associated to the given index.
    template<typename FeatureType>
    FeatureType operator()(const int& elemIdx, const std::vector<FeatureType>& feature_vec) const;

    /// \brief: For general grids, it retunrs the same Entity index.
    ///
    /// \tparam     EntityType  Element type.
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    /// \param [in] element     EntityType object.
    /// \return element index.
    template<typename EntityType, typename GridType = Grid>
    typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int> getOriginIndexFromEntity(const EntityType& elem) const;

    /// \brief: For CpGrid, it returns index of origin cell (parent/equivalent cell when element has no father) in level 0.
    ///
    /// \tparam     EntityType  Element type.
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    /// \param [in] element     EntityType object.
    /// \return element origin index  Index of the origin cell (parent/equivalent cell when element has no father) in level 0.
    template<typename EntityType, typename GridType = Grid>
    typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int> getOriginIndexFromEntity(const EntityType& elem) const;

    /// \brief: For general grids, it retunrs the same element index.
    ///
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    /// \param [in] element index
    /// \return element index
    template<typename GridType>
    typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int> getOriginIndex(const int& elemIdx) const;

    /// \brief: For CpGrid, it returns index of origin cell (parent/equivalent cell when elem has no father) in level 0.
    ///
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    /// \param [in] element index
    /// \return element origin index  Index of the origin cell (parent/equivalent cell when element has no father) in level 0.
    template<typename GridType = Grid>
    typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int> getOriginIndex(const int& elemIdx) const;


protected:
    const GridView& gridView_;
    Dune::MultipleCodimMultipleGeomTypeMapper<GridView> elemMapper_;
}; // end LookUpData class

/// LookUpCartesianData - To search data via CartesianIndex (cartesianMapper)
///
/// Instead of using a specialitation for Dune::CpGrid, we implement std::enable_if
/// to overload methods with different definitions: for Dune:CpGrid and for other
/// Grid types. An auxiliary defualt template parameter (GridType = Grid) is added
/// to deal with the dependent names at template instantiation.
template<typename Grid, typename GridView>
class LookUpCartesianData
{
public:
    /// \brief: Constructor taking a GridView and a CartesianIndexMapper
    ///
    /// \param [in] gridView
    /// \param [in] mapper   Dune::CartesianIndexMapper<Grid>. 
    explicit LookUpCartesianData(const GridView& gridView,
                                 const Dune::CartesianIndexMapper<Grid>& mapper) :
        gridView_(gridView),
        elemMapper_(gridView, Dune::mcmgElementLayout()),
        cartMapper_(&mapper)
    {
    }

    /// \brief: Call operator taking an EntityObject and a FeatureVector.
    ///
    ///         Return feature of the entity, via CARTESIAN INDEX
    ///         For general grids, the feature vector is given for the gridView_.
    ///         For CpGrid, the feature vector is given for level 0.
    ///
    /// \tparam     EntityType  Element type.
    /// \tparam     FeatureType Type of the property of the element, e.g. int, double, float, etc.
    /// \param [in] element     EntityType object.
    /// \param [in] feature_vec Vector with each entry, the feature of an element of
    ///                         the gridView_ [for general grids], or level 0 for CpGrid.
    /// \return feature of the given element.
    template<typename EntityType, typename FeatureType>
    FeatureType operator()(const EntityType& elem,const std::vector<FeatureType>& feature_vec) const;

    /// \brief: Call operator taking an Index and a FeatureVector.
    ///
    ///         Return feature of the entity, via CARTESIAN INDEX
    ///         For general grids, the feature vector is given for the gridView_.
    ///         For CpGrid, the feature vector is given for level 0.
    ///
    /// \tparam     FeatureType     Type of the property of the element, e.g. int, double, float, etc.
    /// \param [in] element index
    /// \param [in] feature_vec     Vector with each entry, the feature of an element of
    ///                             the gridView_ [for general grids], or level 0 for CpGrid.
    /// \return feature of the element of the gridView_ associated to the given index.
    template<typename FeatureType>
    FeatureType operator()(const int& elemIdx, const std::vector<FeatureType>& feature_vec) const;

    /// \brief: For general grids, it retunrs the same Entity index.
    ///
    /// \tparam     EntityType  Element type.
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    /// \param [in] element     EntityType object.
    /// \return element index.
    template<typename EntityType, typename GridType = Grid>
    typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int> getOriginIndexFromEntity(const EntityType& elem) const;

    /// \brief: For CpGrid, it returns index of origin cell (parent/equivalent cell when element has no father) in level 0.
    ///
    /// \tparam     EntityType  Element type.
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    /// \param [in] element     EntityType object.
    /// \return element origin index  Index of the origin cell (parent/equivalent cell when element has no father) in level 0.
    template<typename EntityType, typename GridType = Grid>
    typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int> getOriginIndexFromEntity(const EntityType& elem) const;

    /// \brief: For general grids, it retunrs the same element index.
    ///
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    /// \param [in] element index
    /// \return element index
    template<typename GridType = Grid>
    typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int> getOriginIndex(const int& elemIdx) const;

    /// \brief: For CpGrid, it returns index of origin cell (parent/equivalent cell when elem has no father) in level 0.
    ///
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    /// \param [in] element index
    /// \return element origin index  Index of the origin cell (parent/equivalent cell when element has no father) in level 0.
    template<typename GridType = Grid>
    typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int> getOriginIndex(const int& elemIdx) const;

    /// \brief: It returns the Cartesian index of origin cell (parent/equivalent cell when elem has no father) in level 0.
    ///
    /// \tparam     EntityType
    /// \param [in] element
    /// \return     Cartesian Index of the origin cell (parent/equivalent cell when element has no father) in level 0.
    template<typename EntityType>
    int getCartesianOriginIdxFromEntity(const EntityType& elem) const;

    /// \brief: It returns the Cartesian index of origin cell (parent/equivalent cell when elem has no father) in level 0.
    ///
    /// \param [in] element index
    /// \return     Cartesian Index of the origin cell (parent/equivalent cell when element has no father) in level 0.
    int getCartesianOriginIndex(const int& elemIdx) const;


protected:
    const GridView& gridView_;
    Dune::MultipleCodimMultipleGeomTypeMapper<GridView> elemMapper_;
    const Dune::CartesianIndexMapper<Grid>* cartMapper_;
}; // end LookUpCartesianData class
}
// end namespace Opm



/// LookUpData

template<typename Grid, typename GridView>
template<typename EntityType, typename FeatureType>
FeatureType Opm::LookUpData<Grid,GridView>::operator()(const EntityType& elem, const std::vector<FeatureType>& feature_vec) const
{
    assert( (0 <= this->getOriginIndexFromEntity<EntityType,Grid>(elem)) &&
            (static_cast<int>(feature_vec.size()) > this->getOriginIndexFromEntity<EntityType,Grid>(elem)) );
    return feature_vec[this->getOriginIndexFromEntity<EntityType,Grid>(elem)];
}

template<typename Grid, typename GridView>
template<typename FeatureType>
FeatureType Opm::LookUpData<Grid,GridView>::operator()(const int& elemIdx, const std::vector<FeatureType>& feature_vec) const
{
    assert(0 <= this-> getOriginIndex<Grid>(elemIdx) && static_cast<int>(feature_vec.size()) > this-> getOriginIndex<Grid>(elemIdx));
    return feature_vec[getOriginIndex<Grid>(elemIdx)];
}

template<typename Grid, typename GridView>
template<typename EntityType, typename GridType>
typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpData<Grid,GridView>::getOriginIndexFromEntity(const EntityType& elem) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    return this-> elemMapper_.index(elem);
}

template<typename Grid, typename GridView>
template<typename EntityType, typename GridType>
typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpData<Grid,GridView>::getOriginIndexFromEntity(const EntityType& elem) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    return elem.getOrigin().index();
}

template<typename Grid, typename GridView>
template<typename GridType>
typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpData<Grid,GridView>::getOriginIndex(const int& elemIdx) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    return elemIdx;
}

template<typename Grid, typename GridView>
template<typename GridType>
typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpData<Grid,GridView>::getOriginIndex(const int& elemIdx) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    const auto& elem = Dune::cpgrid::Entity<0>(*(gridView_.grid().current_view_data_), elemIdx, true);
    return elem.getOrigin().index(); // getOrign() returns parent Entity or the equivalent Entity in level 0.
}



/// LookUpCartesianData

template<typename Grid, typename GridView>
template<typename EntityType, typename FeatureType>
FeatureType Opm::LookUpCartesianData<Grid,GridView>::operator()
    (const EntityType& elem, const std::vector<FeatureType>& feature_vec) const
{
    assert(cartMapper_);
    assert( (0 <= this->getOriginIndexFromEntity<EntityType,Grid>(elem)) &&
            (static_cast<int>(feature_vec.size()) > this-> getOriginIndexFromEntity<EntityType,Grid>(elem)) );
    return feature_vec[cartMapper_-> cartesianIndex(this->elemMapper_.index(elem))]; 
}

template<typename Grid, typename GridView>
template<typename FeatureType>
FeatureType Opm::LookUpCartesianData<Grid,GridView>::operator()(const int& elemIdx, const std::vector<FeatureType>& feature_vec) const
{
    assert(cartMapper_);
    assert(0 <= this->getOriginIndex<Grid>(elemIdx) &&
           static_cast<int>(feature_vec.size()) > this-> getOriginIndex<Grid>(elemIdx));
    return feature_vec[cartMapper_-> cartesianIndex(elemIdx)];
}

template<typename Grid, typename GridView>
template<typename EntityType, typename GridType>
typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpCartesianData<Grid,GridView>::getOriginIndexFromEntity(const EntityType& elem) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    return elemMapper_.index(elem);
}

template<typename Grid, typename GridView>
template<typename EntityType, typename GridType>
typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpCartesianData<Grid,GridView>::getOriginIndexFromEntity(const EntityType& elem) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    return elem.getOrigin().index();
}

template<typename Grid, typename GridView>
template<typename GridType>
typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpCartesianData<Grid,GridView>::getOriginIndex(const int& elemIdx) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    return elemIdx;
}

template<typename Grid, typename GridView>
template<typename GridType>
typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpCartesianData<Grid,GridView>::getOriginIndex(const int& elemIdx) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    const auto& elem = Dune::cpgrid::Entity<0>(*(gridView_.grid().current_view_data_), elemIdx, true);
    return elem.getOrigin().index(); // getOrign() returns parent Entity or the equivalent Entity in level 0.
}

template<typename Grid, typename GridView>
template<typename EntityType>
int Opm::LookUpCartesianData<Grid,GridView>::getCartesianOriginIdxFromEntity(const EntityType& elem) const
{
    return cartMapper_-> cartesianIndex(this->elemMapper_.index(elem));
}

template<typename Grid, typename GridView>
int Opm::LookUpCartesianData<Grid,GridView>::getCartesianOriginIndex(const int& elemIdx) const
{
    return cartMapper_-> cartesianIndex(elemIdx);
}

