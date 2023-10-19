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

#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/grid/cpgrid/Entity.hpp>

#include <string>
#include <type_traits>
#include <vector>

namespace Opm
{
class FieldPropsManager;
}

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

/// LookUpData class - To search field properties of leaf grid view elements via element/elementIndex
///
/// Instead of using a specialitation for Dune::CpGrid, we implement std::enable_if
/// to overload methods with different definitions: for Dune:CpGrid and for other
/// Grid types. An auxiliary defualt template parameter (GridType = Grid) is added
/// to deal with the dependent names at template instantiation.
template <typename Grid, typename GridView>
class LookUpData
{
public:
    /// \brief:     Constructor taking a GridView and a bool
    /// \param [in] GridView
    /// \param [in] isFieldPropInLgr   bool: default false (search field property in unrefined grid)
    ///                                      true (search field property in refined grid; LGR-id/level required)
    ///                                      Currently, isFieldPropInLgr_ == false means that all the field
    ///                                      properties are given in the unrefined grid (level 0).
    explicit LookUpData(const  GridView& gridView, bool isFieldPropInLgr = false) :
        gridView_(gridView),
        elemMapper_(gridView, Dune::mcmgElementLayout()),
        isFieldPropInLgr_(isFieldPropInLgr)
    {
    }

    /// \brief: Call operator taking an EntityObject, a FieldPropVector, and a level.
    ///
    ///         Return field property of the entity, via (ACTIVE) INDEX
    ///         For general grids, the field property vector is assumed to be given for the gridView_.
    ///         For CpGrid, the field property vector is assumed to be given for level 0 when isFieldPropLgr_ == false,
    ///         and for certain LGR/level > 0 when isFieldPropLgr_ == true.
    ///
    /// \tparam     EntityType    Element type.
    /// \tparam     FieldPropType Type of the field property of the element, e.g. int, double, float, etc.
    /// \param [in] element       EntityType object.
    /// \param [in] fieldProp     Vector with each entry, the feature of an element of
    ///                           the gridView_ [for general grids], or level 0 for CpGrid.
    /// \return field property of the given leaf grid view element.
    template<typename EntityType, typename FieldPropType>
    FieldPropType operator()(const EntityType& elem, const std::vector<FieldPropType>& fieldProp) const;

    /// \brief: Call operator taking an Index, a FieldPropertyVector, and a level.
    ///
    ///         Return field property of the entity, via (ACTIVE) INDEX
    ///         For general grids, the field property vector is assumed to be given for the gridView_.
    ///         For CpGrid, the field property vector is assumed to be given for level 0 when isFieldPropLgr_ == false,
    ///         and for certain LGR/level > 0 when isFieldPropLgr_ == true.
    ///
    /// \tparam     FieldPropType   Type of the field property of the element, e.g. int, double, float, etc.
    /// \param [in] element index
    /// \param [in] fieldProp       Vector with each entry, the feature of an element of
    ///                             the gridView_ [for general grids], or level 0 for CpGrid.
    /// \return field property of the given leaf grid view element.
    template<typename FieldPropType>
    FieldPropType operator()(const int& elemIdx, const std::vector<FieldPropType>& fieldProp) const;

    /// \brief: fieldPropDouble() taking a FieldPropsManager, std::string, and an element.
    ///                           Return field property of type double of the entity, via (ACTIVE) INDEX
    ///
    /// \tparam     EntityType          Element type.
    /// \param [in] FieldPropsManager
    /// \param [in] propString          String to lookup the field property in the FieldPropsManager
    /// \param [in] element             EntityType object.
    /// \return field property of the given leaf grid view element.
    template<typename EntityType>
    double fieldPropDoubleFromEntity(const FieldPropsManager& fieldPropsManager,
                                     const std::string& propString,
                                     const EntityType& elem) const;


    /// \brief: fieldPropDouble() taking a FieldPropsManager, std::string, and an element index.
    ///                           Return field property of type double of the entity, via (ACTIVE) INDEX
    ///
    /// \tparam     EntityType          Element type.
    /// \param [in] FieldPropsManager
    /// \param [in] propString          String to lookup the field property in the FieldPropsManager
    /// \param [in] elemIdx             Element index
    /// \return field property of the given leaf grid view element.
    double fieldPropDouble(const FieldPropsManager& fieldPropsManager,
                           const std::string& propString,
                           const int& elemIdx) const;

    /// \brief: fieldPropInt() taking a FieldPropsManager, std::string, and an element.
    ///                        Return field property of type int of the entity, via (ACTIVE) INDEX
    ///
    /// \tparam     EntityType          Element type.
    /// \param [in] FieldPropsManager
    /// \param [in] propString          String to lookup the field property in the FieldPropsManager
    /// \param [in] element             EntityType object.
    /// \return field property of the given leaf grid view element.
    template<typename EntityType>
    int fieldPropIntFromEntity(const FieldPropsManager& fieldPropsManager,
                               const std::string& propString,
                               const EntityType& elem) const;

    /// \brief: fieldPropDouble() taking a FieldPropsManager, std::string, and an element index.
    ///                           Return field property of type double of the entity, via (ACTIVE) INDEX
    ///
    /// \tparam     EntityType          Element type.
    /// \param [in] FieldPropsManager
    /// \param [in] propString          String to lookup the field property in the FieldPropsManager
    /// \param [in] elemIdx             Element index
    /// \return field property of the given leaf grid view element.
    int fieldPropInt(const FieldPropsManager& fieldPropsManager,
                     const std::string& propString,
                     const int& elemIdx) const;

    /// \brief: For general grids, it retunrs the same Entity index.
    ///
    /// \tparam     EntityType  Element type.
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    /// \param [in] element     EntityType object.
    /// \return fieldPropIdx    Element Index
    template<typename EntityType, typename GridType = Grid>
    typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int>
    getFieldPropIdxFromEntity(const EntityType& elem) const;

    /// \brief: For CpGrid, it returns index to search for the field propertries.
    ///
    ///         When isFieldPropInLgr_ == false : fieldPropIdx == Index of the origin cell (parent/equivalent cell when element
    ///                                           has nofather) in level 0.
    ///         When isFieldPropInLgr_ == true  : fieldPropIdx == Index of the equivalent cell in the LGR (level>0).
    ///
    /// \tparam     EntityType  Element type.
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    /// \param [in] element     EntityType object.
    /// \return fieldPropIdx
    template<typename EntityType, typename GridType = Grid>
    typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int>
    getFieldPropIdxFromEntity(const EntityType& elem) const;

    /// \brief: For general grids, it retunrs the same element index.
    ///
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    /// \param [in] elemIdx
    /// \return fieldPropIdx    Element index.
    template<typename GridType>
    typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int>
    getFieldPropIdx(const int& elemIdx) const;

    /// \brief: For CpGrid, it returns the index to search for the field properties.
    ///
    ///         When isFieldPropInLgr_ == false : fieldPropIdx == Index of the origin cell (parent/equivalent cell when element
    ///                                           has nofather) in level 0.
    ///         When isFieldPropInLgr_ == true  : fieldPropIdx == Index of the equivalent cell in the LGR (level>0).
    ///
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    /// \param [in] element index
    /// \return fieldPropIdx
    template<typename GridType = Grid>
    typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int>
    getFieldPropIdx(const int& elemIdx) const;


protected:
    const GridView& gridView_;
    Dune::MultipleCodimMultipleGeomTypeMapper<GridView> elemMapper_;
    bool isFieldPropInLgr_;
}; // end LookUpData class

/// LookUpCartesianData - To search field properties of leaf grid view elements via CartesianIndex (cartesianMapper)
///
/// Instead of using a specialitation for Dune::CpGrid, we implement std::enable_if
/// to overload methods with different definitions: for Dune:CpGrid and for other
/// Grid types. An auxiliary defualt template parameter (GridType = Grid) is added
/// to deal with the dependent names at template instantiation.
template<typename Grid, typename GridView>
class LookUpCartesianData
{
public:
    /// \brief: Constructor taking a GridView, a CartesianIndexMapper, and a bool
    ///
    /// \param [in] gridView
    /// \param [in] mapper   Dune::CartesianIndexMapper<Grid>.
    /// \param [in] isFieldPropInLgr bool: default false (search field property in unrefined grid)
    ///                                    true (search field property in refined grid; LGR-id/level required)
    ///                                    Currently, isFieldPropInLgr_ == false means that all the field
    ///                                    properties are given in the unrefined grid (level 0).
    explicit LookUpCartesianData(const GridView& gridView,
                                 const Dune::CartesianIndexMapper<Grid>& mapper,
                                 bool isFieldPropInLgr = false) :
        gridView_(gridView),
        elemMapper_(gridView, Dune::mcmgElementLayout()),
        cartMapper_(&mapper),
        isFieldPropInLgr_(isFieldPropInLgr)
    {
    }

    /// \brief: Call operator taking an EntityObject, a FieldPropVector, and a level.
    ///
    ///         Return field property of the entity, via CARTESIAN INDEX
    ///         For general grids, the field property vector is assumed to be given for the gridView_.
    ///         For CpGrid, the field property vector is assumed to be given for level 0 when isFieldPropLgr_ == false,
    ///         and for certain LGR/level > 0 when isFieldPropLgr_ == true.
    ///
    /// \tparam     EntityType    Element type.
    /// \tparam     FieldPropType Type of the field property of the element, e.g. int, double, float, etc.
    /// \param [in] element       EntityType object.
    /// \param [in] fieldProp     Vector with each entry, the feature of an element of
    ///                           the gridView_ [for general grids], or level 0 for CpGrid.
    /// \return field property of the given leaf grid view element.
    template<typename EntityType, typename FieldPropType>
    FieldPropType operator()(const EntityType& elem,const std::vector<FieldPropType>& fieldProp) const;

    /// \brief: Call operator taking an Index and a FeatureVector. a FieldPropertyVector, and a level.
    ///
    ///         Return field property of the entity, via CARTESIAN INDEX
    ///         For general grids, the field property vector is assumed to be given for the gridView_.
    ///         For CpGrid, the field property vector is assumed to be given for level 0 when isFieldPropLgr_ == false,
    ///         and for certain LGR/level > 0 when isFieldPropLgr_ == true.
    ///
    /// \tparam     FieldPropType   Type of the field property of the element, e.g. int, double, float, etc.
    /// \param [in] element index
    /// \param [in] fieldProp       Vector with each entry, the feature of an element of
    ///                             the gridView_ [for general grids], or level 0 for CpGrid.
    /// \return field property of the given leaf grid view element.
    template<typename FieldPropType>
    FieldPropType operator()(const int& elemIdx, const std::vector<FieldPropType>& fieldProp) const;

    /// \brief: fieldPropDouble() taking a FieldPropsManager, std::string, and an element.
    ///                           Return field property of type double of the entity, via CARTESIAN INDEX
    ///
    /// \tparam     EntityType          Element type.
    /// \param [in] FieldPropsManager
    /// \param [in] propString          String to lookup the field property in the FieldPropsManager
    /// \param [in] element             EntityType object.
    /// \return field property of the given leaf grid view element.
    template<typename EntityType>
    double fieldPropDoubleFromEntity(const FieldPropsManager& fieldPropsManager,
                                     const std::string& propString,
                                     const EntityType& elem) const;

    /// \brief: fieldPropDouble() taking a FieldPropsManager, std::string, and an element index.
    ///                           Return field property of type double of the entity, via (ACTIVE) INDEX
    ///
    /// \tparam     EntityType          Element type.
    /// \param [in] FieldPropsManager
    /// \param [in] propString          String to lookup the field property in the FieldPropsManager
    /// \param [in] elemIdx             Element index
    /// \return field property of the given leaf grid view element.
    double fieldPropDouble(const FieldPropsManager& fieldPropsManager,
                           const std::string& propString,
                           const int& elemIdx) const;

    /// \brief: fieldPropInt() taking a FieldPropsManager, std::string, and an element.
    ///                        Return field property of type int of the entity, via CARTESIAN INDEX
    ///
    /// \tparam     EntityType          Element type.
    /// \param [in] FieldPropsManager
    /// \param [in] propString          String to lookup the field property in the FieldPropsManager
    /// \param [in] element             EntityType object.
    /// \return field property of the given leaf grid view element.
    template<typename EntityType>
    int fieldPropIntFromEntity(const FieldPropsManager& fieldPropsManager,
                               const std::string& propString,
                               const EntityType& elem) const;

    /// \brief: fieldPropDouble() taking a FieldPropsManager, std::string, and an element index.
    ///                           Return field property of type double of the entity, via CARTESIAN INDEX
    ///
    /// \tparam     EntityType          Element type.
    /// \param [in] FieldPropsManager
    /// \param [in] propString          String to lookup the field property in the FieldPropsManager
    /// \param [in] elemIdx             Element index
    /// \return field property of the given leaf grid view element.
    int fieldPropInt(const FieldPropsManager& fieldPropsManager,
                     const std::string& propString,
                     const int& elemIdx) const;

    /// \brief: For general grids, it retunrs the same Entity index.
    ///
    /// \tparam     EntityType  Element type.
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    /// \param [in] element     EntityType object.
    /// \return fieldPropIdx    Element Index
    template<typename EntityType, typename GridType = Grid>
    typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int>
    getFieldPropCartesianIdxFromEntity(const EntityType& elem) const;

    /// \brief: For CpGrid, it returns index to search for the field propertries.
    ///
    ///         When isFieldPropInLgr_ == false : fieldPropIdx == Index of the origin cell (parent/equivalent cell when element
    ///                                           has nofather) in level 0.
    ///         When isFieldPropInLgr_ == true  : fieldPropIdx == Index of the equivalent cell in the LGR (level>0).
    ///
    /// \tparam     EntityType  Element type.
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    /// \param [in] element     EntityType object.
    /// \return fieldPropIdx
    template<typename EntityType, typename GridType = Grid>
    typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int>
    getFieldPropCartesianIdxFromEntity(const EntityType& elem) const;

    /// \brief: For general grids, it retunrs the same element index.
    ///
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    /// \param [in] elemIdx
    /// \return fieldPropIdx    Element index.
    template<typename GridType = Grid>
    typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int>
    getFieldPropCartesianIdx(const int& elemIdx) const;

    /// \brief: For CpGrid, it returns the index to search for the field properties.
    ///
    ///         When isFieldPropInLgr_ == false : fieldPropIdx == Index of the origin cell (parent/equivalent cell when element
    ///                                           has nofather) in level 0.
    ///         When isFieldPropInLgr_ == true  : fieldPropIdx == Index of the equivalent cell in the LGR (level>0).
    ///
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    /// \param [in] element index
    /// \return fieldPropIdx
    template<typename GridType = Grid>
    typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int>
    getFieldPropCartesianIdx(const int& elemIdx) const;

protected:
    const GridView& gridView_;
    Dune::MultipleCodimMultipleGeomTypeMapper<GridView> elemMapper_;
    const Dune::CartesianIndexMapper<Grid>* cartMapper_;
    bool isFieldPropInLgr_;
}; // end LookUpCartesianData class
}
// end namespace Opm



/// LookUpData

template<typename Grid, typename GridView>
template<typename EntityType, typename FieldPropType>
FieldPropType Opm::LookUpData<Grid,GridView>::operator()(const EntityType& elem,
                                                         const std::vector<FieldPropType>& fieldProp) const
{
    const auto& fieldPropIdx = this->getFieldPropIdxFromEntity<EntityType,Grid>(elem);
    assert( (0 <= fieldPropIdx) && (static_cast<int>(fieldProp.size()) > fieldPropIdx));
    return fieldProp[fieldPropIdx];
}

template<typename Grid, typename GridView>
template<typename FieldPropType>
FieldPropType Opm::LookUpData<Grid,GridView>::operator()(const int& elemIdx,
                                                         const std::vector<FieldPropType>& fieldProp) const
{
    const auto& fieldPropIdx = this->getFieldPropIdx<Grid>(elemIdx);
    assert(0 <= fieldPropIdx && static_cast<int>(fieldProp.size()) > fieldPropIdx);
    return fieldProp[fieldPropIdx];
}

template<typename Grid, typename GridView>
template<typename EntityType>
double Opm::LookUpData<Grid,GridView>::fieldPropDoubleFromEntity(const FieldPropsManager& fieldPropsManager,
                                                                 const std::string& propString,
                                                                 const EntityType& elem) const
{
    const auto& fieldPropVec = fieldPropsManager.get_double(propString);
    return this ->operator()(elem,fieldPropVec);
}

template<typename Grid, typename GridView>
double Opm::LookUpData<Grid,GridView>::fieldPropDouble(const FieldPropsManager& fieldPropsManager,
                                                       const std::string& propString,
                                                       const int& elemIdx) const
{
    const auto& fieldPropVec = fieldPropsManager.get_double(propString);
    return this ->operator()(elemIdx,fieldPropVec);
}

template<typename Grid, typename GridView>
template<typename EntityType>
int Opm::LookUpData<Grid,GridView>::fieldPropIntFromEntity(const FieldPropsManager& fieldPropsManager,
                                                           const std::string& propString,
                                                           const EntityType& elem) const
{
    const auto& fieldPropVec = fieldPropsManager.get_int(propString);
    return this ->operator()(elem,fieldPropVec);
}

template<typename Grid, typename GridView>
int Opm::LookUpData<Grid,GridView>::fieldPropInt(const FieldPropsManager& fieldPropsManager,
                                                 const std::string& propString,
                                                 const int& elemIdx) const
{
    const auto& fieldPropVec = fieldPropsManager.get_int(propString);
    return this ->operator()(elemIdx,fieldPropVec);
}

template<typename Grid, typename GridView>
template<typename EntityType, typename GridType>
typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpData<Grid,GridView>::getFieldPropIdxFromEntity(const EntityType& elem) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    assert(elem.level() == 0); // LGRs (level>0) only supported for CpGrid.
    return this-> elemMapper_.index(elem);
}

template<typename Grid, typename GridView>
template<typename EntityType, typename GridType>
typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpData<Grid,GridView>::getFieldPropIdxFromEntity(const EntityType& elem) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    if (isFieldPropInLgr_ && elem.level()) { // level > 0 == true ; level == 0 == false
        // In case some LGRs do not have refined field properties, the next line need to be modified.
        return elem.getLgrElem().index(); // getLgrElem() returns equivalent refined Entity in the LGR the leafGridView cell was born.
    }
    else {
        return elem.getOrigin().index(); // getOrign() returns parent Entity or the equivalent Entity in level 0.
    }
}

template<typename Grid, typename GridView>
template<typename GridType>
typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpData<Grid,GridView>::getFieldPropIdx(const int& elemIdx) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    // Check there are no LGRs. LGRs (level>0) only supported for CpGrid.
    assert(gridView_.grid().maxLevel() == 0);
    return elemIdx;
}

template<typename Grid, typename GridView>
template<typename GridType>
typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpData<Grid,GridView>::getFieldPropIdx(const int& elemIdx) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    const auto& elem = Dune::cpgrid::Entity<0>(*(gridView_.grid().current_view_data_), elemIdx, true);
    if (isFieldPropInLgr_ && elem.level()) { // level > 0 == true ; level == 0 == false
        // In case some LGRs do not have refined field properties, the next line need to be modified.
        return elem.getLgrElem().index(); // getLgrElem() returns equivalent refined Entity in the LGR the leafGridView cell was born.
    }
    else {
        return elem.getOrigin().index(); // getOrign() returns parent Entity or the equivalent Entity in level 0.
    }
}



/// LookUpCartesianData

template<typename Grid, typename GridView>
template<typename EntityType, typename FieldPropType>
FieldPropType Opm::LookUpCartesianData<Grid,GridView>::operator()(const EntityType& elem,
                                                                  const std::vector<FieldPropType>& fieldProp) const
{
    assert(cartMapper_);
    const auto fieldPropCartIdx = this->getFieldPropCartesianIdxFromEntity<EntityType,Grid>(elem);
    assert( (0 <= fieldPropCartIdx) && (static_cast<int>(fieldProp.size()) > fieldPropCartIdx) );
    return fieldProp[fieldPropCartIdx];
}

template<typename Grid, typename GridView>
template<typename FieldPropType>
FieldPropType Opm::LookUpCartesianData<Grid,GridView>::operator()(const int& elemIdx,
                                                                  const std::vector<FieldPropType>& fieldProp) const
{
    assert(cartMapper_);
    const auto fieldPropCartIdx = this->getFieldPropCartesianIdx<Grid>(elemIdx);
    assert(0 <=  fieldPropCartIdx && (static_cast<int>(fieldProp.size()) > fieldPropCartIdx));
    return fieldProp[fieldPropCartIdx];
}



template<typename Grid, typename GridView>
template<typename EntityType>
double Opm::LookUpCartesianData<Grid,GridView>::fieldPropDoubleFromEntity(const FieldPropsManager& fieldPropsManager,
                                                                          const std::string& propString,
                                                                          const EntityType& elem) const
{
    const auto& fieldPropVec = fieldPropsManager.get_double(propString);
    return this ->operator()(elem,fieldPropVec);
}

template<typename Grid, typename GridView>
double Opm::LookUpCartesianData<Grid,GridView>::fieldPropDouble(const FieldPropsManager& fieldPropsManager,
                                                                const std::string& propString,
                                                                const int& elemIdx) const
{
    const auto& fieldPropVec = fieldPropsManager.get_double(propString);
    return this ->operator()(elemIdx,fieldPropVec);
}

template<typename Grid, typename GridView>
template<typename EntityType>
int Opm::LookUpCartesianData<Grid,GridView>::fieldPropIntFromEntity(const FieldPropsManager& fieldPropsManager,
                                                                    const std::string& propString,
                                                                    const EntityType& elem) const
{
    const auto& fieldPropVec = fieldPropsManager.get_int(propString);
    return this ->operator()(elem,fieldPropVec);
}

template<typename Grid, typename GridView>
int Opm::LookUpCartesianData<Grid,GridView>::fieldPropInt(const FieldPropsManager& fieldPropsManager,
                                                          const std::string& propString,
                                                          const int& elemIdx) const
{
    const auto& fieldPropVec = fieldPropsManager.get_int(propString);
    return this ->operator()(elemIdx,fieldPropVec);
}


template<typename Grid, typename GridView>
template<typename EntityType, typename GridType>
typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpCartesianData<Grid,GridView>::getFieldPropCartesianIdxFromEntity(const EntityType& elem) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    // Check there are no LGRs. LGRs (level>0) only supported for CpGrid.
    assert(gridView_.grid().maxLevel() == 0);
    return cartMapper_-> cartesianIndex(this->elemMapper_.index(elem));
}

template<typename Grid, typename GridView>
template<typename EntityType, typename GridType>
typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpCartesianData<Grid,GridView>::getFieldPropCartesianIdxFromEntity(const EntityType& elem) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    if (isFieldPropInLgr_ && elem.level()) { // level == 0 false; level > 0 true
        return elem.getLgrCartesianIdx();
    }
    else {
        return cartMapper_-> cartesianIndex(this->elemMapper_.index(elem));
    }
}

template<typename Grid, typename GridView>
template<typename GridType>
typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpCartesianData<Grid,GridView>::getFieldPropCartesianIdx(const int& elemIdx) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    return  cartMapper_-> cartesianIndex(elemIdx);
}

template<typename Grid, typename GridView>
template<typename GridType>
typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpCartesianData<Grid,GridView>::getFieldPropCartesianIdx(const int& elemIdx) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    const auto& elem = Dune::cpgrid::Entity<0>(*(gridView_.grid().current_view_data_), elemIdx, true);
    return this -> getFieldPropCartesianIdxFromEntity<Dune::cpgrid::Entity<0>,Dune::CpGrid>(elem);
}
