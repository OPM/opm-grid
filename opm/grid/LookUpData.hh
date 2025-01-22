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
#ifndef OPM_LOOKUPDATA_HH
#define OPM_LOOKUPDATA_HH

#include <dune/grid/common/mcmgmapper.hh>

#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/grid/cpgrid/Entity.hpp>

#include <functional>
#include <string>
#include <type_traits>
#include <vector>

namespace Dune
{
class CpGrid;
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

    /// \brief: Get field propertry for an element or index in the leaf grid view, from a vector and element index.
    ///
    ///         For general grids, the field property vector is assumed to be given for the gridView_.
    ///         For CpGrid, the field property vector is assumed to be given for level 0 when isFieldPropLgr_ == false,
    ///         and for certain LGR/level > 0 when isFieldPropLgr_ == true.
    ///
    /// \param [in] elementOrIdx  Element or index of the element in the leaf grid view.
    /// \param [in] fieldProp     Vector (indexable collection) of field properties.
    auto operator()(const auto& elementOrIdx, const auto& fieldProp) const;

    /// \brief: Get field property of type double from field properties manager by name.
    std::vector<double> assignFieldPropsDoubleOnLeaf(const FieldPropsManager& fieldPropsManager,
                                                     const std::string& propString) const;

    /// \brief: Get field property of type int from field properties manager by name.
    template<typename IntType>
    std::vector<IntType> assignFieldPropsIntOnLeaf(const FieldPropsManager& fieldPropsManager,
                                                   const std::string& propString,
                                                   const bool& needsTranslation,
                                                   std::function<void(IntType, int)> valueCheck = [](IntType, int){}) const;

    /// \brief: Get property of type double from field properties manager by name, via element or its index.
    template<typename ElemOrIndex>
    double fieldPropDouble(const FieldPropsManager& fieldPropsManager,
                           const std::string& propString,
                           const ElemOrIndex& elemOrIndex) const;

    /// \brief: Get property of type int from field properties manager by name, via element.
    template<typename ElemOrIndex>
    int fieldPropInt(const FieldPropsManager& fieldPropsManager,
                     const std::string& propString,
                     const ElemOrIndex& elemOrIndex) const;

    /// \brief: Return the same element index for all grids different from CpGrid.
    ///
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    template<typename GridType = Grid>
    auto getFieldPropIdx(const auto& elem) const;

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

    /// \brief: Get field property for an element in the leaf grid view, from a vector, via Cartesian Index.
    ///
    ///         For general grids, the field property vector is assumed to be given for the gridView_.
    ///         For CpGrid, the field property vector is assumed to be given for level 0 when isFieldPropLgr_ == false,
    ///         and for certain LGR/level > 0 when isFieldPropLgr_ == true.
    auto operator()(const auto& elemIdx, const auto& fieldProp) const;


    /// \brief: Get field property of type double from field properties manager by name.
    std::vector<double> assignFieldPropsDoubleOnLeaf(const FieldPropsManager& fieldPropsManager,
                                                     const std::string& propString) const;

    /// \brief: Get field property of type int from field properties manager by name.
    template<typename IntType>
    std::vector<IntType> assignFieldPropsIntOnLeaf(const FieldPropsManager& fieldPropsManager,
                                                   const std::string& propString,
                                                   const bool& needsTranslation,
                                                   std::function<void(IntType, int)> valueCheck = [](IntType, int){}) const;

    /// \brief: Get property of type double from field properties manager by name, via element or its index.
    template<typename ElemOrIndex>
    double fieldPropDouble(const FieldPropsManager& fieldPropsManager,
                           const std::string& propString,
                           const ElemOrIndex& elemOrIndex) const;

    /// \brief: Get property of type int from field properties manager by name, via element or its index.
    template<typename ElemOrIndex>
    int fieldPropInt(const FieldPropsManager& fieldPropsManager,
                     const std::string& propString,
                     const ElemOrIndex& elemOrIndex) const;

    template<typename GridType = Grid>
    auto getFieldPropCartesianIdx(const auto& elemIdx) const;

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
auto Opm::LookUpData<Grid,GridView>::operator()(const auto& elemIdx,
                                                const auto& fieldProp) const
{
    const auto& fieldPropIdx = this->getFieldPropIdx<Grid>(elemIdx);
    assert(0 <= fieldPropIdx && static_cast<int>(fieldProp.size()) > fieldPropIdx);
    return fieldProp[fieldPropIdx];
}

template<typename Grid, typename GridView>
std::vector<double> Opm::LookUpData<Grid,GridView>::assignFieldPropsDoubleOnLeaf(const FieldPropsManager& fieldPropsManager,
                                                                                 const std::string& propString) const
{
    using IndexType = typename Dune::MultipleCodimMultipleGeomTypeMapper<GridView>::Index;

    std::vector<double> fieldPropOnLeaf;
    unsigned int numElements = gridView_.size(0);
    fieldPropOnLeaf.resize(numElements);
    const auto& fieldProp = fieldPropsManager.get_double(propString);
    if ( (propString == "PORV") && (gridView_.grid().maxLevel() > 0)) {
        // PORV poreVolume. LGRs supported (so far) only for CpGrid.
        // For CpGrid with LGRs, poreVolume of a cell on the leaf grid view which has a parent cell on level 0,
        // is computed as  porv[parent] * leafCellVolume / parentCellVolume. In this way, the sum of the pore
        // volume of a parent cell coincides with the sum of the pore volume of its children.
        for (const auto& element : elements(gridView_)) {
            const auto& elemIdx = this-> elemMapper_.index(element);
            const auto& fieldPropIdx = this->getFieldPropIdx<Grid>(elemIdx); // gets parentIdx (or (lgr)levelIdx) for CpGrid with LGRs
            if (element.hasFather()) {
                const auto fatherVolume = element.father().geometry().volume();
                const auto& elemVolume = element.geometry().volume();
                fieldPropOnLeaf[elemIdx] = fieldProp[fieldPropIdx] * elemVolume / fatherVolume;
            }
            else {
                fieldPropOnLeaf[elemIdx] = fieldProp[fieldPropIdx];
            }
        }
    }
    else {
        for (const auto& element : elements(gridView_)) {
            const auto& elemIdx = this-> elemMapper_.index(element);
            const auto& fieldPropIdx = this->getFieldPropIdx<Grid>(elemIdx); // gets parentIdx (or (lgr)levelIdx) for CpGrid with LGRs
            fieldPropOnLeaf[elemIdx] = fieldProp[fieldPropIdx];
        }
    }
    return fieldPropOnLeaf;
}

template<typename Grid, typename GridView>
template<typename IntType>
std::vector<IntType> Opm::LookUpData<Grid,GridView>::assignFieldPropsIntOnLeaf(const FieldPropsManager& fieldPropsManager,
                                                                               const std::string& propString,
                                                                               const bool& needsTranslation,
                                                                               std::function<void(IntType, int)> valueCheck) const
{
    using IndexType = typename Dune::MultipleCodimMultipleGeomTypeMapper<GridView>::Index;
    std::vector<IntType> fieldPropOnLeaf;
    unsigned int numElements = gridView_.size(0);
    fieldPropOnLeaf.resize(numElements);
    const auto& fieldProp = fieldPropsManager.get_int(propString);
    for (const auto& element : elements(gridView_)) {
        const auto& elemIdx = this-> elemMapper_.index(element);
        const auto& fieldPropIdx = this->getFieldPropIdx<Grid>(elemIdx); // gets parentIdx (or (lgr)levelIdx) for CpGrid with LGRs
        fieldPropOnLeaf[elemIdx] = fieldProp[fieldPropIdx] - needsTranslation;
        valueCheck(fieldProp[fieldPropIdx], fieldPropIdx);
    }
    return fieldPropOnLeaf;
}

template<typename Grid, typename GridView>
template<typename ElemOrIndex>
double Opm::LookUpData<Grid,GridView>::fieldPropDouble(const FieldPropsManager& fieldPropsManager,
                                                       const std::string& propString,
                                                       const ElemOrIndex& elemOrIndex) const
{
    const auto& fieldPropVec = fieldPropsManager.get_double(propString);
    return this ->operator()(elemOrIndex,fieldPropVec);
}

template<typename Grid, typename GridView>
template<typename ElemOrIndex>
int Opm::LookUpData<Grid,GridView>::fieldPropInt(const FieldPropsManager& fieldPropsManager,
                                                 const std::string& propString,
                                                 const ElemOrIndex& elemOrIndex) const
{
    const auto& fieldPropVec = fieldPropsManager.get_int(propString);
    return this ->operator()(elemOrIndex,fieldPropVec);
}

template<typename Grid, typename GridView>
template<typename GridType>
auto Opm::LookUpData<Grid,GridView>::getFieldPropIdx(const auto& elementOrIndex) const
{
    using IndexType = std::remove_const_t<std::remove_reference_t<decltype(elementOrIndex)>>;
    constexpr static bool isIntegral = std::is_integral_v<IndexType>;
    if constexpr (std::is_same_v<GridType, Dune::CpGrid>) {
        if constexpr (isIntegral) {
            static_assert(std::is_same_v<Grid,GridType>);
            const auto& elem = Dune::cpgrid::Entity<0>(*(gridView_.grid().currentData().back()), elementOrIndex, true);
            if (isFieldPropInLgr_ && elem.level()) { // level > 0 == true ; level == 0 == false
                // In case some LGRs do not have refined field properties, the next line need to be modified.
                return elem.getLevelElem().index();
            }
            else {
                return elem.getOrigin().index();
            }
        } else {
            static_assert(std::is_same_v<Grid,GridType>);
            static_assert(std::is_same_v<IndexType, Dune::cpgrid::Entity<0>>);
            if (isFieldPropInLgr_ && elementOrIndex.level()) { // level > 0 == true ; level == 0 == false
                // In case some LGRs do not have refined field properties, the next line need to be modified.
                return elementOrIndex.getLevelElem().index();
            }
            else {
                return elementOrIndex.getOrigin().index();
            }
        }
    } else {
        if constexpr (isIntegral) {
            static_assert(std::is_same_v<Grid,GridType>);
            // Check there are no LGRs. LGRs (level>0) only supported for CpGrid.
            assert(gridView_.grid().maxLevel() == 0);
            return elementOrIndex;
        } else {
            static_assert(std::is_same_v<Grid,GridType>);
            assert(elementOrIndex.level() == 0); // LGRs (level>0) only supported for CpGrid.
            return this-> elemMapper_.index(elementOrIndex);
        }
    }
    
}

/// LookUpCartesianData

template<typename Grid, typename GridView>
//template<typename IdType, typename FieldPropType>
//std::enable_if_t<std::is_integral_v<IdType>, FieldPropType>
auto Opm::LookUpCartesianData<Grid,GridView>::operator()(const auto& elementOrIndex,
                                                         const auto& fieldProp) const
{
    using IndexType = std::remove_const_t<std::remove_reference_t<decltype(elementOrIndex)>>;
    constexpr static bool isIntegral = std::is_integral_v<IndexType>;
    if constexpr (isIntegral) {
        assert(cartMapper_);
        const auto fieldPropCartIdx = this->getFieldPropCartesianIdx<Grid>(elementOrIndex);
        assert(0 <=  fieldPropCartIdx && (static_cast<int>(fieldProp.size()) > fieldPropCartIdx));
        return fieldProp[fieldPropCartIdx];
    } else {
        assert(cartMapper_);
        const auto fieldPropCartIdx = this->getFieldPropCartesianIdx<Grid>(elementOrIndex);
        assert( (0 <= fieldPropCartIdx) && (static_cast<int>(fieldProp.size()) > fieldPropCartIdx) );
        return fieldProp[fieldPropCartIdx];
    }    
}

template<typename Grid, typename GridView>
std::vector<double> Opm::LookUpCartesianData<Grid,GridView>::assignFieldPropsDoubleOnLeaf(const FieldPropsManager& fieldPropsManager,
                                                                                          const std::string& propString) const
{
    std::vector<double> fieldPropOnLeaf;
    unsigned int numElements = gridView_.size(0);
    fieldPropOnLeaf.resize(numElements);
    const auto& fieldProp = fieldPropsManager.get_double(propString);
    for (unsigned int elemIdx = 0; elemIdx < numElements; ++elemIdx) {
        const auto fieldPropCartIdx = this->getFieldPropCartesianIdx<Grid>(elemIdx);
        fieldPropOnLeaf[elemIdx] = fieldProp[fieldPropCartIdx];
    }
    return fieldPropOnLeaf;
}

template<typename Grid, typename GridView>
template<typename IntType>
std::vector<IntType> Opm::LookUpCartesianData<Grid,GridView>::assignFieldPropsIntOnLeaf(const FieldPropsManager& fieldPropsManager,
                                                                                        const std::string& propString,
                                                                                        const bool& needsTranslation,
                                                                                        std::function<void(IntType, int)> valueCheck) const
{
    std::vector<IntType> fieldPropOnLeaf;
    unsigned int numElements = gridView_.size(0);
    fieldPropOnLeaf.resize(numElements);
    const auto& fieldProp = fieldPropsManager.get_int(propString);
    for (unsigned elemIdx = 0; elemIdx < numElements; ++elemIdx) {
        const auto fieldPropCartIdx = this->getFieldPropCartesianIdx<Grid>(elemIdx);
        fieldPropOnLeaf[elemIdx] = fieldProp[fieldPropCartIdx] - needsTranslation;
        valueCheck(fieldProp[fieldPropCartIdx], fieldPropCartIdx);
    }
    return fieldPropOnLeaf;
}

template<typename Grid, typename GridView>
template<typename ElemOrIndex>
double Opm::LookUpCartesianData<Grid,GridView>::fieldPropDouble(const FieldPropsManager& fieldPropsManager,
                                                                const std::string& propString,
                                                                const ElemOrIndex& elemOrIndex) const
{
    const auto& fieldPropVec = fieldPropsManager.get_double(propString);
    return this ->operator()(elemOrIndex,fieldPropVec);
}

template<typename Grid, typename GridView>
template<typename ElemOrIndex>
int Opm::LookUpCartesianData<Grid,GridView>::fieldPropInt(const FieldPropsManager& fieldPropsManager,
                                                          const std::string& propString,
                                                          const ElemOrIndex& elemOrIndex) const
{
    const auto& fieldPropVec = fieldPropsManager.get_int(propString);
    return this ->operator()(elemOrIndex,fieldPropVec);
}

template<typename Grid, typename GridView>
template<typename GridType>
auto Opm::LookUpCartesianData<Grid,GridView>::getFieldPropCartesianIdx(const auto& elementOrIndex) const
{
    using IndexType = std::remove_const_t<std::remove_reference_t<decltype(elementOrIndex)>>;
    constexpr static bool isIntegral = std::is_integral_v<IndexType>;
    if constexpr (std::is_same_v<GridType, Dune::CpGrid>) {
        if constexpr (isIntegral) {
            static_assert(std::is_same_v<Grid,GridType>);
            const auto& elem = Dune::cpgrid::Entity<0>(*(gridView_.grid().currentData().back()), elementOrIndex, true);
            return this -> getFieldPropCartesianIdx<Dune::CpGrid>(elem);
        } else {
            static_assert(std::is_same_v<Grid,GridType>);
            if (isFieldPropInLgr_ && elementOrIndex.level()) { // level == 0 false; level > 0 true
                return elementOrIndex.getLevelCartesianIdx();
            }
            else {
                return cartMapper_-> cartesianIndex(this->elemMapper_.index(elementOrIndex));
            }
        }
    } else {
        if constexpr (isIntegral) {
            static_assert(std::is_same_v<Grid,GridType>);
            return  cartMapper_-> cartesianIndex(elementOrIndex);
        } else {
            static_assert(std::is_same_v<Grid,GridType>);
            // Check there are no LGRs. LGRs (level>0) only supported for CpGrid.
            assert(gridView_.grid().maxLevel() == 0);
            return cartMapper_-> cartesianIndex(this->elemMapper_.index(elementOrIndex));            
        }
    }
    
}
#endif
