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
    ///
    /// \tparam ElementOrIndex the type of the element or index passed in.
    /// \tparam FieldProperties the type of the field properties vector. Should be an \c std::vector like type.
    template<class ElementOrIndex, class FieldProperties>
    auto operator()(const ElementOrIndex& elementOrIdx, const FieldProperties& fieldProp) const;

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

    /// @brief calls \c getFieldPropIdx<Grid, ElementType>(elem) 
    template<typename ElementType>
    auto getFieldPropIdx(const ElementType& elem) const;

    /// \brief Return the index used for retrieving field properties, depending on whether
    ///        the grid is a CpGrid or a general grid, and whether the element is in an LGR.
    ///
    /// This method determines the appropriate index under two main conditions:
    ///
    /// - **Non-CpGrid**: Returns the same element index that was passed in. If an entity is
    ///   provided, it uses \c elemMapper_ to retrieve the index. The function asserts
    ///   that \c maxLevel() == 0 (i.e., no LGR support) for non-CpGrid grids.
    ///
    /// - **CpGrid**: Depending on whether \c isFieldPropInLgr_ is \c true and the element
    ///   (or index) is in a refined level (> 0), the returned index is for the equivalent
    ///   LGR cell; otherwise, the index of the origin (level 0) cell is returned.
    ///
    /// \tparam GridType
    ///     Auxiliary type used to specialize the method for CpGrid
    ///     vs. other grids. If \c GridType = \c Dune::CpGrid, local grid refinements (LGR)
    ///     are considered; otherwise, LGR is not supported.
    /// 
    /// \tparam ElementType
    ///     The type of the element or index passed in. 
    ///
    /// \param elementOrIndex
    ///     An integral cell index or a grid entity (e.g., \c Dune::cpgrid::Entity<0>).
    ///
    /// \return
    ///     The integer index to be used for looking up the relevant field properties. For
    ///     non-CpGrid grids, this is the same index or one retrieved from \c elemMapper_.
    ///     For CpGrid grids, it is either the LGR-based index or the origin index,
    ///     depending on \c isFieldPropInLgr_ and whether the cell is refined.
    ///
    /// \note
    ///     If \c GridType is not CpGrid, the method will assert that \c maxLevel() == 0,
    ///     since local grid refinements are only supported for CpGrid.
    template<typename GridType, typename ElementType>
    auto getFieldPropIdx(const ElementType& elem) const;

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
    ///
    /// \tparam ElementOrIndex the type of the element or index passed in.
    /// \tparam FieldProperties the type of the field properties vector. Should be an \c std::vector like type.
    template<class ElementOrIndex, class FieldProperties>
    auto operator()(const ElementOrIndex& elemIdx, const FieldProperties& fieldProp) const;


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

    /// \brief Calls \c getFieldPropCartesianIdx<Grid, ElementType>(elemIdx)
    template<class ElementType>
    auto getFieldPropCartesianIdx(const ElementType& elemIdx) const;

    /// \brief Return the cartesian index used to retrieve field properties, depending on
    ///        whether the grid is a CpGrid or another type of grid (with no LGR).
    ///
    /// This method determines the correct cartesian index under two main conditions:
    ///
    /// - **Non-CpGrid**: If \c GridType is not \c Dune::CpGrid, the cartesian index is
    ///   obtained via \c cartMapper_->cartesianIndex(...). The function asserts that
    ///   \c maxLevel() == 0, i.e., no local grid refinements (LGR). 
    ///
    /// - **CpGrid**: For CpGrid, if \c isFieldPropInLgr_ is true and the element resides
    ///   on a refined level (> 0), then its cartesian index is taken from the LGR-level
    ///   entity. Otherwise, the index is derived from the origin (level 0) entity through
    ///   \c cartMapper_ after mapping the underlying entity index.
    ///
    /// \tparam GridType
    ///     Auxiliary type used to specialize the method for CpGrid
    ///     vs. other grids. If \c GridType = \c Dune::CpGrid, local grid refinements (LGR)
    ///     are considered; otherwise, LGR is not supported.
    ///
    /// \tparam ElementType the type of the element or index passed in.
    ///
    ///
    /// \param elementOrIndex
    ///     An integral cell index or a grid entity (\c EntityType) used to compute the
    ///     cartesian index. For CpGrid, if an integral index is passed, the method creates
    ///     a \c Dune::cpgrid::Entity on the fly. Otherwise, it uses the provided entity.
    ///
    /// \return
    ///     The integer cartesian index to be used for looking up the relevant field
    ///     properties. For non-CpGrid, it is simply the mapped index via
    ///     \c cartMapper_->cartesianIndex(...). For CpGrid, it can be the LGR-based
    ///     cartesian index if \c isFieldPropInLgr_ and the element is refined, or the
    ///     origin-based cartesian index otherwise.
    ///
    /// \note
    ///     For non-CpGrid use-cases, the method asserts that \c maxLevel() == 0, since
    ///     local grid refinements are only supported for CpGrid.
    template<class GridType, class ElementType>
    auto getFieldPropCartesianIdx(const ElementType& elemIdx) const;

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
template<class ElementOrIndex, class FieldProperties>
auto Opm::LookUpData<Grid,GridView>::operator()(const ElementOrIndex& elemIdx,
                                                const FieldProperties& fieldProp) const
{
    const auto& fieldPropIdx = this->getFieldPropIdx(elemIdx);
    assert(0 <= fieldPropIdx && static_cast<int>(fieldProp.size()) > fieldPropIdx);
    return fieldProp[fieldPropIdx];
}

template<typename Grid, typename GridView>
std::vector<double> Opm::LookUpData<Grid,GridView>::assignFieldPropsDoubleOnLeaf(const FieldPropsManager& fieldPropsManager,
                                                                                 const std::string& propString) const
{
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
        const auto& fieldPropIdx = this->getFieldPropIdx(elemIdx); // gets parentIdx (or (lgr)levelIdx) for CpGrid with LGRs
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
template<typename ElementType>
auto Opm::LookUpData<Grid, GridView>::getFieldPropIdx(const ElementType& elementOrIndex) const {
    return this->template getFieldPropIdx<Grid, ElementType>(elementOrIndex);
}

template<typename Grid, typename GridView>
template<typename GridType, typename IndexType>
auto Opm::LookUpData<Grid,GridView>::getFieldPropIdx(const IndexType& elementOrIndex) const
{
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
template<typename IndexType, typename FieldPropType>
auto Opm::LookUpCartesianData<Grid,GridView>::operator()(const IndexType& elementOrIndex,
                                                         const FieldPropType& fieldProp) const
{
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
template<class ElementType>
auto Opm::LookUpCartesianData<Grid, GridView>::getFieldPropCartesianIdx(const ElementType& elemIdx) const {
    return this->getFieldPropCartesianIdx<Grid, ElementType>(elemIdx);
}

template<typename Grid, typename GridView>
template<typename GridType, typename ElementType>
auto Opm::LookUpCartesianData<Grid,GridView>::getFieldPropCartesianIdx(const ElementType& elementOrIndex) const
{
    constexpr static bool isIntegral = std::is_integral_v<ElementType>;
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
