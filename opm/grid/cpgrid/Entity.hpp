//===========================================================================
//
// File: Entity.hpp
//
// Created: Fri May 29 20:26:48 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Brd Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010, 2022 Equinor ASA.

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

#ifndef OPM_ENTITY_HEADER
#define OPM_ENTITY_HEADER

#include <dune/common/version.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/gridenums.hh>

#include "PartitionTypeIndicator.hpp"


namespace Dune
{
    namespace cpgrid
    {

       template<int,int> class Geometry;
       template<int,PartitionIteratorType> class Iterator;
       class IntersectionIterator;
       class HierarchicIterator;
       class CpGridData;
       class LevelGlobalIdSet;


        /// @brief
        /// @todo Doc me!
        /// @tparam
        template <int codim>
        class Entity : public EntityRep<codim>
        {
            friend class LevelGlobalIdSet;
            friend class GlobalIdSet;
            friend class HierarchicIterator;
            friend class CpGridData;

        public:
        /// @brief
        /// @todo Doc me!
            enum { codimension = codim };
            enum { dimension = 3 };
            enum { mydimension = dimension - codimension };
            enum { dimensionworld = 3 };

            // the official DUNE names
            typedef Entity    EntitySeed;

            /// @brief
            /// @todo Doc me!
            /// @tparam
            template <int cd>
            struct Codim
            {
                typedef cpgrid::Entity<cd> Entity;
            };


            typedef cpgrid::Geometry<3-codim,3> Geometry;
            typedef Geometry LocalGeometry;

            typedef cpgrid::IntersectionIterator LeafIntersectionIterator;
            typedef cpgrid::IntersectionIterator LevelIntersectionIterator;
            typedef cpgrid::HierarchicIterator HierarchicIterator;

            typedef double ctype;

            /// Constructor taking a grid and an integer entity representation.
            /// This constructor should probably be removed, since it exposes
            /// details of the implementation of \see EntityRep, see comment in
            /// EntityRep<>::EntityRep(int).
//             Entity(const CpGridData& grid, int entityrep)
//              : EntityRep<codim>(entityrep), pgrid_(&grid)
//          {
//          }

            /// Constructor creating empty entity
            Entity()
                : EntityRep<codim>(), pgrid_( 0 )
            {
            }

            /// Constructor taking a grid and an entity representation.
            Entity(const CpGridData& grid, EntityRep<codim> entityrep)
                : EntityRep<codim>(entityrep), pgrid_(&grid)
            {
            }

            /// Constructor taking a grid, entity index, and orientation.
            Entity(const CpGridData& grid, int index_arg, bool orientation_arg)
                : EntityRep<codim>(index_arg, orientation_arg), pgrid_(&grid)
            {
            }

            /// Equality.
            bool operator==(const Entity& other) const
            {
                return EntityRep<codim>::operator==(other)  &&  pgrid_ == other.pgrid_;
            }

            /// Inequality.
            bool operator!=(const Entity& other) const
            {
                return !operator==(other);
            }

            /// @brief Return an entity seed (light-weight entity).
            ///        EntitySeed objects are used to obtain an Entity back when combined with the corresponding grid.
            ///        For CpGrid, EntitySeed and EntityPtr are the same class.
            EntitySeed seed() const
            {
                return EntitySeed( impl() );
            }

            /// @brief Return the geometry of the entity (does not depend on its orientation).
            const Geometry& geometry() const;

            /// @brief Return the level of the entity in the grid hierarchy. Level = 0 represents the coarsest grid.
            int level() const;

            /// @brief Check if the entity is in the leafview.
            ///
            ///        @TODO: Modify the definition to cover serial and parallel cases.
            ///        Serial: an element is a leaf <-> hbegin and hend return the same iterator
            ///        Parallel: true <-> the element is a leaf entity of the global refinement hierarchy.
            bool isLeaf() const;
            
            /// Refinement is not defined for CpGrid.
            bool isRegular() const
            {
                return true;
            }

            /// @brief For now, the grid is serial and the only partitionType() is InteriorEntity.
            ///        Only needed when distributed_data_ is not empty.
            PartitionType partitionType() const;
            
            /// @brief Return marker object (GeometryType object) representing the reference element of the entity.
            ///        Currently, cube type for all entities (cells and vertices).
            GeometryType type() const
            {
                return Dune::GeometryTypes::cube(3 - codim);
            }
            
            /// @brief Return the number of all subentities of the entity of a given codimension cc.
            unsigned int subEntities ( const unsigned int cc ) const;
            
            /// @brief Obtain subentity.
            ///        Example: If cc = 3 and i = 5, it returns the 5th corner/vertex of the entity.
            template <int cc>
            typename Codim<cc>::Entity subEntity(int i) const;

            /// Start level-iterator for the cell-cell intersections of this entity.
            inline LevelIntersectionIterator ilevelbegin() const;

            /// End level-iterator for the cell-cell intersections of this entity.
            inline LevelIntersectionIterator ilevelend() const;

            /// Start leaf-iterator for the cell-cell intersections of this entity.
            inline LeafIntersectionIterator ileafbegin() const;

            /// End leaf-iterator for the cell-cell intersections of this entity.
            inline LeafIntersectionIterator ileafend() const;

            
            /// @brief Iterator begin over the children. [If requested, also over descendants more than one generation away.]
            HierarchicIterator hbegin(int) const;
            
            /// @brief Iterator end over the children/beyond last child iterator.
            HierarchicIterator hend(int) const; 

            /// \brief Returns true, if the entity has been created during the last call to adapt(). Dummy.
            bool isNew() const
            {
                return false;
            }

            /// \brief Returns true, if entity might disappear during the next call to adapt(). Dummy.
            bool mightVanish() const
            {
                return false;
            }

            /// @brief ONLY FOR CELLS (Entity<0>)
            ///        Check if the entity comes from an LGR, i.e., it has been created via refinement from coarser level.
            ///
            ///        @TODO: When distributed_data_ is not empty, check whether the father element exists on the
            ///        local process, which can be used to test whether it is safe to call father.
            bool hasFather() const;

            /// @brief  ONLY FOR CELLS (Entity<0>). Get the father Entity, in case entity.hasFather() is true.
            ///
            /// @return father-entity
            Entity<0> father() const; 

            /// @brief Return LocalGeometry representing the embedding of the entity inti its father (when hasFather() is true).
            ///        Map from the entity's reference element into the reference element of its father.
            ///        Currently, LGR is built via refinement of a block-shaped patch from the coarse grid. So the LocalGeometry
            ///        of an entity coming from the LGR is one of the refined cells of the unit cube, with suitable amount of cells
            ///        in each direction.
            Dune::cpgrid::Geometry<3,3> geometryInFather() const;
            
            /// Returns true if any of my intersections are on the boundary.
            /// Implementation note:
            /// This is a slow, computed, function. Could be speeded
            /// up by putting boundary info in the CpGrid class.
            bool hasBoundaryIntersections() const;

            // Mimic Dune entity wrapper
            /// @brief Access the actual implementation class behind Entity interface class. 
            const Entity& impl() const
            {
                return *this;
            }

            Entity& impl()
            {
                return *this;
            }

            /// isValid method for EntitySeed
            /// \return return true if seed is pointing to a valid entity
            bool isValid () const;

        protected:
            const CpGridData* pgrid_;
        };

    } // namespace cpgrid
} // namespace Dune

// now we include the Iterators.hh We need to do this here because for hbegin/hend the compiler
// needs to know the size of hierarchicIterator
#include "Iterators.hpp"
#include "Intersection.hpp"
namespace Dune
{
  namespace cpgrid
  {
  template<int codim>
  typename Entity<codim>::LevelIntersectionIterator Entity<codim>::ilevelbegin() const
  {
                static_assert(codim == 0, "");
                return LevelIntersectionIterator(*pgrid_, *this, false);
            }

  template<int codim>
  typename Entity<codim>::LevelIntersectionIterator Entity<codim>::ilevelend() const
  {
      static_assert(codim == 0, "");
      return LevelIntersectionIterator(*pgrid_, *this, true);
  }

  template<int codim>
  typename Entity<codim>::LeafIntersectionIterator Entity<codim>::ileafbegin() const
  {
      static_assert(codim == 0, "");
      return LeafIntersectionIterator(*pgrid_, *this, false);
  }

  template<int codim>
  typename Entity<codim>::LeafIntersectionIterator Entity<codim>::ileafend() const
  {
      static_assert(codim == 0, "");
      return LeafIntersectionIterator(*pgrid_, *this, true);
  }


  template<int codim>
  HierarchicIterator Entity<codim>::hbegin(int maxLevel) const 
  {
      // Creates iterator with first child as target if there is one. Otherwise empty stack and target.
      return HierarchicIterator(*this, maxLevel);
  }

  /// Dummy beyond last child iterator.
  template<int codim>
  HierarchicIterator Entity<codim>::hend(int maxLevel) const
  {
      // Creates iterator with empty stack and target.
      return HierarchicIterator(maxLevel);
  }

    template <int codim>
    PartitionType Entity<codim>::partitionType() const
    {
        return pgrid_->partition_type_indicator_->getPartitionType(*this);
    }
    } // namespace cpgrid
} // namespace Dune


#include <opm/grid/cpgrid/CpGridData.hpp>

namespace Dune {
namespace cpgrid {

template<int codim>
unsigned int Entity<codim>::subEntities ( const unsigned int cc ) const
{
    if (cc == codim) {
            return 1;
    }
    else if ( codim == 0 ){ // Cell/element/Entity<0>
        if ( cc == 3 ) { // Get number of corners of the element.
            return 8;
        }
    }
    return 0;
}

template <int codim>
const typename Entity<codim>::Geometry& Entity<codim>::geometry() const 
{
    return pgrid_->geomVector<codim>()[*this];
}

template <int codim>
template <int cc>
typename Entity<codim>::template Codim<cc>::Entity Entity<codim>::subEntity(int i) const
{
    static_assert(codim == 0, "");
    if (cc == 0) { // Cell/element/Entity<0>
        assert(i == 0);
        typename Codim<cc>::Entity se(*pgrid_, this->index(), this->orientation());
        return se;
    } else if (cc == 3) { // Corner/Entity<3>
        assert(i >= 0 && i < 8);
        int corner_index = pgrid_->cell_to_point_[this->index()][i];
        typename Codim<cc>::Entity se(*pgrid_, corner_index, true);
        return se;
    }
    else {
        OPM_THROW(std::runtime_error,
                  "No subentity exists of codimension " + std::to_string(cc));
    }
}

template <int codim>
bool Entity<codim>::hasBoundaryIntersections() const
{
    // Copied implementation from EntityDefaultImplementation,
    // except for not checking LevelIntersectionIterators.
    typedef LeafIntersectionIterator Iter;
    Iter end = ileafend();
    for (Iter it = ileafbegin(); it != end; ++it) {
        if (it->boundary()) return true;
    }
    return false;
}

template <int codim>
bool Entity<codim>::isValid() const
{
    return pgrid_ ? EntityRep<codim>::index() < pgrid_->size(codim) : false;
}


// level() It simply returns the level of the entity in the grid hierarchy.
template <int codim>
int Entity<codim>::level() const
{
    // if distributed_data_ is not empty, level_data_ptr_ has size 1.
    if ((*(pgrid_ -> level_data_ptr_)).size() == 1){
        return 0; // when there is no refinenment, level_ is not automatically instantiated 
    }
    if (pgrid_ == (*(pgrid_->level_data_ptr_)).back().get()) { // entity on the leafview -> get the level where it was born:
        return pgrid_ -> leaf_to_level_cells_[this-> index()][0]; // leaf_to_level_cells_ leafIdx -> {level/LGR, cell index in LGR}
    }
    else {
        return pgrid_-> level_;
    }
}


// isLeaf() 
// - if distributed_data_ is empty: an element is a leaf <-> hbegin and hend return the same iterator. Then,
//    *cells from level 0 (coarse grid) that are not parents, are Leaf.
//    *cells from any LGR are Leaf, since they do not children (nested refinement not supported). 
// - if distrubuted_data_ is NOT empty: there may be children on a different process. Therefore,
// isLeaf() returns true <-> the element is a leaf entity of the global refinement hierarchy. Equivalently,
// it can be checked whether parent_to_children_cells_ is empty.
template<int codim>
bool Entity<codim>::isLeaf() const
{
    if ((pgrid_ -> parent_to_children_cells_).empty()){ // LGR cells
        return true; 
    }
    else {
        return (std::get<0>((pgrid_ -> parent_to_children_cells_)[this-> index()]) == -1);  // Cells from GLOBAL, not involved in any LGR
    }
}

template<int codim>
bool Entity<codim>::hasFather() const
{
    if (pgrid_ -> child_to_parent_cells_.empty()){
        return false;
    }
    else{
        const auto& [level, parent] = pgrid_-> child_to_parent_cells_[this->index()]; // {-1,-1};
        return  (parent != -1);
    }
}

template<int codim>
Entity<0> Entity<codim>::father() const
{
    if (!(this->hasFather())){
        OPM_THROW(std::logic_error, "Entity has no father.");
    }
    const int& coarse_level = pgrid_ -> child_to_parent_cells_[this->index()][0];
    const int& parent_index = pgrid_ -> child_to_parent_cells_[this->index()][1];
    const auto& coarse_grid = (*(pgrid_ -> level_data_ptr_))[coarse_level].get(); 
    return Entity<0>( *coarse_grid, parent_index, true); 
}

template<int codim>
Dune::cpgrid::Geometry<3,3> Dune::cpgrid::Entity<codim>::geometryInFather() const
{
    if (!(this->hasFather())){
        OPM_THROW(std::logic_error, "Entity has no father.");
    }
    else{
#if 0
        //
        DefaultGeometryPolicy local_geometry;
        std::array<int,8> localEntity_to_point;
        std::array<int,8> allcorners_localEntity;
        EntityVariableBase<cpgrid::Geometry<0,3>>& local_corners = local_geometry.geomVector<codim>();
        local_corners.resize(8);
        // Get IJK index of the entity.
        std::array<int,3> eIJK;
        pgrid_ -> getIJK(this->index(), eIJK);
        // Get dimension of the grid.
        const auto& grid_dim = pgrid_ -> logicalCartesianSize(); // {cells_per_dim[0]*patch_dim[0] ...[1], ...[2]}
        // Get the local coordinates of the entity (in the reference unit cube).
        const std::vector<FieldVector<double, 3>>& local_corners_temp = {
            // corner '0'
            { double(eIJK[0])/grid_dim[0], double(eIJK[1])/grid_dim[1], double(eIJK[2])/grid_dim[2] },
            // corner '1'
            { (double(eIJK[0])+1)/grid_dim[0], double(eIJK[1])/grid_dim[1], double(eIJK[2])/grid_dim[2] },
            // corner '2'
            { double(eIJK[0])/grid_dim[0], (double(eIJK[1])+1)/grid_dim[1], double(eIJK[2])/grid_dim[2] },
            // corner '3'
            { (double(eIJK[0])+1)/grid_dim[0], (double(eIJK[1])+1)/grid_dim[1], double(eIJK[2])/grid_dim[2] },
            // corner '4'
            { double(eIJK[0])/grid_dim[0], double(eIJK[1])/grid_dim[1], (double(eIJK[2])+1)/grid_dim[2] },
            // corner '5'
            { (double(eIJK[0])+1)/grid_dim[0], double(eIJK[1])/grid_dim[1], (double(eIJK[2])+1)/grid_dim[2] },
            // corner '6'
            { double(eIJK[0])/grid_dim[0], (double(eIJK[1])+1)/grid_dim[1], (double(eIJK[2])+1)/grid_dim[2] },
            // corner '7'
            { (double(eIJK[0])+1)/grid_dim[0], (double(eIJK[1])+1)/grid_dim[1], (double(eIJK[2])+1)/grid_dim[2] }};
        // Compute the center of the 'local-entity'.
        Dune::FieldVector<double, 3> local_center = {0., 0.,0.};
        for (int corn = 0; corn < 8; ++corn) {
            local_corners[corn] = local_corners_temp[corn];
            local_center += local_corners[corn].center()/8.;
        }
        // Compute the volume of the 'local-entity'.
        double local_volume = double(1)/(grid_dim[0]*grid_dim[1]*grid_dim[2]);
        // Indices of 'all the corners', in this case, 0-7 (required to construct a Geometry<3,3> object).
        allcorners_localEntity= {0,1,2,3,4,5,6,7};
        // Create a pointer to the first element of "cellfiedPatch_to_point" (required to construct a Geometry<3,3> object).
        const int* localEntity_indices_storage_ptr = &allcorners_localEntity[0];
        // Construct (and return) the Geometry<3,3> of the 'cellified patch'.
        return Dune::cpgrid::Geometry<3,3>(local_center, local_volume,
                                           local_geometry.geomVector<codim>(), localEntity_indices_storage_ptr);
#endif
        OPM_THROW(std::logic_error, "geometryInFather() not implemented");
        return {};
    }
}

} // namespace cpgrid
} // namespace Dune


#endif // OPM_ENTITY_HEADER
