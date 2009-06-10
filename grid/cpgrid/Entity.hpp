//===========================================================================
//
// File: Entity.hpp
//
// Created: Fri May 29 20:26:48 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
Copyright 2009 SINTEF ICT, Applied Mathematics.
Copyright 2009 Statoil ASA.

This file is part of The Open Reservoir Simulator Project (OpenRS).

OpenRS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OpenRS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPENRS_ENTITY_HEADER
#define OPENRS_ENTITY_HEADER

#include <boost/static_assert.hpp>
#include <dune/common/geometrytype.hh>
#include <dune/grid/common/gridenums.hh>
#include "EntityRep.hpp"

namespace Dune
{
    namespace cpgrid
    {




	template <int codim, class GridType>
	class Entity : public EntityRep<codim>
	{
	public:

	    enum { codimension = codim };
	    enum { dimension = 3 };
	    enum { mydimension = dimension - codimension };
	    enum { dimensionworld = 3 };

	    template <int cd>
	    struct Codim
	    {
		typedef typename GridType::template Codim<cd>::EntityPointer EntityPointer;
	    };

	    typedef typename GridType::template Codim<codim>::Geometry Geometry;

	    typedef double ctype;

	    /// Constructor taking a grid and an integer entity representation.
	    /// This constructor should probably be removed, since it exposes
	    /// details of the implementation of \see EntityRep, see comment in
	    /// EntityRep<>::EntityRep(int).
	    Entity(const GridType& grid, int entityrep)
		: EntityRep<codim>(entityrep), pgrid_(&grid)
	    {
	    }

	    /// Constructor taking a grid and an entity representation.
	    Entity(const GridType& grid, EntityRep<codim> entityrep)
		: EntityRep<codim>(entityrep), pgrid_(&grid)
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

	    /// Returns the geometry of the entity (does not depend on its orientation).
	    const Geometry& geometry() const
	    {
		return (*pgrid_).template geomVector<codim>()[*this];
	    }

	    /// We do not support refinement, so level() is always 0.
	    int level() const
	    {
		return 0;
	    }

	    /// The entity is always on the leaf grid, since we have no refinement.
	    bool isLeaf() const
	    {
		return true;
	    }

	    /// For now, the grid is serial and the only partitionType() is InteriorEntity.
	    PartitionType partitionType() const
	    {
		return InteriorEntity;
	    }

	    /// Using a singular as GeometryType for our entities.
	    GeometryType type() const
	    {
		GeometryType t;
		t.makeSingular(3 - codim);
		return t;
	    }

	    /// The count of subentities of codimension cc
	    template <int cc>
	    int count() const
	    {
		BOOST_STATIC_ASSERT(codim == 0);
		if (cc == 1) {
		    return pgrid_->cell_to_face_[*this].size();
		} else {
		    return 0;
		}
	    }

	    /// Obtain subentity.
	    template <int cd>
	    typename Codim<cd>::EntityPointer subEntity(int i) const
	    {
	    }

	    /// Start iterator for the cell-cell intersections of this entity.
	    typename GridType::Traits::LevelIntersectionIterator ilevelbegin() const
	    {
		BOOST_STATIC_ASSERT(codim == 0);
		return typename GridType::Traits::LevelIntersectionIterator(*pgrid_, *this, false);
	    }

	    /// End iterator for the cell-cell intersections of this entity.
	    typename GridType::Traits::LevelIntersectionIterator ilevelend() const
	    {
		BOOST_STATIC_ASSERT(codim == 0);
		return typename GridType::Traits::LevelIntersectionIterator(*pgrid_, *this, true);
	    }

	    /// Start iterator for the cell-cell intersections of this entity.
	    typename GridType::Traits::LeafIntersectionIterator ileafbegin() const
	    {
		BOOST_STATIC_ASSERT(codim == 0);
		return typename GridType::Traits::LeafIntersectionIterator(*pgrid_, *this, false);
	    }

	    /// End iterator for the cell-cell intersections of this entity.
	    typename GridType::Traits::LeafIntersectionIterator ileafend() const
	    {
		BOOST_STATIC_ASSERT(codim == 0);
		return typename GridType::Traits::LeafIntersectionIterator(*pgrid_, *this, true);
	    }

	    /// Dummy first son iterator.
	    typename GridType::Traits::HierarchicIterator hbegin(int) const
	    {
		return typename GridType::Traits::HierarchicIterator();
	    }

	    /// Dummy beyond last son iterator.
	    typename GridType::Traits::HierarchicIterator hend(int) const
	    {
		return typename GridType::Traits::HierarchicIterator();
	    }


	protected:
	    const GridType* pgrid_;
	};




	/// \brief Class representing a pointer to an entity.
	/// Implementation note:
	/// Since our entities are quite lightweight, we have chosen
	/// to implement EntityPointer by inheritance from
	/// Entity. Thus all dereferencing operators return the object
	/// itself as an Entity.
	template <int codim, class GridType>
	class EntityPointer : public Entity<codim, GridType>
	{
	public:
	    /// Constructor taking a grid and entity representation.
	    EntityPointer(const GridType& grid, int entityrep)
		: Entity<codim, GridType>(grid, entityrep)
	    {
	    }

	    /// Construction from entity.
	    explicit EntityPointer(const Entity<codim, GridType>& e)
		: Entity<codim, GridType>(e)
	    {
	    }

	    /// Member by pointer operator.
	    Entity<codim, GridType>* operator->()
	    {
		return this;

	    /// Const member by pointer operator.
	    }
	    const Entity<codim, GridType>* operator->() const
	    {
		return this;
	    }

	    /// Dereferencing operator.
	    Entity<codim, GridType>& operator*()
	    {
		return *this;
	    }

	    /// Const dereferencing operator.
	    const Entity<codim, GridType>& operator*() const
	    {
		return *this;
	    }

	    /// Minimizes memory usage.
	    /// Nothing to do in our case.
	    void compactify()
	    {
	    }
	};




    } // namespace cpgrid
} // namespace Dune


#endif // OPENRS_ENTITY_HEADER
