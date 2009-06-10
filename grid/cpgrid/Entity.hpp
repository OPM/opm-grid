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
	    /// Constructor taking a grid and an integer entity representation.
	    /// This constructor should probably be removed, since it exposes
	    /// details of the implementation of \see EntityRep, see comment in
	    /// EntityRep<>::EntityRep(int).
	    Entity(const GridType& grid, int entityrep)
		: EntityRep<codim>(entityrep), grid_(grid)
	    {
	    }

	    /// Constructor taking a grid and an entity representation.
	    Entity(const GridType& grid, EntityRep<codim> entityrep)
		: EntityRep<codim>(entityrep), grid_(grid)
	    {
	    }

	    /// Equality.
	    bool operator==(const Entity& other) const
	    {
		return EntityRep<codim>::operator==(other)  &&  &grid_ == &other.grid_;
	    }

	    /// Inequality.
	    bool operator!=(const Entity& other) const
	    {
		return !operator==(other);
	    }

	    /// Returns the geometry of the entity (does not depend on its orientation).
	    typedef typename GridType::template Codim<codim>::Geometry Geometry;
	    const Geometry& geometry() const
	    {
		return grid_.template geomVector<codim>()[*this];
	    }

	    /// We do not support refinement, so level() is always 0.
	    int level() const
	    {
		return 0;
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
		    return grid_.cell_to_face_[*this].size();
		} else {
		    return 0;
		}
	    }

	    /// Start iterator for the cell-cell intersections of this entity.
	    typename GridType::Traits::LeafIntersectionIterator ileafbegin() const
	    {
		BOOST_STATIC_ASSERT(codim == 0);
		return typename GridType::Traits::LeafIntersectionIterator(grid_, *this, false);
	    }

	    /// End iterator for the cell-cell intersections of this entity.
	    typename GridType::Traits::LeafIntersectionIterator ileafend() const
	    {
		BOOST_STATIC_ASSERT(codim == 0);
		return typename GridType::Traits::LeafIntersectionIterator(grid_, *this, true);
	    }

	protected:
	    const GridType& grid_;
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
	};




    } // namespace cpgrid
} // namespace Dune


#endif // OPENRS_ENTITY_HEADER
