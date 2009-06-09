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
	    Entity(const GridType& grid, int entityrep)
		: EntityRep<codim>(entityrep), grid_(grid)
	    {
	    }
	    Entity(const GridType& grid, EntityRep<codim> entityrep)
		: EntityRep<codim>(entityrep), grid_(grid)
	    {
	    }

	    bool operator!=(const Entity& other) const
	    {
		return EntityRep<codim>::operator!=(other)  ||  &grid_ != &other.grid_;
	    }

	    typedef typename GridType::template Codim<codim>::Geometry Geometry;
	    const Geometry& geometry() const
	    {
		return grid_.template geomVector<codim>()[*this];
	    }

	    int level() const
	    {
		return 0;
	    }

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

	    typename GridType::Traits::LeafIntersectionIterator ileafbegin() const
	    {
		return typename GridType::Traits::LeafIntersectionIterator(grid_, *this, false);
	    }

	    typename GridType::Traits::LeafIntersectionIterator ileafend() const
	    {
		return typename GridType::Traits::LeafIntersectionIterator(grid_, *this, true);
	    }

	protected:
	    const GridType& grid_;
	};





	template <int codim, class GridType>
	class EntityPointer : public Entity<codim, GridType>
	{
	public:
	    EntityPointer(const GridType& grid, int index)
		: Entity<codim, GridType>(grid, index)
	    {
	    }

	    Entity<codim, GridType>* operator->()
	    {
		return this;
	    }
	    const Entity<codim, GridType>* operator->() const
	    {
		return this;
	    }

	    Entity<codim, GridType>& operator*()
	    {
		return *this;
	    }
	    const Entity<codim, GridType>& operator*() const
	    {
		return *this;
	    }
	};




    } // namespace cpgrid
} // namespace Dune


#endif // OPENRS_ENTITY_HEADER
