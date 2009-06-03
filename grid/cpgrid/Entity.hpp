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


namespace Dune
{
    namespace cpgrid
    {


	template <int cd, class GridType>
	class Entity
	{
	public:
	    Entity(const GridType& grid, int index)
		: grid_(grid), index_(index)
	    {
	    }

	    bool operator!=(const Entity& other) const
	    {
		return index_ != other.index_  ||  &grid_ != &other.grid_;
	    }

	    typedef typename GridType::template Codim<cd>::Geometry Geometry;
	    const Geometry& geometry() const
	    {
		return grid_.template geomVector<cd>()[index_];
	    }

	    /// Using a hexahedron as GeometryType for our entities.
	    GeometryType type() const
	    {
		return GeometryType(3);
	    }

	    /// The index of an entity. Not a Dune interface function.
	    int index() const
	    {
		return index_;
	    }

	protected:
	    const GridType& grid_;
	    int index_;

	};

	template <int cd, class GridType>
	class EntityPointer : public Entity<cd, GridType>
	{
	public:
	    EntityPointer(const GridType& grid, int index)
		: Entity<cd, GridType>(grid, index)
	    {
	    }

	    Entity<cd, GridType>* operator->()
	    {
		return this;
	    }
	    const Entity<cd, GridType>* operator->() const
	    {
		return this;
	    }

	    Entity<cd, GridType>& operator*()
	    {
		return *this;
	    }
	    const Entity<cd, GridType>& operator*() const
	    {
		return *this;
	    }
	};


    } // namespace cpgrid
} // namespace Dune


#endif // OPENRS_ENTITY_HEADER
