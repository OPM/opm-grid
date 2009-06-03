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
#include "../common/SparseTable.hpp"

namespace Dune
{
    namespace cpgrid
    {


	// Forward declarations.
	template <typename T, class EntityType> class EntityVariable;
	template <typename T, class EntityType> class SignedEntityVariable;



	template <int codim, class GridType>
	class Entity
	{
	public:
	    Entity(const GridType& grid, int entityrep)
		: grid_(grid), entityrep_(entityrep)
	    {
	    }

	    bool operator!=(const Entity& other) const
	    {
		return entityrep_ != other.entityrep_  ||  &grid_ != &other.grid_;
	    }

	    typedef typename GridType::template Codim<codim>::Geometry Geometry;
	    const Geometry& geometry() const
	    {
		return grid_.template geomVector<codim>()[*this];
	    }

	    /// Using a hexahedron as GeometryType for our entities.
	    GeometryType type() const
	    {
		return GeometryType(3);
	    }

	    /// The count of subentities of codimension cc
	    template <int cc>
	    int count() const
	    {
		BOOST_STATIC_ASSERT(codim == 0);
		BOOST_STATIC_ASSERT(cc == 1);
		return grid_.cell_to_face_[*this].size();
	    }

	    /// The (positive) index of an entity.
	    int index() const
	    {
		return entityrep_ < 0 ? ~entityrep_ : entityrep_;
	    }

	protected:
	    /// Returns true if the entity has positive orientation.
	    bool pos_orient() const
	    {
		return entityrep_ >= 0;
	    }

	    const GridType& grid_;
	    int entityrep_;

	    template <typename T, class ET>
	    friend class EntityVariable;
	    template <typename T, class ET>
	    friend class SignedEntityVariable;
	};




	template <typename T>
	class EntityVariableBase : private std::vector<T>
	{
	public:
	    typedef std::vector<T> V;
	    using V::size;
	    using V::assign;
	protected:
	    const T& get(int i) const
	    {
		return V::operator[](i);
	    }
	};

	template <typename T, class EntityType>
	class EntityVariable : public EntityVariableBase<T>
	{
	public:
	    const T& operator[](const EntityType& e) const
	    {
		return get(e.index());
	    }
	};

	template <typename T, class EntityType>
	class SignedEntityVariable : public EntityVariableBase<T>
	{
	public:
	    const T operator[](const EntityType& e) const
	    {
		return e.entityrep_ < 0 ?
		    -get(~e.entityrep_)
		    : get(e.entityrep_);
	    }
	};

	template <class GridType, int codim_to>
	class OrientedEntityRange : private SparseTable<int>::row_type
	{
	public:
	    typedef SparseTable<int>::row_type R;
	    typedef typename GridType::Traits::template Codim<codim_to>::Entity EntityToType;
	    OrientedEntityRange(const GridType& grid, bool pos_orient, const R& r)
		: R(r), grid_(grid), pos_orient_(pos_orient)
	    {
	    }
	    using R::size;
	    using R::empty;
	    EntityToType operator[](int subindex) const
	    {
		int erep = R::operator[](subindex);
		return EntityToType(grid_, pos_orient_ ? erep : ~erep);
	    }
	private:
	    const GridType& grid_;
	    bool pos_orient_;
	};

	template <class GridType, int codim_from, int codim_to>
	class OrientedEntityTable : private SparseTable<int>
	{
	public:
	    typedef typename GridType::Traits::template Codim<codim_from>::Entity EntityFromType;
	    typedef typename GridType::Traits::template Codim<codim_to>::Entity EntityToType;
	    using SparseTable<int>::empty;
	    using SparseTable<int>::size;
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
