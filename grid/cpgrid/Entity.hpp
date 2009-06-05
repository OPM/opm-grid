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
#include "../common/SparseTable.hpp"
#include <map>

namespace Dune
{
    namespace cpgrid
    {


	// Forward declarations.
	template <typename T, int codim> class EntityVariable;
	template <typename T, int codim> class SignedEntityVariable;


	template <int codim>
	class EntityRep
	{
	public:
	    explicit EntityRep(int erep)
		: entityrep_(erep)
	    {
	    }
	    /// The (positive) index of an entity. Not a Dune interface method.
	    int index() const
	    {
		return entityrep_ < 0 ? ~entityrep_ : entityrep_;
	    }

	    /// Returns true if the entity has positive orientation. Not a Dune interface method.
	    bool orientation() const
	    {
		return entityrep_ >= 0;
	    }

	    /// Ordering relation used for maps etc. Sorting on index and then orientation,
	    /// with positive orientations first.
	    bool operator<(const EntityRep& other) const
	    {
		int i1 = index();
		int i2 = other.index();
		if (i1 < i2) return true;
		if (orientation() && !other.orientation()) return true;
		return false;
	    }

	    /// Equality.
	    bool operator==(const EntityRep& other) const
	    {
		return entityrep_ == other.entityrep_;
	    }

	    /// Inequality.
	    bool operator!=(const EntityRep& other) const
	    {
		return entityrep_ != other.entityrep_;
	    }

	protected:
	    int entityrep_;

	    template <typename T, int cd>
	    friend class EntityVariable;
	    template <typename T, int cd>
	    friend class SignedEntityVariable;
	    template <int cf, int ct>
	    friend class OrientedEntityTable;
	};


	template <int codim, class GridType>
	class Entity : public EntityRep<codim>
	{
	public:
	    Entity(const GridType& grid, int entityrep)
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




	template <typename T>
	class EntityVariableBase : private std::vector<T>
	{
	public:
	    typedef std::vector<T> V;
	    using V::empty;
	    using V::size;
	    using V::assign;
	protected:
	    const T& get(int i) const
	    {
		return V::operator[](i);
	    }
	};

	template <typename T, int codim>
	class EntityVariable : public EntityVariableBase<T>
	{
	public:
	    const T& operator[](const EntityRep<codim>& e) const
	    {
		return get(e.index());
	    }
	};

	template <typename T, int codim>
	class SignedEntityVariable : public EntityVariableBase<T>
	{
	public:
	    const T operator[](const EntityRep<codim>& e) const
	    {
		return e.entityrep_ < 0 ?
		    -get(~e.entityrep_)
		    : get(e.entityrep_);
	    }
	};

	template <int codim_to>
	class OrientedEntityRange : private SparseTable<int>::row_type
	{
	public:
	    typedef SparseTable<int>::row_type R;
	    typedef EntityRep<codim_to> ToType;
	    OrientedEntityRange(const R& r, bool orientation)
		: R(r), orientation_(orientation)
	    {
	    }
	    using R::size;
	    using R::empty;
	    ToType operator[](int subindex) const
	    {
		int erep = R::operator[](subindex);
		return ToType(orientation_ ? erep : ~erep);
	    }
	private:
	    bool orientation_;
	};


	template <int codim_from, int codim_to>
	class OrientedEntityTable : private SparseTable<int>
	{
	public:
	    typedef EntityRep<codim_from> FromType;
	    typedef OrientedEntityRange<codim_to> row_type;

	    OrientedEntityTable()
	    {
	    }

	    template <typename DataIter, typename IntegerIter>
	    OrientedEntityTable(DataIter data_beg, DataIter data_end,
				IntegerIter rowsize_beg, IntegerIter rowsize_end)
		: SparseTable<int>(data_beg, data_end, rowsize_beg, rowsize_end)
	    {
	    }

	    using SparseTable<int>::empty;
	    using SparseTable<int>::size;
	    row_type operator[](const FromType& e) const
	    {
		return row_type(SparseTable<int>::operator[](e.index()), e.orientation());
	    }

	    bool operator==(const OrientedEntityTable& other) const
	    {
		return SparseTable<int>::operator==(other);
	    }

	    void printRelationMatrix()
	    {
	    }

	    void makeInverseRelation(OrientedEntityTable<codim_to, codim_from>& inv) const
	    {
		typedef std::multimap<int, int> RelationMap;
		RelationMap rm;
		for (int i = 0; i < size(); ++i) {
		    EntityRep<codim_from> from_ent(i);
		    row_type r = operator[](from_ent);
		    for (int j = 0; j < r.size(); ++j) {
			EntityRep<codim_to> to_ent = r[j];
			// Making sure all the keys we insert are positive.
			int from = to_ent.orientation() ? i : ~i;
			rm.insert(std::make_pair(to_ent.index(), from));
		    }
		}
		ASSERT(int(rm.size()) == SparseTable<int>::dataSize());
		std::vector<int> new_data;
		new_data.reserve(rm.size());
		int last = (--rm.end())->first; // The last key.
		std::vector<int> new_sizes(last + 1, 0);
		for (typename RelationMap::iterator it = rm.begin(); it != rm.end(); ++it) {
		    new_data.push_back(it->second);
		    ++new_sizes[it->first];
		}
		ASSERT(new_data.size() == rm.size());
		inv = OrientedEntityTable<codim_to, codim_from>(new_data.begin(),
								new_data.end(),
								new_sizes.begin(),
								new_sizes.end());
	    }

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
