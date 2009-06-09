//===========================================================================
//
// File: EntityRep.hpp
//
// Created: Tue Jun  9 11:11:24 2009
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

#ifndef OPENRS_ENTITYREP_HEADER
#define OPENRS_ENTITYREP_HEADER



#include <boost/static_assert.hpp>
#include "../common/SparseTable.hpp"
#include <map>
#include <boost/algorithm/minmax_element.hpp>


namespace Dune
{
    namespace cpgrid
    {


	/// \brief Represents an entity of a given codim, with positive or negative orientation.
	/// This class is not a part of the Dune interface, but of our implementation.
	/// Since this class has a few friends, and for aid in debugging, we document its
	/// interior representation here:
	/// The interior representation consists of an integer entityrep_
	/// which, if positive or zero, indicates the index of the entity.
	/// In that case, the entity's orientation is positive.
	/// If entityrep_ is negative, the orientation is negative, and the index
	/// is given by ~entityrep_ (we cannot use -entityrep_, since 0 is a valid index).
	/// We may consider changing this representation to using something like a
	/// std::pair<int, bool> instead.
	template <int codim>
	class EntityRep
	{
	public:
	    /// Constructor taking an integer representation directly.
	    /// This is one of the few places where the private representation is exposed,
	    /// the others being in the classes that are friends of this one. These places
	    /// need to be modified if we change the representation.
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

	    /// Equality operator.
	    bool operator==(const EntityRep& other) const
	    {
		return entityrep_ == other.entityrep_;
	    }

	    /// Inequality operator.
	    bool operator!=(const EntityRep& other) const
	    {
		return !operator==(other);
	    }

	protected:
	    // Interior representation is documented in class main comment.
	    int entityrep_;
	    // These friends need not be forward declared since they are templates.
	    template <typename T, int cd>
	    friend class EntityVariable;
	    template <typename T, int cd>
	    friend class SignedEntityVariable;
	    template <int cf, int ct>
	    friend class OrientedEntityTable;
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
	    typedef EntityRep<codim_to> ToType;
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

	    void printRelationMatrix(std::ostream& os) const
	    {
		int columns = numberOfColumns();
		for (int i = 0; i < size(); ++i) {
		    FromType from_ent(i);
		    row_type r  = operator[](from_ent);
		    int cur_col = 0;
		    int next_ent = 0;
		    ToType to_ent = r[next_ent];
		    int next_print = to_ent.index();
		    while (cur_col < columns) {
			if (cur_col == next_print) {
			    if (to_ent.orientation()) {
				os << "  1";
			    } else {
				os << " -1";
			    }
			    ++next_ent;
			    if (next_ent >= r.size()) {
				next_print = columns;
			    } else {
				to_ent = r[next_ent];
				next_print = to_ent.index();
			    }
			} else {
			    os << "  0";
			}
			++cur_col;
		    }
		    os << '\n';
		}
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

	private:
	    int numberOfColumns() const
	    {
		int maxind = 0;
		for (int i = 0; i < size(); ++i) {
		    FromType from_ent(i);
		    row_type r  = operator[](from_ent);
		    for (int j = 0; j < r.size(); ++j) {
			maxind = std::max(maxind, r[j].index());
		    }
		}
		return maxind + 1;
	    }
	};


    } // namespace cpgrid
} // namespace Dune




#endif // OPENRS_ENTITYREP_HEADER
