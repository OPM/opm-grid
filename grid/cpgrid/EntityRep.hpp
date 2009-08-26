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



//#include <boost/static_assert.hpp>
//#include <dune/common/SparseTable.hpp>
#include <dune/common/ErrorMacros.hpp>
#include <climits>
//#include <boost/algorithm/minmax_element.hpp>
#include <vector>

/// The namespace Dune is the main namespace for all Dune code.
namespace Dune
{
    namespace cpgrid
    {

	/// @brief Represents an entity of a given codim, with positive or negative orientation.
	///
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
	/// @tparam codim Codimension

	template <int codim>
	class EntityRep
	{
	public:
	    /// Default constructor.
	    EntityRep()
		: entityrep_(0)
	    {
	    }
	    /// @brief Constructor taking an integer representation directly.
	    ///
	    /// This is one of the few places where the private representation is exposed,
	    /// the others being in the classes that inherit this one. These places
	    /// need to be modified if we change the representation, then we should remove
	    /// this constructor.
	    /// @param erep Entity representation.
	    explicit EntityRep(int erep)
		: entityrep_(erep)
	    {
	    }
	    /// @brief Constructor taking an entity index and an orientation.
	    /// @param index Entity index
	    /// @param orientation True if the entity's orientations is positive.
	    EntityRep(int index, bool orientation)
		: entityrep_(orientation ? index : ~index)
	    {
		ASSERT(index >= 0);
	    }
	    /// @brief Set entity value.
	    /// @param index Entity index
	    /// @param orientation True if the entity's orientations is positive.
	    void setValue(int index, bool orientation)
	    {
		ASSERT(index >= 0);
		entityrep_ = orientation ? index : ~index;
	    }
	    /// @brief The (positive) index of an entity. Not a Dune interface method.
	    /// @return the (positive) index of an entity.
	    int index() const
	    {
		return entityrep_ < 0 ? ~entityrep_ : entityrep_;
	    }

	    /// @brief Returns true if the entity has positive orientation.
	    /// Not a Dune interface method.
	    ///
	    /// @return true if the entity has positive orientation.
	    bool orientation() const
	    {
		return entityrep_ >= 0;
	    }

	    /// @brief Returns an EntityRep with opposite orientation.
	    /// @return an EntityRep with opposite orientation.
	    EntityRep opposite() const
	    {
		return EntityRep(~entityrep_);
	    }

	    /// @brief Ordering relation used for maps etc.
	    ///
	    /// Sorting on index and then orientation, with positive orientations first.
	    /// @param other The other entity representation.
	    /// @return true if \b this element is less than the \b other.
	    bool operator<(const EntityRep& other) const
	    {
		int i1 = index();
		int i2 = other.index();
		if (i1 < i2) return true;
		if (orientation() && !other.orientation()) return true;
		return false;
	    }

	    /// @brief Equality operator.
	    /// @param other The other entity representation.
	    /// @return true if \b this and the \b other element are equal.
	    bool operator==(const EntityRep& other) const
	    {
		return entityrep_ == other.entityrep_;
	    }

	    /// @brief Inequality operator.
	    /// @param other The other entity representation.
	    /// @return true if \b this and the \b other element are \b not equal.
	    bool operator!=(const EntityRep& other) const
	    {
		return !operator==(other);
	    }

	    enum { InvalidIndex = INT_MAX };

	protected:
	    // Interior representation is documented in class main comment.
	    int entityrep_;
	};



	/// @brief Base class for EntityVariable and SignedEntityVariable.
	/// Forwards a restricted subset of the std::vector interface.
	/// @tparam T A value type for the variable,
	/// such as double for pressure etc.
	template <typename T>
	class EntityVariableBase : private std::vector<T>
	{
	public:
	    typedef std::vector<T> V;
	    using V::empty;
	    using V::size;
	    using V::assign;
	    /// Default constructor.
	    EntityVariableBase()
	    {
	    }

	protected:
	    const T& get(int i) const
	    {
		return V::operator[](i);
	    }
	};




	/// @brief A class design to hold a variable with a value for
	/// each entity of the given codimension, where the variable
	/// is \b not changing in sign with orientation. Examples include
	/// pressures and positions.
	/// @tparam T A value type for the variable,
	///           such as double for pressure etc.
	/// @tparam codim Codimension.
	template <typename T, int codim>
	class EntityVariable : public EntityVariableBase<T>
	{
	public:
	    /// Default constructor.
	    EntityVariable()
	    {
	    }
	    /// @brief Random access to the variable through an EntityRep.
	    /// @param e Entity representation.
	    /// @return a const reference to the varable, at e.
	    const T& operator[](const EntityRep<codim>& e) const
	    {
		return get(e.index());
	    }
	};





	/// @brief A class design to hold a variable with a value for
	/// each entity of the given codimension, where the variable
	/// \b is changing in sign with orientation. An example is
	/// velocity fields.
	/// @tparam T A value type for the variable,
	///           such as double for pressure etc.
	/// @tparam codim Codimension.
	template <typename T, int codim>
	class SignedEntityVariable : public EntityVariableBase<T>
	{
	public:
	    /// Default constructor.
	    SignedEntityVariable()
	    {
	    }
	    /// @brief Random access to the variable through an EntityRep.
	    /// Note that this operator always returns a copy, not a
	    /// reference, since we may need to flip the sign.
	    const T operator[](const EntityRep<codim>& e) const
	    {
		return e.orientation() ? get(e.index()) : -get(e.index());
	    }
	};


    } // namespace cpgrid
} // namespace Dune




#endif // OPENRS_ENTITYREP_HEADER
