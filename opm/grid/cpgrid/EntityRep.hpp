//===========================================================================
//
// File: EntityRep.hpp
//
// Created: Tue Jun  9 11:11:24 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            B�rd Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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

#ifndef OPM_ENTITYREP_HEADER
#define OPM_ENTITYREP_HEADER


// -------------------------------------------------------------------
// -> Layering violation. --------------------------------------------
//
//    We need a unary operator-() for class Dune::FieldVector<K,n>
//    within method Dune::SignedEntityVariable<T,codim>::operator[](),
//    but Dune::FieldVector<K,n> does not provide such an operator.
//

// Warning suppression for Dune includes.
#include <opm/grid/utility/platform_dependent/disable_warnings.h>
#include <dune/common/fvector.hh>
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>

namespace Dune
{
    template<typename K, int n>
    FieldVector<K,n>
    operator- (const FieldVector<K,n>& v)
    {
        // Assume 'K' supports a single parameter constructor.  The
        // assumption holds for all standard C++ built-in arithmetic
        // types such as 'int', 'float', and 'complex<double>'.
        //
        return FieldVector<K,n>(K(0)) - v;
    }
}
//
// <- Layering violation. --------------------------------------------
// -------------------------------------------------------------------


//#include <opm/core/utility/SparseTable.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <climits>
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
            enum{ codimension=codim};

            /// Default constructor.
            EntityRep()
                : entityrep_(0)
            {
            }
            /// @brief Constructor taking an entity index and an orientation.
            /// @param index_arg Entity index
            /// @param orientation_arg True if the entity's orientation is positive.
            EntityRep(int index_arg, bool orientation_arg)
                : entityrep_(orientation_arg ? index_arg : ~index_arg)
            {
                assert(index_arg >= 0);
            }
            /// @brief Set entity value.
            /// @param index_arg Entity index
            /// @param orientation_arg True if the entity's orientation is positive.
            void setValue(int index_arg, bool orientation_arg)
            {
                assert(index_arg >= 0);
                entityrep_ = orientation_arg ? index_arg : ~index_arg;
            }
            /// @brief The (positive) index of an entity. Not a Dune interface method.
            /// @return the (positive) index of an entity.
            int index() const
            {
                return entityrep_ < 0 ? ~entityrep_ : entityrep_;
            }

            /// @brief The signed index that also tells us the orientation
            int signedIndex() const
            {
                return entityrep_;
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

            /// @brief Increments the entityrep's index() by one.
            void increment()
            {
                if (entityrep_ < 0) {
                    --entityrep_;
                } else {
                    ++entityrep_;
                }
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

        private:
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
            friend class CpGridData;
        public:
            typedef std::vector<T> V;
            typedef typename std::vector<T>::iterator iterator;
            typedef typename std::vector<T>::const_iterator const_iterator;

            using V::empty;
            using V::size;
            using V::assign;
            using V::begin;
            using V::end;
            using typename V::value_type;
            using V::reserve;
            using V::push_back;
            using V::data;
            using V::operator[];
            using V::resize;

            /// Default constructor.
            EntityVariableBase()
            {
            }

            const T& get(int i) const
            {
                return V::operator[](i);
            }

            T& get(int i)
            {
                return V::operator[](i);
            }

            void swap(EntityVariableBase& other)
            {
                V::swap(static_cast<V&>(other));
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
                return EntityVariableBase<T>::get(e.index());
            }
            /// @brief Random access to the variable through an EntityRep.
            /// @param e Entity representation.
            /// @return a mutable reference to the varable, at e.
            T& operator[](const EntityRep<codim>& e)
            {
                return EntityVariableBase<T>::get(e.index());
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
                return e.orientation() ?
                    EntityVariableBase<T>::get(e.index()) :
                    -EntityVariableBase<T>::get(e.index());
            }
        };


    } // namespace cpgrid
} // namespace Dune




#endif // OPM_ENTITYREP_HEADER
