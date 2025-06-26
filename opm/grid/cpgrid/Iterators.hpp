//===========================================================================
//
// File: Iterators.hpp
//
// Created: Fri May 29 23:29:09 2009
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

#ifndef OPM_ITERATORS_HEADER
#define OPM_ITERATORS_HEADER

#include <dune/grid/common/gridenums.hh>
#include "PartitionIteratorRule.hpp"
#include <opm/common/ErrorMacros.hpp>
#include "CpGridData.hpp"


#include <stack>

namespace Dune
{
    namespace cpgrid
    {
        class CpGridData;


        /// Iterator intended to be used as LeafIterator and LevelIterator
        /// (no difference due to no adaptivity) for CpGrid.
        /// This could have been a random access iterator, perhaps we will
        /// use a facade to do this later.
        template<int cd, PartitionIteratorType pitype>
        class Iterator : public Entity<cd>
        {
        public:
            using Reference = const Entity<cd>&;
            /// @brief
            /// @todo Doc me!
            /// @param
            Iterator(const CpGridData& grid, int index, bool orientation);

            Iterator() = default;

            /// Increment operator.
            /// Implementation note: This class is a friend of
            /// \see EntityRep (which is a private base class of
            /// Entity) in order to actually access the private
            /// variable entityrep_. We may want to change EntityRep,
            /// then this must change, too.
            Iterator& operator++()
            {
                EntityRep<cd>::increment();
                if(rule_.fullSet || rule_.emptySet)
                    return *this;
                while(this->index()<noEntities_ && rule_.isInvalid(*this))
                    EntityRep<cd>::increment();
                return *this;
            }

            Iterator operator++(int)
            {
                Iterator tmp(*this);
                ++(*this);
                return tmp;
            }

            /// Const member by pointer operator.
            const Entity<cd>* operator->() const
            {
                assert(Entity<cd>::isValid());
                return (this);
            }

            /// Const dereferencing operator.
            const Entity<cd>& operator*() const
            {
                assert(Entity<cd>::isValid());
                return (*this);
            }

        private:
            /// \brief The number of Entities with codim cd.   (no = number)
            int noEntities_;
            PartitionIteratorRule<pitype> rule_;
        };




        /// Only needs to provide interface for doing nothing.
    class HierarchicIterator
        {

        public:
            using Reference = const Entity<0>&;
            /// @brief
            /// @todo Doc me!
            /// @param
            explicit HierarchicIterator(const CpGridData& grid)
                : virtualEntity_(grid, EntityRep<0>::InvalidIndex, true )
            {
            }

            explicit HierarchicIterator() = default;

            // Constructor with Entity<0> target and maxLevel (begin iterator).
            HierarchicIterator(Entity<0> target, int maxLevel)
                : virtualEntity_(target), maxLevel_(maxLevel)
            {
                // Load sons of target onto the iterator stack
                stackChildren_(target);

                // Set entity target to the next child if exists
                resetEntity_();
            }


            // Constructor without valid element (end iterator).
            explicit HierarchicIterator(int maxLevel)
                : maxLevel_(maxLevel)
            {
                resetEntity_();
            }

            /// Equality.
            bool operator==(const HierarchicIterator& other) const
            {
                return virtualEntity_ == other.virtualEntity_;
            }

            /// Inequality.
            bool operator!=(const HierarchicIterator& other) const
            {
                return !this->operator==(other);
            }

            /// @brief
            /// @todo Doc me!
            /// @param
            HierarchicIterator& operator++()
            {
                if (elemStack_.empty()){
                    return *this;
                }
                // Reference to the top element of elemStack_
                auto target = elemStack_.top();
                // Remove the element on top of elemStack_
                elemStack_.pop();
                // Load sons of previous target onto elemStack_
                stackChildren_(target);
                // Set entity target to the next stacked element if exists
                resetEntity_();
                return *this;
            }

            /// @brief
            /// @todo Doc me!
            /// @param
            HierarchicIterator operator++(int)
            {
                if (elemStack_.empty()){
                    return *this;
                }
                auto copy = *this;
                // Reference to the top element of elemStack_
                auto target = elemStack_.top();
                // Remove the element on top of elemStack_
                elemStack_.pop();
                // Load sons of previous target onto elemStack_
                stackChildren_(target);
                // Set entity target to the next stacked element if exists
                resetEntity_();
                return copy;
            }

            /// Const member by pointer operator.
            const Entity<0>* operator->() const
            {
                assert(this -> virtualEntity_.isValid());
                return &virtualEntity_;
            }

            /// Const dereferencing operator.
            const Entity<0>& operator*() const
            {
                assert(this-> virtualEntity_.isValid());
                return virtualEntity_;
            }

        private:
            void stackChildren_(const Entity<0>& target);

            void resetEntity_();

            Entity<0> virtualEntity_;

            //! max level to iterate over
            int maxLevel_;

            // For depth-first search
            std::stack<Entity<0>> elemStack_;

        }; // end class HierarchicIterator

    } // namespace cpgrid
} // namespace Dune

namespace std
{
    template< int codim, Dune::PartitionIteratorType pitype >
    struct iterator_traits< Dune::cpgrid::Iterator< codim, pitype > >
    {
        typedef Dune::cpgrid::Iterator< codim, pitype >     Iterator;
        typedef ptrdiff_t                                   difference_type;
        typedef typename Iterator::Entity                   value_type;
        typedef value_type*                                 pointer;
        typedef value_type&                                 reference;
        typedef forward_iterator_tag                        iterator_category;
    };

    template <>
    struct iterator_traits< Dune::cpgrid::HierarchicIterator >
    {
        typedef ptrdiff_t                                   difference_type;
        typedef Dune::cpgrid::Entity<0>                     value_type;
        typedef value_type*                                 pointer;
        typedef value_type&                                 reference;
        typedef forward_iterator_tag                        iterator_category;
    };

} // namespace std


#include <opm/grid/cpgrid/CpGridData.hpp>
#include "Entity.hpp"

namespace Dune {
namespace cpgrid {

template<int cd, PartitionIteratorType pitype>
Iterator<cd, pitype>::Iterator(const CpGridData& grid, int index, bool orientation)
    : Entity<cd>(grid,
                 // If the partition is empty, goto to end iterator!
                 EntityRep<cd>(PartitionIteratorRule<pitype>::emptySet?grid.size(cd):index,
                               orientation)),
      noEntities_(grid.size(cd))
{
    if(rule_.fullSet || rule_.emptySet)
        return;

    while(this->index()<noEntities_ && rule_.isInvalid(*this))
        EntityRep<cd>::increment();
}
}}


#endif // OPM_ITERATORS_HEADER
