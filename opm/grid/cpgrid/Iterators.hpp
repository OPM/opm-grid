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

#ifndef OPM_ITERATORS_HEADER
#define OPM_ITERATORS_HEADER

#include <dune/grid/common/gridenums.hh>
#include "PartitionIteratorRule.hpp"
#include <opm/grid/utility/ErrorMacros.hpp>

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
        class Iterator : public EntityPointer<cd>
        {
        public:
            /// @brief
            /// @todo Doc me!
            /// @param
            Iterator(const CpGridData& grid, int index, bool orientation);

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
        private:
            /// \brief The number of Entities with codim cd.
            int noEntities_;
            PartitionIteratorRule<pitype> rule_;
        };




        /// Only needs to provide interface for doing nothing.
        class HierarchicIterator : public EntityPointer<0>
        {
        public:
            /// @brief
            /// @todo Doc me!
            /// @param
            HierarchicIterator(const CpGridData& grid)
                : EntityPointer<0>(grid, EntityRep<0>::InvalidIndex, true )
            {
            }

            /// @brief
            /// @todo Doc me!
            /// @param
            HierarchicIterator& operator++()
            {
                OPM_THROW(std::runtime_error, "Calling operator++() on HierarchicIterator for CpGrid, which has no refinement.");
                return *this;
            }
        };





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
        typedef Dune::cpgrid::HierarchicIterator::Entity    value_type;
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
    : EntityPointer<cd>(grid,
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
