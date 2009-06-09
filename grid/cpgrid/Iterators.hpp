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

#ifndef OPENRS_ITERATORS_HEADER
#define OPENRS_ITERATORS_HEADER

#include <dune/grid/common/gridenums.hh>
#include "Entity.hpp"
#include "../common/ErrorMacros.hpp"

namespace Dune
{
    namespace cpgrid
    {




	template<int cd, PartitionIteratorType pitype, class GridType>
	class Iterator : public EntityPointer<cd, GridType>
	{
	public:
	    Iterator(const GridType& grid, int index)
		: EntityPointer<cd, GridType>(grid, index)
	    {
	    }
	    Iterator& operator++()
	    {
		ASSERT((Entity<cd, GridType>::entityrep_) >= 0);
		++Entity<cd, GridType>::entityrep_;
		return *this;
	    }
	};





	class HierarchicIterator
	{
	};





    } // namespace cpgrid
} // namespace Dune

#endif // OPENRS_ITERATORS_HEADER
