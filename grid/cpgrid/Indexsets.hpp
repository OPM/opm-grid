//===========================================================================
//
// File: Indexsets.hpp
//
// Created: Fri May 29 23:30:01 2009
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

#ifndef OPENRS_INDEXSETS_HEADER
#define OPENRS_INDEXSETS_HEADER

#include <dune/common/geometrytype.hh>
#include "../common/ErrorMacros.hpp"
namespace Dune
{
    namespace cpgrid
    {

	template <class GridType>
	class IndexSet
	{
	public:
	    // typedef GridType::GeometryType GeometryType;
	    IndexSet(const GridType& grid)
		: grid_(grid)
	    {
		GeometryType t;
		t.makeSingular(3);
		gt_.push_back(t);
	    }

	    const std::vector<GeometryType>& geomTypes(int /*codim*/) const
	    {
		return gt_;
	    }

	    int size(GeometryType type) const
	    {
		if (!type.isSingular()) {
		    THROW("IndexSet::size(GeometryType) not implemented for its dim != 3 types.");
		}
		return grid_.size(0);
	    }

	    int size(int codim) const
	    {
		return grid_.size(codim);
	    }

	    template<int cd>
	    int index(const typename GridType::template Codim<cd>::Entity& e) const 
	    {
		return e.index(); 
	    }

	    template<class EntityType>
	    int index(const EntityType& e) const 
	    {
		return e.index();
	    }

// 	    template<class EntityType>
// 	    int subIndex(const EntityType& e, int i) const 
// 	    {
// 		return grid_.cell_to_face_[e][i].index();
// 	    }

	    template <int cc>
	    int subIndex(const typename GridType::Traits::template Codim<0>::Entity& e, int i) const 
	    {
		if (cc == 1) {
		    return grid_.cell_to_face_[e][i].index();
		} else {
		    THROW("No entities defined for codimension 2 and up. Cannot evaluate subIndex().");
		}
	    }

	private:
	    const GridType& grid_;
	    std::vector<GeometryType> gt_;
	};


	template <class GridType>
	class IdSet : private IndexSet<GridType>
	{
	private:
	    typedef IndexSet<GridType> super_t;
	public:
	    IdSet(const GridType& grid)
		: super_t(grid)
	    {
	    }

	    template<int cd>
	    int id(const typename GridType::template Codim<cd>::Entity& e) const 
	    {
		return super_t::index<cd>(e);
	    }

	    template<class EntityType>
	    int id(const EntityType& e) const 
	    {
		return super_t::index(e);
	    }

	};


    } // namespace cpgrid
} // namespace Dune

#endif // OPENRS_INDEXSETS_HEADER
