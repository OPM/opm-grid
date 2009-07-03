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
	    typedef int IndexType;

	    IndexSet(const GridType& grid)
		: grid_(grid)
	    {
		GeometryType t;
		t.makeSingular(3);
		geom_types_.push_back(t);
	    }

	    const std::vector<GeometryType>& geomTypes(int codim) const
	    {
		return codim == 0 ? geom_types_ : empty_;
	    }

	    int size(GeometryType type) const
	    {
		return grid_.size(type);
	    }

	    int size(int codim) const
	    {
		return grid_.size(codim);
	    }

	    template<int cd>
	    IndexType index(const typename GridType::template Codim<cd>::Entity& e) const 
	    {
		return e.index(); 
	    }

	    template<class EntityType>
	    IndexType index(const EntityType& e) const 
	    {
		return e.index();
	    }

// 	    template<class EntityType>
// 	    int subIndex(const EntityType& e, int i) const 
// 	    {
// 		return grid_.cell_to_face_[e][i].index();
// 	    }

	    template <int cc>
	    IndexType subIndex(const typename GridType::template Codim<0>::Entity& e, int i) const 
	    {
		return index(e.subEntity<cc>(i));
	    }

	    IndexType subIndex(const typename GridType::template Codim<0>::Entity& e, int i, unsigned int cc) const 
	    {
		switch(cc) {
		case 0: return index(e.template subEntity<0>(i));
		case 1: return index(e.template subEntity<1>(i));
		case 2: return index(e.template subEntity<2>(i));
		case 3: return index(e.template subEntity<3>(i));
		default: THROW("Codimension " << cc << " not supported.");
		}
	    }

	    template <class EntityType>
	    bool contains(const EntityType& e) const
	    {
		return EntityType::codimension == 0;
	    }

	private:
	    const GridType& grid_;
	    std::vector<GeometryType> geom_types_;
	    std::vector<GeometryType> empty_;
	};


	template <class GridType>
	class IdSet
	{
	public:
	    typedef int IdType;

	    IdSet(const GridType& grid)
		: grid_(grid)
	    {
		cumul_sizes[0] = 0;
		cumul_sizes[1] = cumul_sizes[0] + grid.size(0);
		cumul_sizes[2] = cumul_sizes[1] + grid.size(1);
		cumul_sizes[3] = cumul_sizes[2] + grid.size(2);
	    }

	    template<int cc>
	    IdType id(const typename GridType::template Codim<cc>::Entity& e) const 
	    {
		return id(e);
	    }

	    template<class EntityType>
	    IdType id(const EntityType& e) const 
	    {
		return cumul_sizes[EntityType::codimension] + e.index();
	    }

	    template<int cc>
	    IdType subId(const typename GridType::template Codim<0>::Entity& e, int i) const 
	    {
		return id(e.subEntity<cc>(i));
	    }

	    IdType subId(const typename GridType::template Codim<0>::Entity& e, int i, int cc) const
	    {
		switch (cc) {
		case 0: return id(e.template subEntity<0>(i));
		case 1: return id(e.template subEntity<1>(i));
		case 2: return id(e.template subEntity<2>(i));
		case 3: return id(e.template subEntity<3>(i));
		default: THROW("Cannot get subId of codimension " << cc);
		}
		return -1;
	    }
	private:
	    const GridType& grid_;
	    int cumul_sizes[4];
	};


    } // namespace cpgrid
} // namespace Dune

#endif // OPENRS_INDEXSETS_HEADER
