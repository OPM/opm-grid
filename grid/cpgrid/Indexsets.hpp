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

namespace Dune
{
    namespace cpgrid
    {

	template <class GridType>
	class IndexSet
	{
	public:
	    // typedef GridType::GeometryType GeometryType;
	    IndexSet(const GridType* grid)
		: grid_(grid)
	    {
	    }

	    const std::vector<GeometryType>& geomTypes(int /*codim*/) const
	    {
		return gt_;
	    }

	    int size (GeometryType type) const
	    {
		return 0;
	    }

	    int size (int codim) const
	    {
		return grid_->size(codim);
	    }

	private:
	    const GridType* grid_;
	    std::vector<GeometryType> gt_;
	};


	template <class GridType>
	class IdSet : public IndexSet<GridType>
	{
	public:
	    IdSet(const GridType* grid)
		: IndexSet<GridType>(grid)
	    {
	    }
	};


    } // namespace cpgrid
} // namespace Dune

#endif // OPENRS_INDEXSETS_HEADER
