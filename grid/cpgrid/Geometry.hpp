//===========================================================================
//
// File: Geometry.hpp
//
// Created: Fri May 29 23:29:24 2009
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

#ifndef OPENRS_GEOMETRY_HEADER
#define OPENRS_GEOMETRY_HEADER

namespace Dune
{
    namespace cpgrid
    {

	template <int dim, int dimworld>
	class Geometry
	{
	public:
	    typedef double ct;
	    typedef FieldVector<ct, dimworld> WorldPointType;
	    typedef FieldVector<ct, dim> LocalPointType;

	    Geometry(const WorldPointType& pos, ct vol)
		: pos_(pos), vol_(vol)
	    {
	    }

	    const WorldPointType& global(const LocalPointType&) const
	    {
		return pos_;
	    }

	    double integrationElement(const LocalPointType&) const
	    {
		return vol_;
	    }

	    GeometryType type() const
	    {
		GeometryType t;
		t.makeSingular(dim);
		return t;
	    }

	    /// The corner method requires the points, which we may not necessarily want to provide.
	    /// We will need it for visualization purposes though. For now returning a default point.
	    WorldPointType corner(int i) const
	    {
		return pos_;
	    }

	    ct volume() const
	    {
		return vol_;
	    }

	private:
	    WorldPointType pos_;
	    double vol_;
	};


    } // namespace cpgrid
} // namespace Dune

#endif // OPENRS_GEOMETRY_HEADER
