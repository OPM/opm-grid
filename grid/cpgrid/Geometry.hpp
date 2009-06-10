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

	    enum { dimension = 3 };
	    enum { mydimension = dim };
	    enum { coorddimension = 3 };
	    enum { dimensionworld = 3 };

	    typedef double ctype;

	    typedef FieldVector<ctype, dimworld> WorldPointType;
	    typedef FieldVector<ctype, dim> LocalPointType;

	    Geometry(const WorldPointType& pos, ctype vol)
		: pos_(pos), vol_(vol)
	    {
	    }

	    const WorldPointType& global(const LocalPointType&) const
	    {
		return pos_;
	    }

	    LocalPointType local(const LocalPointType&) const
	    {
		LocalPointType dummy(0.0);
		return dummy;
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

	    /// The number of corners of this convex polytope.
	    /// Returning 1 as long as we are using the singular geometry/refelem approach.
	    int corners() const
	    {
		return 1;
	    }

	    /// The corner method requires the points, which we may not necessarily want to provide.
	    /// We will need it for visualization purposes though. For now returning a single "corner".
	    WorldPointType corner(int i) const
	    {
		return pos_;
	    }

	    ctype volume() const
	    {
		return vol_;
	    }


	    const FieldMatrix<ctype, mydimension, coorddimension>&
	    jacobianTransposed(const LocalPointType& local) const
	    {
		THROW("Meaningless to call jacobianTransposed() on singular geometries.");
		static FieldMatrix<ctype, mydimension, coorddimension> dummy;
		return dummy;
	    }


	    const FieldMatrix<ctype, coorddimension, mydimension>&
	    jacobianInverseTransposed(const LocalPointType& local) const
	    {
		THROW("Meaningless to call jacobianInverseTransposed() on singular geometries.");
		static FieldMatrix<ctype, coorddimension, mydimension> dummy;
		return dummy;
	    }

	private:
	    WorldPointType pos_;
	    double vol_;
	};


    } // namespace cpgrid
} // namespace Dune

#endif // OPENRS_GEOMETRY_HEADER
