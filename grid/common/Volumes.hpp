//===========================================================================
//
// File: Volumes.hpp
//
// Created: Mon Jun 22 15:46:32 2009
//
// Author(s): Jan B Thomassen <jbt@sintef.no>
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

#ifndef OPENRS_VOLUMES_HEADER
#define OPENRS_VOLUMES_HEADER

#include <numeric>

#include <dune/common/misc.hh>
#include <dune/common/fvector.hh>

namespace Dune
{


    template <typename T>
    FieldVector<T, 3> cross(const FieldVector<T, 3>& a, const FieldVector<T, 3>& b)
    {
	FieldVector<T, 3> res;
	res[0] = a[1]*b[2] - a[2]*b[1];
	res[1] = a[2]*b[0] - a[0]*b[2];
	res[2] = a[0]*b[1] - a[1]*b[0];
	return res;
    }

    template <class Vector>
    typename Vector::value_type inner(const Vector& a, const Vector& b)
    {
	return std::inner_product(a.begin(), a.end(), b.begin(), typename Vector::value_type());
    }

    /// Calculates the determinant of a 2 x 2 matrix, represented in memory as an
    /// array of two-dimensional points.  Same function also exists for 3 x 3
    /// matrices.
    template<typename T, template <typename, int> class Point>
    inline T determinantOf(const Point<T, 2>* a) 
    {
	return a[0][0] * a[1][1] - a[1][0] * a[0][1];
    };


    /// Calculates the determinant of a 3 x 3 matrix, represented in memory as an
    /// array of three-dimensional points.  Same function also exists for 2 x 2
    /// matrices.
    template<typename T, template <typename, int> class Point>
    inline T determinantOf(const Point<T, 3>* a) 
    {
	return 
	    a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2]) -
	    a[0][1] * (a[1][0] * a[2][2] - a[2][0] * a[1][2]) +
	    a[0][2] * (a[1][0] * a[2][1] - a[2][0] * a[1][1]);
    };


    /// Computes the volume of a simplex consisting of (Dim+1) vertices embedded
    /// in Euclidean space of dimension (Dim)
    template<typename T, template <typename, int> class Point, int Dim>
    inline T simplex_volume(const Point<T, Dim>* a)
    {
	Point<T, Dim> tmp[Dim];
	for (int i = 0; i < Dim; ++i) {
	    tmp[i] = a[i+1] - a[i];
	}
	return determinantOf(tmp) / double(Factorial<Dim>::factorial);
	// determinant / factorial
    }


    /// Computes the area of a 2-dimensional triangle.  Input is an array of
    /// corner points.  Same function also exists for 3-dimensional triangles.
    template <typename T, template <typename, int> class Point>
    inline T area(const Point<T, 2>* c)
    { return simplex_volume(c); }


    /// Computes the area of a 3-dimensional triangle.  Input is an array of
    /// corner points.  Same function also exists for 2-dimensional triangles.
    template <typename T, template <typename, int> class Point>
    inline T area(const Point<T, 3>* c)
    {
	// Using the one-half cross product rule
	Point<T, 3> d0 = c[1] - c[0];
	Point<T, 3> d1 = c[2] - c[0];
	Point<T, 3> crossprod = cross(d0,d1);
	return 0.5 * crossprod.two_norm();
    }


    /// Computes the volume of a 3D simplex (embedded i 3D space).
    template <typename T, template <typename, int> class Point>
    inline T volume(const Point<T, 3>* c)
    { return simplex_volume(c); }


    /// Computes the signed area of a triangle embedded in 3D space. Input is an
    /// array of corner points and a normal to determine the sign.
    template <typename T, template <typename, int> class Point>
    T signed_area(const Point<T, 3>* c, const Point<T, 3>& normal)
    {
	// Using the one-half cross product rule
	Point<T, 3> d0 = c[1] - c[0];
	Point<T, 3> d1 = c[2] - c[0];
	Point<T, 3> crossprod = cross(d0, d1);
	if (inner(crossprod, normal) > 0) {
	    return 0.5 * crossprod.two_norm();
	} else {
	    return -0.5 * crossprod.two_norm();
	}
    }


} // namespace Dune



#endif // OPENRS_VOLUMES_HEADER
