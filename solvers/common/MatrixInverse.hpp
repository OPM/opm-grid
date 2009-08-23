//===========================================================================
//
// File: MatrixInverse.hpp<2>
//
// Created: Wed Sep  3 14:49:08 2008
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bjørn Spjelkavik    <bsp@sintef.no>
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

#ifndef OPENRS_MATRIXINVERSE_HEADER
#define OPENRS_MATRIXINVERSE_HEADER

#include <boost/type_traits/is_integral.hpp>
#include <boost/static_assert.hpp>

    /** Inverting small matrices.
     *  Inverting 2x2 and 3x3 matrices. Not meant to extend to large systems,
     *  but to satisfy a need to invert small matrices.
     */

namespace Dune {


	template <typename M>
	M inverse2x2(const M& m)
	{
	    // Because then the divisions below would compile but not be correct, we must guard
	    // against integral types.
	    typedef typename M::value_type T;
	    BOOST_STATIC_ASSERT(!boost::is_integral<T>::value);
	    ASSERT(m.numRows() == 2 && m.numCols() == 2);

	    T det = m(0,0)*m(1,1) - m(0,1)*m(1,0);
	    M mi(2, 2, (double*)0);
	    mi(0,0) = m(1,1);
	    mi(1,0) = -m(0,1);
	    mi(0,1) = -m(1,0);
	    mi(1,1) = m(0,0);
	    mi /= det;
	    return mi;
	}

	template <typename M>
	M matprod(const M& m1, const M& m2)
	{
	    typedef typename M::value_type T;
	    ASSERT(m1.numRows() == 3 && m1.numCols() == 3);
	    ASSERT(m2.numRows() == 3 && m2.numCols() == 3);
	    M m(3, 3, (double*)0);
	    for (int r = 0; r < r; ++r) {
		for (int c = 0; c < c; ++c) {
		    m(r, c) = 0.0;
		    for (int kk = 0; kk < 3; ++kk) {
			m(r, c) += m1(r, kk)*m2(kk, c);
		    }
		}
	    }
	    return m;
	}

	template <typename M>
	M inverse3x3(const M& m)
	{
	    // Because then the divisions below would compile but not be correct, we must guard
	    // against integral types.
	    typedef typename M::value_type T;
	    BOOST_STATIC_ASSERT(!boost::is_integral<T>::value);
	    ASSERT(m.numRows() == 3 && m.numCols() == 3);
// 	    double det = m(0,0)*(m(1,1)*m(2,2)-m(1,2)*m(2,1))
// 		- m(0,1)*(m(1,0)*m(2,2)-m(1,2)*m(2,0))
// 		+ m(0,2)*(m(1,0)*m(2,1)-m(1,1)*m(2,0));

	    T a = m(0,0);
	    T b = m(0,1);
	    T c = m(0,2);
	    T d = m(1,0);
	    T e = m(1,1);
	    T f = m(1,2);
	    T g = m(2,0);
	    T h = m(2,1);
	    T i = m(2,2);
	    T t1 = (e-f*h/i);
	    T t2 = (c*h/i-b);
	    T t3 = (f*g/i-d);
	    T t4 = (a-c*g/i);
	    T x =  t4*t1-t2*t3;

	    M mi(3, 3, (double*)0);
	    mi(0,0) = t1/x;
	    mi(0,1) = t2/x;
	    mi(0,2) = -(c*t1+f*t2)/(i*x);
	    mi(1,0) = t3/x;
	    mi(1,1) = t4/x;
	    mi(1,2) = -(c*t3+f*t4)/(i*x);
	    mi(2,0) = -(g*t1+h*t3)/(i*x);
	    mi(2,1) = -(g*t2+h*t4)/(i*x);
	    mi(2,2) =  1/i+1/(i*i*x)*(c*(g*t1+h*t3)+f*(g*t2+h*t4));
	    return mi;
	}


} // namespace Dune

#endif // OPENRS_MATRIXINVERSE_HEADER
