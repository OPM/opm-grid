//===========================================================================
//
// File: GeometryHelpers.hpp
//
// Created: Mon Jun 22 13:43:23 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Halvor M Nilsen     <hnil@sintef.no>
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

#ifndef OPENRS_GEOMETRYHELPERS_HEADER
#define OPENRS_GEOMETRYHELPERS_HEADER

#include "ErrorMacros.hpp"
#include "Volumes.hpp"

namespace Dune
{

    namespace GeometryHelpers
    {

	template <class Point, template <class> class Vector>
	Point average(const Vector<Point>& points)
	{
	    int num_points = points.size();
	    ASSERT(num_points > 0);
	    Point pt = points[0];
	    for (int i = 1; i < num_points; ++i) {
		pt += points[i];
	    }
	    pt /= double(num_points);
	    return pt;
	}


	template <class Point, template <class> class Vector>
	double polygonArea(const Vector<Point>& points,
			   const Point& centroid)
	{
	    double tot_area = 0.0;
	    int num_points = points.size();
	    for (int i = 0; i < num_points; ++i) {
		Point tri[3] = { centroid, points[i], points[(i+1)%num_points] };
		tot_area += area(tri);
	    }
	    return tot_area;
	}



	template <class Point, template <class> class Vector>
	Point polygonCentroid(const Vector<Point>& points,
			      const Point& inpoint)
	{
	    double tot_area = 0.0;
	    Point tot_centroid(0.0);
	    int num_points = points.size();
	    for (int i = 0; i < num_points; ++i) {
		Point tri[3] = { inpoint, points[i], points[(i+1)%num_points] };
		double tri_area = area(tri);
		Point tri_w_mid = (tri[0] + tri[1] + tri[2]);
		tri_w_mid *= tri_area/3.0;
		tot_area += tri_area;
		tot_centroid += tri_w_mid;
	    }
	    tot_centroid /= tot_area;	
	    return tot_centroid;
	}



	template <class Point, template <class> class Vector>
	Point polygonNormal(const Vector<Point>& points,
			    const Point& centroid)
	{
	    Point tot_normal(0.0);
	    int num_points = points.size();
	    for (int i = 0; i < num_points; ++i) {
		Point tri[3] = { centroid, points[i], points[(i+1)%num_points] };
		Point d0 = tri[1] - tri[0];
		Point d1 = tri[2] - tri[0];
		Point w_normal = cross(d0, d1);
		w_normal *= area(tri);
		tot_normal += w_normal;
	    }
	    tot_normal /= tot_normal.two_norm();	
	    return tot_normal;
	}



	template <class Point, template <class> class Vector>
	double polygonCellVolume(const Vector<Point>& points,
				 const Point& face_centroid,
				 const Point& cell_centroid)
	{
	    double tot_volume = 0.0;
	    int num_points = points.size();
	    for (int i = 0; i < num_points; ++i) {
		Point tet[4] = { cell_centroid, face_centroid, points[i], points[(i+1)%num_points] };
		double small_volume = simplex_volume(tet);
		small_volume *= -1.0;
		ASSERT(small_volume > 0);
		tot_volume += small_volume;
	    }
	    ASSERT(tot_volume>0);
	    return tot_volume;
	}



	template <class Point, template <class> class Vector>
	Point polygonCellCentroid(const Vector<Point>& points,
				  const Point& face_centroid,
				  const Point& cell_centroid)
	{
	    Point centroid(0.0);
	    double tot_volume = 0.0;
	    int num_points = points.size();
	    for (int i = 0; i < num_points; ++i) {
		Point tet[4] = { cell_centroid, face_centroid, points[i], points[(i+1)%num_points] };
		double small_volume = simplex_volume(tet);
		small_volume *= -1.0;
		ASSERT(small_volume > 0);
		Point small_centroid = tet[0];
		for(int i = 1; i < 4; ++i){
		    small_centroid += tet[i];
		}
		centroid += small_centroid*small_volume/4.0;
		ASSERT(small_volume >0);
		tot_volume += small_volume;
	    }
	    centroid /= tot_volume;
	    ASSERT(tot_volume>0);
	    return centroid;
	
	}


    } // namespace GeometryHelpers

} // namespace Dune


#endif // OPENRS_GEOMETRYHELPERS_HEADER
