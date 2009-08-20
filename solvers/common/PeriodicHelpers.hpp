//===========================================================================
//
// File: PeriodicHelpers.hpp
//
// Created: Fri Aug 14 10:59:57 2009
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

#ifndef OPENRS_PERIODICHELPERS_HEADER
#define OPENRS_PERIODICHELPERS_HEADER

#include "BoundaryConditions.hpp"
#include <dune/common/array.hh>
#include <dune/common/fvector.hh>
#include <algorithm>
#include <boost/array.hpp>

namespace Dune
{


    /// @brief
    /// @todo Doc me!
    struct BoundaryFaceInfo
    {
	/// @brief
	/// @todo Doc me!
	int face_index;
	/// @brief
	/// @todo Doc me!
	int bid;
	/// @brief
	/// @todo Doc me!
	int canon_pos;
	/// @brief
	/// @todo Doc me!
	int partner_face_index;
	/// @brief
	/// @todo Doc me!
	int partner_bid;
	/// @brief
	/// @todo Doc me!
	double area;
	/// @brief
	/// @todo Doc me!
	FieldVector<double,3> centroid;

	/// @brief
	/// @todo Doc me!
	/// @param
	/// @return
	bool operator<(const BoundaryFaceInfo& other) const
	{
	    return cmpval() < other.cmpval();
	}

	/// @brief
	/// @todo Doc me!
	/// @return
	double cmpval() const
	{
	    return centroid[(canon_pos/2 + 1)%3] + M_PI*centroid[(canon_pos/2 + 2)%3];
	}
    };


    /// @brief
    /// @todo Doc me!
    /// @param
    /// @return
    bool match(std::vector<BoundaryFaceInfo>& bfaces, int face, int lower, int upper)
    {
	const double area_tol = 1e-6;
	const double centroid_tol = 1e-6;
	int cp = bfaces[face].canon_pos;
	int target_cp = (cp%2 == 0) ? cp + 1 : cp - 1;
	FieldVector<double, 3> cent_this = bfaces[face].centroid;
	for (int j = lower; j < upper; ++j) {
	    if (bfaces[j].canon_pos == target_cp) {
		if (fabs(bfaces[face].area - bfaces[j].area) <= area_tol) {
		    FieldVector<double, 3> cent_other = bfaces[j].centroid;
		    cent_other -= cent_this;
		    double dist = cent_other.two_norm();
		    if (dist <= centroid_tol) {
			bfaces[face].partner_face_index = bfaces[j].face_index;
			bfaces[face].partner_bid = bfaces[j].bid;
			bfaces[j].partner_face_index = bfaces[face].face_index;
			bfaces[j].partner_bid = bfaces[face].bid;
			break;
		    }
		}
	    }
	}
	return (bfaces[face].partner_face_index != -1);
    }


    /// @brief Makes a boundary condition object representing
    /// periodic boundary conditions in any cartesian directions.
    /// The grid interface needs to export boundary ids that are
    /// unique for each boundary face for this to be possible.
    /// @tparam BCs the boundary condition class
    /// @tparam GridInterface grid interface class
    template <class BCs, class GridInterface>
    void createPeriodic(BCs& fbcs,
			const GridInterface& g,
			const boost::array<FlowBC, 2*GridInterface::Dimension> conditions,
			double spatial_tolerance = 1e-6)
    {
	// Check the conditions given.
	for (int i = 0; i < GridInterface::Dimension; ++i) {
	    if (conditions[2*i].isPeriodic()) {
		ASSERT(conditions[2*i + 1].isPeriodic());
		ASSERT(conditions[2*i].pressureDifference() == -conditions[2*i + 1].pressureDifference());
	    }
	}

	// Pick out all boundary faces, simultaneously find the
	// bounding box of their centroids, and the max id.
	typedef typename GridInterface::CellIterator CI;
	typedef typename CI::FaceIterator FI;
	typedef typename GridInterface::Vector Vector;
	std::vector<FI> bface_iters;
	Vector low(1e100);
	Vector hi(-1e100);
	int max_bid = 0;
	for (CI c = g.cellbegin(); c != g.cellend(); ++c) {
	    for (FI f = c->facebegin(); f != c->faceend(); ++f) {
		if (f->boundaryId()) {
		    bface_iters.push_back(f);
		    Vector fcent = f->centroid();
		    for (int dd = 0; dd < GridInterface::Dimension; ++dd) {
			low[dd] = std::min(low[dd], fcent[dd]);
			hi[dd] = std::max(hi[dd], fcent[dd]);
		    }
		    max_bid = std::max(max_bid, f->boundaryId());
		}
	    }
	}
	int num_bdy = bface_iters.size();
	if (max_bid != num_bdy) {
	    THROW("createPeriodic() assumes that every boundary face has a unique boundary id. That seems to be violated.");
	}

	// Store boundary face info in a suitable structure. Also find side total volumes.
	array<double, 2*GridInterface::Dimension> side_areas;
	std::fill(side_areas.begin(), side_areas.end(), 0.0);
	std::vector<BoundaryFaceInfo> bfinfo;
	bfinfo.reserve(num_bdy);
	for (int i = 0; i < num_bdy; ++i) {
	    BoundaryFaceInfo bf;
	    bf.face_index = i;
	    bf.bid = bface_iters[i]->boundaryId();
	    bf.canon_pos = -1;
	    bf.partner_face_index = -1;
	    bf.partner_bid = 0;
	    bf.area = bface_iters[i]->area();
	    bf.centroid = bface_iters[i]->centroid();
	    for (int dd = 0; dd < GridInterface::Dimension; ++dd) {
		double coord = bf.centroid[dd];
		if (fabs(coord - low[dd]) <= spatial_tolerance) {
		    bf.canon_pos = 2*dd;
		    break;
		} else if (fabs(coord - hi[dd]) <= spatial_tolerance) {
		    bf.canon_pos = 2*dd + 1;
		    break;
		}
	    }
	    if (bf.canon_pos == -1) {
		std::cerr << "Centroid: " << bf.centroid << "\n";
		std::cerr << "Bounding box min: " << low << "\n";
		std::cerr << "Bounding box max: " << hi << "\n";
		THROW("Boundary face centroid not on bounding box. Maybe the grid is not an axis-aligned shoe-box?");
	    }
	    side_areas[bf.canon_pos] += bf.area;
	    bf.centroid[bf.canon_pos/2] = 0.0;
	    bfinfo.push_back(bf);
	}

	std::sort(bfinfo.begin(), bfinfo.end());

	// Identify partners.
	for (int i = 0; i < num_bdy; ++i) {
	    if (bfinfo[i].partner_face_index != -1) {
		continue;
	    }
	    if (!conditions[bfinfo[i].canon_pos].isPeriodic()) {
		continue;
	    }
	    int lower = std::max(0, i - 10);
	    int upper = std::min(num_bdy, i + 10);
	    bool ok = match(bfinfo, i, lower, upper);
	    if (!ok) {
		// We have not found a partner.
		ok = match(bfinfo, i, 0, num_bdy);
		if (!ok) {
		    MESSAGE("Warning: No partner found for boundary id " << bfinfo[i].bid);
		    // THROW("No partner found.");
		}
	    }
	}

	// Resize the conditions object. We clear it first to make sure it's all defaults.
	typedef FlowBC FBC;
	fbcs.clear();
	fbcs.resize(max_bid + 1);

	// Now we have all the info, and the object to store it in. So we store it...
	for (int i = 0; i < num_bdy; ++i) {
	    FBC master_bc = conditions[bfinfo[i].canon_pos];
	    FBC face_bc;
	    if (master_bc.isPeriodic()) {
		face_bc = FBC(FBC::Periodic, master_bc.pressureDifference());
		int bid1 = bfinfo[i].bid;
		int bid2 = bfinfo[i].partner_bid;
		if (bid1 < bid2) {
		    fbcs.setPeriodicPartners(bid1, bid2);
		}
	    } else if (master_bc.isNeumann()) {
		face_bc = FBC(FBC::Neumann, master_bc.outflux()*bfinfo[i].area/side_areas[bfinfo[i].canon_pos]);
	    } else if (master_bc.isDirichlet()) {
		face_bc = FBC(FBC::Dirichlet, master_bc.pressure());
	    } else {
		THROW("Boundary condition is unknown.");
	    }
	    fbcs.flowCond(bfinfo[i].bid) = face_bc;
	}

	// Check that all boundary faces were visited.
// 	for (int i = 1; i <= max_bid; ++i) {
// 	    if (fbcs[i].isDirichlet()) {
// 		if (fbcs[i].pressure() == insane_pressure) {
// 		    THROW("Not all boundary faces have been visited.");
// 		}
// 	    }
// 	}
    }



} // namespace Dune

#endif // OPENRS_PERIODICHELPERS_HEADER
