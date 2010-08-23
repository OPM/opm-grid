//===========================================================================
//
// File: EclipseGridInspector.C
//
// Created: Mon Jun  2 12:17:51 2008
//
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//
// $Date$
//
// $Revision$
//
// Revision: $Id: EclipseGridInspector.C,v 1.2 2008/08/18 14:16:13 atgeirr Exp $
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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

#if HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdexcept>
#include <numeric>
#include <cmath>
#include <cfloat>
#include <algorithm>
#include "EclipseGridInspector.hpp"
#include "EclipseGridParser.hpp"
#include "SpecialEclipseFields.hpp"

namespace Dune
{

EclipseGridInspector::EclipseGridInspector(const EclipseGridParser& parser)
    : parser_(parser)
{
    std::vector<std::string> keywords;
    keywords.push_back("COORD");
    keywords.push_back("ZCORN");

    if (!parser_.hasFields(keywords)) {
	THROW("Needed field is missing in file");
    }

    if (parser_.hasField("SPECGRID")) {
        const SPECGRID& sgr = parser.getSPECGRID();
        logical_gridsize_[0] = sgr.dimensions[0];
        logical_gridsize_[1] = sgr.dimensions[1];
        logical_gridsize_[2] = sgr.dimensions[2];
    } else if (parser_.hasField("DIMENS")) {
        const std::vector<int>& dim = parser.getIntegerValue("DIMENS");
        logical_gridsize_[0] = dim[0];
        logical_gridsize_[1] = dim[1];
        logical_gridsize_[2] = dim[2];
    } else {
        THROW("Found neither SPECGRID nor DIMENS in file. At least one is needed.");
    }

}

double EclipseGridInspector::cellVolumeVerticalPillars(int i, int j, int k) const
{
    // Checking parameters and obtaining values from parser.
    checkLogicalCoords(i, j, k);
    const std::vector<double>& pillc = parser_.getFloatingPointValue("COORD");
    int num_pillars = (logical_gridsize_[0] + 1)*(logical_gridsize_[1] + 1);
    if (6*num_pillars != int(pillc.size())) {
	throw std::runtime_error("Wrong size of COORD field.");
    }
    const std::vector<double>& z = parser_.getFloatingPointValue("ZCORN");
    int num_cells = logical_gridsize_[0]*logical_gridsize_[1]*logical_gridsize_[2];
    if (8*num_cells != int(z.size())) {
	throw std::runtime_error("Wrong size of ZCORN field");
    }

    // Computing the base area as half the 2d cross product of the diagonals.
    int numxpill = logical_gridsize_[0] + 1;
    int pix = i + j*numxpill;
    double px[4] = { pillc[6*pix],
		     pillc[6*(pix + 1)],
		     pillc[6*(pix + numxpill)],
		     pillc[6*(pix + numxpill + 1)] };
    double py[4] = { pillc[6*pix + 1],
		     pillc[6*(pix + 1) + 1],
		     pillc[6*(pix + numxpill) + 1],
		     pillc[6*(pix + numxpill + 1) + 1] };
    double diag1[2] = { px[3] - px[0], py[3] - py[0] };
    double diag2[2] = { px[2] - px[1], py[2] - py[1] };
    double area = 0.5*(diag1[0]*diag2[1] - diag1[1]*diag2[0]);

    // Computing the average of the z-differences along each pillar.
    int delta[3] = { 1,
		     2*logical_gridsize_[0],
		     4*logical_gridsize_[0]*logical_gridsize_[1] };
    int ix = 2*(i*delta[0] + j*delta[1] + k*delta[2]);
    double cellz[8] = { z[ix], z[ix + delta[0]],
			z[ix + delta[1]], z[ix + delta[1] + delta[0]],
			z[ix + delta[2]], z[ix + delta[2] + delta[0]],
			z[ix + delta[2] + delta[1]], z[ix + delta[2] + delta[1] + delta[0]] };
    double diffz[4] = { cellz[4] - cellz[0],
			cellz[5] - cellz[1],
			cellz[6] - cellz[2],
			cellz[7] - cellz[3] };
    double averzdiff = 0.25*std::accumulate(diffz, diffz + 4, 0.0);
    return averzdiff*area;
}


double EclipseGridInspector::cellVolumeVerticalPillars(int cell_idx) const
{
    int i, j, k;
    int horIdx = (cell_idx+1) -
        int(std::floor(((double)(cell_idx+1))/
                       ((double)(logical_gridsize_[0] * logical_gridsize_[1])))) *
        logical_gridsize_[0]*logical_gridsize_[1]; // index in the corresponding horizon
    if (horIdx == 0) {
        horIdx = logical_gridsize_[0] * logical_gridsize_[1];
    }
    i = horIdx - int(std::floor(((double)horIdx)/((double)logical_gridsize_[0]))) * logical_gridsize_[0];
    if (i == 0) {
        i = logical_gridsize_[1];
    }
    j = (horIdx-i)/logical_gridsize_[0] + 1;
    k = ((cell_idx+1)-logical_gridsize_[0]*(j-1)-1)/(logical_gridsize_[0]*logical_gridsize_[1]) + 1;
    return cellVolumeVerticalPillars(i-1, j-1, k-1);
}

void EclipseGridInspector::checkLogicalCoords(int i, int j, int k) const
{
    if (i < 0 || i >= logical_gridsize_[0])
	throw std::runtime_error("First coordinate out of bounds");
    if (j < 0 || j >= logical_gridsize_[1])
	throw std::runtime_error("Second coordinate out of bounds");
    if (k < 0 || k >= logical_gridsize_[2])
	throw std::runtime_error("Third coordinate out of bounds");
}


boost::array<double, 6> EclipseGridInspector::getGridLimits() const
{
    if (! (parser_.hasField("COORD") && parser_.hasField("ZCORN") && parser_.hasField("SPECGRID")) ) {
        throw std::runtime_error("EclipseGridInspector: Grid does not have SPECGRID, COORD, and ZCORN, can't find dimensions.");
    }

    std::vector<int> specgrid = parser_.getIntegerValue("SPECGRID");
    std::vector<double> coord = parser_.getFloatingPointValue("COORD");
    std::vector<double> zcorn = parser_.getFloatingPointValue("ZCORN");

    double xmin = +DBL_MAX;
    double xmax = -DBL_MAX;
    double ymin = +DBL_MAX;
    double ymax = -DBL_MAX;


    int pillars = (specgrid[0]+1) * (specgrid[1]+1);

    for (int pillarindex = 0; pillarindex < pillars; ++pillarindex) {
        if        (coord[pillarindex * 6 + 0] > xmax)
            xmax = coord[pillarindex * 6 + 0];
        if        (coord[pillarindex * 6 + 0] < xmin)
            xmin = coord[pillarindex * 6 + 0];
        if        (coord[pillarindex * 6 + 1] > ymax)
            ymax = coord[pillarindex * 6 + 1];
        if        (coord[pillarindex * 6 + 1] < ymin)
            ymin = coord[pillarindex * 6 + 1];
        if        (coord[pillarindex * 6 + 3] > xmax)
            xmax = coord[pillarindex * 6 + 3];
        if        (coord[pillarindex * 6 + 3] < xmin)
            xmin = coord[pillarindex * 6 + 3];
        if        (coord[pillarindex * 6 + 4] > ymax)
            ymax = coord[pillarindex * 6 + 4];
        if        (coord[pillarindex * 6 + 4] < ymin)
            ymin = coord[pillarindex * 6 + 4];
    }

    boost::array<double, 6> gridlimits = {{ xmin, xmax, ymin, ymax,
                                            *min_element(zcorn.begin(), zcorn.end()),
                                            *max_element(zcorn.begin(), zcorn.end()) }};
    return gridlimits;
}



boost::array<int, 3> EclipseGridInspector::gridSize() const
{
    boost::array<int, 3> retval = {{ logical_gridsize_[0],
				     logical_gridsize_[1],
				     logical_gridsize_[2] }};
    return retval;
}


boost::array<double, 8> EclipseGridInspector::cellZvals(int i, int j, int k) const
{
    // Get the zcorn field.
    const std::vector<double>& z = parser_.getFloatingPointValue("ZCORN");
    int num_cells = logical_gridsize_[0]*logical_gridsize_[1]*logical_gridsize_[2];
    if (8*num_cells != int(z.size())) {
	throw std::runtime_error("Wrong size of ZCORN field");
    }

    // Make the coordinate array.
    int delta[3] = { 1,
		     2*logical_gridsize_[0],
		     4*logical_gridsize_[0]*logical_gridsize_[1] };
    int ix = 2*(i*delta[0] + j*delta[1] + k*delta[2]);
    boost::array<double, 8> cellz = {{ z[ix], z[ix + delta[0]],
				       z[ix + delta[1]], z[ix + delta[1] + delta[0]],
				       z[ix + delta[2]], z[ix + delta[2] + delta[0]],
				       z[ix + delta[2] + delta[1]], z[ix + delta[2] + delta[1] + delta[0]] }};
    return cellz;
}


} // namespace Dune
