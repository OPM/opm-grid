/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <config.h>
#include <opm/grid/FaceQuadrature.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/grid/UnstructuredGrid.h>
#include <cmath>

namespace {

/// Calculates the cross product of two 3-vectors, represented as
/// 3-element arrays. Calculates res = a X b. The res array must
/// already be allocated with room for 3 elements.
inline void cross(const double* a, const double* b, double* res)
{
    res[0] = a[1]*b[2] - a[2]*b[1];
    res[1] = a[2]*b[0] - a[0]*b[2];
    res[2] = a[0]*b[1] - a[1]*b[0];
}

/// Calculates the area of a triangle consisting of 3 vertices
/// with 3-dimensional coordinates
inline double triangleArea3d(const double* p0,
                             const double* p1,
                             const double* p2)
{
    double a[3] = { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };
    double b[3] = { p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2] };
    double cr[3];
    cross(a, b, cr);
    return 0.5*std::sqrt(cr[0]*cr[0] + cr[1]*cr[1] + cr[2]*cr[2]);
}

} // anonymous namespace

namespace Opm {

FaceQuadrature::FaceQuadrature(const UnstructuredGrid& grid,
                               const int face,
                               const int degree)
    : grid_(grid), face_(face), degree_(degree)
{
    if (grid_.dimensions > 3) {
        OPM_THROW(std::runtime_error, "FaceQuadrature only implemented for up to 3 dimensions.");
    }
    if (degree_ > 2) {
        OPM_THROW(std::runtime_error, "FaceQuadrature exact for polynomial degrees > 2 not implemented.");
    }
}

int FaceQuadrature::numQuadPts() const
{
    if (degree_ < 2 || grid_.dimensions < 2) {
        return 1;
    }
    // Degree 2 case.
    if (grid_.dimensions == 2) {
        return 3;
    } else {
        return 2 * (grid_.face_nodepos[face_ + 1] - grid_.face_nodepos[face_]);
    }
}

void FaceQuadrature::quadPtCoord(const int index, double* coord) const
{
    const int dim = grid_.dimensions;
    const double* fc = grid_.face_centroids + dim*face_;
    if (degree_ < 2 || dim < 2) {
        std::copy(fc, fc + dim, coord);
        return;
    }
    // Degree 2 case.
    const int nn = grid_.face_nodepos[face_ + 1] - grid_.face_nodepos[face_];
    const int* fnodes = grid_.face_nodes + grid_.face_nodepos[face_];
    const double* nc = grid_.node_coordinates;
    if (dim == 2) {
        assert(nn == 2);
        const double* pa[3] = { nc + dim*fnodes[0], fc, nc + dim*fnodes[1] };
        std::copy(pa[index], pa[index] + dim, coord);
        return;
    }
    assert(dim == 3);
    if (index < nn) {
        // Boundary edge midpoint.
        const int node0 = fnodes[index];
        const int node1 = fnodes[(index + 1)%nn];
        for (int dd = 0; dd < dim; ++dd) {
            coord[dd] = 0.5*(nc[dim*node0 + dd] + nc[dim*node1 + dd]);
        }
    } else {
        // Interiour edge midpoint.
        // Recall that index is now in [nn, 2*nn).
        const int node = fnodes[index - nn];
        for (int dd = 0; dd < dim; ++dd) {
            coord[dd] = 0.5*(nc[dim*node + dd] + fc[dd]);
        }
    }
}

double FaceQuadrature::quadPtWeight(const int index) const
{
    if (degree_ < 2) {
        return grid_.face_areas[face_];
    }
    // Degree 2 case.
    const int dim = grid_.dimensions;
    if (dim == 2) {
        const double simpsonw[3] = { 1.0/6.0, 4.0/6.0, 1.0/6.0 };
        return grid_.face_areas[face_]*simpsonw[index];
    }
    assert(dim == 3);
    const double* fc = grid_.face_centroids + dim*face_;
    const int nn = grid_.face_nodepos[face_ + 1] - grid_.face_nodepos[face_];
    const int* fnodes = grid_.face_nodes + grid_.face_nodepos[face_];
    const double* nc = grid_.node_coordinates;
    if (index < nn) {
        // Boundary edge midpoint.
        const int node0 = fnodes[index];
        const int node1 = fnodes[(index + 1)%nn];
        const double area = triangleArea3d(nc + dim*node1, nc + dim*node0, fc);
        return area / 3.0;
    } else {
        // Interiour edge midpoint.
        // Recall that index is now in [nn, 2*nn).
        const int node0 = fnodes[(index - 1) % nn];
        const int node1 = fnodes[index - nn];
        const int node2 = fnodes[(index + 1) % nn];
        const double area0 = triangleArea3d(nc + dim*node1, nc + dim*node0, fc);
        const double area1 = triangleArea3d(nc + dim*node2, nc + dim*node1, fc);
        return (area0 + area1) / 3.0;
    }
}

} // namespace Opm
