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

#ifndef OPM_CELLQUADRATURE_HEADER_INCLUDED
#define OPM_CELLQUADRATURE_HEADER_INCLUDED

struct UnstructuredGrid;

namespace Opm
{
    /// A class providing numerical quadrature for cells.
    /// In general: \int_{cell} g(x) dx = \sum_{i=0}^{n-1} w_i g(x_i).
    /// Note that this class does multiply weights by cell volume,
    /// so weights always sum to cell volume.
    ///
    /// Degree 1 method:
    ///     Midpoint (centroid) method.
    ///         n = 1, w_0 = cell volume, x_0 = cell centroid
    ///
    /// Degree 2 method for 2d (but see the note):
    ///    Based on subdivision of the cell into triangles,
    ///    with the centroid as a common vertex, and the triangle
    ///    edge midpoint rule.
    ///    Triangle i consists of the centroid C, nodes N_i and N_{i+1}.
    ///    Its area is A_i.
    ///        n = 2 * nn  (nn = num nodes in face)
    ///        For i = 0..(nn-1):
    ///        w_i      = 1/3 A_i.
    ///        w_{nn+i} = 1/3 A_{i-1} + 1/3 A_i
    ///        x_i      = (N_i + N_{i+1})/2
    ///        x_{nn+i} = (C + N_i)/2
    ///    All N and A indices are interpreted cyclic, modulus nn.
    ///    Note: for simplicity of implementation, we currently use
    ///        n = 3 * nn
    ///        For i = 0..(nn-1):
    ///        w_{3*i + {0,1,2}} = 1/3 A_i
    ///        x_{3*i}           = (N_i + N_{i+1})/2
    ///        x_{3*i + {1,2}}   = (C + N_{i,i+1})/2
    ///    This is simpler, because we can implement it easily
    ///    based on iteration over faces without requiring any
    ///    particular (cyclic) ordering.
    ///
    /// Degree 2 method for 3d:
    ///    Based on subdivision of each cell face into triangles
    ///    with the face centroid as a common vertex, and then
    ///    subdividing the cell into tetrahedra with the cell
    ///    centroid as a common vertex. Then apply the tetrahedron
    ///    rule with the following 4 nodes (uniform weights):
    ///        a = 0.138196601125010515179541316563436
    ///        x_i has all barycentric coordinates = a, except for
    ///            the i'th coordinate which is = 1 - 3a.
    ///    This rule is from http://nines.cs.kuleuven.be/ecf,
    ///    it is the second degree 2 4-point rule for tets,
    ///    referenced to Stroud(1971).
    ///    The tetrahedra are numbered T_{i,j}, and are given by the
    ///    cell centroid C, the face centroid FC_i, and two nodes
    ///    of face i: FN_{i,j}, FN_{i,j+1}.
    class CellQuadrature
    {
    public:
        CellQuadrature(const UnstructuredGrid& grid,
                       const int cell,
                       const int degree);

        int numQuadPts() const;
        void quadPtCoord(const int index, double* coord) const;
        double quadPtWeight(const int index) const;

    private:
        const UnstructuredGrid& grid_;
        const int cell_;
        const int degree_;
    };

} // namespace Opm

#endif // OPM_CELLQUADRATURE_HEADER_INCLUDED
