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

#ifndef OPM_FACEQUADRATURE_HEADER_INCLUDED
#define OPM_FACEQUADRATURE_HEADER_INCLUDED

struct UnstructuredGrid;

namespace Opm
{

    /// A class providing numerical quadrature for faces.
    /// In general: \int_{face} g(x) dx = \sum_{i=0}^{n-1} w_i g(x_i).
    /// Note that this class does multiply weights by face area,
    /// so weights always sum to face area.
    ///
    /// Degree 1 method:
    ///     Midpoint (centroid) method.
    ///         n = 1, w_0 = face area, x_0 = face centroid
    ///
    /// Degree 2 method for 2d:
    ///    Simpson's method (actually this is degree 3).
    ///
    /// Degree 2 method for 3d:
    ///    Based on subdivision of the face into triangles,
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
    class FaceQuadrature
    {
    public:
        FaceQuadrature(const UnstructuredGrid& grid,
                       const int face,
                       const int degree);

        int numQuadPts() const;
        void quadPtCoord(const int index, double* coord) const;
        double quadPtWeight(const int index) const;

    private:
        const UnstructuredGrid& grid_;
        const int face_;
        const int degree_;
    };

} // namespace Opm

#endif // OPM_FACEQUADRATURE_HEADER_INCLUDED
