/*
  Copyright 2026 Equinor ASA

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
#ifndef OPM_TESTS_CORNERPOINTTESTGRID_HEADER_INCLUDED
#define OPM_TESTS_CORNERPOINTTESTGRID_HEADER_INCLUDED

#include <array>
#include <cstddef>
#include <vector>

namespace Opm::TestGrids
{

// A corner-point grid with vertical pillars at the integer positions x = 0..nx, y = 0..ny.
// The depth of the corner (ix, jy, zs) of the cell (i, j, k), where ix and jy tell which of
// the pillars of the cell and zs whether it is the top or the bottom, is given by depth().
struct CornerPointGrid
{
    std::array<int, 3> dims;
    std::vector<double> coord;
    std::vector<double> zcorn;
    // For each cell its eight corners, in the order of the Dune reference cube: ix + 2 * jy + 4 * zs.
    std::vector<std::array<double, 3>> expectedCorners;

    template <class Depth>
    CornerPointGrid(const int nx, const int ny, const int nz, Depth depth)
        : dims{nx, ny, nz}
    {
        for (int py = 0; py <= ny; ++py) {
            for (int px = 0; px <= nx; ++px) {
                for (const double z : {-10.0, 10.0}) {
                    coord.insert(coord.end(), {double(px), double(py), z});
                }
            }
        }

        // The order of ZCORN: layers from top to bottom, within each the rows of corners in y,
        // and within each row the corners in x. Each cell has two corners in each direction.
        zcorn.resize(8 * nx * ny * nz);
        expectedCorners.resize(8 * nx * ny * nz);
        for (int k = 0; k < nz; ++k) {
            for (int zs = 0; zs < 2; ++zs) {
                for (int j = 0; j < ny; ++j) {
                    for (int jy = 0; jy < 2; ++jy) {
                        for (int i = 0; i < nx; ++i) {
                            for (int ix = 0; ix < 2; ++ix) {
                                const double z = depth(i, j, k, zs, ix, jy);
                                const std::size_t index = ((2 * k + zs) * 2 * ny + 2 * j + jy) * 2 * nx + 2 * i + ix;
                                zcorn[index] = z;
                                const std::size_t cell = i + nx * (j + ny * k);
                                expectedCorners[8 * cell + ix + 2 * jy + 4 * zs] = {double(i + ix), double(j + jy), z};
                            }
                        }
                    }
                }
            }
        }
    }
};

// Layers of unit thickness, with the corners moved down by shift(i, j, k, ix, jy). The top of the
// layer k is at the depth k and the bottom at k + 1, before the shift. Without a shift the grid
// is regular, and seen along y a grid with two columns and two layers looks like
//
//        x=0       x=1       x=2
//  z=0   +---------+---------+
//        |         |         |
//  z=1   +---------+---------+
//        |         |         |
//  z=2   +---------+---------+
//
// A shift that depends on i moves a whole column, which makes a fault between columns. A shift that
// depends on ix or jy as well tilts the layers of the cell.
template <class Shift>
auto layersWithShift(Shift shift)
{
    return [shift](int i, int j, int k, int zs, int ix, int jy) { return k + zs + shift(i, j, k, ix, jy); };
}


} // namespace Opm::TestGrids

#endif // OPM_TESTS_CORNERPOINTTESTGRID_HEADER_INCLUDED
