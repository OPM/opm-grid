/*
  Copyright 2026 Equinor ASA.

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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <opm/grid/cpgrid/refinement/conforming/GrdeclRefinement.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/utility/numeric/CornerPointRefinement.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm
{
namespace Refinement
{

std::array<int,3> RefinedBlockGrdecl::intDims() const
{
    return { static_cast<int>(dims[0]), static_cast<int>(dims[1]), static_cast<int>(dims[2]) };
}

RefinedBlockGrdecl refineBlock(const std::array<int,3>& parentDims,
                               const double* coord,
                               const double* zcorn,
                               const int* actnum,
                               const BlockRefinement& request)
{
    for (int dim = 0; dim < 3; ++dim) {
        if (request.startIJK[dim] < 0 || request.endIJK[dim] > parentDims[dim]
            || request.startIJK[dim] >= request.endIJK[dim]) {
            OPM_THROW(std::invalid_argument, "Refinement box '" + request.name
                                        + "' does not fit inside the parent dimensions.");
        }
        if (request.subdivision[dim].empty() && request.cellsPerDim[dim] < 1) {
            OPM_THROW(std::invalid_argument, "Refinement '" + request.name
                                        + "' has non-positive subdivisions.");
        }
    }
    // The parent dimensions as sizes, for the array lengths and indices below.
    const std::size_t nx = parentDims[0];
    const std::size_t ny = parentDims[1];
    const std::size_t nz = parentDims[2];

    const std::array<Opm::Refinement::AxisSubdivision,3> subs = { axisSubdivision(request, 0),
                                                axisSubdivision(request, 1),
                                                axisSubdivision(request, 2) };

    // The number of refined cells in each direction is the number of refined columns.
    RefinedBlockGrdecl out;
    out.dims = { subs[0].size(), subs[1].size(), subs[2].size() };

    // --- COORD and ZCORN: by the functions the LGR grids opm-common builds use,
    // so every block that shares a pillar or a face computes it alike.
    Opm::CornerPointRefinement::Box box { parentDims, request.startIJK, {} };
    for (int dim = 0; dim < 3; ++dim) {
        box.axes[dim] = { subs[dim].parentOffset, subs[dim].fracLo, subs[dim].fracHi };
    }
    out.coord = Opm::CornerPointRefinement::refinedCoord(box, { coord, 6 * (nx + 1) * (ny + 1) });
    out.zcorn = Opm::CornerPointRefinement::refinedZcorn(box, { zcorn, 8 * nx * ny * nz });

    // --- ACTNUM: children inherit the parent flag -------------------------
    // Child (refined) cell (ci, cj, ck) lies in parent cell (pi, pj, pk).
    out.actnum.assign(out.dims[0] * out.dims[1] * out.dims[2], 1);
    if (actnum) {
        for (std::size_t ck = 0; ck < out.dims[2]; ++ck) {
            const std::size_t pk = request.startIJK[2] + subs[2].parentOffset[ck];
            for (std::size_t cj = 0; cj < out.dims[1]; ++cj) {
                const std::size_t pj = request.startIJK[1] + subs[1].parentOffset[cj];
                for (std::size_t ci = 0; ci < out.dims[0]; ++ci) {
                    const std::size_t pi = request.startIJK[0] + subs[0].parentOffset[ci];
                    const std::size_t parentIdx = pi + nx*pj + nx*ny*pk;
                    const std::size_t childIdx = ci + out.dims[0]*cj + out.dims[0]*out.dims[1]*ck;
                    out.actnum[childIdx] = actnum[parentIdx];
                }
            }
        }
    }

    // ... and the block's own MINPV takes out what it takes out.
    if (! request.minpvRemoved.empty()) {
        if (request.minpvRemoved.size() != out.actnum.size()) {
            OPM_THROW(std::invalid_argument, "Refinement '" + request.name
                                        + "' has a MINPV removal mask of the wrong size.");
        }
        for (std::size_t cell = 0; cell < out.actnum.size(); ++cell) {
            if (request.minpvRemoved[cell]) {
                out.actnum[cell] = 0;
            }
        }
    }

    return out;
}

} // namespace Refinement
} // namespace Opm
