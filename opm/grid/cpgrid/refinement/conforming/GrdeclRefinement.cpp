/*
  Copyright 2026 SINTEF Digital.

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

RefinedBlockGrdecl refineBlock(const std::array<int,3>& parentDims,
                               const double* coord,
                               const double* zcorn,
                               const int* actnum,
                               const BlockRefinement& request)
{
    const auto& [nx, ny, nz] = parentDims;
    for (int c = 0; c < 3; ++c) {
        if (request.startIJK[c] < 0 || request.endIJK[c] > parentDims[c]
            || request.startIJK[c] >= request.endIJK[c]) {
            OPM_THROW(std::invalid_argument, "Refinement box '" + request.name
                                        + "' does not fit inside the parent dimensions.");
        }
        if (request.subdivision[c].empty() && request.cellsPerDim[c] < 1) {
            OPM_THROW(std::invalid_argument, "Refinement '" + request.name
                                        + "' has non-positive subdivisions.");
        }
    }

    const std::array<Opm::Refinement::AxisSubdivision,3> subs = { axisSubdivision(request, 0),
                                                axisSubdivision(request, 1),
                                                axisSubdivision(request, 2) };

    RefinedBlockGrdecl out;
    out.dims = refinedDims(request);

    // --- COORD and ZCORN: by the functions the LGR grids opm-common builds use,
    // so every block that shares a pillar or a face computes it alike.
    Opm::CornerPointRefinement::Box box { parentDims, request.startIJK, {} };
    for (int c = 0; c < 3; ++c) {
        box.axes[c] = { subs[c].parentOffset, subs[c].fracLo, subs[c].fracHi };
    }
    out.coord = Opm::CornerPointRefinement::refinedCoord(
        box, { coord, 6 * static_cast<std::size_t>(nx + 1) * (ny + 1) });
    out.zcorn = Opm::CornerPointRefinement::refinedZcorn(
        box, { zcorn, 8 * static_cast<std::size_t>(nx) * ny * nz });

    // --- ACTNUM: children inherit the parent flag -------------------------
    out.actnum.assign(static_cast<std::size_t>(out.dims[0]) * out.dims[1] * out.dims[2], 1);
    if (actnum) {
        for (int kr = 0; kr < out.dims[2]; ++kr) {
            const int ck = request.startIJK[2] + subs[2].parentOffset[kr];
            for (int jr = 0; jr < out.dims[1]; ++jr) {
                const int cj = request.startIJK[1] + subs[1].parentOffset[jr];
                for (int ir = 0; ir < out.dims[0]; ++ir) {
                    const int ci = request.startIJK[0] + subs[0].parentOffset[ir];
                    const std::size_t parentIdx = static_cast<std::size_t>(ci)
                        + static_cast<std::size_t>(nx)*cj
                        + static_cast<std::size_t>(nx)*ny*ck;
                    const std::size_t childIdx = static_cast<std::size_t>(ir)
                        + static_cast<std::size_t>(out.dims[0])*jr
                        + static_cast<std::size_t>(out.dims[0])*out.dims[1]*kr;
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
