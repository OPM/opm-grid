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

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace
{

// Lateral position of a refined pillar/corner line within the parent block:
// which parent column it belongs to and the fractional coordinate in it.
// The last refined line is assigned to the last column with fraction 1 so
// that every (column, fraction) pair is well defined.
struct LateralPos
{
    int cell;       // parent cell index in this direction (absolute)
    double frac;    // in [0, 1]
};

// A refined line is the low edge of its own column, except the last, which
// closes the previous one.
LateralPos lateralPos(int refinedIdx, const Opm::Refinement::AxisSubdivision& sub, int start)
{
    const int last = static_cast<int>(sub.size()) - 1;
    if (refinedIdx > last) {
        return { start + sub.parentOffset[last], sub.fracHi[last] };
    }
    return { start + sub.parentOffset[refinedIdx], sub.fracLo[refinedIdx] };
}

} // anonymous namespace

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
            throw std::invalid_argument("Refinement box '" + request.name
                                        + "' does not fit inside the parent dimensions.");
        }
        if (request.subdivision[c].empty() && request.cellsPerDim[c] < 1) {
            throw std::invalid_argument("Refinement '" + request.name
                                        + "' has non-positive subdivisions.");
        }
    }

    const std::array<Opm::Refinement::AxisSubdivision,3> subs = { axisSubdivision(request, 0),
                                                axisSubdivision(request, 1),
                                                axisSubdivision(request, 2) };

    RefinedBlockGrdecl out;
    out.dims = refinedDims(request);

    // --- COORD: sub-pillars ---------------------------------------------
    // Inside a parent column the line passes through the bilinear interpolation of the
    // parent's corners on the top and bottom of the column's first non-collapsed layer, so
    // boxes stacked in a column share pillars (pillarsFromBoxLayer: the box's own first
    // layer, as the reference). Parent pillars are copied; an all-collapsed column falls
    // back to endpoint-wise interpolation of the four parent pillars.
    const auto parentPillar = [&](int i, int j) {
        return coord + 6*(static_cast<std::size_t>(j)*(nx + 1) + i);
    };
    const auto parentCorner = [&](int i, int j, int k, int di, int dj, int dk) {
        const double* pil = parentPillar(i + di, j + dj);
        const double z = zcorn[static_cast<std::size_t>(2*i + di)
                               + 2*static_cast<std::size_t>(nx)*(2*j + dj)
                               + 4*static_cast<std::size_t>(nx)*ny*(2*k + dk)];
        const double t = (pil[5] != pil[2]) ? (z - pil[2]) / (pil[5] - pil[2]) : 0.0;
        return std::array<double,3>{ pil[0] + t*(pil[3] - pil[0]), pil[1] + t*(pil[4] - pil[1]), z };
    };

    // A line on a parent boundary inside the box belongs to the lower column, as in the
    // reference; across a fault the two columns' corners differ.
    const auto pillarPos = [&](int refinedIdx, int dim) {
        auto pos = lateralPos(refinedIdx, subs[dim], request.startIJK[dim]);
        if ((pos.frac == 0.0) && (pos.cell > request.startIJK[dim])) {
            pos = { pos.cell - 1, 1.0 };
        }
        return pos;
    };

    out.coord.resize(6 * static_cast<std::size_t>(out.dims[0] + 1) * (out.dims[1] + 1));
    for (int jr = 0; jr <= out.dims[1]; ++jr) {
        const auto [cj, b] = pillarPos(jr, 1);
        for (int ir = 0; ir <= out.dims[0]; ++ir) {
            const auto [ci, a] = pillarPos(ir, 0);

            const double* p00 = parentPillar(ci,     cj);
            const double* p10 = parentPillar(ci + 1, cj);
            const double* p01 = parentPillar(ci,     cj + 1);
            const double* p11 = parentPillar(ci + 1, cj + 1);

            double* sub = &out.coord[6*(static_cast<std::size_t>(jr)*(out.dims[0] + 1) + ir)];
            for (int comp = 0; comp < 6; ++comp) {
                sub[comp] = (1.0 - a)*(1.0 - b)*p00[comp] + a*(1.0 - b)*p10[comp]
                          + (1.0 - a)*b*p01[comp] + a*b*p11[comp];
            }
            if ((a == 0.0 || a == 1.0) && (b == 0.0 || b == 1.0)) {
                continue;
            }
            const int kFirst = request.pillarsFromBoxLayer ? request.startIJK[2] : 0;
            const int kEnd = request.pillarsFromBoxLayer ? request.endIJK[2] : nz;
            for (int ck = kFirst; ck < kEnd; ++ck) {
                std::array<std::array<double,3>,2> ends{};
                for (int dk = 0; dk < 2; ++dk) {
                    for (int dj = 0; dj < 2; ++dj) {
                        for (int di = 0; di < 2; ++di) {
                            const double w = (di ? a : 1.0 - a) * (dj ? b : 1.0 - b);
                            const auto x = parentCorner(ci, cj, ck, di, dj, dk);
                            for (int c = 0; c < 3; ++c) {
                                ends[dk][c] += w*x[c];
                            }
                        }
                    }
                }
                const double dz = ends[1][2] - ends[0][2];
                if (std::abs(dz) > 1.0e-6) {
                    // Keep the endpoint depths, so only the line changes.
                    for (double* end : {sub, sub + 3}) {
                        const double t = (end[2] - ends[0][2]) / dz;
                        end[0] = ends[0][0] + t*(ends[1][0] - ends[0][0]);
                        end[1] = ends[0][1] + t*(ends[1][1] - ends[0][1]);
                    }
                    break;
                }
            }
        }
    }

    // --- ZCORN: trilinear resampling per parent cell ----------------------
    const auto parentZ = [&](int i_, int j_, int k_) {
        return zcorn[static_cast<std::size_t>(i_)
                     + 2*static_cast<std::size_t>(nx)*j_
                     + 4*static_cast<std::size_t>(nx)*ny*k_];
    };
    const auto refinedZIndex = [&](int i_, int j_, int k_) {
        return static_cast<std::size_t>(i_)
               + 2*static_cast<std::size_t>(out.dims[0])*j_
               + 4*static_cast<std::size_t>(out.dims[0])*out.dims[1]*k_;
    };

    out.zcorn.resize(8 * static_cast<std::size_t>(out.dims[0]) * out.dims[1] * out.dims[2]);
    for (int kr = 0; kr < out.dims[2]; ++kr) {
        const int ck = request.startIJK[2] + subs[2].parentOffset[kr];
        for (int jr = 0; jr < out.dims[1]; ++jr) {
            const int cj = request.startIJK[1] + subs[1].parentOffset[jr];
            for (int ir = 0; ir < out.dims[0]; ++ir) {
                const int ci = request.startIJK[0] + subs[0].parentOffset[ir];

                for (int dk = 0; dk < 2; ++dk) {
                    const double c = dk ? subs[2].fracHi[kr] : subs[2].fracLo[kr];
                    for (int dj = 0; dj < 2; ++dj) {
                        const double b = dj ? subs[1].fracHi[jr] : subs[1].fracLo[jr];
                        for (int di = 0; di < 2; ++di) {
                            const double a = di ? subs[0].fracHi[ir] : subs[0].fracLo[ir];

                            // Where the parent is collapsed (as MINPV leaves it), so are
                            // the children; rounding would otherwise invert some.
                            double parentDz = 0.0;
                            for (int pj = 0; pj < 2; ++pj) {
                                const double wj = (pj == 0) ? (1.0 - b) : b;
                                for (int pi = 0; pi < 2; ++pi) {
                                    const double wi = (pi == 0) ? (1.0 - a) : a;
                                    parentDz += wi*wj*(parentZ(2*ci + pi, 2*cj + pj, 2*ck + 1)
                                                       - parentZ(2*ci + pi, 2*cj + pj, 2*ck));
                                }
                            }
                            const double cc = (parentDz == 0.0) ? 0.0 : c;

                            double z = 0.0;
                            for (int pk = 0; pk < 2; ++pk) {
                                const double wk = (pk == 0) ? (1.0 - cc) : cc;
                                for (int pj = 0; pj < 2; ++pj) {
                                    const double wj = (pj == 0) ? (1.0 - b) : b;
                                    for (int pi = 0; pi < 2; ++pi) {
                                        const double wi = (pi == 0) ? (1.0 - a) : a;
                                        z += wi*wj*wk*parentZ(2*ci + pi, 2*cj + pj, 2*ck + pk);
                                    }
                                }
                            }
                            out.zcorn[refinedZIndex(2*ir + di, 2*jr + dj, 2*kr + dk)] = z;
                        }
                    }
                }
            }
        }
    }

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
            throw std::invalid_argument("Refinement '" + request.name
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
