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

#include <opm/grid/cpgrid/refinement/conforming/EdgeConformal.hpp>

#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/Geometry.hpp>
#include <opm/grid/cpgrid/refinement/GridStateWriter.hpp>
#include <opm/common/utility/SparseTable.hpp>

#include <dune/common/fvector.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <tuple>
#include <vector>

namespace Opm
{
namespace Refinement
{

void edgeConformalizeLeaf(Dune::cpgrid::CpGridData& leaf)
{
    using Coord = Dune::FieldVector<double,3>;

    const auto& cornerGeom =
        GridStateWriter::geometry(leaf).geomVector(std::integral_constant<int,3>());
    const int numCorners = static_cast<int>(cornerGeom.size());
    if (numCorners == 0) {
        return;
    }

    std::vector<Coord> P(numCorners);
    Coord lo = cornerGeom.get(0).center();
    Coord hi = lo;
    for (int i = 0; i < numCorners; ++i) {
        P[i] = cornerGeom.get(i).center();
        for (int c = 0; c < 3; ++c) {
            lo[c] = std::min(lo[c], P[i][c]);
            hi[c] = std::max(hi[c], P[i][c]);
        }
    }
    double diag = 0.0;
    for (int c = 0; c < 3; ++c) {
        diag += (hi[c] - lo[c]) * (hi[c] - lo[c]);
    }
    diag = std::sqrt(diag);
    if (diag <= 0.0) {
        return;
    }

    // Spatial hash so each edge only tests nearby nodes.
    const double h = diag / 256.0 + 1e-30;
    const double tol = 1e-9 * diag;
    using Bucket = std::tuple<long long, long long, long long>;
    const auto bucketOf = [&](const Coord& p) {
        return Bucket{ std::llround(std::floor((p[0] - lo[0]) / h)),
                       std::llround(std::floor((p[1] - lo[1]) / h)),
                       std::llround(std::floor((p[2] - lo[2]) / h)) };
    };
    std::map<Bucket, std::vector<int>> hash;
    for (int i = 0; i < numCorners; ++i) {
        hash[bucketOf(P[i])].push_back(i);
    }

    auto& faceToPoint = GridStateWriter::faceToPoint(leaf);
    Opm::SparseTable<int> rebuilt;
    std::vector<int> newRow;
    std::vector<std::pair<double,int>> mids;

    for (int face = 0; face < faceToPoint.size(); ++face) {
        auto row = faceToPoint[face];
        const int n = row.size();
        newRow.clear();
        for (int e = 0; e < n; ++e) {
            const int a = row[e];
            const int b = row[(e + 1) % n];
            newRow.push_back(a);

            Coord d = P[b];
            d -= P[a];
            const double len2 = d.two_norm2();
            if (len2 <= tol * tol) {
                continue;
            }

            // Candidate nodes: all buckets overlapping the edge's bounding box.
            const long long x0 = std::llround(std::floor((std::min(P[a][0], P[b][0]) - lo[0] - tol) / h));
            const long long x1 = std::llround(std::floor((std::max(P[a][0], P[b][0]) - lo[0] + tol) / h));
            const long long y0 = std::llround(std::floor((std::min(P[a][1], P[b][1]) - lo[1] - tol) / h));
            const long long y1 = std::llround(std::floor((std::max(P[a][1], P[b][1]) - lo[1] + tol) / h));
            const long long z0 = std::llround(std::floor((std::min(P[a][2], P[b][2]) - lo[2] - tol) / h));
            const long long z1 = std::llround(std::floor((std::max(P[a][2], P[b][2]) - lo[2] + tol) / h));

            mids.clear();
            for (long long bx = x0; bx <= x1; ++bx) {
                for (long long by = y0; by <= y1; ++by) {
                    for (long long bz = z0; bz <= z1; ++bz) {
                        const auto it = hash.find({bx, by, bz});
                        if (it == hash.end()) {
                            continue;
                        }
                        for (const int c : it->second) {
                            if (c == a || c == b) {
                                continue;
                            }
                            Coord ac = P[c];
                            ac -= P[a];
                            const double t = (ac * d) / len2;
                            if (t <= 0.0 || t >= 1.0) {
                                continue; // not strictly between the endpoints
                            }
                            // Distance from the node to the line a-b.
                            Coord closest = d;
                            closest *= t;
                            closest += P[a];
                            Coord diff = P[c];
                            diff -= closest;
                            if (diff.two_norm() < tol) {
                                mids.emplace_back(t, c);
                            }
                        }
                    }
                }
            }
            std::sort(mids.begin(), mids.end());
            for (const auto& [t, c] : mids) {
                newRow.push_back(c);
            }
        }
        rebuilt.appendRow(newRow.begin(), newRow.end());
    }

    faceToPoint.swap(rebuilt);
}

} // namespace Refinement
} // namespace Opm
