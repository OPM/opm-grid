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

#ifndef OPM_CONFORMINGTESTHELPERS_HEADER_INCLUDED
#define OPM_CONFORMINGTESTHELPERS_HEADER_INCLUDED

#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/refinement/conforming/ConformingBlockBuilder.hpp>
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>
#include <opm/grid/cpgrid/refinement/RefinementRequest.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <dune/common/parallel/mpihelper.hh>

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// Helpers shared by the tests of the conforming refinement: small corner-point grids
// built from a pillar and a depth function, refinement requests, and checks of the
// refined leaf.
//
// A test file includes this header, says `using namespace Opm::ConformingTest;` and
// keeps its own `BOOST_GLOBAL_FIXTURE(MPIFixture);`.

namespace Opm::ConformingTest
{

struct MPIFixture
{
    MPIFixture()
    {
        auto& argv = boost::unit_test::framework::master_test_suite().argv;
        auto& argc = boost::unit_test::framework::master_test_suite().argc;
        Dune::MPIHelper::instance(argc, argv);
    }
};

/// A corner-point grid description that owns its arrays. An empty actnum means
/// every cell is active.
struct TestGrdecl
{
    std::array<int,3> dims;
    std::vector<double> coord;
    std::vector<double> zcorn;
    std::vector<int> actnum;

    /// The corner-point processor's input; it points into this object.
    grdecl raw() const
    {
        grdecl g;
        g.dims[0] = dims[0]; g.dims[1] = dims[1]; g.dims[2] = dims[2];
        g.coord = coord.data();
        g.zcorn = zcorn.data();
        g.actnum = actnum.empty() ? nullptr : actnum.data();
        return g;
    }
};

// ZCORN runs over a doubled lattice: index i_ = 2*i + side in each direction. These
// give the cell i, the side (0 = low, 1 = high) and the pillar or layer node i + side.
inline int cellOf(int doubled) { return doubled / 2; }
inline int sideOf(int doubled) { return doubled % 2; }
inline int nodeOf(int doubled) { return cellOf(doubled) + sideOf(doubled); }

/// Pillar (i, j) as COORD's top and bottom points: x, y, z, x, y, z.
using Pillar = std::function<std::array<double,6>(int, int)>;
/// Depth of the corner (i_, j_, k_) on the doubled lattice.
using Depth = std::function<double(int, int, int)>;

/// Vertical unit-spaced pillars.
inline std::array<double,6> vertical(int i, int j)
{
    return { double(i), double(j), 0.0, double(i), double(j), 100.0 };
}

/// J runs toward -y: a left-handed grid.
inline std::array<double,6> leftHanded(int i, int j)
{
    return { double(i), -double(j), 0.0, double(i), -double(j), 100.0 };
}

/// Flat layers 2 thick.
inline double flat(int, int, int k_) { return 2.0*nodeOf(k_); }

inline TestGrdecl makeGrid(const std::array<int,3>& dims, const Pillar& pillar, const Depth& depth)
{
    TestGrdecl g;
    g.dims = dims;
    const auto& [nx, ny, nz] = dims;
    g.coord.resize(6*(nx + 1)*(ny + 1));
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            const auto p = pillar(i, j);
            std::copy(p.begin(), p.end(), &g.coord[6*(static_cast<std::size_t>(j)*(nx + 1) + i)]);
        }
    }
    g.zcorn.resize(8*static_cast<std::size_t>(nx)*ny*nz);
    for (int k_ = 0; k_ < 2*nz; ++k_) {
        for (int j_ = 0; j_ < 2*ny; ++j_) {
            for (int i_ = 0; i_ < 2*nx; ++i_) {
                g.zcorn[static_cast<std::size_t>(i_) + 2*static_cast<std::size_t>(nx)*j_
                        + 4*static_cast<std::size_t>(nx)*ny*k_] = depth(i_, j_, k_);
            }
        }
    }
    return g;
}

/// A grid on vertical pillars.
inline TestGrdecl makeGrid(const std::array<int,3>& dims, const Depth& depth)
{
    return makeGrid(dims, vertical, depth);
}

/// A uniform refinement request for the cells startIJK <= ijk < endIJK of the grid
/// parent.
inline Opm::Refinement::BlockRefinement box(const std::string& name, std::array<int,3> cellsPerDim,
                                            std::array<int,3> start, std::array<int,3> end,
                                            const std::string& parent = "GLOBAL")
{
    return {
        .name = name,
        .parentGridName = parent,
        .cellsPerDim = cellsPerDim,
        .startIJK = start,
        .endIJK = end
    };
}

/// Selects the conforming refinement, with a builder for the description g.
inline void useConformingBuilder(Dune::CpGrid& grid, const TestGrdecl& g, bool edgeConformal = false)
{
    grid.setLgrBackend(Opm::Refinement::Backend::Conforming);
    grid.setRefinementBuilder(std::make_shared<Opm::Refinement::ConformingBlockBuilder>(
        g.dims, g.coord, g.zcorn, g.actnum, edgeConformal));
}

/// Refines a grid already processed from g, by the conforming refinement.
inline void refine(Dune::CpGrid& grid, const TestGrdecl& g,
                   const std::vector<Opm::Refinement::BlockRefinement>& requests,
                   bool edgeConformal = false)
{
    useConformingBuilder(grid, g, edgeConformal);
    grid.addLgrsUpdateLeafView(requests);
}

inline double totalVolume(const Dune::CpGrid& grid)
{
    double volume = 0.0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        volume += element.geometry().volume();
    }
    return volume;
}

/// Largest |sum of area * outward normal| of a leaf cell, relative to its surface.
/// The largest closure error of a leaf cell, |sum of its outward face area
/// vectors|, in absolute terms (worstClosure divides it by the cell's surface).
inline double maxClosure(const Dune::CpGrid& grid)
{
    double worst = 0.0;
    for (const auto& e : Dune::elements(grid.leafGridView())) {
        Dune::FieldVector<double,3> closure(0.0);
        for (const auto& is : Dune::intersections(grid.leafGridView(), e)) {
            auto n = is.centerUnitOuterNormal();
            n *= is.geometry().volume();
            closure += n;
        }
        worst = std::max(worst, closure.two_norm());
    }
    return worst;
}

inline double worstClosure(const Dune::CpGrid& grid)
{
    double worst = 0.0;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        Dune::FieldVector<double,3> sum(0.0);
        double area = 0.0;
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            auto n = is.centerUnitOuterNormal();
            n *= is.geometry().volume();
            sum += n;
            area += is.geometry().volume();
        }
        worst = std::max(worst, sum.two_norm() / area);
    }
    return worst;
}

/// Faces between each pair of leaf cells, each face counted once.
inline std::map<std::pair<int,int>, int> facesPerPair(const Dune::CpGrid& grid)
{
    std::map<std::pair<int,int>, int> pairs;
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (is.neighbor() && (is.inside().index() < is.outside().index())) {
                ++pairs[{is.inside().index(), is.outside().index()}];
            }
        }
    }
    return pairs;
}

/// The most faces any pair of leaf cells has.
inline int maxFacesPerPair(const Dune::CpGrid& grid)
{
    int most = 0;
    for (const auto& [cells, count] : facesPerPair(grid)) {
        most = std::max(most, count);
    }
    return most;
}

} // namespace Opm::ConformingTest

#endif // OPM_CONFORMINGTESTHELPERS_HEADER_INCLUDED
