/*
  Copyright 2026 SINTEF Digital, Mathematics and Cybernetics.

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

#define NVERBOSE

#define BOOST_TEST_MODULE PartitionOfUnityTest
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <opm/grid/CpGrid.hpp>

#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>

#include <array>
#include <vector>

namespace
{

/// Build a faulted 4x4x2 corner-point grid: the right half (i >= 2) is
/// shifted down by half a cell.  Face processing along the fault then
/// produces nodes that are NOT among the eight canonical corners of the
/// neighbouring cells (hanging nodes) -- the situation in which point
/// (codim 3) communication interfaces built from cell_to_point_ miss
/// nodes entirely.
Opm::EclipseGrid faultedGrid(const int nx, const int ny, const int nz)
{
    const double shift = 0.5;

    // Straight vertical pillars on the unit lattice.
    std::vector<double> coord;
    coord.reserve((nx + 1) * (ny + 1) * 6);
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            coord.push_back(i);      // x top
            coord.push_back(j);      // y top
            coord.push_back(0.0);    // z top
            coord.push_back(i);      // x bottom
            coord.push_back(j);      // y bottom
            coord.push_back(nz + 1.0); // z bottom
        }
    }

    // ZCORN in Eclipse ordering: per layer, top surface then bottom surface,
    // each surface row-wise with two entries per cell in each direction.
    std::vector<double> zcorn;
    zcorn.reserve(8 * nx * ny * nz);
    for (int k = 0; k < nz; ++k) {
        for (int face = 0; face < 2; ++face) {          // 0 = top, 1 = bottom
            for (int j = 0; j < ny; ++j) {
                for (int jj = 0; jj < 2; ++jj) {
                    (void) jj;
                    for (int i = 0; i < nx; ++i) {
                        const double dz = (i >= nx / 2) ? shift : 0.0;
                        const double z = k + face + dz;
                        zcorn.push_back(z);
                        zcorn.push_back(z);
                    }
                }
            }
        }
    }

    return Opm::EclipseGrid(std::array<int, 3>{nx, ny, nz}, coord, zcorn);
}

/// Codim-3 data handle: senders write their rank for every vertex in the
/// communication interface; receivers count the messages per vertex.
class VertexReachedHandle
{
public:
    using DataType = int;

    VertexReachedHandle(const Dune::CpGrid& grid,
                        std::vector<int>& received,
                        const int rank)
        : grid_(grid), received_(received), rank_(rank)
    {}

    bool fixedSize(int /*dim*/, int /*codim*/)
    {
        return true;
    }
    bool contains(int dim, int codim)
    {
        return dim == 3 && codim == 3;
    }
    template<class T>
    std::size_t size(const T&)
    {
        return 1;
    }
    template<class B, class T>
    void gather(B& buffer, const T&)
    {
        buffer.write(rank_);
    }
    template<class B, class T>
    void scatter(B& buffer, const T& t, std::size_t)
    {
        int sender = -1;
        buffer.read(sender);
        ++received_[grid_.leafIndexSet().index(t)];
    }

private:
    const Dune::CpGrid& grid_;
    std::vector<int>& received_;
    int rank_;
};

} // anonymous namespace

// Every vertex on the processor boundary (partition type != interior) must be
// reached by codim-3 communication.  On a faulted grid the fault-face hanging
// nodes are not canonical corners of the neighbouring cells; interfaces built
// from the canonical corners only never include them, so they receive no
// message and their ownership cannot be made consistent across ranks (the
// "Owner is not a partition of unity" failure seen in vertex-based mechanics).
BOOST_AUTO_TEST_CASE(FaultedGridVerticesReached)
{
    const auto& helper = Dune::MPIHelper::instance(
        boost::unit_test::framework::master_test_suite().argc,
        boost::unit_test::framework::master_test_suite().argv);
    const int size = helper.size();
    const int rank = helper.rank();

    if (size < 2) {
        return; // needs at least two processes to have processor boundaries
    }

    const int nx = 4, ny = 4, nz = 2;
    const auto eclGrid = faultedGrid(nx, ny, nz);

    Dune::CpGrid grid;
    grid.processEclipseFormat(&eclGrid, nullptr,
                              /* periodic_extension = */ false,
                              /* turn_normals = */ false,
                              /* clip_z = */ false,
                              /* pinchActive = */ false,
                              /* edge_conformal = */ false);

    // Partition in vertical slabs of columns so the fault plane (i == nx/2)
    // coincides with a processor boundary.
    std::vector<int> parts(grid.size(0));
    const auto& gv = grid.leafGridView();
    for (const auto& element : elements(gv)) {
        const int cell = gv.indexSet().index(element);
        const auto center = element.geometry().center();
        const int i = static_cast<int>(center[0]); // unit lattice
        const int j = static_cast<int>(center[1]);
        if (size >= 4) {
            parts[cell] = (i >= nx / 2) + 2 * (j >= ny / 2);
        } else {
            parts[cell] = (i >= nx / 2);
        }
    }

    grid.loadBalance(parts, /* ownersFirst = */ false,
                     /* addCornerCells = */ true, /* overlapLayers = */ 1);

    const int numVertices = grid.leafIndexSet().size(3);
    std::vector<int> received(numVertices, 0);

    VertexReachedHandle handle(grid, received, rank);
    grid.communicate(handle, Dune::All_All_Interface,
                     Dune::ForwardCommunication);

    // Every non-interior vertex sits on a processor boundary and is shared
    // with at least one other rank, so it must have received at least one
    // message.  Hanging nodes on the fault violate this if the interface was
    // built from the canonical corners only.
    int unreached = 0;
    for (const auto& vertex : vertices(grid.leafGridView())) {
        const int idx = grid.leafIndexSet().index(vertex);
        if (vertex.partitionType() != Dune::InteriorEntity
            && received[idx] == 0) {
            ++unreached;
        }
    }

    BOOST_CHECK_MESSAGE(unreached == 0,
                        "rank " << rank << ": " << unreached
                        << " non-interior vertices were not reached by codim-3"
                        " communication (hanging nodes missing from the point"
                        " interface)");
}

bool
init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
