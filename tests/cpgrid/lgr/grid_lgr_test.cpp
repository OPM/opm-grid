/*
  Copyright 2022-2023, 2025 Equinor ASA.

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
#include "config.h"

#define BOOST_TEST_MODULE AddLgrsUpdateLeafViewInAllActiveCartesianGridTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>

#include <tests/cpgrid/LgrChecks.hpp>

struct Fixture
{
    Fixture()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(m_argc, m_argv);
        Opm::OpmLog::setupSimpleDefaultLogging();
    }
};

BOOST_GLOBAL_FIXTURE(Fixture);

BOOST_AUTO_TEST_CASE(refineCellBlockWithDifferentCellSizes)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {2.,1.,4.});
    // Refine a few cells, each into 2x2x2 children (2 subdivisions per x-,y-,z- direction).
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}},
                               /* startIJK_vec = */ {{1,0,1}},
                               /* endIJK_vec = */ {{3,2,3}},
                               /* lgr_name_vec = */ {"LGR1"});

    Opm::checkGridWithLgrs(grid, /* cells_per_dim_vec = */ {{2,2,2}}, /* lgr_name_vec = */ {"LGR1"}, /* isGlobalRefined = */ false);

    /** To-do: move out patchesShareFace. */ 
    BOOST_CHECK_EQUAL(grid.currentData()[0]->patchesShareFace({{1,0,1}}, {{3,2,3}}), false);
}

BOOST_AUTO_TEST_CASE(refineDisjointCellBlocksWithDifferentSubdivisionsPerDirection)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.,1.,1.});
    // Refine 3 cell-blocks, 2x2x2, 3x3x3, and 4x4x4 subdivisions per x-,y-,z- direction, respectively.
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}, {3,3,3}, {4,4,4}},
                               /* startIJK_vec = */{{0,0,0}, {0,0,2}, {3,2,2}},
                               /* endIJK_vec = */  {{2,1,1}, {1,1,3}, {4,3,3}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2", "LGR3"});

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{2,2,2}, {3,3,3}, {4,4,4}},
                           /* lgr_name_vec = */  {"LGR1", "LGR2", "LGR3"},
                           /* isGlobalRefined = */ false);

    BOOST_CHECK_EQUAL(grid.currentData()[0]->patchesShareFace({{0,0,0}, {0,0,2}, {3,2,2}}, {{2,1,1}, {1,1,3}, {4,3,3}}), false);
}

BOOST_AUTO_TEST_CASE(parentCellBlocksShareAcorner)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.,1.,1.});
    // Refine 3 cell-blocks, 2x2x2, 3x3x3, and 4x4x4 subdivisions per x-,y-,z- direction, respectively.
    // LGR1 parent cell and LGR2 parent cell share one corner.
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}, {3,3,3}, {4,4,4}},
                               /* startIJK_vec = */ {{0,0,0}, {1,1,1}, {3,2,2}},
                               /* endIJK_vec = */   {{1,1,1}, {2,2,2}, {4,3,3}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2", "LGR3"});

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{2,2,2}, {3,3,3}, {4,4,4}},
                           /* lgr_name_vec = */  {"LGR1", "LGR2", "LGR3"},
                           /* isGlobalRefined = */ false);

    BOOST_CHECK_EQUAL(grid.currentData()[0]->patchesShareFace({{0,0,0}, {1,1,1}, {3,2,2}}, {{1,1,1}, {2,2,2}, {4,3,3}}), false);
}

BOOST_AUTO_TEST_CASE(parentCellBlocksShareAnEdge)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.,1.,1.});
    // Refine 3 cell-blocks, 2x2x2 subdivisions per x-,y-,z- direction.
    // LGR1 parent cell (idx = 1) and LGR2 parent cell (idx = 14) share an edge.
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}, {2,2,2}, {2,2,2}},
                               /* startIJK_vec = */  {{0,0,0}, {2,0,1}, {3,2,2}},
                               /* endIJK_vec = */  {{2,1,1}, {3,1,2}, {4,3,3}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2", "LGR3"});

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{2,2,2}, {2,2,2}, {2,2,2}},
                           /* lgr_name_vec = */  {"LGR1", "LGR2", "LGR3"},
                           /* isGlobalRefined = */ false);

    BOOST_CHECK_EQUAL(grid.currentData()[0]->patchesShareFace({{0,0,0}, {2,0,1}, {3,2,2}}, {{2,1,1}, {3,1,2}, {4,3,3}}), false);
}

BOOST_AUTO_TEST_CASE(parentCellBloackShareAFaceWithCompatibleSubdivisions)
{   /** To-do: add/modify LGRs to have all possible sharing faces: I, J, and K */
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.,1.,1.});
    // Refine 3 cell-blocks, 2x2x2 subdivisions per x-,y-,z- direction.
    // LGR1 parent cell (idx = 1) and LGR2 parent cell (idx = 2) share an I_FACE.
    // Refinement takes place since LGR1 and LGR2 (y,z)-subdivisions are both (3,2).
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{5,3,2}, {4,3,2}, {2,2,2}},
                               /* startIJK_vec = */  {{0,0,0}, {2,0,0}, {3,2,2}},
                               /* endIJK_vec = */   {{2,1,1}, {3,1,1}, {4,3,3}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2", "LGR3"});

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{5,3,2}, {4,3,2}, {2,2,2}},
                           /* lgr_name_vec = */  {"LGR1", "LGR2", "LGR3"},
                           /* isGlobalRefined = */ false);
    BOOST_CHECK_EQUAL(grid.currentData()[0]->patchesShareFace( {{0,0,0}, {2,0,0}, {3,2,2}},  {{2,1,1}, {3,1,1}, {4,3,3}}), true);
}

BOOST_AUTO_TEST_CASE(throwIfUncompatibleNumberOfSubdivisionsInSharedFaces)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    const std::vector<std::array<int,3>> cells_per_dim_vec = {{4,2,5}, {3,2,2}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};

    // Parent cells LGR1 and LGR2 share I_FACE with uncompatible number of subdivisions.
    // LGR1 (y,z)-subdivisions: (2,5) != (2,2) LGR2 (y,z)-subdivisions.
    const std::vector<std::array<int,3>> startIJK_iFace = {{0,0,0}, {1,0,0}};
    const std::vector<std::array<int,3>> endIJK_iFace = {{1,1,1}, {2,1,1}};
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_iFace, endIJK_iFace, lgr_name_vec),
                      std::logic_error);

    // Parent cells LGR1 and LGR2 share J_FACE with uncompatible number of subdivisions.
    // LGR1 (x,z)-subdivisions: (4,5) != (3,2) LGR2 (x,z)-subdivisions.
    const std::vector<std::array<int,3>> startIJK_jFace = {{0,0,0}, {0,1,0}};
    const std::vector<std::array<int,3>> endIJK_jFace = {{1,1,1}, {1,2,1}};
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_jFace, endIJK_jFace, lgr_name_vec),
                      std::logic_error);

    // Parent cells LGR1 and LGR2 share K_FACE with uncompatible number of subdivisions.
    // LGR1 (x,y)-subdivisions: (4,2) != (3,2) LGR2 (x,y)-subdivisions.
    const std::vector<std::array<int,3>> startIJK_kFace = {{0,0,0}, {0,0,1}};
    const std::vector<std::array<int,3>> endIJK_kFace = {{1,1,1}, {1,1,2}};
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_kFace, endIJK_kFace, lgr_name_vec),
                      std::logic_error);
}

BOOST_AUTO_TEST_CASE(refineEntireGrid)
{
    // Create a 4x3x3 grid with sizes 4x3x3
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    // Refine each cell into 2x2x2 children (2 subdivisions per x-,y-,z- direction).
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}},
                               /* startIJK_vec = */ {{0,0,0}},
                               /* endIJK_vec = */ { {4,3,3}},
                               /* lgr_name_vec = */ {"LGR1"});

    // Create a 8x6x6 grid with sizes 4x3x3
    Dune::CpGrid fine_grid;
    fine_grid.createCartesian(/* grid_dim = */ {8,6,6}, /* cell_sizes = */ {0.5, 0.5, 0.5});

    Opm::checkLeafGridGeometryEquality(grid, fine_grid);
}


BOOST_AUTO_TEST_CASE(doNothingIfSubdivisionsInAllDirectionsIsOne)
{
    // Create a 4x3x3 grid with sizes 4x3x3
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1., 1., 1.});
    // Refine each cell into 1 child does nothing to the grid. (1 subdivision per x-,y-,z- direction).
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{1,1,1}},
                               /* startIJK_vec = */ {{0,0,0}},
                               /* endIJK_vec = */ { {4,3,3}},
                               /* lgr_name_vec = */ {"LGR1"});

    // Create a 4x3x3 grid with sizes 4x3x3
    Dune::CpGrid fine_grid;
    fine_grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.,1.,1.});

    Opm::checkLeafGridGeometryEquality(grid, fine_grid);
}
