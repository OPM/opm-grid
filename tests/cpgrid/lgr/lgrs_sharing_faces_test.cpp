/*
  Copyright 2025 Equinor ASA.

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

#define BOOST_TEST_MODULE LgrsSharingFacesTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>

#include <tests/cpgrid/lgr/LgrChecks.hpp>

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

BOOST_AUTO_TEST_CASE(lgrsSharingIFaces)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {2,1,1}, /* cell_sizes = */ {1.,1.,1.});

    // I-face coincides between cells (0–1).
    // Ensure NY and NZ (y- and z-direction subdivisions) are identical for LGR1 and LGR2.
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{3,2,4}, {5,2,4}},
                               /* startIJK_vec = */ {{0,0,0}, {1,0,0}},
                               /* endIJK_vec = */ {{1,1,1}, {2,1,1}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2"});

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{3,2,4}, {5,2,4}},
                           /* lgr_name_vec = */ {"LGR1", "LGR2"},
                           /* isGlobalRefined = */ false);
}

BOOST_AUTO_TEST_CASE(lgrSharingJFaces)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {1,2,1}, /* cell_sizes = */ {1.,1.,1.});

    // J-face coincides between cells (0–1).
    // Ensure NX and NZ (x- and z-direction subdivisions) are identical for LGR1 and LGR2.
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{4,3,2}, {4,5,2}},
                               /* startIJK_vec = */ {{0,0,0}, {0,1,0}},
                               /* endIJK_vec = */ {{1,1,1}, {1,2,1}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2"});

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{4,3,2}, {4,5,2}},
                           /* lgr_name_vec = */ {"LGR1", "LGR2"},
                           /* isGlobalRefined = */ false);
}


BOOST_AUTO_TEST_CASE(lgrSharingKFaces)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {1,1,2}, /* cell_sizes = */ {1.,1.,1.});

    // K-face coincides between cells (0–1).
    // Ensure NX and NY (x- and y-direction subdivisions) are identical for LGR1 and LGR2.
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{4,2,3}, {4,2,5}},
                               /* startIJK_vec = */ {{0,0,0}, {0,0,1}},
                               /* endIJK_vec = */ {{1,1,1}, {1,1,2}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2"});

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{4,2,3}, {4,2,5}},
                           /* lgr_name_vec = */ {"LGR1", "LGR2"},
                           /* isGlobalRefined = */ false);
}

BOOST_AUTO_TEST_CASE(lgrsSharingJFaces)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.,1.,1.});
    //                          LGR1 parent cells    LGR2 parent cells   LGR3 parent cells   LGR4 parent cells
    // k = 2   32 33 34 35 |                                  33 34
    //         28 29 30 31 |          29 30                                                 28
    //         24 25 26 27 |          25 26
    // --------------------|    ------------------  ------------------  ------------------  ------------------
    // k = 1   20 21 22 23 |                                  21 22
    //         16 17 18 19 |          17 18                                                 16
    //         12 13 14 15 |          13 14
    // --------------------|    ------------------  ------------------  ------------------  ------------------
    // k = 0    8  9 10 11 |                                   9 10
    //          4  5  6  7 |                                                                 4
    //          0  1  2  3 |                                                 1  2
    //---------------------|
    // I-faces coincide between cells (16–17) and (28–29).
    // Ensure NY and NZ (y- and z-direction subdivisions) are identical for LGR1 and LGR4.
    // J-faces coincide between cells (17–21), (18–22), (29–33), and (30–34).
    // Ensure NX and NZ (x- and z-direction subdivisions) are identical for LGR1 and LGR2.
    // K-faces coincide between cells (1–13) and (2–14).
    // Ensure NX and NY (x- and y-direction subdivisions) are identical for LGR1 and LGR3.
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{4,1,2}, {4,3,2}, {4,1,3}, {5,1,2}},
                               /* startIJK_vec = */ {{1,0,1}, {1,2,0}, {1,0,0}, {0,1,0}},
                               /* endIJK_vec = */ {{3,2,3}, {3,3,3}, {3,1,1}, {1,2,3}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2", "LGR3", "LGR4"});

    Opm::checkGridWithLgrs(grid,
                           /* cells_per_dim_vec = */ {{4,1,2}, {4,3,2}, {4,1,3}, {5,1,2}},
                           /* lgr_name_vec = */ {"LGR1", "LGR2", "LGR3", "LGR4"},
                           /* isGlobalRefined = */ false);
}
