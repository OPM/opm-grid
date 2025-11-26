/*
  Copyright 2023, 2025 Equinor ASA.

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

#define BOOST_TEST_MODULE CuboidShapeTest
#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/LgrHelpers.hpp>

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

bool checkCuboidShape(const Dune::cpgrid::Entity<0>& element, const Dune::CpGrid& grid)
{
    const auto& leafGrid = grid.currentLeafData();
    const auto cellToPoint = leafGrid.cellToPoint(element.index());

    // Select four corners to approximate cuboid volume
    //    2 ---- 3                          6 ---- 7
    //   /      /   bottom face            /      /    top face
    //  0 ---- 1                          4 ---- 5
    std::vector<Dune::cpgrid::Geometry<0,3>::GlobalCoordinate> aFewCorners;
    aFewCorners.resize(4); // {'0', '1', '3', '5'}
    aFewCorners[0] = (*(leafGrid.getGeometry().geomVector(std::integral_constant<int,3>()))).get(cellToPoint[0]).center();
    aFewCorners[1] = (*(leafGrid.getGeometry().geomVector(std::integral_constant<int,3>()))).get(cellToPoint[1]).center();
    aFewCorners[2] = (*(leafGrid.getGeometry().geomVector(std::integral_constant<int,3>()))).get(cellToPoint[3]).center();
    aFewCorners[3] = (*(leafGrid.getGeometry().geomVector(std::integral_constant<int,3>()))).get(cellToPoint[5]).center();

    auto distance = [](const auto& p1, const auto& p2) {
        double dx = p2[0] - p1[0];
        double dy = p2[1] - p1[1];
        double dz = p2[2] - p1[2];
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    };

    double length  = distance(aFewCorners[0], aFewCorners[1]);
    double breadth = distance(aFewCorners[1], aFewCorners[2]);
    double height  = distance(aFewCorners[1], aFewCorners[3]);

    const double cuboidVolume = length*breadth*height;
    const auto actualVolume =  (*(leafGrid.getGeometry().geomVector(std::integral_constant<int,0>())))[Dune::cpgrid::EntityRep<0>(element.index(), true)].volume();

    return (std::abs(cuboidVolume - actualVolume) <  1e-6);
}

BOOST_AUTO_TEST_CASE(refineCellWithFewerThanEightCornersThrows)
{
    const std::string deckString =
        R"(RUNSPEC
        DIMENS
        1  1  1/
        GRID
        -- COORD: 4 pillars × 6 values (x1 y1 z1  x2 y2 z2)
        COORD
        0 0 0  0 0 1   -- Pillar 1: bottom at (0,0,0), top at (0,0,1)
        1 0 0  1 0 0   -- Pillar 2: bottom at (1,0,0) = same at top
        0 1 0  0 1 1   -- Pillar 3: bottom at (0,1,0), top at (0,1,1)
        1 1 0  1 1 0   -- Pillar 4: bottom at (1,1,0) = same at top
        /
        -- ZCORN: eight Z values: top 4, then bottom 4
        ZCORN
        6*0 -- top end of the pillars
        2*1 -- bottom end of the pillars
        /
        ACTNUM
        1*1
        /
        )";

    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);
    Opm::EclipseGrid eclGrid(deck);

    Dune::CpGrid grid;
    grid.processEclipseFormat(&eclGrid, nullptr, false, false, false);

    for (const auto& element : Dune::elements(grid.leafGridView())) {
        // Check the single-cell from level zero is not a cuboid
        BOOST_CHECK(!checkCuboidShape(element, grid));
        // Check the single-cell from level zero has fewer than 8 corners
        BOOST_CHECK_THROW(Opm::Lgr::containsEightDifferentCorners( grid.currentLeafData().cellToPoint(element.index())),
                          std::logic_error);

        // Mark element for refinment, to check (below) that calling adapt() throws.
        grid.mark(1,element);
    }
    // Throw for all four refinement methods.
    BOOST_CHECK_THROW(grid.addLgrsUpdateLeafView({{2,2,3}}, // cells_per_dim_vec
                                                 {{0,0,0}}, // startIJK_vec
                                                 {{1,1,1}}, // endIJK_vec
                                                 {"LGR1"}), // lgr_name_vec
                      std::logic_error);

    BOOST_CHECK_THROW(grid.globalRefine(2), std::logic_error);

    BOOST_CHECK_THROW(grid.autoRefine({3,5,7}), // refinement factors in (x,y,z) directions
                      std::logic_error);

    // Single-cell has been marked for refinement, even though refinement won't take place (<8 vertices)
    grid.preAdapt();
    BOOST_CHECK_THROW(grid.adapt(), std::logic_error);
}

// Create a test grid from a simple deckstring.
// Refine it with one of:
// 0-> addLgrsUpdateLeafView
// 1-> globalRefine
// 2-> autoRefine
// 3-> adapt
void createAndRefineTestGrid(const std::string& deckString,
                             int selectMethod)
{
    // Create the grid from string
    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);
    Opm::EclipseGrid eclGrid(deck);
    Dune::CpGrid grid;
    grid.processEclipseFormat(&eclGrid, nullptr, false, false, false);

    for (const auto& element : Dune::elements(grid.leafGridView())) {
        // Check the cell(s) from level zero is(are) not a cuboid
        BOOST_CHECK(!checkCuboidShape(element, grid));

        // Check the cell(s) from level zero has(have) 8 corners
        Opm::Lgr::containsEightDifferentCorners(grid.currentLeafData().cellToPoint(element.index()));

        // If selectMethod == 3-> adapt -> mark at least one element
        if ((selectMethod == 3) && (element.index()==0)){
            grid.mark(1,element);
        }
    }

    if (selectMethod == 0) {
        grid.addLgrsUpdateLeafView({{3,2,4}}, // cells_per_dim_vec
                                   {{0,0,0}}, // startIJK_vec
                                   {{1,1,1}}, // endIJK_vec
                                   {"LGR1"}); // lgr_name_vec
    }
    else if (selectMethod == 1){
        grid.globalRefine(2);
    }
    else if (selectMethod == 2) {
        grid.autoRefine(std::array{3,5,7}); // refinement factors in (x,y,z) directions
    }
    else if (selectMethod == 3) {
        grid.preAdapt();
        grid.adapt();
        grid.postAdapt();
    }

    for (const auto& element : Dune::elements(grid.leafGridView())) {
        // Check all leaf cells are also not cuboids, after refinement of any kind
        BOOST_CHECK(!checkCuboidShape(element, grid));

        // Check the cell has 8 different corners
        Opm::Lgr::containsEightDifferentCorners(grid.currentLeafData().cellToPoint(element.index()));
    }
}

BOOST_AUTO_TEST_CASE(singleHexahedronNotCuboidCanBeRefined)
{
    const std::string singleHexa = R"(
RUNSPEC
DIMENS
 1 1 1 /
GRID
-- COORD: 4 pillars × 6 values (x1 y1 z1  x2 y2 z2)
COORD
 0 0 0    0 0 2    -- Pillar 1: bottom at (0,0,0), top at (0,0,2)
 2 0 0    3 1 3    -- Pillar 2: bottom at (2,0,0), top at (3,1,3)
 0 2 0   -1 3 1    -- Pillar 3: bottom at (0,2,0), top at (-1,3,1)
 2 2 0    4 4 4    -- Pillar 4: bottom at (2,2,0), top at (4,4,4)
/
-- ZCORN: eight Z values: top 4, then bottom 4
ZCORN
 0 0 0 0 -- top end of the pillars
 2 3 1 4 -- bottom end of the pillars
/
ACTNUM
 1
/
)";

    createAndRefineTestGrid(singleHexa, 0); // 0-> refinement via addLgrsUpdateLeafView
    createAndRefineTestGrid(singleHexa, 1); // 1-> refinement via globalRefine
    createAndRefineTestGrid(singleHexa, 2); // 2-> refinement via autoRefine
    createAndRefineTestGrid(singleHexa, 3); // 3-> refinement via adpat
}

BOOST_AUTO_TEST_CASE(twoHexahedronNotCuboidCanBeRefined)
{

    const std::string twoHexaZ = R"(
RUNSPEC
DIMENS
 1 1 2 /
GRID
-- COORD: 4 pillars × 6 values (x1 y1 z1  x2 y2 z2)
COORD
 0 0 0    0 0 6     -- Pillar 1 bottom at (0,0,0), top at (0,0,6)
 2 0 0    3 1 7     -- Pillar 2 bottom at (2,0,0), top at (3,1,7)
 0 2 0   -1 3 5     -- Pillar 3 bottom at (0,2,0), top at (-1,3,5)
 2 2 0    4 4 8     -- Pillar 4 bottom at (2,2,0), top at (4,4,8)
/
ZCORN
-- ZCORN: 2 cells stacked vertically -> 16 values
-- Cell 2 (above cell 1)
2 3 1 4   -- top end of the pillars
5 7 4 8   -- bottom end of the pillars = same as cell 1 top
-- Cell 1 (bottom)
0 0 0 0   -- top end of the pillars
2 3 1 4   -- bottom end of the pillars
/
ACTNUM
 2*1
/
)";

    createAndRefineTestGrid(twoHexaZ,0); // 0-> refinement via addLgrsUpdateLeafView
    createAndRefineTestGrid(twoHexaZ,1); // 1-> refinement via globalRefine
    createAndRefineTestGrid(twoHexaZ,2); // 2-> refinement via addLgrsUpdateLeafView
    createAndRefineTestGrid(twoHexaZ,3); // 3-> refinement via adapt
}
