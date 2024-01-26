//===========================================================================
//
// File: lookupdataCpGrid_test.cpp
//
// Created: Thurs 25.05.2023 16:05:00
//
// Author(s): Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================
/*
  Copyright 2023 Equinor ASA.

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

#define BOOST_TEST_MODULE LookUpDataCpGridTests
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/LookUpData.hh>

#include <dune/grid/common/mcmgmapper.hh>


#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Deck/DeckSection.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckKeyword.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>

struct Fixture
{
    Fixture()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(m_argc, m_argv);
        Opm::OpmLog::setupSimpleDefaultLogging();
    }

    static int rank()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        return Dune::MPIHelper::instance(m_argc, m_argv).rank();
    }
};

BOOST_GLOBAL_FIXTURE(Fixture);

void lookup_check(const Dune::CpGrid& grid)
{
    const auto& data = grid.data_;

    std::vector<int> fake_feature(data[0]->size(0), 0);
    std::iota(fake_feature.begin(), fake_feature.end(), 3);

    std::vector<double> fake_feature_double(data[0]->size(0), 0.);
    std::iota(fake_feature_double.begin(), fake_feature_double.end(), .5);

    std::vector<std::vector<int>> fakeLgrFeatures;
    fakeLgrFeatures.resize(data.size()-1);
    // Creating fake field properties for each LGR
    if (data.size()>1) {
        for (long unsigned int lgr = 1; lgr < data.size(); ++lgr)
        {
            std::vector<int> fake_feature_lgr(data[lgr]->size(0), lgr);
            fakeLgrFeatures[lgr-1] = fake_feature_lgr;
        }
    }


    // LookUpData
    const auto& leaf_view = grid.leafGridView();
    const Opm::LookUpData<Dune::CpGrid, Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>> lookUpData(leaf_view);
    // LookUpData with LGR field properties
    const Opm::LookUpData<Dune::CpGrid, Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>> lookUpDataLGR(leaf_view, true);
    // LookUpCartesianData
    const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapper(grid);
    const Opm::LookUpCartesianData<Dune::CpGrid, Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>
        lookUpCartesianData(leaf_view, cartMapper);
    // LookUpCartesianData with LGR field properties
    const Opm::LookUpCartesianData<Dune::CpGrid, Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>
        lookUpCartesianDataLGR(leaf_view, cartMapper, true);

    const auto& level0_view = grid.levelGridView(0);
    const Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LeafGridView> leafMapper(leaf_view, Dune::mcmgElementLayout());
    const Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LevelGridView> level0Mapper(level0_view, Dune::mcmgElementLayout());

    const auto& leaf_idSet = (*data.back()).localIdSet();
    const auto& level0_idSet = (*data[0]).localIdSet();

    for (const auto& elem : elements(leaf_view)) {
        // Search via Entity/Element
        const auto featureInElem = lookUpData(elem, fake_feature);
        const auto featureInElemDouble = lookUpData(elem, fake_feature_double);
        const auto featureInElemCartesian = lookUpCartesianData(elem, fake_feature);
        const auto featureInElemDoubleCartesian = lookUpCartesianData(elem,fake_feature_double);
        BOOST_CHECK(featureInElem == lookUpData.getFieldPropIdx(elem) +3);
        BOOST_CHECK(featureInElemDouble == lookUpData.getFieldPropIdx(elem) +.5);
        BOOST_CHECK(featureInElemCartesian == lookUpData.getFieldPropIdx(elem) +3);
        BOOST_CHECK(featureInElemDoubleCartesian == lookUpData.getFieldPropIdx(elem) +.5);
        BOOST_CHECK(elem.getOrigin().index() == lookUpData.getFieldPropIdx(elem));
        // Search via INDEX
        const auto featureInElemIDX = lookUpData(elem.index(), fake_feature);
        const auto featureInElemDoubleIDX = lookUpData(elem.index(), fake_feature_double);
        const auto featureInElemCartesianIDX = lookUpCartesianData(elem.index(), fake_feature);
        const auto featureInElemDoubleCartesianIDX = lookUpCartesianData(elem.index(),fake_feature_double);
        BOOST_CHECK(featureInElemIDX == lookUpData.getFieldPropIdx(elem.index()) +3);
        BOOST_CHECK(featureInElemDoubleIDX == lookUpData.getFieldPropIdx(elem.index()) +.5);
        BOOST_CHECK(featureInElemCartesianIDX == lookUpData.getFieldPropIdx(elem.index()) +3);
        BOOST_CHECK(featureInElemDoubleCartesianIDX == lookUpData.getFieldPropIdx(elem.index()) +.5);
        BOOST_CHECK(elem.getOrigin().index() == lookUpData.getFieldPropIdx(elem.index()));
        BOOST_CHECK(featureInElemIDX == featureInElem);
        BOOST_CHECK(featureInElemDoubleIDX == featureInElemDouble);
        BOOST_CHECK(featureInElemCartesianIDX == featureInElemCartesian);
        BOOST_CHECK(featureInElemDoubleCartesianIDX == featureInElemDoubleCartesian);
        // Extra checks related to Cartesian Index
        const auto cartIdx = cartMapper.cartesianIndex(elem.index());
        BOOST_CHECK(cartIdx == lookUpCartesianData.getFieldPropCartesianIdx(elem));
        BOOST_CHECK(cartIdx == lookUpCartesianData.getFieldPropCartesianIdx(elem.index()));
        // Extra checks related to Cartesian Coordinate
        std::array<int,3> ijk;
        cartMapper.cartesianCoordinate(elem.index(), ijk); // this ijk corresponds to the parent/equivalent cell in level 0.
        std::array<int,3> ijkLevel0;
        cartMapper.cartesianCoordinateLevel(elem.getOrigin().index(), ijkLevel0, 0);
        BOOST_CHECK(ijk == ijkLevel0);
        // Throw for level < 0 or level > maxLevel()
        std::array<int,3> ijkThrow;
        BOOST_CHECK_THROW(cartMapper.cartesianCoordinateLevel(elem.index(), ijkThrow, -3), std::invalid_argument);
        BOOST_CHECK_THROW(cartMapper.cartesianCoordinateLevel(elem.index(), ijkThrow, grid.maxLevel() + 1), std::invalid_argument);
        // Checks related to LGR field properties
        if (elem.level())
        {
            const auto featureInLGR = lookUpDataLGR(elem, fakeLgrFeatures[elem.level()-1]);
            const auto featureInLGR_Cartesian = lookUpCartesianDataLGR(elem, fakeLgrFeatures[elem.level()-1]);
            const auto featureInLGR_FromIdx = lookUpDataLGR(elem.index(), fakeLgrFeatures[elem.level()-1]);
            const auto featureInLGR_Cartesian_FromIdx = lookUpCartesianDataLGR(elem.index(), fakeLgrFeatures[elem.level()-1]);
            BOOST_CHECK(featureInLGR == elem.level());
            BOOST_CHECK(featureInLGR == featureInLGR_Cartesian);
            BOOST_CHECK(featureInLGR == featureInLGR_FromIdx);
            BOOST_CHECK(featureInLGR == featureInLGR_Cartesian_FromIdx);
            // Checks for CartesianCoordinateLevel
            const auto idxOnLevel = elem.getLevelElem().index(); // getLevelElemt throws when entity does not belong to the leafGridView
            std::array<int,3> ijkLevelGrid;
            (*grid.data_[elem.level()]).getIJK(idxOnLevel, ijkLevelGrid);
            std::array<int,3> ijkLevel;
            cartMapper.cartesianCoordinateLevel(idxOnLevel, ijkLevel, elem.level());
            BOOST_CHECK( ijkLevelGrid == ijkLevel);
        }
        // Extra checks related to ElemMapper
        BOOST_CHECK(featureInElem == level0Mapper.index(elem.getOrigin()) +3);
        BOOST_CHECK(featureInElem == fake_feature[lookUpData.getFieldPropIdx(elem)]);
        if (elem.hasFather()) { // leaf_cell has a father!
            const auto id = leaf_idSet.id(elem);
            const auto parent_id = level0_idSet.id(elem.father());
            BOOST_CHECK(elem.index() == id);
            BOOST_CHECK(elem.index() == leafMapper.index(elem));
            BOOST_CHECK(elem.father().index() == featureInElem -3);
            BOOST_CHECK(elem.father().index() == parent_id);
            BOOST_CHECK(elem.father().index() == level0Mapper.index(elem.father()));
        }
    }
}

BOOST_AUTO_TEST_CASE(one_lgr_grid)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);
    // Add LGRs and update LeafGridView
    const std::array<int, 3> cells_per_dim = {2,2,2};
    const std::array<int, 3> startIJK = {1,0,1};
    const std::array<int, 3> endIJK = {3,2,3};  // patch_dim = {3-1, 2-0, 3-1} ={2,2,2}
    const std::string lgr_name = {"LGR1"};
    grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

    lookup_check(grid);
}

BOOST_AUTO_TEST_CASE(single_cell_lgr_grid)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);
    // Add LGRs and update LeafGridView
    const std::array<int, 3> cells_per_dim = {2,2,2};
    const std::array<int, 3> startIJK = {1,0,1};
    const std::array<int, 3> endIJK = {2,1,2};  // patch_dim = {2-1, 1-0, 2-1} ={1,1,1} -> Single Cell!
    const std::string lgr_name = {"LGR1"};
    grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

    lookup_check(grid);
}

BOOST_AUTO_TEST_CASE(lgrs_grid_A)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);
    // Add LGRs and update LeafGridView
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}, {4,4,4}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,2}, {3,2,2}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,1,1}, {1,1,3}, {4,3,3}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3"};
    grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

    lookup_check(grid);
}

BOOST_AUTO_TEST_CASE(lgrs_grid_B)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);
    // Add LGRs and update LeafGridView
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {3,2,0}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,2,1}, {4,3,3}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

    lookup_check(grid);
}


BOOST_AUTO_TEST_CASE(lgrs_grid_C)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {5,4,4};
    grid.createCartesian(grid_dim, cell_sizes);
    // Add LGRs and update LeafGridView
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,3,4}, {3,2,4}, {4,3,2}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {4,0,0}, {4,3,3}};
    const std::vector<std::array<int,3>> endIJK_vec = {{3,2,2}, {5,2,1}, {5,4,4}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3"};
    grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

    lookup_check(grid);
}

BOOST_AUTO_TEST_CASE(no_lgrs_grid)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);
    lookup_check(grid);
}


void fieldProp_check(const Dune::CpGrid& grid, Opm::EclipseGrid eclGrid, std::string deck_string)
{
    Opm::Deck deck = Opm::Parser{}.parseString(deck_string);
    Opm::FieldPropsManager fpm(deck, Opm::Phases{true, true, true}, eclGrid, Opm::TableManager());
    const auto& poro = fpm.get_double("PORO");
    const auto& eqlnum =  fpm.get_int("EQLNUM");

    // LookUpData
    const auto& leaf_view = grid.leafGridView();
    const Opm::LookUpData<Dune::CpGrid,Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>> lookUpData(leaf_view);
    // LookUpCartesianData
    const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapper(grid);
    const Opm::LookUpCartesianData<Dune::CpGrid, Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>
        lookUpCartesianData(leaf_view, cartMapper);
    // Element mapper
    const Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>
        mapper(leaf_view, Dune::mcmgElementLayout());

    const auto& poroOnLeaf = lookUpData.assignFieldPropsDoubleOnLeaf(fpm, "PORO");
    const auto& poroOnLeafCart = lookUpCartesianData.assignFieldPropsDoubleOnLeaf(fpm, "PORO");

    const auto& eqlnumOnLeaf = lookUpData.assignFieldPropsIntOnLeaf<int>(fpm, "EQLNUM", true);
    const auto& eqlnumOnLeafCart = lookUpCartesianData.assignFieldPropsIntOnLeaf<int>(fpm, "EQLNUM", true);

    for (const auto& elem : elements(leaf_view))
    {
        const auto elemIdx = mapper.index(elem);
        const auto elemOriginIdx = elem.getOrigin().index();
        // PORO
        BOOST_CHECK_EQUAL(poro[elemOriginIdx], lookUpData.fieldPropDouble<Dune::cpgrid::Entity<0>>(fpm, "PORO", elem));
        BOOST_CHECK_EQUAL(poro[elemOriginIdx], lookUpData.fieldPropDouble<int>(fpm, "PORO", elemIdx));
        BOOST_CHECK_EQUAL(poro[elemOriginIdx], lookUpCartesianData.fieldPropDouble<Dune::cpgrid::Entity<0>>(fpm, "PORO", elem));
        BOOST_CHECK_EQUAL(poro[elemOriginIdx], lookUpCartesianData.fieldPropDouble(fpm, "PORO", elemIdx));
        BOOST_CHECK_EQUAL(poro[elemOriginIdx], poroOnLeaf[elemIdx]);
        BOOST_CHECK_EQUAL(poro[elemOriginIdx], poroOnLeafCart[elemIdx]);
        // EQLNUM
        BOOST_CHECK_EQUAL(eqlnum[elemOriginIdx], lookUpData.fieldPropInt<Dune::cpgrid::Entity<0>>(fpm, "EQLNUM", elem));
        BOOST_CHECK_EQUAL(eqlnum[elemOriginIdx], lookUpData.fieldPropInt<int>(fpm, "EQLNUM", elemIdx));
        BOOST_CHECK_EQUAL(eqlnum[elemOriginIdx], lookUpCartesianData.fieldPropInt<Dune::cpgrid::Entity<0>>(fpm, "EQLNUM", elem));
        BOOST_CHECK_EQUAL(eqlnum[elemOriginIdx], lookUpCartesianData.fieldPropInt(fpm, "EQLNUM", elemIdx));
        BOOST_CHECK_EQUAL(eqlnum[elemOriginIdx]-true, eqlnumOnLeaf[elemIdx]);
        BOOST_CHECK_EQUAL(eqlnum[elemOriginIdx]-true, eqlnumOnLeafCart[elemIdx]);
    }
}


BOOST_AUTO_TEST_CASE(fieldProp) {
  const std::string deckString =
        R"( RUNSPEC
        DIMENS
        1  1  5 /
        GRID
        COORD
        0 0 0
        0 0 1
        1 0 0
        1 0 1
        0 1 0
        0 1 1
        1 1 0
        1 1 1
        /
        ZCORN
        4*0
        8*1
        8*2
        8*3
        8*4
        4*5
        /
        ACTNUM
        5*1
        /
        PORO
        1.0 2.0 3.0 4.0 5.0
        /
        REGIONS

        EQLNUM
        1 2 3 4 5
        /)";


    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);

    Dune::CpGrid grid;
    Opm::EclipseGrid eclGrid(deck);

    grid.processEclipseFormat(&eclGrid, nullptr, false, false, false);

    fieldProp_check(grid, eclGrid, deckString);
}


BOOST_AUTO_TEST_CASE(fieldPropLgr) {
  const std::string deckString = R"( RUNSPEC
DIMENS
1 1 5
/
GRID
DX
-- There are in total 1 cell with length 1ft in x-direction
1*1
/
DY
-- There are in total 1 cell with length 1ft in y-direction
1*1
/
DZ
-- The layers are 2 ft thick, in each layer there is 1 cell
1*2 1*2 1*2 1*2 1*2
/
TOPS
40*1
/
ACTNUM
5*1
/
PORO
1.0 2.0 3.0 4.0 5.0
/
REGIONS

EQLNUM
1 2 3 4 5
/)";

    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);

    Dune::CpGrid grid;
    Opm::EclipseGrid eclGrid(deck);

    grid.processEclipseFormat(&eclGrid, nullptr, false, false, false);

    // Add LGRs and update LeafGridView
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {2,2,2}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,3}};
    const std::vector<std::array<int,3>> endIJK_vec = {{1,1,1}, {1,1,4}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

    fieldProp_check(grid, eclGrid, deckString);
}

