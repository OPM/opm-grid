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

#define BOOST_TEST_MODULE LookUpDataCpGridTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>
#include <opm/grid/LookUpData.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
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
};

BOOST_GLOBAL_FIXTURE(Fixture);

using LeafMapper = Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LeafGridView>;
using LevelMapper = Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LevelGridView>;

using LookUpData =  const Opm::LookUpData<Dune::CpGrid, Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>;
using LookUpCartesianData = const Opm::LookUpCartesianData<Dune::CpGrid, Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>;

void lookup_check(const Dune::CpGrid& grid)
{
    const auto& level0_view = grid.levelGridView(0);
    // Create a vector representing a field property of type int,
    // in (level zero) GLOBAL grid.
    std::vector<int> fake_feature(level0_view.size(0), 0);
    std::iota(fake_feature.begin(), fake_feature.end(), 3);

    // Create a vector representing a field property of type double,
    // in (level zero) GLOBAL grid.
    std::vector<double> fake_feature_double(level0_view.size(0), 0.);
    std::iota(fake_feature_double.begin(), fake_feature_double.end(), .5);

    // For each refined level grid (LGR), create a vector representing
    // a field property of type int.
    std::vector<std::vector<int>> fakeLgrFeatures;
    fakeLgrFeatures.resize(grid.maxLevel());
    if (grid.maxLevel()>0) {
        for (int level = 1; level <= grid.maxLevel(); ++level) {
            std::vector<int> fake_feature_lgr(grid.levelGridView(level).size(0), level);
            fakeLgrFeatures[level-1] = fake_feature_lgr;
        }
    }

    const auto& leaf_view = grid.leafGridView();

    const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapper(grid);
    const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);

    const LeafMapper leafMapper(leaf_view, Dune::mcmgElementLayout());
    const LevelMapper level0Mapper(level0_view, Dune::mcmgElementLayout());

    const LookUpData lookUpData(leaf_view);
    const LookUpData lookUpDataLGR(leaf_view, /* isFieldPropInLgr = */ true);

    const LookUpCartesianData lookUpCartesianData(leaf_view, cartMapper);
    const LookUpCartesianData lookUpCartesianDataLGR(leaf_view, cartMapper, /* isFieldPropInLgr = */ true);

    const auto& data = grid.currentData();
    const auto& level0_localIdSet = (*data[0]).localIdSet();

    for (const auto& elem : elements(leaf_view)) {
        // Search via Entity/Element
        const auto featureInElem = lookUpData(elem, fake_feature);
        BOOST_CHECK(featureInElem == lookUpData.getFieldPropIdx(elem) +3);

        const auto featureInElemDouble = lookUpData(elem, fake_feature_double);
        BOOST_CHECK(featureInElemDouble == lookUpData.getFieldPropIdx(elem) +.5);

        const auto featureInElemCartesian = lookUpCartesianData(elem, fake_feature);
        BOOST_CHECK(featureInElemCartesian == lookUpData.getFieldPropIdx(elem) +3);

        const auto featureInElemDoubleCartesian = lookUpCartesianData(elem,fake_feature_double);
        BOOST_CHECK(featureInElemDoubleCartesian == lookUpData.getFieldPropIdx(elem) +.5);

        BOOST_CHECK(elem.getOrigin().index() == lookUpData.getFieldPropIdx(elem));

        // Search via INDEX
        const auto featureInElemIDX = lookUpData(elem.index(), fake_feature);
        BOOST_CHECK(featureInElemIDX == lookUpData.getFieldPropIdx(elem.index()) +3);

        const auto featureInElemDoubleIDX = lookUpData(elem.index(), fake_feature_double);
        BOOST_CHECK(featureInElemDoubleIDX == lookUpData.getFieldPropIdx(elem.index()) +.5);

        const auto featureInElemCartesianIDX = lookUpCartesianData(elem.index(), fake_feature);
        BOOST_CHECK(featureInElemCartesianIDX == lookUpData.getFieldPropIdx(elem.index()) +3);

        const auto featureInElemDoubleCartesianIDX = lookUpCartesianData(elem.index(),fake_feature_double);
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
        levelCartMapp.cartesianCoordinate(elem.getOrigin().index(), ijkLevel0, 0);
        BOOST_CHECK(ijk == ijkLevel0);

        // Throw for level < 0 or level > maxLevel()
        std::array<int,3> ijkThrow;
        BOOST_CHECK_THROW(levelCartMapp.cartesianCoordinate(elem.index(), ijkThrow, -3), std::invalid_argument);
        BOOST_CHECK_THROW(levelCartMapp.cartesianCoordinate(elem.index(), ijkThrow, grid.maxLevel() + 1), std::invalid_argument);

        // Checks related to LGR field properties
        if (elem.level())
        {
            const auto featureInLGR = lookUpDataLGR(elem, fakeLgrFeatures[elem.level()-1]); // shifted since level zero is excluded
            const auto featureInLGR_Cartesian = lookUpCartesianDataLGR(elem, fakeLgrFeatures[elem.level()-1]);
            const auto featureInLGR_FromIdx = lookUpDataLGR(elem.index(), fakeLgrFeatures[elem.level()-1]);
            const auto featureInLGR_Cartesian_FromIdx = lookUpCartesianDataLGR(elem.index(), fakeLgrFeatures[elem.level()-1]);
            BOOST_CHECK(featureInLGR == elem.level());
            BOOST_CHECK(featureInLGR == featureInLGR_Cartesian);
            BOOST_CHECK(featureInLGR == featureInLGR_FromIdx);
            BOOST_CHECK(featureInLGR == featureInLGR_Cartesian_FromIdx);

            // Checks for CartesianCoordinateLevel
            const auto idxOnLevel = elem.getLevelElem().index();
            std::array<int,3> ijkLevelGrid;
            (*data[elem.level()]).getIJK(idxOnLevel, ijkLevelGrid);
            std::array<int,3> ijkLevel;
            levelCartMapp.cartesianCoordinate(idxOnLevel, ijkLevel, elem.level());
            BOOST_CHECK( ijkLevelGrid == ijkLevel);
        }
        // Extra checks related to ElemMapper
        BOOST_CHECK(featureInElem == level0Mapper.index(elem.getOrigin()) +3);
        BOOST_CHECK(featureInElem == fake_feature[lookUpData.getFieldPropIdx(elem)]);
        if (elem.hasFather()) {
            const auto parent_localId = level0_localIdSet.id(elem.father());
            BOOST_CHECK(elem.index() == leafMapper.index(elem));
            BOOST_CHECK(elem.father().index() == featureInElem -3);
            BOOST_CHECK(elem.father().index() == parent_localId);
            BOOST_CHECK(elem.father().index() == level0Mapper.index(elem.father()));
        }
    }
}

BOOST_AUTO_TEST_CASE(grid_without_lgrs)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    lookup_check(grid);
}


BOOST_AUTO_TEST_CASE(grid_with_one_lgr)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    grid.addLgrsUpdateLeafView( /* cells_per_dim_vec = */ {{2,2,2}},
                                /* startIJK_vec = */{{1,0,1}},
                                /* endIJK_vec = */ {{3,2,3}},
                                /* lgr_name_vec = */ {"LGR1"});
    lookup_check(grid);
}

BOOST_AUTO_TEST_CASE(grid_with_multiple_lgrs)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    grid.addLgrsUpdateLeafView( /* cells_per_dim_vec = */ {{2,2,2}, {3,3,3}, {4,4,4}},
                                /* startIJK_vec = */ {{0,0,0}, {0,0,2}, {3,2,2}},
                                /* endIJK_vec = */ {{2,1,1}, {1,1,3}, {4,3,3}},
                                /* lgr_name_vec = */ {"LGR1", "LGR2", "LGR3"});
    lookup_check(grid);
}

void fieldProp_check(const Dune::CpGrid& grid, Opm::EclipseGrid eclGrid, const std::string& deck_string)
{
    Opm::Deck deck = Opm::Parser{}.parseString(deck_string);
    Opm::FieldPropsManager fpm(deck, Opm::Phases{true, true, true}, eclGrid, Opm::TableManager());
    const auto& poro = fpm.get_double("PORO");
    const auto& eqlnum =  fpm.get_int("EQLNUM");
    const auto& porv = fpm.get_double("PORV");

    const auto& leaf_view = grid.leafGridView();
    const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapper(grid);
    const LeafMapper  mapper(leaf_view, Dune::mcmgElementLayout());

    const LookUpData lookUpData(leaf_view);
    const LookUpCartesianData lookUpCartesianData(leaf_view, cartMapper);

    const auto& poroOnLeaf = lookUpData.assignFieldPropsDoubleOnLeaf(fpm, "PORO");
    const auto& poroOnLeafCart = lookUpCartesianData.assignFieldPropsDoubleOnLeaf(fpm, "PORO");

    const auto& eqlnumOnLeaf = lookUpData.assignFieldPropsIntOnLeaf<int>(fpm, "EQLNUM", /* needsTranslation = */ true);
    const auto& eqlnumOnLeafCart = lookUpCartesianData.assignFieldPropsIntOnLeaf<int>(fpm, "EQLNUM",  /* needsTranslation = */ true);

    const auto& porvOnLeaf = lookUpData.assignFieldPropsDoubleOnLeaf(fpm, "PORV");

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

        // PORV
        if (elem.hasFather()) {
            // Pore volume of the father must be equalivalent to the sum of the pore volume of its chidren.
            // Remark: not optimal, repeted computation for children with the same parent cell.
            const auto& parentIdx = elem.father().index();
            const auto& parentPorv = porv[parentIdx];
            // Remark: children_list indices are the indices on the LGR - Not on the leaf grid View.
            const auto& [lgr, children_list] = (*grid.currentData()[0]).getParentToChildren()[parentIdx];
            // Get child indices on the leaf grid view, get their porv value, sum them up, and compare
            // the sum with the pore volume of their parent.
            double sumChildrenPorv = 0.;
            for (const auto& child : children_list) {
                const auto& childIdxOnLeaf = (*grid.currentData()[lgr]).getLeafIdxFromLevelIdx(child);
                sumChildrenPorv += porvOnLeaf[childIdxOnLeaf];
            }
            BOOST_CHECK_CLOSE(parentPorv, sumChildrenPorv, 1e-6);
        }
        else {
            BOOST_CHECK_EQUAL(porv[elemOriginIdx], lookUpData.fieldPropDouble<Dune::cpgrid::Entity<0>>(fpm, "PORV", elem));
            BOOST_CHECK_EQUAL(porv[elemOriginIdx], lookUpData.fieldPropDouble<int>(fpm, "PORV", elemIdx));
            BOOST_CHECK_EQUAL(porv[elemOriginIdx], porvOnLeaf[elemIdx]);
        }
    }
}

BOOST_AUTO_TEST_CASE(poro_porv_eqlnum_grid_withOUT_lgrs) {
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
        PORV
        .5 1.0 1.5 2.0 2.5
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


BOOST_AUTO_TEST_CASE(poro_porv_eqlnum_grid_with_lgrs) {
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
PORV
0.5 1.0 1.5 2. 2.5
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

    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}, {2,2,2}},
                               /* startIJK_vec = */ {{0,0,0}, {0,0,3}},
                               /* endIJK_vec = */ {{1,1,1}, {1,1,4}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2"});

    fieldProp_check(grid, eclGrid, deckString);
}
