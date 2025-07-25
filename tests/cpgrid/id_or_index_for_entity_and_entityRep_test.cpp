//===========================================================================
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

#define BOOST_TEST_MODULE TemplateOverloadsForIdAndIndexTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>

#include <unordered_set>

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

template<int codim, typename EntityRange, typename IndexSet>
void checkEntityRepAndEntityIndices(const EntityRange& entities, const IndexSet& indexSet)
{
    for (const auto& entity : entities) {
        const auto entityRep = Dune::cpgrid::EntityRep<codim>(entity.index(), true);
        BOOST_CHECK_EQUAL(indexSet.index(entityRep), indexSet.index(entity));
    }
}

template<int codim, typename EntityRange, typename IdSet>
bool checkEntityRepAndEntityIds(const EntityRange& entities, const IdSet& idSet)
{
    for (const auto& entity : entities) {
        const auto entityRep = Dune::cpgrid::EntityRep<codim>(entity.index(), true);
        if (idSet.id(entityRep) != idSet.id(entity)){
            return false;
        }
    }
    return true;
}

void entityRepIndexAndEntityIndexCoincide(const Dune::CpGrid& grid)
{
    const auto& leafIndexSet = grid.leafIndexSet();
    const auto leafView = grid.leafGridView();
    checkEntityRepAndEntityIndices<0>(Dune::elements(leafView), leafIndexSet);
    checkEntityRepAndEntityIndices<3>(Dune::vertices(leafView), leafIndexSet);

    for (int level = 0; level <= grid.maxLevel(); ++level) {
        const auto& levelIndexSet = grid.levelIndexSet(level);
        const auto levelView = grid.levelGridView(level);
        checkEntityRepAndEntityIndices<0>(Dune::elements(levelView), levelIndexSet);
        checkEntityRepAndEntityIndices<3>(Dune::vertices(levelView), levelIndexSet);
    }
}

void globalIdsEntityRepAndEntityInLevelZeroGrid(const Dune::CpGrid& grid, bool coincide)
{
    const auto levelZeroView = grid.levelGridView(0);
    const auto& levelZeroGlobalIdSet = grid.currentData()[0]->globalIdSet();
    BOOST_CHECK_EQUAL( checkEntityRepAndEntityIds<0>(Dune::elements(levelZeroView), levelZeroGlobalIdSet), coincide);
    BOOST_CHECK_EQUAL( checkEntityRepAndEntityIds<3>(Dune::vertices(levelZeroView), levelZeroGlobalIdSet), coincide);
}

void localIdsEntityRepAndEntityInLevelZeroGrid(const Dune::CpGrid& grid, bool coincide)
{
    const auto levelZeroView = grid.levelGridView(0);
    const auto& levelZeroLocalIdSet = grid.currentData()[0]->localIdSet();
    BOOST_CHECK_EQUAL( checkEntityRepAndEntityIds<0>(Dune::elements(levelZeroView), levelZeroLocalIdSet), coincide);
    BOOST_CHECK_EQUAL( checkEntityRepAndEntityIds<3>(Dune::vertices(levelZeroView), levelZeroLocalIdSet), coincide);
}

void localIdsEntityRepAndEntityInRefinedLevelGrids(const Dune::CpGrid& grid, bool coincide)
{
    for (int level = 1; level <= grid.maxLevel(); ++level) {
        const auto levelView = grid.levelGridView(level);
        if (levelView.size(0)>0) {
            const auto& levelLocalIdSet = grid.currentData()[level]->localIdSet();
            BOOST_CHECK_EQUAL(checkEntityRepAndEntityIds<0>(Dune::elements(levelView), levelLocalIdSet), coincide);
            BOOST_CHECK_EQUAL(checkEntityRepAndEntityIds<3>(Dune::vertices(levelView), levelLocalIdSet), coincide);
        }
    }
}

void localIdsEntityRepAndEntityInLeafGrid(const Dune::CpGrid& grid, bool coincide)
{
    const auto leafView = grid.leafGridView();
    const auto& leafLocalIdSet = grid.currentData().back()->localIdSet();
    BOOST_CHECK_EQUAL(checkEntityRepAndEntityIds<0>(Dune::elements(leafView), leafLocalIdSet), coincide);
    BOOST_CHECK_EQUAL(checkEntityRepAndEntityIds<3>(Dune::vertices(leafView), leafLocalIdSet), coincide);
}

void globalIdsEntityRepAndEntityInRefinedLevelGrids(const Dune::CpGrid& grid, bool coincide)
{
    for (int level = 1; level <= grid.maxLevel(); ++level) {
        const auto levelView = grid.levelGridView(level);
        if (levelView.size(0)>0) {
            const auto& levelGlobalIdSet = grid.currentData()[level]->globalIdSet();
            BOOST_CHECK_EQUAL(checkEntityRepAndEntityIds<0>(Dune::elements(levelView), levelGlobalIdSet), coincide);
            BOOST_CHECK_EQUAL(checkEntityRepAndEntityIds<3>(Dune::vertices(levelView), levelGlobalIdSet), coincide);
        }
    }
}

void globalIdsEntityRepAndEntityInLeafGrid(const Dune::CpGrid& grid, bool coincide)
{
    const auto leafView = grid.leafGridView();
    if (leafView.size(0)>0) {
        const auto& leafGlobalIdSet = grid.currentData().back()->globalIdSet();
        BOOST_CHECK_EQUAL(checkEntityRepAndEntityIds<0>(Dune::elements(leafView), leafGlobalIdSet), coincide);
        BOOST_CHECK_EQUAL(checkEntityRepAndEntityIds<3>(Dune::vertices(leafView), leafGlobalIdSet), coincide);
    }
}

void serialBeforeRefinementChecks(const Dune::CpGrid& grid)
{
    /**  BEFORE refinement via addLgrsUpdateLeafView, adapt, or globalRefine(1) in serial */
    // All indices and local/global ids for Entity<codim> and EntityRep<codim> codim = 0 and 3 coincide.
    entityRepIndexAndEntityIndexCoincide(grid);

    localIdsEntityRepAndEntityInLevelZeroGrid(grid, true /* coincide */);
    localIdsEntityRepAndEntityInRefinedLevelGrids(grid, true /* coincide */);
    localIdsEntityRepAndEntityInLeafGrid(grid, true /* coincide */);

    globalIdsEntityRepAndEntityInLevelZeroGrid(grid, true /* coincide */);
    globalIdsEntityRepAndEntityInRefinedLevelGrids(grid, true /* coincide */);
    globalIdsEntityRepAndEntityInLeafGrid(grid, true /* coincide */);
}

void parallelBeforeRefinementChecks(const Dune::CpGrid& grid)
{
    /**  BEFORE refinement via addLgrsUpdateLeafView, adapt, or globalRefine(1) in parallel */
    // All indices and local ids for Entity<codim> and EntityRep<codim> codim = 0 and 3 coincide.
    entityRepIndexAndEntityIndexCoincide(grid);

    localIdsEntityRepAndEntityInLevelZeroGrid(grid, true /* coincide */);
    localIdsEntityRepAndEntityInRefinedLevelGrids(grid, true /* coincide */);
    localIdsEntityRepAndEntityInLeafGrid(grid, true /* coincide */);

    globalIdsEntityRepAndEntityInRefinedLevelGrids(grid, true /* coincide */);
    /** In constrast to serial, in level zero grid id(Entity) != id(EntityRep) */
    globalIdsEntityRepAndEntityInLevelZeroGrid(grid, false /* coincide */);
    /** In contrast to the serial, in leaf grid (here, level zero) id(Entity) != id(EntityRep) */
    globalIdsEntityRepAndEntityInLeafGrid(grid, false /* coincide */);
}


void serialAfterRefinementChecks(const Dune::CpGrid& grid)
{
    /**  After refinement via addLgrsUpdateLeafView, adapt, or globalRefine(1) in serial */
    entityRepIndexAndEntityIndexCoincide(grid);

    localIdsEntityRepAndEntityInLevelZeroGrid(grid, true /* coincide */);
    localIdsEntityRepAndEntityInRefinedLevelGrids(grid, false /* DO NOT coincide */);
    localIdsEntityRepAndEntityInLeafGrid(grid, false /* DO NOT coincide */);

    globalIdsEntityRepAndEntityInLevelZeroGrid(grid, true /* coincide */);
    globalIdsEntityRepAndEntityInRefinedLevelGrids(grid, false /* DO NOT coincide */);
    globalIdsEntityRepAndEntityInLeafGrid(grid, false /* DO NOT coincide */);
}

void parallelAfterRefinementChecks(const Dune::CpGrid& grid)
{
    /**  After refinement via addLgrsUpdateLeafView, adapt, or globalRefine(1) in parallel */
    entityRepIndexAndEntityIndexCoincide(grid);

    localIdsEntityRepAndEntityInLevelZeroGrid(grid, true /* coincide */);
    localIdsEntityRepAndEntityInRefinedLevelGrids(grid, false /* DO NOT coincide */);
    localIdsEntityRepAndEntityInLeafGrid(grid, false /* DO NOT coincide */);

    globalIdsEntityRepAndEntityInLevelZeroGrid(grid, false /* DO NOT coincide */);
    globalIdsEntityRepAndEntityInRefinedLevelGrids(grid, false /* DO NOT coincide */);
    globalIdsEntityRepAndEntityInLeafGrid(grid, false /* DO NOT coincide */);
}

BOOST_AUTO_TEST_CASE(idEntityRep_and_idEntity_differ_in_refinedLevelAndLeafGrids_viaAddLgrs)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    serialBeforeRefinementChecks(grid);

    if (grid.comm().size()>1) {
        grid.loadBalance();
        parallelBeforeRefinementChecks(grid);
    }

    grid.addLgrsUpdateLeafView( /* cells_per_dim = */ {{2,2,2}},
                                /* startIJK = */ {{1,1,1}},
                                /* endIJK = */  {{3,2,2}}, // block cell indices = {17, 18}
                                /* lgr_name = */  {"LGR1"});

    if (grid.comm().size() == 1) { // Serial
        serialAfterRefinementChecks(grid);
    }
    else { // Parallel
        // BEFORE applying the same refinement to the global view.
        parallelAfterRefinementChecks(grid);

        grid.switchToGlobalView();
        // Synchronizing cell ids requires that both the global
        // and distributed views undergo the same refinement process.
        grid.addLgrsUpdateLeafView( /* cells_per_dim = */ {{2,2,2}},
                                    /* startIJK = */ {{1,1,1}},
                                    /* endIJK = */  {{3,2,2}}, // block cell global ids = {17, 18}
                                    /* lgr_name = */  {"LGR1"});
        grid.switchToDistributedView();

        // BEFORE synchronizing cell ids
        parallelAfterRefinementChecks(grid);

        grid.syncDistributedGlobalCellIds();

        // AFTER synchronizing cell ids
        parallelAfterRefinementChecks(grid);
    }
}

BOOST_AUTO_TEST_CASE(idEntityRep_and_idEntity_differ_in_refinedLevelAndLeafGrids_viaAdapt)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    serialBeforeRefinementChecks(grid);

    if (grid.comm().size()>1) {
        grid.loadBalance();
        parallelBeforeRefinementChecks(grid);
    }

    std::unordered_set<int> markedCells = {17,18}; // parent cell global ids
    // Mark selected elements for refinement. In this moment, level zero and leaf grids coincide.
    for (const auto& element : elements(grid.leafGridView())) {
        const auto& id = grid.globalIdSet().id(element);
        if (markedCells.count(id) > 0) {
            grid.mark(1, element);
        }
    }
    grid.preAdapt();
    grid.adapt(); // Default subdivisions per cell 2x2x2 in x-,y-, and z-direction.
    grid.postAdapt();

    if (grid.comm().size() == 1) { // Serial
        serialAfterRefinementChecks(grid);
    }
    else { // Parallel
        // BEFORE applying the same refinement to the global view
        parallelAfterRefinementChecks(grid);

        grid.switchToGlobalView();
        // Synchronizing cell ids requires that both the global
        // and distributed views undergo the same refinement process.
        for (const auto& element : elements(grid.levelGridView(0))) {
            const auto& id = grid.globalIdSet().id(element);
            if (markedCells.count(id) > 0) {
                grid.mark(1, element);
            }
        }
        grid.preAdapt();
        grid.adapt(); // Default subdivisions per cell 2x2x2 in x-,y-, and z-direction.
        grid.postAdapt();
        grid.switchToDistributedView();

        // BEFORE synchronizing cell ids
        parallelAfterRefinementChecks(grid);

        //grid.syncDistributedGlobalCellIds();

        // AFTER synchronozing cell ids
        parallelAfterRefinementChecks(grid);
    }
}

BOOST_AUTO_TEST_CASE(idEntityRep_and_idEntity_differ_in_refinedLevelAndLeafGrids_viaGlobalRefine)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});
    serialBeforeRefinementChecks(grid);

    if (grid.comm().size()>1) {
        grid.loadBalance();
        parallelBeforeRefinementChecks(grid);
    }

    grid.globalRefine(1);

    if (grid.comm().size() == 1) { // Serial
        serialAfterRefinementChecks(grid);
    }
    else { // Parallel
        // BEFORE applying the same refinement to the global view.
        parallelAfterRefinementChecks(grid);

        grid.switchToGlobalView();
        // Synchronizing cell ids requires that both the global
        // and distributed views undergo the same refinement process.
        grid.globalRefine(1);
        grid.switchToDistributedView();

        // BEFORE synchronizing cell ids
        parallelAfterRefinementChecks(grid);

        grid.syncDistributedGlobalCellIds();

        // AFTER synchronizing cell ids
        parallelAfterRefinementChecks(grid);
    }
}
