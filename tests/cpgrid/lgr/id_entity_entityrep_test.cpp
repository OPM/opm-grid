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

// This fine attempts to verify that the ids and indices returned by Entity and EntityRep are the same
// where they are supposed to be the same, and different where they are supposed to be different.
// This is a test to make sure that modofications to the overloads of id() and index() functions
// in CpGridData and IndexSet does not break the current behaviour.

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

template<typename GridView, typename EntityToIndex, typename EntityRepToIndex>
bool checkEntityIndex(const GridView& gridView, EntityToIndex&& entityToIndex, EntityRepToIndex&& entityRepToIndex, bool coincide)
{
    std::size_t comp_count = 0;
    std::size_t ref_count = 0;
    Dune::Hybrid::forEach(std::make_index_sequence<GridView::dimension + 1>{}, [&](auto codim){
        if constexpr (Dune::Capabilities::hasEntity<typename GridView::Grid, codim>::v)
            for (const Dune::cpgrid::Entity<codim>& entity : entities(gridView, Dune::Codim<codim>())) {
                const auto& entityRep = static_cast<const Dune::cpgrid::EntityRep<codim>&>(entity);
                auto entityIndex = entityToIndex(entity);
                auto entityRepIndex = entityRepToIndex(entityRep);
                comp_count += (entityIndex == entityRepIndex);
                ref_count++;
                if (coincide)
                    BOOST_WARN_EQUAL(entityIndex, entityRepIndex);
                else
                    BOOST_WARN_NE(entityIndex, entityRepIndex);
            }
    });
    BOOST_TEST_MESSAGE("       └─ [" << gridView.comm().rank()
                                     << "]: Compared " + std::to_string(ref_count) + " entities, "
                           + (coincide ? std::to_string(comp_count) + " matched."
                                       : std::to_string(ref_count - comp_count) + " differed."));
    comp_count = gridView.comm().sum(comp_count);
    ref_count = gridView.comm().sum(ref_count);
    bool succed = coincide ? (comp_count == ref_count) : (comp_count < ref_count);

    return succed;
}

std::string coincideString(bool coincide)
{
    return coincide ? "coincide (all ids/indices are equal)" : "differ (at least one id/index is different)";
}

template<typename IdSet, typename EntityRep>
auto idRep(const IdSet& idSet, const EntityRep& entityRep)
{
    return idSet.id(entityRep);
}

void entityRepIndexAndEntityIndexCoincide(const Dune::CpGrid& grid)
{
    if (grid.comm().rank() == 0)
        BOOST_TEST_MESSAGE("    └─ Comparing Entity vs EntityRep indices - lvl leaf     | Expected to " + coincideString(true));
    const auto& leafIndexSet = grid.leafIndexSet();
    const auto leafView = grid.leafGridView();
    BOOST_CHECK(
        checkEntityIndex(leafView,
            [&](const auto& e){ return leafIndexSet.index(e); },
            [&](const auto& er){ return leafIndexSet.index(er); },
            true /* mustCoincide */
        ) || grid.comm().rank() != 0
    );

    for (int level = 0; level <= grid.maxLevel(); ++level) {
        if (grid.comm().rank() == 0)
            BOOST_TEST_MESSAGE("    └─ Comparing Entity vs EntityRep indices - lvl " + std::to_string(level) + "        | Expected to " + coincideString(true));
        const auto& levelIndexSet = grid.levelIndexSet(level);
        const auto levelView = grid.levelGridView(level);
        BOOST_CHECK(
            checkEntityIndex(levelView,
                [&](const auto& e){ return levelIndexSet.index(e); },
                [&](const auto& er){ return levelIndexSet.index(er); },
                true /* mustCoincide */
            ) || grid.comm().rank() != 0
        );
    }
}

void globalIdsEntityRepAndEntityInLevelZeroGrid(const Dune::CpGrid& grid, bool coincide)
{
    if (grid.comm().rank() == 0)
        BOOST_TEST_MESSAGE("    └─ Comparing Entity vs EntityRep global ids - lvl 0     | Expected to " + coincideString(coincide));
    const auto levelZeroView = grid.levelGridView(0);
    const auto& levelZeroGlobalIdSet = grid.currentData()[0]->globalIdSet();
    const auto& globalIdSet = levelZeroGlobalIdSet; // TODO -> grid.globalIdSet();
    BOOST_CHECK(
        checkEntityIndex(levelZeroView,
            [&](const auto& e){ return globalIdSet.template id<std::decay_t<decltype(e)>::codimension>(e); },
            [&](const auto& er){ return idRep(levelZeroGlobalIdSet, er); },
            coincide /* mustCoincide */
        ) || grid.comm().rank() != 0
    );
}

void localIdsEntityRepAndEntityInLevelZeroGrid(const Dune::CpGrid& grid, bool coincide)
{
    if (grid.comm().rank() == 0)
        BOOST_TEST_MESSAGE("    └─ Comparing Entity vs EntityRep local ids - lvl 0      | Expected to " + coincideString(coincide));
    const auto levelZeroView = grid.levelGridView(0);
    const auto& levelZeroLocalIdSet = grid.currentData()[0]->localIdSet();
    const auto& localIdSet = levelZeroLocalIdSet; // TODO -> grid.localIdSet();
    BOOST_CHECK(
        checkEntityIndex(levelZeroView,
            [&](const auto& e){ return localIdSet.id(e); },
            [&](const auto& er){ return idRep(levelZeroLocalIdSet, er); },
            coincide /* mustCoincide */
        ) || grid.comm().rank() != 0
    );
}

void localIdsEntityRepAndEntityInRefinedLevelGrids(const Dune::CpGrid& grid, bool coincide)
{
    for (int level = 1; level <= grid.maxLevel(); ++level) {
        if (grid.comm().rank() == 0)
            BOOST_TEST_MESSAGE("    └─ Comparing Entity vs EntityRep local ids - lvl " + std::to_string(level) + "      | Expected to " + coincideString(coincide));
        const auto levelView = grid.levelGridView(level);
        const auto& levelLocalIdSet = grid.currentData()[level]->localIdSet();
        const auto& localIdSet = levelLocalIdSet; // TODO -> grid.localIdSet();
        BOOST_CHECK(
            checkEntityIndex(levelView,
                [&](const auto& e){ return localIdSet.id(e); },
                [&](const auto& er){ return idRep(levelLocalIdSet, er); },
                coincide /* mustCoincide */
            ) || grid.comm().rank() != 0
        );
    }
}

void localIdsEntityRepAndEntityInLeafGrid(const Dune::CpGrid& grid, bool coincide)
{
    if (grid.comm().rank() == 0)
        BOOST_TEST_MESSAGE("    └─ Comparing Entity vs EntityRep local ids - lvl leaf   | Expected to " + coincideString(coincide));
    const auto leafView = grid.leafGridView();
    const auto& leafLocalIdSet = grid.currentData().back()->localIdSet();
    const auto& localIdSet = leafLocalIdSet; // TODO -> grid.localIdSet();
    BOOST_CHECK(
        checkEntityIndex(leafView,
            [&](const auto& e){ return localIdSet.id(e); },
            [&](const auto& er){ return idRep(leafLocalIdSet, er); },
            coincide /* mustCoincide */
        ) || grid.comm().rank() != 0
    );
}

void globalIdsEntityRepAndEntityInRefinedLevelGrids(const Dune::CpGrid& grid, bool coincide)
{
    for (int level = 1; level <= grid.maxLevel(); ++level) {
        if (grid.comm().rank() == 0)
            BOOST_TEST_MESSAGE("    └─ Comparing Entity vs EntityRep global ids - lvl " + std::to_string(level) + "     | Expected to " + coincideString(coincide));
        const auto levelView = grid.levelGridView(level);
        const auto& levelGlobalIdSet = grid.currentData()[level]->globalIdSet();
        const auto& globalIdSet = grid.currentData()[level]->globalIdSet(); // TODO: -> grid.globalIdSet();
        BOOST_CHECK(
            checkEntityIndex(levelView,
                [&](const auto& e){ return globalIdSet.id(e); },
                [&](const auto& er){ return idRep(levelGlobalIdSet, er); },
                coincide /* mustCoincide */
            ) || grid.comm().rank() != 0
        );
    }
}

void globalIdsEntityRepAndEntityInLeafGrid(const Dune::CpGrid& grid, bool coincide)
{
    if (grid.comm().rank() == 0)
        BOOST_TEST_MESSAGE("    └─ Comparing Entity vs EntityRep global ids - lvl leaf  | Expected to " + coincideString(coincide));
    const auto leafView = grid.leafGridView();
    const auto& leafGlobalIdSet = grid.currentData().back()->globalIdSet();
    const auto& globalIdSet = leafGlobalIdSet; // TODO -> grid.globalIdSet();
    BOOST_CHECK(
        checkEntityIndex(leafView,
            [&](const auto& e){ return globalIdSet.id(e); },
            [&](const auto& er){ return idRep(leafGlobalIdSet, er); },
            coincide /* mustCoincide */
        ) || grid.comm().rank() != 0
    );
}

void checkEntityIds(const Dune::CpGrid& grid)
{
    bool grid_is_distributed = grid.size(0) != grid.comm().sum(grid.size(0));
    bool grid_is_refined = grid.maxLevel() > 0;
    std::string grid_type = grid_is_distributed ? "distributed" : "serial";
    grid_type += grid_is_refined ? " refined" : " coarse";
    if (grid.comm().rank() == 0)
        BOOST_TEST_MESSAGE("  └─ Checking ids of Entity and EntityRep in " + grid_type + " grid");

    // All test below compare ids between Entity and EntityRep.
    // EntityRep id represent the level zero grid ids, but their indices always coincide.

    // Indices ALWAYS coincide
    entityRepIndexAndEntityIndexCoincide(grid);

    // local ids on level zero grid view ALWAYS coincide
    localIdsEntityRepAndEntityInLevelZeroGrid(grid, true);

    // local ids on level zero grid view NEVER coincide
    localIdsEntityRepAndEntityInRefinedLevelGrids(grid, false);

    // TODO: this test currently fails
    // local ids on leaf grid view coincide ONLY if the grid is NOT refined
    localIdsEntityRepAndEntityInLeafGrid(grid, !grid_is_refined);

    // global ids on level zero grid view ALWAYS coincide
    globalIdsEntityRepAndEntityInLevelZeroGrid(grid, true);

    // global ids on refined level grids view NEVER coincide
    if (grid.comm().size() == 0) // TODO: this test currently fails in parallel
        globalIdsEntityRepAndEntityInRefinedLevelGrids(grid, false);

    // global ids on leaf grid view coincide ONLY if the grid is NOT distributed and NOT refined
    if (grid.comm().size() == 0) // TODO: this test currently fails in parallel
        globalIdsEntityRepAndEntityInLeafGrid(grid, !grid_is_distributed and !grid_is_refined);
}

BOOST_AUTO_TEST_CASE(idChecks_viaAddLgrs)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    if (grid.comm().rank() == 0)
        BOOST_TEST_MESSAGE("└─ Checking grid refinement via addLgrsUpdateLeafView");

    checkEntityIds(grid);

    if (grid.comm().size()>1) {
        grid.loadBalance();
        checkEntityIds(grid);
    }

    grid.addLgrsUpdateLeafView( /* cells_per_dim = */ {{2,2,2}},
                                /* startIJK = */ {{1,1,1}},
                                /* endIJK = */  {{3,2,2}}, // block cell indices = {17, 18}
                                /* lgr_name = */  {"LGR1"});

    if (grid.comm().size() == 1) { // Serial
        checkEntityIds(grid);
    }
    else { // Parallel
        // BEFORE applying the same refinement to the global view.
        checkEntityIds(grid);

        grid.switchToGlobalView();
        // Synchronizing cell ids requires that both the global
        // and distributed views undergo the same refinement process.
        grid.addLgrsUpdateLeafView( /* cells_per_dim = */ {{2,2,2}},
                                    /* startIJK = */ {{1,1,1}},
                                    /* endIJK = */  {{3,2,2}}, // block cell global ids = {17, 18}
                                    /* lgr_name = */  {"LGR1"});
        grid.switchToDistributedView();

        grid.syncDistributedGlobalCellIds();

        // AFTER synchronizing cell ids
        checkEntityIds(grid);
    }
}

BOOST_AUTO_TEST_CASE(idChecks_viaAdapt)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    if (grid.comm().rank() == 0)
        BOOST_TEST_MESSAGE("└─ Checking grid refinement via adapt");

    checkEntityIds(grid);

    if (grid.comm().size()>1) {
        grid.loadBalance();
        checkEntityIds(grid);
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
        checkEntityIds(grid);
    }
    else { // Parallel
        // BEFORE applying the same refinement to the global view
        checkEntityIds(grid);

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

        grid.syncDistributedGlobalCellIds();

        // AFTER synchronozing cell ids
        checkEntityIds(grid);
    }
}

BOOST_AUTO_TEST_CASE(idChecks_viaGlobalRefine)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    if (grid.comm().rank() == 0)
        BOOST_TEST_MESSAGE("└─ Checking grid refinement via globalRefine");

    checkEntityIds(grid);

    if (grid.comm().size()>1) {
        grid.loadBalance();
        checkEntityIds(grid);
    }

    grid.globalRefine(1);

    if (grid.comm().size() == 1) { // Serial
        checkEntityIds(grid);
    }
    else { // Parallel
        // BEFORE applying the same refinement to the global view.
        checkEntityIds(grid);

        grid.switchToGlobalView();
        // Synchronizing cell ids requires that both the global
        // and distributed views undergo the same refinement process.
        grid.globalRefine(1);
        grid.switchToDistributedView();

        grid.syncDistributedGlobalCellIds();

        // AFTER synchronizing cell ids
        checkEntityIds(grid);
    }
}
