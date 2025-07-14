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


BOOST_AUTO_TEST_CASE(idCodim_and_idEntityRep_differ_in_refinedLevelAndLeafGrids)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    grid.addLgrsUpdateLeafView( /* cells_per_dim = */ {{2,2,2}},
                                /* startIJK = */ {{1,1,1}},
                                /* endIJK = */  {{3,2,2}}, // block cell indices = {17, 18}
                                /* lgr_name = */  {"LGR1"});

    const auto& leafIndexSet = grid.leafIndexSet();
    const auto& leafGlobalIdSet = grid.currentData().back()->globalIdSet();

    for (const auto& element : Dune::elements(grid.leafGridView())) {
        // This should be replaced by finding in "leaf_data" ->geomVector<0>()[cpgrid::EntityRep<0>(element.index(), true)]?
        const auto entityRep = Dune::cpgrid::EntityRep<0>(element.index(), true);

        // Consistency between template index<EntityType>(EntityType) with EntityType = EntityRep<0> and template index<codim>(Entity<0>)
        BOOST_CHECK_EQUAL( leafIndexSet.template index(entityRep), leafIndexSet.template index<0>(element));

        // Consistency between template id<EntityType>(EntityType) with EntityType = EntityRep<0> and template id<codim>(Entity<0>)
        // BOOST_CHECK_EQUAL( leafGlobalIdSet.template id(entityRep), leafGlobalIdSet.template id<0>(element)); /** FAILS!!! */

        // Extra check: calling id from the leaf grid view or the grid is equivalent.
        // Note: grid.template id<EntityRep<0>> does not exist.
        BOOST_CHECK_EQUAL( grid.globalIdSet().template id<0>(element), leafGlobalIdSet.template id<0>(element));


        // For intersections, only the overload with EntityRep<codim> is implemented (there is no index<1> nor id<1>; codim=1)
        for (const auto& intersection : Dune::intersections(grid.leafGridView(), element)) {

            const auto intersectionRep = Dune::cpgrid::EntityRep<1>(intersection.id(), true);

            leafIndexSet.template index(intersectionRep);
            // leafIndexSet.template index(intersection); this does not compile

            leafGlobalIdSet.template id(intersectionRep);
            //leafGlobalIdSet.template id(intersection); this does not compile
        }

    }

    for (const auto& vertex : Dune::vertices(grid.leafGridView())) {
        // This should be replaced by finding in "leaf_data" ->geomVector<3>()[cpgrid::EntityRep<3>(vertex.index(), true)]?
        const auto entityRep = Dune::cpgrid::EntityRep<3>(vertex.index(), true);
        // Consistency between template index<EntityType>(EntityType) with EntityType = EntityRep<3> and template index<codim>(Entity<3>)
        BOOST_CHECK_EQUAL( leafIndexSet.template index(entityRep), leafIndexSet.template index<3>(vertex));

        // Consistency between template id<EntityType>(EntityType) with EntityType = EntityRep<0> and template id<codim>(Entity<3>)
        // BOOST_CHECK_EQUAL( leafGlobalIdSet.template id(entityRep), leafGlobalIdSet.template id<3>(vertex)); /** FAILS!!! */

        // Extra check: calling id from the leaf grid view or the grid is equivalent.
        // Note: grid.template id<EntityRep<3>> does not exist.
        BOOST_CHECK_EQUAL( grid.globalIdSet().template id<3>(vertex), leafGlobalIdSet.template id<3>(vertex));
    }


    for (int level = 1; level <= grid.maxLevel(); ++level) {

        const auto& levelIndexSet = grid.levelIndexSet(level);
        const auto& levelGlobalIdSet = grid.currentData()[level]->globalIdSet();

        for (const auto& element : Dune::elements(grid.levelGridView(level))) {
            // This should be replaced by finding in "level_data" ->geomVector<0>()[cpgrid::EntityRep<0>(element.index(), true)]?
            const auto entityRep = Dune::cpgrid::EntityRep<0>(element.index(), true);

            // Consistency between template index<EntityType>(EntityType) with EntityType = EntityRep<0> and template index<codim>(Entity<0>)
            BOOST_CHECK_EQUAL( levelIndexSet.template index(entityRep), levelIndexSet.template index<0>(element));

            // Consistency between template id<EntityType>(EntityType) with EntityType = EntityRep<0> and template id<codim>(Entity<0>)
            // BOOST_CHECK_EQUAL( levelGlobalIdSet.template id(entityRep), levelGlobalIdSet.template id<0>(element)); /** FAILS!!! */

            // Extra check: calling id from the level grid or the grid is equivalent.
            // Note: grid.template id<EntityRep<0>> does not exist.
            BOOST_CHECK_EQUAL( grid.globalIdSet().template id<0>(element), levelGlobalIdSet.template id<0>(element));

            // For intersections, only the overload with EntityRep<codim> is implemented (there is no index<1> nor id<1>; codim=1)
            for (const auto& intersection : Dune::intersections(grid.levelGridView(level), element)) {

                const auto intersectionRep = Dune::cpgrid::EntityRep<1>(intersection.id(), true);

                levelIndexSet.template index(intersectionRep);
                // levelIndexSet.template index(intersection); this does not compile

                levelGlobalIdSet.template id(intersectionRep);
                //levelGlobalIdSet.template id(intersection); this does not compile
            }
        }

        for (const auto& vertex : Dune::vertices(grid.levelGridView(level))) {
            // This should be replaced by finding in "level_data" ->geomVector<3>()[cpgrid::EntityRep<3>(vertex.index(), true)]?
            const auto entityRep = Dune::cpgrid::EntityRep<3>(vertex.index(), true);

            // Consistency between template index<EntityType>(EntityType) with EntityType = EntityRep<3> and template index<codim>(Entity<3>)
            BOOST_CHECK_EQUAL( levelIndexSet.template index(entityRep), levelIndexSet.template index<3>(vertex));

            // Consistency between template id<EntityType>(EntityType) with EntityType = EntityRep<0> and template id<codim>(Entity<3>)
            // BOOST_CHECK_EQUAL( levelGlobalIdSet.template id(entityRep), levelGlobalIdSet.template id<3>(vertex)); /** FAILS!!! */

            // Extra check: calling id from the level grid  or the grid is equivalent.
            // Note: grid.template id<EntityRep<3>> does not exist.
            BOOST_CHECK_EQUAL( grid.globalIdSet().template id<3>(vertex), levelGlobalIdSet.template id<3>(vertex));
        }
    }

    // For level zero, there is no difference
    const auto& levelZeroIndexSet = grid.levelIndexSet(0);
    const auto& levelZeroGlobalIdSet = grid.currentData()[0]->globalIdSet();

    for (const auto& element : Dune::elements(grid.levelGridView(0))) {
        // This should be replaced by finding in "level_data" ->geomVector<0>()[cpgrid::EntityRep<0>(element.index(), true)]?
        const auto entityRep = Dune::cpgrid::EntityRep<0>(element.index(), true);

        // Consistency between template index<EntityType>(EntityType) with EntityType = EntityRep<0> and template index<codim>(Entity<0>)
        BOOST_CHECK_EQUAL( levelZeroIndexSet.template index(entityRep), levelZeroIndexSet.template index<0>(element));

        // Consistency between template id<EntityType>(EntityType) with EntityType = EntityRep<0> and template id<codim>(Entity<0>)
        BOOST_CHECK_EQUAL( levelZeroGlobalIdSet.template id(entityRep), levelZeroGlobalIdSet.template id<0>(element));

        // Extra check: calling id from the level zero grid or the grid is equivalent.
        // Note: grid.template id<EntityRep<0>> does not exist.
        BOOST_CHECK_EQUAL( grid.globalIdSet().template id<0>(element), levelZeroGlobalIdSet.template id<0>(element));

        // For intersections, only the overload with EntityRep<codim> is implemented,
        // there is no index<1> nor id<1>; codim=1, or index<Dune::cpgrid::Intersection> or id<Dune::cpgrid::Intersection>
        for (const auto& intersection : Dune::intersections(grid.levelGridView(0), element)) {

            const auto intersectionRep = Dune::cpgrid::EntityRep<1>(intersection.id(), true);

            levelZeroIndexSet.template index(intersectionRep);
            // levelZeroIndexSet.template index(intersection); this does not compile

            levelZeroGlobalIdSet.template id(intersectionRep);
            //levelZeroGlobalIdSet.template id(intersection); this does not compile
        }
    }

    for (const auto& vertex : Dune::vertices(grid.levelGridView(0))) {
        // This should be replaced by finding in "level_data" ->geomVector<3>()[cpgrid::EntityRep<3>(vertex.index(), true)]?
        const auto entityRep = Dune::cpgrid::EntityRep<3>(vertex.index(), true);

        // Consistency between template index<EntityType>(EntityType) with EntityType = EntityRep<3> and template index<codim>(Entity<3>)
        BOOST_CHECK_EQUAL( levelZeroIndexSet.template index(entityRep), levelZeroIndexSet.template index<3>(vertex));

        // Consistency between template id<EntityType>(EntityType) with EntityType = EntityRep<0> and template id<codim>(Entity<3>)
        BOOST_CHECK_EQUAL( levelZeroGlobalIdSet.template id(entityRep), levelZeroGlobalIdSet.template id<3>(vertex));

        // Extra check: calling id from the level zero grid  or the grid is equivalent.
        // Note: grid.template id<EntityRep<3>> does not exist.
        BOOST_CHECK_EQUAL( grid.globalIdSet().template id<3>(vertex), levelZeroGlobalIdSet.template id<3>(vertex));
    }
}
