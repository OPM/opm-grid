// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>

#include <dune/common/version.hh>

#define BOOST_TEST_MODULE DistributeGridWithLgrsAndWellsTest
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/GraphOfGrid.hpp>
#include <opm/grid/GraphOfGridWrappers.hpp>
#include <opm/grid/common/WellConnections.hpp>
#include <opm/grid/utility/OpmWellType.hpp>

#include <tests/cpgrid/lgr/LgrChecks.hpp>

#include <opm/input/eclipse/EclipseState/Aquifer/NumericalAquifer/NumericalAquifers.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/ScheduleTypes.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQActive.hpp>
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/input/eclipse/Python/Python.hpp>

#include <array>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
// create Wells, we only use well name and cell locations
auto createConnection (int i, int j, int k)
{
    return Opm::Connection(i,j,k, 0, 0, Opm::Connection::State::OPEN,
                           Opm::Connection::Direction::Z,
                           Opm::Connection::CTFKind::DeckValue, 0,
                           5.,Opm::Connection::CTFProperties(),0,false);
}
auto createWell (const std::string& name)
{
    using namespace Opm;
    return Dune::cpgrid::OpmWellType(name,name,0,0,0,0,0.,WellType(),
                                     Well::ProducerCMode(),Connection::Order(),UnitSystem(),
                                     0.,0.,false,false,0,Well::GasInflowEquation());
};
} // end anonymous namespace

BOOST_AUTO_TEST_CASE(add_wells_and_loadBalance_level_zero_of_cartesian_cpgrid_with_lgrs)
{
    // Create a grid with lgrs
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dims = */ {1,2,4}, /* size = */ {1.,1.,1.});

    grid.addLgrsUpdateLeafView( /* cells_per_dim_vec = */ {{2,2,1}},
                                /* startIJK_vec = */ {{0,0,1}},
                                /* endIJK_vec = */ {{1,2,3}},
                                /* lgr_name_vec = */ {"LGR1"});

    // level zero grid view
    //
    // k = 0   1   | k = 1  [3]  | k = 2   [5]   | k = 3   7   |
    //         0   |        [2]  |         [4]   |         6   |

    // Add wells in level zero grid
    auto wellCon = std::make_shared<Opm::WellConnections>();
    wellCon->add(createConnection(0,1,0)); // {level 0, cell idx 1} -> ijk = {0,1,0}
    wellCon->add(createConnection(0,1,1)); // {level 0, cell idx 3} -> ijk = {0,1,1}

    std::vector<Dune::cpgrid::OpmWellType> wells;
    wells.push_back(createWell("well1"));
    wells[0].updateConnections(wellCon,true);

    wellCon = std::make_shared<Opm::WellConnections>(); // reset
    wellCon->add(createConnection(0,0,2)); // {level 0, cell idx 4} -> ijk = {0,0,2}
    wellCon->add(createConnection(0,0,3)); // {level 0, cell idx 6} -> ijk = {0,0,3}

    wells.push_back(createWell("well2"));
    wells[1].updateConnections(wellCon,true);

    wells.push_back(createWell("well3"));

    std::unordered_map<std::string, std::set<int>> futureConnections;
    futureConnections.emplace("well3",std::set<int>{0,2});
    Dune::cpgrid::WellConnections wellConnections(wells, futureConnections, grid);

    if (grid.comm().size()>1) {

        grid.loadBalance(&wells,
                         /* possibleFutureConnections = */ {},
                         /* transmissibilities = */ nullptr,
                         /* overlapLayers = */ 1,
                         /* partitionMethod = */ Dune::PartitionMethod::zoltanGoG,
                         /* level = */ 0);

        grid.addLgrsUpdateLeafView( /* cells_per_dim_vec = */ {{2,2,1}},
                                    /* startIJK_vec = */ {{0,0,1}},
                                    /* endIJK_vec = */ {{1,2,3}},
                                    /* lgr_name_vec = */ {"LGR1"});

        grid.syncDistributedGlobalCellIds();

        const auto& data = grid.currentData();
        Opm::checkCellGlobalIdUniquenessForInteriorCells(grid, data);
        Opm::checkConsecutiveChildGlobalIdsPerParent(grid);

    }
}

BOOST_AUTO_TEST_CASE(extract_wells_from_deck_and_loadBalance_level_zero_of_cpgrid_with_lgrs_created_from_deck){
    const std::string deck_string = R"(
RUNSPEC
DIMENS
   10 10 3 /
START
   1 'JAN' 2015 /
WELLDIMS
   2 1 1 2 /
GRID
CARFIN
-- NAME I1-I2 J1-J2 K1-K2 NX NY NZ
'LGR1'  5  6  5  6  1  3  6  6  9 /
ENDFIN
CARFIN
-- NAME I1-I2 J1-J2 K1-K2 NX NY NZ
'LGR2'  7  8  7  8  1  3  6  6  9 /
ENDFIN
DX
	300*1000 /
DY
	300*1000 /
DZ
	100*20 100*30 100*50 /
TOPS
	100*8325 /
PORO
   	300*0.3 /
PERMX
	100*500 100*50 100*200 /
PERMY
	100*500 100*50 100*200 /
PERMZ
	100*500 100*50 100*200 /
SCHEDULE
WELSPECS
	'PROD'	'G1'	10	10	8400	'OIL' /
	'INJ'	'G1'	1	1	8335	'GAS' /
/
COMPDAT
	'PROD'	10	10	3	3	'OPEN'	1*	1*	0.5 /
	'INJ'	1	1	1	1	'OPEN'	1*	1*	0.5 /
/
END
)";

    Opm::Parser parser;
    const auto deck = parser.parseString(deck_string);
    Opm::EclipseState ecl_state(deck);
    Opm::EclipseGrid ecl_grid = ecl_state.getInputGrid();

    // Create CpGrid from deck and add LGRs
    Dune::CpGrid grid;
    grid.processEclipseFormat(&ecl_grid, &ecl_state, false, false, false);
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{3,3,3}, {3,3,3}},
                               /* startIJK_vec = */ {{4,4,0}, {6,6,0}},
                               /* endIJK_vec = */ {{6,6,3}, {8,8,3}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2"});


    // Extract wells from SCHEDULE
    const Opm::TableManager table ( deck );
    const Opm::FieldPropsManager fp( deck, Opm::Phases{true, true, true}, ecl_grid, table);
    const Opm::Runspec runspec (deck);
    const Opm::Schedule schedule { deck, ecl_grid, fp, Opm::NumericalAquifers{}, runspec, std::make_shared<Opm::Python>()};
    auto wells = schedule.getWellsatEnd();

    if (grid.comm().size()>1) {

        grid.loadBalance(&wells, /* possibleFutureConnections = */ {},
                         /* transmissibilities = */ nullptr,
                         /* overlapLayers = */ 1,
                         /* partitionMethod = */ Dune::PartitionMethod::zoltanGoG,
                         /* level = */ 0);

        grid.addLgrsUpdateLeafView( /* cells_per_dim_vec = */ {{3,3,3}, {3,3,3}},
                                    /* startIJK_vec = */ {{4,4,0}, {6,6,0}},
                                    /* endIJK_vec = */ {{6,6,3}, {8,8,3}},
                                    /* lgr_name_vec = */ {"LGR1", "LGR2"});


        grid.syncDistributedGlobalCellIds();

        const auto& data = grid.currentData();
        Opm::checkCellGlobalIdUniquenessForInteriorCells(grid, data);
        Opm::checkConsecutiveChildGlobalIdsPerParent(grid);
    }

}

BOOST_AUTO_TEST_CASE(add_wells_and_loadBalance_level_zero_of_global_refined_cpgrid)
{
    // Create grid and refine it globally.
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {2.0, 2.0, 2.0});
    grid.globalRefine(1);

    // Add wells in level zero grid
    std::vector<Dune::cpgrid::OpmWellType> wells;
    auto wellCon = std::make_shared<Opm::WellConnections>();

    wellCon->add(createConnection(0,0,0)); // {level 0, cell idx 0}
    wellCon->add(createConnection(0,1,0)); // {level 0, cell idx 1}
    wellCon->add(createConnection(0,1,1)); // {level 0, cell idx 3}
    wells.push_back(createWell("first"));
    wells[0].updateConnections(wellCon,true);

    wellCon = std::make_shared<Opm::WellConnections>(); // reset
    wellCon->add(createConnection(0,0,2)); // {level 0, cell idx 4}
    wellCon->add(createConnection(0,1,2)); // {level 0, cell idx 5}
    wells.push_back(createWell("second"));
    wells[1].updateConnections(wellCon,true);

    wells.push_back(createWell("third"));
    std::unordered_map<std::string, std::set<int>> futureConnections;
    futureConnections.emplace("third", std::set<int>{6,7});
    Dune::cpgrid::WellConnections wellConnections(wells, futureConnections, grid);

    if (grid.comm().size()>1) {

        grid.loadBalance(&wells, futureConnections,
                         /* transmissibilities = */ nullptr,
                         /* overlapLayers = */ 1,
                         /* partitionMethod = */ Dune::PartitionMethod::zoltanGoG,
                         /* level = */ 0);

        grid.globalRefine(1);

        grid.syncDistributedGlobalCellIds();

        const auto& data = grid.currentData();
        Opm::checkCellGlobalIdUniquenessForInteriorCells(grid, data);
        Opm::checkConsecutiveChildGlobalIdsPerParent(grid);

    }
}

bool
init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    boost::unit_test::unit_test_main(&init_unit_test_func,
                                     argc, argv);
}
