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
#include <tests/cpgrid/LgrChecks.hpp>
#include <opm/grid/common/WellConnections.hpp>

#include <opm/grid/utility/OpmWellType.hpp>
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

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
    // create a grid with lgrs and wells
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

BOOST_AUTO_TEST_CASE(add_wells_and_loadBalance_level_zero_of_cpgrid_with_lgrs_created_from_deck){
    const std::string deck_string = R"(
RUNSPEC
DIMENS
  3 3 4 /
GRID
CARFIN
-- NAME  I1-I2 J1-J2 K1-K2 NX NY NZ
'LGR1'   2  3  2  3  1  2  4  4  4/
ENDFIN
DX
  36*1000 /
DY
	36*1000 /
DZ
	36*20 /
TOPS
	36*8325 /
 ACTNUM
1 1 1
1 1 1
1 1 1

1 1 1
1 1 1
1 1 1

1 1 1
1 1 1
1 1 1

1 1 1
1 1 1
1 1 1
 /
PORO
  36*0.15 /
PERMX
  36*1 /
COPY
  PERMX PERMZ /
  PERMX PERMY /
/
EDIT
OIL
GAS
TITLE
The title
START
16 JUN 1988 /
PROPS
REGIONS
SOLUTION
SCHEDULE
)";

    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deck_string,
                              /* cells_per_dim_vec = */ {{2,2,2}},
                              /* startIJK_vec = */ {{1,1,0}},
                              /* endIJK_vec = */ {{3,3,2}},
                              /* lgr_name_vec = */ {"LGR1"});

    auto wellCon = std::make_shared<Opm::WellConnections>();
    wellCon->add(createConnection(0,1,1)); // {level 0, cell idx 12} -> ijk = {0,1,1}
    wellCon->add(createConnection(0,1,2)); // {level 0, cell idx 21} -> ijk = {0,1,2}
    wellCon->add(createConnection(0,1,3)); // {level 0, cell idx 30} -> ijk = {0,1,3}

    std::vector<Dune::cpgrid::OpmWellType> wells;
    wells.push_back(createWell("first"));
    wells[0].updateConnections(wellCon,true);

    wellCon = std::make_shared<Opm::WellConnections>(); // reset
    wellCon->add(createConnection(2,1,1)); // {level 0, cell idx 14} -> ijk = {2,1,1}
    wellCon->add(createConnection(2,1,2)); // {level 0, cell idx 23} -> ijk = {2,1,2}
    wellCon->add(createConnection(2,1,3)); // {level 0, cell idx 32} -> ijk = {2,1,3}

    wells.push_back(createWell("second"));
    wells[1].updateConnections(wellCon,true);

    wells.push_back(createWell("third"));

    std::unordered_map<std::string, std::set<int>> futureConnections;
    futureConnections.emplace("third",std::set<int>{6,7});
    Dune::cpgrid::WellConnections wellConnections(wells, futureConnections, grid);


    if (grid.comm().size()>1) {

        grid.loadBalance(&wells, /* possibleFutureConnections = */ {},
                         /* transmissibilities = */ nullptr,
                         /* overlapLayers = */ 1,
                         /* partitionMethod = */ Dune::PartitionMethod::zoltanGoG,
                         /* level = */ 0);

        grid.addLgrsUpdateLeafView( /* cells_per_dim_vec = */ {{2,2,2}},
                                    /* startIJK_vec = */ {{1,1,0}},
                                    /* endIJK_vec = */ {{3,3,2}},
                                    /* lgr_name_vec = */ {"LGR1"});


        grid.syncDistributedGlobalCellIds();

        const auto& data = grid.currentData();
        Opm::checkCellGlobalIdUniquenessForInteriorCells(grid, data);
        Opm::checkConsecutiveChildGlobalIdsPerParent(grid);

    }

}

BOOST_AUTO_TEST_CASE(add_wells_and_loadBalance_level_zero_of_global_refined_cpgrid)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {2.0, 2.0, 2.0});
    grid.globalRefine(1);

    const std::array<int,3> expected_nxnynz = {8, 6, 6};
    const auto& actual_nxnynz = grid.logicalCartesianSize(); // For GR, logicalCartesianSize is well defined!
    BOOST_CHECK_EQUAL( expected_nxnynz[0], actual_nxnynz[0]);
    BOOST_CHECK_EQUAL( expected_nxnynz[1], actual_nxnynz[1]);
    BOOST_CHECK_EQUAL( expected_nxnynz[2], actual_nxnynz[2]);

    auto wellCon = std::make_shared<Opm::WellConnections>();
    wellCon->add(createConnection(0,0,0)); // {level 0, cell idx 0}
    wellCon->add(createConnection(0,1,0)); // {level 0, cell idx 1}
    wellCon->add(createConnection(0,1,1)); // {level 0, cell idx 3}

    std::vector<Dune::cpgrid::OpmWellType> wells;
    wells.push_back(createWell("first"));
    wells[0].updateConnections(wellCon,true);

    wellCon = std::make_shared<Opm::WellConnections>(); // reset
    wellCon->add(createConnection(0,0,2)); // {level 0, cell idx 4}
    wellCon->add(createConnection(0,1,2)); // {level 0, cell idx 5}
    wells.push_back(createWell("second"));
    wells[1].updateConnections(wellCon,true);

    wells.push_back(createWell("third"));

    std::unordered_map<std::string, std::set<int>> futureConnections;
    futureConnections.emplace("third",std::set<int>{6,7});
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
