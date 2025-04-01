// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2024 Equinor ASA.

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
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/grid/CpGrid.hpp>

#include <opm/grid/GraphOfGrid.hpp>
#include <opm/grid/GraphOfGridWrappers.hpp>

#include <opm/grid/utility/OpmWellType.hpp>
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <algorithm>
#include <string>


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

BOOST_AUTO_TEST_CASE(loadBalance_grid_with_lgrs_and_coarseCellWells)
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
    //

    // leaf grid view (LGR1 parent cell indices = {2,3,4,5}, LGR1 dimension 2x4x2)
    //
    // k = 0   1   | k = 1  [6 7 8 9]  | k = 2   [14 15 16 17]   | k = 3   19  |
    //         0   |        [2 3 4 5]  |         [10 11 12 13]   |         18  |

    // coarse cells leaf indices 0 and 1 -> well1
    // (equivalent coarse cells in level zero with indices 0 and 1 -> well1)
    auto wellCon = std::make_shared<Opm::WellConnections>();
    wellCon->add(createConnection(0,0,0)); // {level 0, cell idx 0} -> ijk = {0,0,0}
    wellCon->add(createConnection(0,1,0)); // {level 0, cell idx 1} -> ijk = {0,1,0}
    std::vector<Dune::cpgrid::OpmWellType> wells;
    wells.push_back(createWell("well1"));
    wells[0].updateConnections(wellCon,true);


    // coarse cells leaf indices 18 and 19 -> well2
    // (equivalent coarse cells in level zero with indices 6 and 7 -> well2)
    wellCon = std::make_shared<Opm::WellConnections>(); // reset
    wellCon->add(createConnection(0,0,3)); // {level 0, cell idx 6} -> ijk = {0,0,3}
    wellCon->add(createConnection(0,1,3)); // {level 0, cell idx 7} -> ijk = {0,1,3}

    wells.push_back(createWell("well2"));
    wells[1].updateConnections(wellCon,true);

    if (grid.comm().size()>1) {

        grid.loadBalance(&wells,
                         /* possibleFutureConnections = */ {},
                         /* transmissibilities = */ nullptr,
                         /* overlapLayers = */ 1,
                         /* partitionMethod = */ Dune::PartitionMethod::zoltanGoG,
                         /* level = */ 0);

    }
}

BOOST_AUTO_TEST_CASE(loadBalance_grid_with_lgrs_and_refinedCellWells)
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
    //

    // leaf grid view (LGR1 parent cell indices = {2,3,4,5}, LGR1 dimension 2x4x2)
    //
    // k = 0   1   | k = 1  [6 7 8 9]  | k = 2   [14 15 16 17]   | k = 3   19  |
    //         0   |        [2 3 4 5]  |         [10 11 12 13]   |         18  |

    // refined cells leaf indices [2 3 4 5] and [10 11 12 13] -> well1
    // (parent cells in level zero with indices 2 and 4 -> well1)
    auto wellCon = std::make_shared<Opm::WellConnections>();
    wellCon->add(createConnection(0,0,1)); // {level 0, cell idx 2} -> ijk = {0,0,1}
    wellCon->add(createConnection(0,0,2)); // {level 0, cell idx 4} -> ijk = {0,0,2}


    std::vector<Dune::cpgrid::OpmWellType> wells;
    wells.push_back(createWell("well1"));
    wells[0].updateConnections(wellCon,true);


    // refined cells leaf indices [6 7 8 9] and [14 15 16 17] -> well2
    // (parent cells level zero with indices 3 and 5 -> well2)
    wellCon = std::make_shared<Opm::WellConnections>(); // reset
    wellCon->add(createConnection(0,1,1)); // {level 0, cell idx 3} -> ijk = {0,1,1}
    wellCon->add(createConnection(0,1,2)); // {level 0, cell idx 5} -> ijk = {0,1,2}

    wells.push_back(createWell("well2"));
    wells[1].updateConnections(wellCon,true);

    if (grid.comm().size()>1) {

        grid.loadBalance(&wells,
                         /* possibleFutureConnections = */ {},
                         /* transmissibilities = */ nullptr,
                         /* overlapLayers = */ 1,
                         /* partitionMethod = */ Dune::PartitionMethod::zoltanGoG,
                         /* level = */ 0);

    }
}



BOOST_AUTO_TEST_CASE(loadBalance_grid_with_lgrs_and_coarseAndRefinedCellWells)
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
    //

    // leaf grid view (LGR1 parent cell indices = {2,3,4,5}, LGR1 dimension 2x4x2)
    //
    // k = 0   1   | k = 1  [6 7 8 9]  | k = 2   [14 15 16 17]   | k = 3   19  |
    //         0   |        [2 3 4 5]  |         [10 11 12 13]   |         18  |

    // coarse cell leaf index 1 and refined cells leaf indices [6 7 8 9] -> well1
    // (equivalent cell in level zero with index 1 and parent cell with index 3 -> well1)
    auto wellCon = std::make_shared<Opm::WellConnections>();
    wellCon->add(createConnection(0,1,0)); // {level 0, cell idx 1} -> ijk = {0,1,0}
    wellCon->add(createConnection(0,1,1)); // {level 0, cell idx 3} -> ijk = {0,1,1}


    std::vector<Dune::cpgrid::OpmWellType> wells;
    wells.push_back(createWell("well1"));
    wells[0].updateConnections(wellCon,true);


    // refined cells leaf indices [10 11 12 13] and coarse cell index 18 -> well2
    // (parent cell level zero with indices 4 and equivalent cell index 6 -> well2)
    wellCon = std::make_shared<Opm::WellConnections>(); // reset
    wellCon->add(createConnection(0,0,2)); // {level 0, cell idx 4} -> ijk = {0,0,2}
    wellCon->add(createConnection(0,0,3)); // {level 0, cell idx 6} -> ijk = {0,0,3}

    wells.push_back(createWell("well2"));
    wells[1].updateConnections(wellCon,true);

    if (grid.comm().size()>1) {

        grid.loadBalance(&wells,
                         /* possibleFutureConnections = */ {},
                         /* transmissibilities = */ nullptr,
                         /* overlapLayers = */ 1,
                         /* partitionMethod = */ Dune::PartitionMethod::zoltanGoG,
                         /* level = */ 0);

    }
}


BOOST_AUTO_TEST_CASE(loadBalance_grid_with_lgrs_and_oneCellWells)
{
    // create a grid with lgrs and wells
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dims = */ {1,2,4}, /* size = */ {1.,1.,1.});

    grid.addLgrsUpdateLeafView( /* cells_per_dim_vec = */ {{2,2,1}},
                                    /* startIJK_vec = */ {{0,0,1}},
                                    /* endIJK_vec = */ {{1,2,3}},
                                    /* lgr_name_vec = */ {"LGR1"});
     
     // refined cells leaf indices [6 7 8 9] (children of {level 0, cell idx 3})-> well1
    auto wellCon = std::make_shared<Opm::WellConnections>();
    wellCon->add(createConnection(0,1,1)); // {level 0, cell idx 3} -> ijk = {0,1,1}

    std::vector<Dune::cpgrid::OpmWellType> wells;
    wells.push_back(createWell("first"));
    wells[0].updateConnections(wellCon,true);

    wellCon = std::make_shared<Opm::WellConnections>(); // reset
    wellCon->add(createConnection(0,0,3)); // {level 0, cell idx 6} -> ijk = {0,0,3}

    wells.push_back(createWell("second"));
    wells[1].updateConnections(wellCon,true);


    if (grid.comm().size()>1) {

        grid.loadBalance(&wells, /* possibleFutureConnections = */ {},
                         /* transmissibilities = */ nullptr,
                         /* overlapLayers = */ 1,
                         /* partitionMethod = */ Dune::PartitionMethod::zoltanGoG,
                         /* level = */ 0);

    }
}




BOOST_AUTO_TEST_CASE(futureConnection_in_grid_with_lgrs, *boost::unit_test::disabled())
{
    // create a grid with lgrs and wells
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dims = */ {1,2,4}, /* size = */ {1.,1.,1.});

    grid.addLgrsUpdateLeafView( /* cells_per_dim_vec = */ {{2,2,1}},
                                    /* startIJK_vec = */ {{0,0,1}},
                                    /* endIJK_vec = */ {{1,2,3}},
                                    /* lgr_name_vec = */ {"LGR1"});

    auto wellCon = std::make_shared<Opm::WellConnections>();
    wellCon->add(createConnection(0,0,0));
    wellCon->add(createConnection(0,1,0));
    wellCon->add(createConnection(0,1,1));
    std::vector<Dune::cpgrid::OpmWellType> wells;
    wells.push_back(createWell("first"));
    wells[0].updateConnections(wellCon,true);

    wellCon = std::make_shared<Opm::WellConnections>(); // reset
    wellCon->add(createConnection(0,0,2));
    wellCon->add(createConnection(0,1,2));
    wells.push_back(createWell("second"));
    wells[1].updateConnections(wellCon,true);

    wells.push_back(createWell("third"));

  
    std::unordered_map<std::string, std::set<int>> futureConnections;
    futureConnections.emplace("third",std::set<int>{6,7});
    Dune::cpgrid::WellConnections wellConnections(wells,futureConnections,grid);

    for (const auto& [key, value] : futureConnections)
    {
        std::cout << key << std::endl;
        for (const auto& el : value)
        {
            std::cout<< el << std::endl;
        }
    }
    
 

     if (grid.comm().size()>1) {

        grid.loadBalance(&wells, futureConnections,
                         /* transmissibilities = */ nullptr,
                         /* overlapLayers = */ 1,
                         /* partitionMethod = */ Dune::PartitionMethod::zoltanGoG,
                         /* level = */ 0);

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
