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

#define BOOST_TEST_MODULE LGRTests
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/LookUpDataCpGrid.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <sstream>
#include <iostream>

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
    const auto& leaf_view = grid.leafGridView();

    Dune::LookUpData<Dune::CpGrid> lookUpData(grid);

    const auto& level0_view = grid.levelGridView(0);
    Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LeafGridView> leafMapper(leaf_view, Dune::mcmgElementLayout());
    Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LevelGridView> level0Mapper(level0_view, Dune::mcmgElementLayout());

    const auto& leaf_idSet = (*data.back()).local_id_set_;
    const auto& level0_idSet = (*data[0]).local_id_set_;

    for (const auto& elem : elements(leaf_view)) {
        auto featureInElem = lookUpData(elem, fake_feature);
        BOOST_CHECK(featureInElem == level0Mapper.index(elem.getOrigin()) +3);
        if (elem.hasFather()) { // leaf_cell has a father!
            const auto& id = (*leaf_idSet).id(elem);
            const auto& parent_id = (*level0_idSet).id(elem.father());
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
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);

    const std::array<int, 3> cells_per_dim = {2,2,2};
    const std::array<int, 3> startIJK = {1,0,1};
    const std::array<int, 3> endIJK = {3,2,3};  // patch_dim = {3-1, 2-0, 3-1} ={2,2,2}
    const std::string lgr_name = {"LGR1"};
    coarse_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

    lookup_check(coarse_grid);
}

BOOST_AUTO_TEST_CASE(single_cell_lgr_grid)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);

    const std::array<int, 3> cells_per_dim = {2,2,2};
    const std::array<int, 3> startIJK = {1,0,1};
    const std::array<int, 3> endIJK = {2,1,2};  // patch_dim = {2-1, 1-0, 2-1} ={1,1,1} -> Single Cell!
    const std::string lgr_name = {"LGR1"};
    coarse_grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

    lookup_check(coarse_grid);
}

BOOST_AUTO_TEST_CASE(lgrs_grid_A)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);

    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}, {4,4,4}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,2}, {3,2,2}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,1,1}, {1,1,3}, {4,3,3}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3"};
    coarse_grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

    lookup_check(coarse_grid);
}

BOOST_AUTO_TEST_CASE(lgrs_grid_B)
{
    // Create a grid
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);

    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {3,2,0}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,2,1}, {4,3,3}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    coarse_grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

    lookup_check(coarse_grid);
}
