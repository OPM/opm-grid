/*
  Copyright 2022-2023 Equinor ASA.

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
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/DefaultGeometryPolicy.hpp>
#include <opm/grid/cpgrid/EntityRep.hpp>
#include <opm/grid/cpgrid/Geometry.hpp>

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

void check_refinedPatch_grid(const std::array<int,3>& cells_per_dim,
                             const std::array<int,3>& start_ijk,
                             const std::array<int,3>& end_ijk,
                             const Dune::cpgrid::EntityVariable<Dune::cpgrid::Geometry<3, 3>,0>& refined_cells,
                             const Dune::cpgrid::EntityVariable<Dune::cpgrid::Geometry<2,3>,1>& refined_faces,
                             const Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& refined_corners)
{
    const std::array<int,3> patch_dim = {end_ijk[0]-start_ijk[0], end_ijk[1]-start_ijk[1], end_ijk[2]-start_ijk[2]};
    if ((patch_dim[0] == 0) || (patch_dim[1] == 0) || (patch_dim[2] == 0)) {
                    OPM_THROW(std::logic_error, "Empty patch. Cannot convert patch into cell.");
    }
    // Check amount of refined faces.
    int count_faces = (cells_per_dim[0]*patch_dim[0]*cells_per_dim[1]*patch_dim[1]*((cells_per_dim[2]*patch_dim[2])+1)) // 'bottom/top faces'
        +  (((cells_per_dim[0]*patch_dim[0])+1)*cells_per_dim[1]*patch_dim[1]*cells_per_dim[2]*patch_dim[2]) // 'front/back faces'
        + (cells_per_dim[0]*patch_dim[0]*((cells_per_dim[1]*patch_dim[1]) +1)*cells_per_dim[2]*patch_dim[2]);  // 'left/right faces'
    BOOST_CHECK_EQUAL(refined_faces.size(), count_faces);
    // Check amount of refined corners.
    int count_corners = ((cells_per_dim[0]*patch_dim[0])+1)*((cells_per_dim[1]*patch_dim[1])+1)*((cells_per_dim[2]*patch_dim[2])+1);
    BOOST_CHECK_EQUAL(refined_corners.size(), count_corners);

    int count_cells = cells_per_dim[0]*patch_dim[0]*cells_per_dim[1]*patch_dim[1]*cells_per_dim[2]*patch_dim[2];
    BOOST_CHECK_EQUAL(refined_cells.size(), count_cells);
}


void refinePatch_and_check(Dune::CpGrid& coarse_grid,
                           const std::array<int, 3>& cells_per_dim,
                           const std::array<int,3>& start_ijk,
                           const std::array<int,3>& end_ijk)
{
    
    // Call getLeafView2LevelsPatch()
    auto& data = coarse_grid.data_;
    coarse_grid.getLeafView2LevelsPatch(cells_per_dim, start_ijk, end_ijk);
    BOOST_CHECK(data.size()==3);
    check_refinedPatch_grid(cells_per_dim, start_ijk, end_ijk,
                            (*data[1]).geometry_.template geomVector<0>(),
                            (*data[1]).geometry_.template geomVector<1>(),
                            (*data[1]).geometry_.template geomVector<3>());
    
    /*  cpgrid::OrientedEntityTable<1,0> face_to_cell_computed;
    cell_to_face.makeInverseRelation(face_to_cell_computed);
    BOOST_CHECK(face_to_cell_computed == face_to_cell); */
} 

BOOST_AUTO_TEST_CASE(refine_patch)
{
     // Create a grid
    Dune::CpGrid coarse_grid;
    std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    std::array<int, 3> grid_dim = {4,3,3};
    std::array<int, 3> cells_per_dim_patch = {2,2,2};   
    std::array<int, 3> start_ijk = {1,0,1};
    std::array<int, 3> end_ijk = {3,2,3};  // then patch_dim = {3-1, 2-0, 3-1} ={2,2,2}
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    refinePatch_and_check(coarse_grid,
                          cells_per_dim_patch,
                          start_ijk,
                          end_ijk);
}

