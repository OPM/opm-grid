/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.
  Copyright 2022 Equinor ASA.

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

#define BOOST_TEST_MODULE GeometryTests
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
using namespace Dune;

class Null;

BOOST_AUTO_TEST_CASE(vertexgeom)
{
    typedef cpgrid::Geometry<0, 3> Geometry;
    // Default construction.
    Geometry g_default;

    // Construction from point.
    Geometry::GlobalCoordinate c(3.0);
    Geometry g(c);

    // Verification of properties.
    BOOST_CHECK(g.type().isVertex());
    BOOST_CHECK(g.affine());
    BOOST_CHECK_EQUAL(g.corners(), 1);
    BOOST_CHECK_EQUAL(g.corner(0), c);
    Geometry::LocalCoordinate lc(0.0);
    BOOST_CHECK_EQUAL(g.global(lc), c);
    // BOOST_CHECK_THROW(g.local(c), std::exception);
    BOOST_CHECK_EQUAL(g.integrationElement(lc), 1.0);
    BOOST_CHECK_EQUAL(g.volume(), 1.0);
    BOOST_CHECK_EQUAL(g.center(), c);
    // BOOST_CHECK_THROW(g.jacobianTransposed(lc), std::exception);
    // BOOST_CHECK_THROW(g.jacobianInverseTransposed(lc), std::exception);
}


BOOST_AUTO_TEST_CASE(intersectiongeom)
{
    typedef cpgrid::Geometry<2, 3> Geometry;
    // Default construction.
    Geometry g_default;

    // Construction from point and volume.
    Geometry::GlobalCoordinate c(3.0);
    Geometry::ctype v = 8.0;
    Geometry g(c, v);

    // Verification of properties.
    BOOST_CHECK(g.type().isNone());
    BOOST_CHECK(g.affine());
    BOOST_CHECK_EQUAL(g.corners(), 0); // "corners()" returns '0' (None Type -> singular geometry -> no corners)
    // BOOST_CHECK_THROW(g.corner(0), std::exception);
    Geometry::LocalCoordinate lc(0.0);
    BOOST_CHECK_THROW(g.global(lc), std::exception);
    BOOST_CHECK_THROW(g.local(c), std::exception);
    BOOST_CHECK_EQUAL(g.integrationElement(lc), v);
    BOOST_CHECK_EQUAL(g.volume(), v);
    BOOST_CHECK_EQUAL(g.center(), c);
    BOOST_CHECK_THROW(g.jacobianTransposed(lc), std::exception);
    BOOST_CHECK_THROW(g.jacobianInverseTransposed(lc), std::exception);

}


BOOST_AUTO_TEST_CASE(cellgeom)
{
    typedef cpgrid::Geometry<3, 3> Geometry;
    // Default construction.
    Geometry g_default;

    // Construction from point and volume.
    // This is a dangerous constructor kept for backwards compatibility,
    // these checks may be removed if constructor removed.
    // This possibly dangerous constructor is available for the benefit of
    // CpGrid::readSintefLegacyFormat().
    Geometry::GlobalCoordinate c(3.0);
    Geometry::ctype v = 8.0;
    Geometry g_dangerous(c, v);

    // Construction from point, volume, points and pointindices.
    // Construction from
    // 1. center
    // 2. volume
    // 3. all corners (a pointer of type 'const cpgrid::Geometry<0, 3>')
    // 4. corner indices (a pointer to the first element of "global_refined_cell8corners_indices_storage")
    // First a unit cube, i.e. the mapping represented is the identity.
    typedef Geometry::GlobalCoordinate GC;
    c = GC(0.5);
    v = 1.0;
    GC corners[8];
    GC cor;
    for (int k = 0; k < 2; ++k) {
        cor[2] = k;
        for (int j = 0; j < 2; ++j) {
            cor[1] = j;
            for (int i = 0; i < 2; ++i) {
                cor[0] = i;
                corners[4*k + 2*j + i] = cor;
            }
        }
    }
//     for (int i = 0; i < 8; ++i) {
//         std::cout << corners[i] << std::endl;
//     }
    cpgrid::EntityVariable<cpgrid::Geometry<0, 3>, 3> pg;
    pg.reserve(8);
    for (const auto& crn : corners)
    {
        pg.push_back(cpgrid::Geometry<0, 3>(crn));
    }

    int cor_idx[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
    Geometry g(c, v, pg, cor_idx);

    // Verification of properties.
    BOOST_CHECK(g.type().isCube());
    BOOST_CHECK(!g.affine());
    BOOST_CHECK_EQUAL(g.corners(), 8);
    for (int i = 0; i < 8; ++i) {
        BOOST_CHECK_EQUAL(g.corner(i), corners[i]);   // "corner(i)" returns 'pg[cor_idx_[i]].center()'
    }
    BOOST_CHECK_EQUAL(g.volume(), v);
    BOOST_CHECK_EQUAL(g.center(), c);

    // Verification of properties that depend on the mapping.'
    // Construction refined (local) corners.
    typedef Geometry::LocalCoordinate LC;
    const int N = 5;  //  4 in each direction
    const int num_pts = N*N*N;  // total amount of corners
    LC testpts[num_pts] = { LC(0.0) };
    LC pt(0.0);
    for (int k = 0; k < N; ++k) {
        pt[2] = double(k)/double(N-1);
        for (int j = 0; j < N; ++j) {
            pt[1] = double(j)/double(N-1);
            for (int i = 0; i < N; ++i) {
                pt[0] = double(i)/double(N-1);
                testpts[k*N*N + j*N + i] = pt;
//                 std::cout << pt << std::endl;
            }
        }
    }
    Geometry::JacobianTransposed id(0.0);
    id[0][0] = id[1][1] = id[2][2] = 1.0;
    for (int i = 0; i < num_pts; ++i) { // all refined corners
        BOOST_CHECK_EQUAL(g.global(testpts[i]), testpts[i]);
        BOOST_CHECK_EQUAL(g.local(g.global(testpts[i])), testpts[i]);
        BOOST_CHECK_EQUAL(g.integrationElement(testpts[i]), 1.0);   // 'MatrixHelperType::template sqrtDetAAT<3, 3>(testpts[i])'
        BOOST_CHECK_EQUAL(g.jacobianTransposed(testpts[i]), id);
        BOOST_CHECK_EQUAL(g.jacobianInverseTransposed(testpts[i]), id);
    }

    // Next testcase: a degenerate hexahedron, wedge shaped.
    typedef Geometry::GlobalCoordinate GC;
    c = GC(1.0/3.0); c[2] = 0.5;
    v = 0.5;
    corners[5][2] = 0.0;
    corners[7][2] = 0.0;

    cpgrid::EntityVariable<cpgrid::Geometry<0, 3>, 3> pg1;
    pg1.reserve(8);
    for (const auto& crn : corners)
    {
        pg1.push_back(cpgrid::Geometry<0, 3>(crn));
    }
    g = Geometry(c, v, pg1, cor_idx);

    // Verification of properties.
    BOOST_CHECK(g.type().isCube());
    BOOST_CHECK(!g.affine());
    BOOST_CHECK_EQUAL(g.corners(), 8);
    for (int i = 0; i < 8; ++i) {
        BOOST_CHECK_EQUAL(g.corner(i), corners[i]);
    }
    BOOST_CHECK_EQUAL(g.volume(), v);
    BOOST_CHECK_EQUAL(g.center(), c);

    struct Wedge
    {
        static GC global(const LC& lc)
        {
            GC gc(0.0);
            gc[0] = lc[0];
            gc[1] = lc[1];
            gc[2] = (1.0 - lc[0])*lc[2];
            return gc;
        }
        static double integrationElement(const LC& lc)
        {
            return 1.0 - lc[0];
        }
        static Geometry::JacobianTransposed jacobianTransposed(const LC& lc)
        {
            Geometry::JacobianTransposed Jt(0.0);
            Jt[0][0] = 1.0;
            Jt[0][2] = -lc[2];
            Jt[1][1] = 1.0;
            Jt[2][2] = 1.0 - lc[0];
            return Jt;
        }
    };

    // Verification of properties that depend on the mapping.
    const double tolerance = 1e-14;
    for (int i = 0; i < num_pts; ++i) {
        GC gl = g.global(testpts[i]);
        BOOST_CHECK_EQUAL(gl, Wedge::global(testpts[i]));
        BOOST_CHECK_EQUAL(g.integrationElement(testpts[i]), Wedge::integrationElement(testpts[i]));
        Geometry::JacobianTransposed Jt = Wedge::jacobianTransposed(testpts[i]);
        BOOST_CHECK_EQUAL(g.jacobianTransposed(testpts[i]), Jt);
        if (testpts[i][0] < 1.0) {
            // Only do this test if we are away from the degeneracy.
            LC diff = g.local(gl);
            diff -= testpts[i];
            BOOST_CHECK_SMALL(diff.two_norm(), tolerance);
            Geometry::Jacobian Jit = Jt; // This implicitly assumes that the Jacobian is square.
            Jit.invert();
            BOOST_CHECK_EQUAL(g.jacobianInverseTransposed(testpts[i]), Jit);
        }
    }
}


#define CHECK_COORDINATES(c1, c2)                                        \
    for (int c = 0; c < 3; c++) {                                        \
        BOOST_TEST(c1[c] == c2[c], boost::test_tools::tolerance(1e-12)); \
    }

void
check_refined_grid(const cpgrid::Geometry<3, 3>& parent,
                   const cpgrid::EntityVariable<cpgrid::Geometry<3, 3>,0>& refined,
                   const cpgrid::EntityVariable<cpgrid::Geometry<2,3>,1>& refined_faces,
                   const cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>& refined_corners,
                   const std::array<int, 3>& cells_per_dim)
{
    // Check amount of refined faces.
    int count_faces = (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1)) // 'bottom/top faces'
                    + (cells_per_dim[0]*(cells_per_dim[1]+1)*cells_per_dim[2]) // 'front/back faces'
                    + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2]);  // 'left/right faces'
    BOOST_CHECK_EQUAL(refined_faces.size(), count_faces);
    // Check amount of refined corners.
    int count_corners = (cells_per_dim[0]+1)*(cells_per_dim[1]+1)*(cells_per_dim[2]+1);
    BOOST_CHECK_EQUAL(refined_corners.size(), count_corners);


    using Geometry = cpgrid::Geometry<3, 3>;
    using GlobalCoordinate = Geometry::GlobalCoordinate;

    int count = cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2];
    BOOST_CHECK_EQUAL(refined.size(), count);


    // Check that the corresponding refined corners match the
    // the parent cell corners. Notice that 'each parent corner' coincide with
    // one refined corner of one of (at most 8 different) refined cells.
    // (l *(cells_per_dim[2]-1)* slice) + (m*(cells_per_dim[1]-1) * cells_per_dim[0]) + (n*cells_per_dim[0]-1)
    auto& cell_corn0 = refined.get(0); // l=0,m=0,n=0
    auto& cell_corn1 = refined.get(cells_per_dim[0]-1); // l=0,m=0,n=1
    auto& cell_corn2 = refined.get((cells_per_dim[1]-1) * cells_per_dim[0]); // l=0,m=1,n=0
    auto& cell_corn3 = refined.get(((cells_per_dim[1]-1) * cells_per_dim[0]) + (cells_per_dim[0]-1)); // l=0,m=1,n=1
    auto& cell_corn4 = refined.get((cells_per_dim[2]-1)*cells_per_dim[0] * cells_per_dim[1]); // l=1,m=0,n=0
    auto& cell_corn5 = refined.get(((cells_per_dim[2]-1)*cells_per_dim[0] * cells_per_dim[1]) +(cells_per_dim[0]-1)); // l=1,m=0,n=1
    auto& cell_corn6 = refined.get(((cells_per_dim[2]-1)*  cells_per_dim[0] * cells_per_dim[1])
                                   +((cells_per_dim[1]-1) * cells_per_dim[0])); // l=1,m=1,n=0
    auto& cell_corn7 = refined.get(((cells_per_dim[2]-1)*  cells_per_dim[0] * cells_per_dim[1])
                                   +((cells_per_dim[1]-1) * cells_per_dim[0]) + (cells_per_dim[0]-1)); // l=1,m=1,n=1
    CHECK_COORDINATES(cell_corn0.corner(0), parent.corner(0));
    CHECK_COORDINATES(cell_corn1.corner(1), parent.corner(1));
    CHECK_COORDINATES(cell_corn2.corner(2), parent.corner(2));
    CHECK_COORDINATES(cell_corn3.corner(3), parent.corner(3));
    CHECK_COORDINATES(cell_corn4.corner(4), parent.corner(4));
    CHECK_COORDINATES(cell_corn5.corner(5), parent.corner(5));
    CHECK_COORDINATES(cell_corn6.corner(6), parent.corner(6));
    CHECK_COORDINATES(cell_corn7.corner(7), parent.corner(7));

    // Make sure the corners of neighboring cells overlap.
    for (int k = 0; k < cells_per_dim[2]; k++) {
        int slice = cells_per_dim[1] * cells_per_dim[0];
        for (int j = 0; j < cells_per_dim[1]; j++) {
            for (int i = 0; i < cells_per_dim[0]; i++) {
                auto& r0 = refined.get(k * slice + cells_per_dim[0] * j + i);
                if (i < cells_per_dim[0] - 1) {
                    auto& r1 = refined.get(k * slice + cells_per_dim[0] * j + i + 1);
                    CHECK_COORDINATES(r0.corner(1), r1.corner(0));
                    CHECK_COORDINATES(r0.corner(3), r1.corner(2));
                    CHECK_COORDINATES(r0.corner(5), r1.corner(4));
                    CHECK_COORDINATES(r0.corner(7), r1.corner(6));
                }
                if (j < cells_per_dim[1] - 1) {
                    auto& r1 = refined.get(k * slice + cells_per_dim[0] * (j + 1) + i);
                    CHECK_COORDINATES(r0.corner(2), r1.corner(0));
                    CHECK_COORDINATES(r0.corner(3), r1.corner(1));
                    CHECK_COORDINATES(r0.corner(6), r1.corner(4));
                    CHECK_COORDINATES(r0.corner(7), r1.corner(5));
                }
                if (k < cells_per_dim[2] - 1) {
                    auto& r1 = refined.get((k + 1) * slice + cells_per_dim[0] * j + i);
                    CHECK_COORDINATES(r0.corner(4), r1.corner(0));
                    CHECK_COORDINATES(r0.corner(5), r1.corner(1));
                    CHECK_COORDINATES(r0.corner(6), r1.corner(2));
                    CHECK_COORDINATES(r0.corner(7), r1.corner(3));
                }
            }
        }
    }

    // Check the centers of the cells.
    for (auto r : refined) {
        GlobalCoordinate center = {0.0, 0.0, 0.0};
        for (int h = 0; h < 8; h++) {
            for (int c = 0; c < 3; c++) {
                center[c] += r.corner(h)[c] / 8.0;
            }
        }
        CHECK_COORDINATES(r.center(), center);
    }

    // Check that the weighted mean of all centers equals the parent center
    GlobalCoordinate center = {0.0, 0.0, 0.0};
    for (auto r : refined) {
        for (int c = 0; c < 3; c++) {
            center[c] += r.center()[c] * r.volume()
                / parent.volume();
        }
    }
    CHECK_COORDINATES(parent.center(), center);

    // Check that mean of all corners equals the center of the parent.
    center = {0.0, 0.0, 0.0};
    for (auto r : refined) {
        for (int h = 0; h < 8; h++) {
            for (int c = 0; c < 3; c++) {
                center[c] += r.corner(h)[c] / count / 8;
            }
        }
    }
    CHECK_COORDINATES(parent.center(), center);

    // Check the total volume against the parent volume.
    Geometry::ctype volume = 0.0;
    for (auto r : refined) {
        volume += r.volume();
    }
    BOOST_CHECK_CLOSE(volume, parent.volume(), 1e-6);
}

void refine_and_check(const cpgrid::Geometry<3, 3>& parent_geometry,
                      const std::array<int, 3>& cells,
                      bool is_simple = false)
{
    using cpgrid::DefaultGeometryPolicy;
    CpGrid refined_grid;
    auto& child_view_data = *refined_grid.current_view_data_;
    cpgrid::OrientedEntityTable<0, 1>& cell_to_face = child_view_data.cell_to_face_;
    Opm::SparseTable<int>& face_to_point = child_view_data.face_to_point_;
    DefaultGeometryPolicy& geometries = child_view_data.geometry_;
    std::vector<std::array<int, 8>>& cell_to_point = child_view_data.cell_to_point_;
    cpgrid::OrientedEntityTable<1,0>& face_to_cell = child_view_data.face_to_cell_;
    cpgrid::EntityVariable<enum face_tag, 1>& face_tags = child_view_data.face_tag_;
    cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>, 1>& face_normals = child_view_data.face_normals_;

    parent_geometry.refine(cells, geometries, cell_to_point,
                           cell_to_face, face_to_point, face_to_cell,
                           face_tags, face_normals);
    check_refined_grid(parent_geometry, geometries.template geomVector<0>(),
                       geometries.template geomVector<1>(),
                       geometries.template geomVector<3>(), cells);
    cpgrid::OrientedEntityTable<1,0> face_to_cell_computed;
    cell_to_face.makeInverseRelation(face_to_cell_computed);
    BOOST_CHECK(face_to_cell_computed == face_to_cell);

    if (is_simple)
    {
        // Create a grid that is equivalent to the refinement
        Dune::CpGrid equivalent_refined_grid;
        std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
        std::size_t i = 0;
        for (auto& size: cell_sizes)
            size /= cells[i++];

        equivalent_refined_grid.createCartesian(cells, cell_sizes);

        // Check that the sizes match
        BOOST_CHECK(equivalent_refined_grid.size(0) == refined_grid.size(0));
        BOOST_CHECK(equivalent_refined_grid.size(1) == refined_grid.size(1));
        BOOST_CHECK(equivalent_refined_grid.size(3) == refined_grid.size(3));

        // Check that the points (ordering/coordinates) matches
        auto equiv_point_iter = equivalent_refined_grid.current_view_data_
            ->geometry_.geomVector<3>().begin();
        for(const auto& point: geometries.geomVector<3>())
        {
            CHECK_COORDINATES(point.center(), equiv_point_iter->center());
            ++equiv_point_iter;
        }

        // Check that cell to point matches (Implicitly checks for correct cell ordering)
        BOOST_CHECK(refined_grid.current_view_data_->cell_to_point_ ==
                    equivalent_refined_grid.current_view_data_->cell_to_point_);

        // Check that the cell centers and volume match
        auto equiv_cell_iter = equivalent_refined_grid.current_view_data_
            ->geometry_.geomVector<0>().begin();
        for(const auto& cell: geometries.geomVector<0>())
        {
            CHECK_COORDINATES(cell.center(), equiv_cell_iter->center());
            BOOST_CHECK_CLOSE(cell.volume(), equiv_cell_iter->volume(), 1e-6);
            ++equiv_cell_iter;
        }

        const auto& grid_view = refined_grid.leafGridView();
        const auto& equiv_grid_view = equivalent_refined_grid.leafGridView();

        auto equiv_element_iter = equiv_grid_view.begin<0>();
        for(const auto& element: elements(grid_view))
        {
            for(const auto& intersection: intersections(grid_view, element))
            {
                // find matching intersection (needed as ordering is allowed to be different
                bool matching_intersection_found = false;
                for(auto& intersection_match: intersections(equiv_grid_view, *equiv_element_iter))
                {
                    if(intersection_match.indexInInside() == intersection.indexInInside())
                    {
                        BOOST_CHECK(intersection_match.neighbor() == intersection.neighbor());

                        if(intersection.neighbor())
                        {
                            BOOST_CHECK(intersection_match.indexInOutside() == intersection.indexInOutside());
                        }

                        CHECK_COORDINATES(intersection_match.centerUnitOuterNormal(), intersection.centerUnitOuterNormal());
                        const auto& geom_match = intersection_match.geometry();
                        BOOST_TEST(0 == 1e-11, boost::test_tools::tolerance(1e-8));
                        const auto& geom =  intersection.geometry();
                        BOOST_CHECK_CLOSE(geom_match.volume(), geom.volume(), 1e-6);
                        CHECK_COORDINATES(geom_match.center(), geom.center());
                        BOOST_CHECK(geom_match.corners() == geom.corners());

                        decltype(geom.corner(0)) sum_match{}, sum{};

                        for(int cor = 0; cor < geom.corners(); ++cor)
                        {
                            sum += geom.corner(cor);
                            sum_match += geom_match.corner(1);
                        }
                        CHECK_COORDINATES(sum, sum_match);
                        matching_intersection_found = true;
                        break;
                    }
                }
                BOOST_CHECK(matching_intersection_found);
            }
            ++equiv_element_iter;
        }
    }
    // Create a grid that is equivalent to the refinement
    /*  Dune::CpGrid coarse_grid;
    std::array<double, 3> cell_sizes_new = {1.0, 1.0, 1.0};
    std::array<int, 3> coarse_grid_dim = {4,3,3};
    std::array<int, 3> cells_per_dim_patch = {2,2,2};
    std::array<int, 3> start_ijk = {1,0,1};
    std::array<int, 3> end_ijk = {3,2,3};  // then patch_dim = {3-1,2-0,3-1} ={2,2,2}
    coarse_grid.createCartesian(coarse_grid_dim, cell_sizes_new);
    // Call refinedBlockPatch()
    coarse_grid.current_view_data_->refineBlockPatch(cells_per_dim_patch, start_ijk, end_ijk);
    // Create a pointer pointing at the CpGridData object coarse_grid.current_view_data_.
    std::shared_ptr<Dune::cpgrid::CpGridData> coarse_grid_ptr =  std::make_shared<Dune::cpgrid::CpGridData>();
    *coarse_grid.current_view_data_ = *coarse_grid_ptr;
    // Create a vector of shared pointers of CpGridData type.
    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>> data;
    // Add coarse_grid_ptr to data.
    data.push_back(coarse_grid_ptr);
    // Call getLeafView2Levels()
    coarse_grid.getLeafView2Levels(data, cells_per_dim_patch, start_ijk, end_ijk);
    // Call addLevel()
    const int level_to_refine = 0;
    std::vector<std::array<int,2>> future_leaf_corners;
    std::vector<std::array<int,2>> future_leaf_faces;
    std::vector<std::array<int,2>> future_leaf_cells;
    coarse_grid.addLevel(data, level_to_refine, cells_per_dim_patch, start_ijk, end_ijk,
    future_leaf_corners, future_leaf_faces, future_leaf_cells);*/
}

BOOST_AUTO_TEST_CASE(refine_simple_cube)
{
    using Geometry = cpgrid::Geometry<3, 3>;
    using GlobalCoordinate = Geometry::GlobalCoordinate;

    const GlobalCoordinate corners[8] = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {1.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 0.0, 1.0},
        {0.0, 1.0, 1.0},
        {1.0, 1.0, 1.0},
    };
    const GlobalCoordinate c = {0.5, 0.5, 0.5};
    const Geometry::ctype v = 1.0;

    cpgrid::EntityVariable<cpgrid::Geometry<0, 3>, 3> pg;
    for (const auto& crn : corners) {
        pg.push_back(cpgrid::Geometry<0, 3>(crn));
    }

    int cor_idx[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    Geometry g(c, v, pg, cor_idx);

    refine_and_check(g, {1, 1, 1}, true);
    refine_and_check(g, {2, 3, 4}, true);
}


BOOST_AUTO_TEST_CASE(refine_distorted_cube)
{
    using Geometry = cpgrid::Geometry<3, 3>;
    using GlobalCoordinate = Geometry::GlobalCoordinate;

    // Distorted cube:
    const GlobalCoordinate corners[8] = {
        {0.1, 0.2, 0.3},
        {1.2, 0.3, 0.4},
        {0.3, 1.4, 0.5},
        {1.4, 1.5, 0.6},
        {0.5, 0.6, 1.7},
        {1.6, 0.7, 1.8},
        {0.7, 1.8, 1.9},
        {1.8, 1.9, 2.0},
    };

    // Arbitrary volume:
    const Geometry::ctype v = 123.0;

    // Calculate the centroid:
    GlobalCoordinate center = {0.0, 0.0, 0.0};
    for (int h = 0; h < 8; h++) {
        for (int c = 0; c < 3; c++) {
            center[c] += corners[h][c] / 8.0;
        }
    }

    cpgrid::EntityVariable<cpgrid::Geometry<0, 3>, 3> pg;
    for (const auto& crn : corners) {
        pg.push_back(cpgrid::Geometry<0, 3>(crn));
    }

    int cor_idx[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    Geometry g(center, v, pg, cor_idx);
    refine_and_check(g, {1, 1, 1});
    refine_and_check(g, {2, 3, 4});

}

void refinePatch_and_check(const std::array<int,3>&,
                           const std::array<int,3>&,
                           const std::array<int,3>&)
{
    // Create a grid that is equivalent to the refinement
    Dune::CpGrid coarse_grid;
    std::array<double, 3> cell_sizes_new = {1.0, 1.0, 1.0};
    std::array<int, 3> coarse_grid_dim = {4,3,3};
    std::array<int, 3> cells_per_dim_patch = {2,2,2};
    std::array<int, 3> start_ijk = {1,0,1};
    std::array<int, 3> end_ijk = {3,2,3};  // then patch_dim = {3-1,2-0,3-1} ={2,2,2}
    coarse_grid.createCartesian(coarse_grid_dim, cell_sizes_new);
    // Call refinedBlockPatch()
    coarse_grid.current_view_data_->refineBlockPatch(cells_per_dim_patch, start_ijk, end_ijk);

}
BOOST_AUTO_TEST_CASE(refine_patch)
{
    refinePatch_and_check({}, {}, {});
}
