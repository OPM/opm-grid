#include <config.h>

#define NVERBOSE // to suppress our messages when throwing

#define BOOST_TEST_MODULE Process_Grdecl
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

#include <opm/grid/cpgpreprocess/preprocess.h>
#include <opm/grid/cpgpreprocess/make_edge_conformal.hpp>

#include <array>
#include <cstddef>
#include <initializer_list>
#include <optional>
#include <vector>

namespace {
    class TestGrid
    {
    public:
        explicit TestGrid(const std::array<int,3>& dims)
            : dims_ { dims }
        {}

        ~TestGrid()
        {
            if (this->g_.has_value()) {
                free_processed_grid(&*this->g_);
            }
        }

        TestGrid& coord(const std::vector<double>& c)
        {
            this->coord_ = c;
            return *this;
        }

        TestGrid& zcorn(const std::vector<double>& z)
        {
            this->zcorn_ = z;
            return *this;
        }

        TestGrid& actnum(const std::vector<int>& a)
        {
            this->actnum_ = a;
            return *this;
        }

        TestGrid& ztol(const double tol)
        {
            this->ztol_ = tol;
            return *this;
        }

        TestGrid& pinchActive(const bool active)
        {
            this->pinch_active_ = static_cast<int>(active);
            return *this;
        }

        TestGrid& edgeConformal(const bool conformal)
        {
            this->edge_conformal_ = static_cast<int>(conformal);
            return *this;
        }

        TestGrid& process()
        {
            auto input = grdecl{};

            std::copy(this->dims_.begin(), this->dims_.end(), input.dims);

            input.coord  = this->coord_.data();
            input.zcorn  = this->zcorn_.data();
            input.actnum = this->actnum_.data();

            this->g_.emplace(processed_grid{});

            this->status_ =
                process_grdecl(this->pinch_active_,
                               this->edge_conformal_,
                               this->ztol_,
                               &input,
                               /* is_aquifer_cell = */ nullptr,
                               &*this->g_);

            if ((this->status_ != 0) && (this->edge_conformal_ != 0)) {
                add_cell_face_mapping(&*this->g_);
                this->status_ = make_edge_conformal(&*this->g_);
            }

            return *this;
        }

        int status() const { return this->status_; }

        const processed_grid& grid() const { return *this->g_; }

    private:
        std::array<int,3> dims_{};
        std::vector<double> coord_{};
        std::vector<double> zcorn_{};
        std::vector<int> actnum_{};

        double ztol_{0.0};
        int pinch_active_{0};
        int edge_conformal_{0};
        int status_{};

        std::optional<processed_grid> g_{};
    };

    TestGrid zeroThicknessMiddle()
    {
        return TestGrid {{ 1, 1, 3 }}
            .coord({
                    0.0, 0.0, 0.0,   0.0, 0.0, 2.0,
                    1.0, 0.0, 0.0,   1.0, 0.0, 2.0,
                    0.0, 1.0, 0.0,   0.0, 1.0, 2.0,
                    1.0, 1.0, 0.0,   1.0, 1.0, 2.0,
                })

            .zcorn({
                    // Top cell--thickness 1
                    0.0, 0.0,
                    0.0, 0.0,
                    1.0, 1.0,
                    1.0, 1.0,

                    // Middle cell--thickness 0
                    1.0, 1.0,
                    1.0, 1.0,
                    1.0, 1.0,
                    1.0, 1.0,

                    // Bottom cell--thickness 1
                    1.0, 1.0,
                    1.0, 1.0,
                    2.0, 2.0,
                    2.0, 2.0,
                })

            .actnum({ 1, 0, 1, });
    }

    TestGrid unitThicknessMiddle()
    {
        return TestGrid {{ 1, 1, 3 }}
            .coord({
                    0.0, 0.0, 0.0,   0.0, 0.0, 2.0,
                    1.0, 0.0, 0.0,   1.0, 0.0, 2.0,
                    0.0, 1.0, 0.0,   0.0, 1.0, 2.0,
                    1.0, 1.0, 0.0,   1.0, 1.0, 2.0,
                })

            .zcorn({
                    // Top cell--thickness 1
                    0.0, 0.0,
                    0.0, 0.0,
                    1.0, 1.0,
                    1.0, 1.0,

                    // Middle cell--thickness 1
                    1.0, 1.0,
                    1.0, 1.0,
                    2.0, 2.0,
                    2.0, 2.0,

                    // Bottom cell--thickness 1
                    2.0, 2.0,
                    2.0, 2.0,
                    3.0, 3.0,
                    3.0, 3.0,
                })

            .actnum({ 1, 0, 1, });
    }

    TestGrid singleFaultConnection()
    {
        return TestGrid {{ 2, 1, 1 }}
            .coord({
                    0.0, 0.0, 0.0,   0.0, 0.0, 2.0,
                    1.0, 0.0, 0.0,   1.0, 0.0, 2.0,
                    2.0, 0.0, 0.0,   2.0, 0.0, 2.0,
                    0.0, 1.0, 0.0,   0.0, 1.0, 2.0,
                    1.0, 1.0, 0.0,   1.0, 1.0, 2.0,
                    2.0, 1.0, 0.0,   2.0, 1.0, 2.0,
                })

            .zcorn({
                    // Left cell, back row, top surface
                    1.0, 1.0,

                    // Right cell, back row, top surface
                    1.5, 1.5,

                    // Left cell, front row, top surface
                    1.0, 1.0,

                    // Right cell, front row, top surface
                    0.5, 1.0,

                    // ----------------------------------

                    // Left cell, back row, bottom surface
                    2.0, 2.0,

                    // Right cell, back row, bottom surface
                    2.0, 2.0,

                    // Left cell, front row, bottom surface
                    2.0, 2.0,

                    // Right cell, front row, bottom surface
                    2.0, 2.0,
                })

            .actnum({ 1, 1, });
    }
} // Anonymous namespace

BOOST_AUTO_TEST_SUITE(Regular_Processing)

BOOST_AUTO_TEST_CASE(Zero_Thickness_Middle_No_Pinch)
{
    auto testCase = zeroThicknessMiddle()
        .pinchActive(false)
        .ztol(0.0);

    testCase.process();
    BOOST_REQUIRE_EQUAL(testCase.status(), 1);

    const auto& out = testCase.grid();

    BOOST_CHECK_EQUAL(out.dimensions[0], 1);
    BOOST_CHECK_EQUAL(out.dimensions[1], 1);
    BOOST_CHECK_EQUAL(out.dimensions[2], 3);

    BOOST_CHECK_EQUAL(out.number_of_nodes, 12);
    BOOST_CHECK_EQUAL(out.number_of_faces, 12);
    BOOST_CHECK_EQUAL(out.number_of_cells,  2);

    // Local-to-global cell mapping
    {
        const auto expect = std::vector { 0, 2 };
        BOOST_CHECK_EQUAL_COLLECTIONS(out.local_cell_index,
                                      out.local_cell_index + out.number_of_cells,
                                      expect.begin(), expect.end());
    }

    // Face-to-cell mapping
    {
        const auto expect = std::vector {
            -1, 0,              // I-min (top)
            -1, 1,              // I-min (bot)
            0, -1,              // I-max (top)
            1, -1,              // I-max (bot)

            // -------------------------------

            -1, 0,              // J-min (top)
            -1, 1,              // J-min (bot)
            0, -1,              // J-max (top)
            1, -1,              // J-max (bot)

            // -------------------------------

            -1, 0,              // K-min (top)
            0, -1,              // K-max (top)
            -1, 1,              // K-min (bot)
            1, -1,              // K-max (bot)
        };

        BOOST_CHECK_EQUAL_COLLECTIONS(out.face_neighbors,
                                      out.face_neighbors + 2*out.number_of_faces,
                                      expect.begin(), expect.end());
    }

    // Face taxonomy
    {
        const auto expect = std::vector {
            face_tag::I_FACE, face_tag::I_FACE, face_tag::I_FACE, face_tag::I_FACE,
            face_tag::J_FACE, face_tag::J_FACE, face_tag::J_FACE, face_tag::J_FACE,
            face_tag::K_FACE, face_tag::K_FACE, face_tag::K_FACE, face_tag::K_FACE,
        };

        BOOST_CHECK_EQUAL_COLLECTIONS(out.face_tag,
                                      out.face_tag + out.number_of_faces,
                                      expect.begin(), expect.end());
    }
}

BOOST_AUTO_TEST_CASE(Zero_Thickness_Middle_With_Pinch)
{
    auto testCase = zeroThicknessMiddle()
        .pinchActive(true)
        .ztol(0.0);

    testCase.process();
    BOOST_REQUIRE_EQUAL(testCase.status(), 1);

    const auto& out = testCase.grid();

    BOOST_CHECK_EQUAL(out.dimensions[0], 1);
    BOOST_CHECK_EQUAL(out.dimensions[1], 1);
    BOOST_CHECK_EQUAL(out.dimensions[2], 3);

    BOOST_CHECK_EQUAL(out.number_of_nodes, 12);
    BOOST_CHECK_EQUAL(out.number_of_faces, 11);
    BOOST_CHECK_EQUAL(out.number_of_cells,  2);

    // Local-to-global cell mapping
    {
        const auto expect = std::vector { 0, 2 };
        BOOST_CHECK_EQUAL_COLLECTIONS(out.local_cell_index,
                                      out.local_cell_index + out.number_of_cells,
                                      expect.begin(), expect.end());
    }

    // Face-to-cell mapping
    {
        const auto expect = std::vector {
            -1, 0,              // I-min (top)
            -1, 1,              // I-min (bot)
            0, -1,              // I-max (top)
            1, -1,              // I-max (bot)

            // -------------------------------

            -1, 0,              // J-min (top)
            -1, 1,              // J-min (bot)
            0, -1,              // J-max (top)
            1, -1,              // J-max (bot)

            // -------------------------------

            -1, 0,              // K-min (top)
            0,  1,              // K-max (top)/K-min (bot)
            1, -1,              // K-max (bot)
        };

        BOOST_CHECK_EQUAL_COLLECTIONS(out.face_neighbors,
                                      out.face_neighbors + 2*out.number_of_faces,
                                      expect.begin(), expect.end());
    }

    // Face taxonomy
    {
        const auto expect = std::vector {
            face_tag::I_FACE, face_tag::I_FACE, face_tag::I_FACE, face_tag::I_FACE,
            face_tag::J_FACE, face_tag::J_FACE, face_tag::J_FACE, face_tag::J_FACE,
            face_tag::K_FACE, face_tag::K_FACE, face_tag::K_FACE,
        };

        BOOST_CHECK_EQUAL_COLLECTIONS(out.face_tag,
                                      out.face_tag + out.number_of_faces,
                                      expect.begin(), expect.end());
    }
}

BOOST_AUTO_TEST_CASE(Unit_Thickness_Middle_No_Pinch)
{
    auto testCase = unitThicknessMiddle()
        .pinchActive(false)
        .ztol(0.0);

    testCase.process();
    BOOST_REQUIRE_EQUAL(testCase.status(), 1);

    const auto& out = testCase.grid();

    BOOST_CHECK_EQUAL(out.dimensions[0], 1);
    BOOST_CHECK_EQUAL(out.dimensions[1], 1);
    BOOST_CHECK_EQUAL(out.dimensions[2], 3);

    BOOST_CHECK_EQUAL(out.number_of_nodes, 16);
    BOOST_CHECK_EQUAL(out.number_of_faces, 12);
    BOOST_CHECK_EQUAL(out.number_of_cells,  2);

    // Local-to-global cell mapping
    {
        const auto expect = std::vector { 0, 2 };
        BOOST_CHECK_EQUAL_COLLECTIONS(out.local_cell_index,
                                      out.local_cell_index + out.number_of_cells,
                                      expect.begin(), expect.end());
    }

    // Face-to-cell mapping
    {
        const auto expect = std::vector {
            -1, 0,              // I-min (top)
            -1, 1,              // I-min (bot)
            0, -1,              // I-max (top)
            1, -1,              // I-max (bot)

            // -------------------------------

            -1, 0,              // J-min (top)
            -1, 1,              // J-min (bot)
            0, -1,              // J-max (top)
            1, -1,              // J-max (bot)

            // -------------------------------

            -1, 0,              // K-min (top)
            0, -1,              // K-max (top)
            -1, 1,              // K-min (bot)
            1, -1,              // K-max (bot)
        };

        BOOST_CHECK_EQUAL_COLLECTIONS(out.face_neighbors,
                                      out.face_neighbors + 2*out.number_of_faces,
                                      expect.begin(), expect.end());
    }

    // Face taxonomy
    {
        const auto expect = std::vector {
            face_tag::I_FACE, face_tag::I_FACE, face_tag::I_FACE, face_tag::I_FACE,
            face_tag::J_FACE, face_tag::J_FACE, face_tag::J_FACE, face_tag::J_FACE,
            face_tag::K_FACE, face_tag::K_FACE, face_tag::K_FACE, face_tag::K_FACE,
        };

        BOOST_CHECK_EQUAL_COLLECTIONS(out.face_tag,
                                      out.face_tag + out.number_of_faces,
                                      expect.begin(), expect.end());
    }
}

BOOST_AUTO_TEST_CASE(Unit_Thickness_Middle_With_Pinch)
{
    auto testCase = unitThicknessMiddle()
        .pinchActive(true)
        .ztol(0.0);

    testCase.process();
    BOOST_REQUIRE_EQUAL(testCase.status(), 1);

    const auto& out = testCase.grid();

    BOOST_CHECK_EQUAL(out.dimensions[0], 1);
    BOOST_CHECK_EQUAL(out.dimensions[1], 1);
    BOOST_CHECK_EQUAL(out.dimensions[2], 3);

    BOOST_CHECK_EQUAL(out.number_of_nodes, 16);
    BOOST_CHECK_EQUAL(out.number_of_faces, 12);
    BOOST_CHECK_EQUAL(out.number_of_cells,  2);

    // Local-to-global cell mapping
    {
        const auto expect = std::vector { 0, 2 };
        BOOST_CHECK_EQUAL_COLLECTIONS(out.local_cell_index,
                                      out.local_cell_index + out.number_of_cells,
                                      expect.begin(), expect.end());
    }

    // Face-to-cell mapping
    {
        const auto expect = std::vector {
            -1, 0,              // I-min (top)
            -1, 1,              // I-min (bot)
            0, -1,              // I-max (top)
            1, -1,              // I-max (bot)

            // -------------------------------

            -1, 0,              // J-min (top)
            -1, 1,              // J-min (bot)
            0, -1,              // J-max (top)
            1, -1,              // J-max (bot)

            // -------------------------------

            -1, 0,              // K-min (top)
            0, -1,              // K-max (top)
            -1, 1,              // K-min (bot)
            1, -1,              // K-max (bot)
        };

        BOOST_CHECK_EQUAL_COLLECTIONS(out.face_neighbors,
                                      out.face_neighbors + 2*out.number_of_faces,
                                      expect.begin(), expect.end());
    }

    // Face taxonomy
    {
        const auto expect = std::vector {
            face_tag::I_FACE, face_tag::I_FACE, face_tag::I_FACE, face_tag::I_FACE,
            face_tag::J_FACE, face_tag::J_FACE, face_tag::J_FACE, face_tag::J_FACE,
            face_tag::K_FACE, face_tag::K_FACE, face_tag::K_FACE, face_tag::K_FACE,
        };

        BOOST_CHECK_EQUAL_COLLECTIONS(out.face_tag,
                                      out.face_tag + out.number_of_faces,
                                      expect.begin(), expect.end());
    }
}

BOOST_AUTO_TEST_CASE(Single_Fault_Connection)
{
    auto testCase = singleFaultConnection()
        .pinchActive(true)
        .ztol(0.0);

    testCase.process();
    BOOST_REQUIRE_EQUAL(testCase.status(), 1);

    const auto& out = testCase.grid();

    BOOST_CHECK_EQUAL(out.dimensions[0], 2);
    BOOST_CHECK_EQUAL(out.dimensions[1], 1);
    BOOST_CHECK_EQUAL(out.dimensions[2], 1);

    BOOST_CHECK_EQUAL(out.number_of_nodes, 15);
    BOOST_CHECK_EQUAL(out.number_of_faces, 13);
    BOOST_CHECK_EQUAL(out.number_of_cells,  2);

    // Face-to-cell mapping
    {
        const auto expect = std::vector {
            -1, 0,              // 0, I-min (left)
            -1, 1,              // 1, I-min (right)
            0, -1,              // 2, I-max (1, left)
            0,  1,              // 3, I-max (2, left)
            1, -1,              // 4, I-max (right)

            // -------------------------------

            -1, 0,              // 5, J-min (left)
            -1, 1,              // 6, J-min (right)
            0, -1,              // 7, J-max (left)
            1, -1,              // 8, J-max (right)

            // -------------------------------

            -1, 0,              //  9, K-min (left)
            0, -1,              // 10, K-max (left)
            -1, 1,              // 11, K-min (right)
            1, -1,              // 12, K-max (right)
        };

        BOOST_CHECK_EQUAL_COLLECTIONS(out.face_neighbors,
                                      out.face_neighbors + 2*out.number_of_faces,
                                      expect.begin(), expect.end());
    }

    // Face taxonomy
    {
        const auto expect = std::vector {
            face_tag::I_FACE, face_tag::I_FACE, face_tag::I_FACE, face_tag::I_FACE, face_tag::I_FACE,
            face_tag::J_FACE, face_tag::J_FACE, face_tag::J_FACE, face_tag::J_FACE,
            face_tag::K_FACE, face_tag::K_FACE, face_tag::K_FACE, face_tag::K_FACE,
        };

        BOOST_CHECK_EQUAL_COLLECTIONS(out.face_tag,
                                      out.face_tag + out.number_of_faces,
                                      expect.begin(), expect.end());
    }

    // Vertex coordinates
    {
        const auto expect = std::array {
            // Back left pillar
            std::array { 0.0, 0.0, 1.0 }, // 0
            std::array { 0.0, 0.0, 2.0 }, // 1

            // Back centre pillar
            std::array { 1.0, 0.0, 1.0 }, // 2
            std::array { 1.0, 0.0, 1.5 }, // 3
            std::array { 1.0, 0.0, 2.0 }, // 4

            // Back right pillar
            std::array { 2.0, 0.0, 1.5 }, // 5
            std::array { 2.0, 0.0, 2.0 }, // 6

            // -------------------------------

            // Front left pillar
            std::array { 0.0, 1.0, 1.0 }, // 7
            std::array { 0.0, 1.0, 2.0 }, // 8

            // Front centre pillar
            std::array { 1.0, 1.0, 0.5 }, //  9
            std::array { 1.0, 1.0, 1.0 }, // 10
            std::array { 1.0, 1.0, 2.0 }, // 11

            // Front right pillar
            std::array { 2.0, 1.0, 1.0 }, // 12
            std::array { 2.0, 1.0, 2.0 }, // 13

            // --------------------------------

            // Fault node
            std::array { 1.0, 0.5, 1.0 }, // 14
        };

        BOOST_REQUIRE_EQUAL(expect.size(), static_cast<std::size_t>(out.number_of_nodes));

        auto vertexId = 0*expect.size();
        const auto* coord = out.node_coordinates;
        for (const auto& vertex : expect) {
            auto coordCompId = 0*vertex.size();
            for (const auto& component : vertex) {
                BOOST_TEST_MESSAGE("Vertex " << vertexId <<
                                   ", component " << coordCompId);
                BOOST_CHECK_CLOSE(*coord++, component, 1.0e-8);

                ++coordCompId;
            }

            ++vertexId;
        }
    }

    // Face vertices (offsets)
    {
        const auto expect = std::array {
            // I faces (0..4)
            0, 4, 7, 10, 15,

            // J faces (5..9)
            19, 23, 27, 31,

            // K faces (10..13)
            35, 39, 43, 47,

            // End
            51,
        };

        BOOST_CHECK_EQUAL_COLLECTIONS(out.face_node_ptr, out.face_node_ptr + out.number_of_faces + 1,
                                      expect.begin(), expect.end());
    }

    // Face vertices (vertex IDs)
    {
        const auto expect = std::array {
            0, 7, 8, 1,         // 0, I-min (left)
            14, 9, 10,          // 1, I-min (right)
            2, 14, 3,           // 2, I-max (1, left)
            3, 14, 10, 11, 4,   // 3, I-max (2, left)
            5, 12, 13, 6,       // 4, I-max (right)

            // --------------------------------------

            2, 0, 1, 4,         // 5, J-min (left)
            5, 3, 4, 6,         // 6, J-max (left)
            10, 7, 8, 11,       // 7, J-min (right)
            12, 9, 11, 13,      // 8, J-max (right)

            // --------------------------------------

            0, 2, 10, 7,        //  9, K-min (left)
            1, 4, 11, 8,        // 10, K-max (left)
            3, 5, 12, 9,        // 11, K-min (right)
            4, 6, 13, 11,       // 12, K-max (right)
        };

        BOOST_REQUIRE_EQUAL(expect.size(), static_cast<std::size_t>(out.face_node_ptr[out.number_of_faces]));

        BOOST_CHECK_EQUAL_COLLECTIONS(out.face_nodes, out.face_nodes + out.face_node_ptr[out.number_of_faces],
                                      expect.begin(), expect.end());
    }
}

BOOST_AUTO_TEST_SUITE_END()     // Regular_Processing

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Edge_Conformal_Processing)

BOOST_AUTO_TEST_CASE(Single_Fault_Connection)
{
    auto testCase = singleFaultConnection()
        .edgeConformal(true)
        .pinchActive(true)
        .ztol(0.0);

    testCase.process();
    BOOST_REQUIRE_EQUAL(testCase.status(), 1);

    const auto& out = testCase.grid();

    BOOST_CHECK_EQUAL(out.dimensions[0], 2);
    BOOST_CHECK_EQUAL(out.dimensions[1], 1);
    BOOST_CHECK_EQUAL(out.dimensions[2], 1);

    BOOST_CHECK_EQUAL(out.number_of_nodes, 15);
    BOOST_CHECK_EQUAL(out.number_of_faces, 13);
    BOOST_CHECK_EQUAL(out.number_of_cells,  2);

    // Face-to-cell mapping
    {
        const auto expect = std::vector {
            -1, 0,              // 0, I-min (left)
            -1, 1,              // 1, I-min (right)
            0, -1,              // 2, I-max (1, left)
            0,  1,              // 3, I-max (2, left)
            1, -1,              // 4, I-max (right)

            // -------------------------------

            -1, 0,              // 5, J-min (left)
            -1, 1,              // 6, J-min (right)
            0, -1,              // 7, J-max (left)
            1, -1,              // 8, J-max (right)

            // -------------------------------

            -1, 0,              //  9, K-min (left)
            0, -1,              // 10, K-max (left)
            -1, 1,              // 11, K-min (right)
            1, -1,              // 12, K-max (right)
        };

        BOOST_CHECK_EQUAL_COLLECTIONS(out.face_neighbors,
                                      out.face_neighbors + 2*out.number_of_faces,
                                      expect.begin(), expect.end());
    }

    // Face taxonomy
    {
        const auto expect = std::vector {
            face_tag::I_FACE, face_tag::I_FACE, face_tag::I_FACE, face_tag::I_FACE, face_tag::I_FACE,
            face_tag::J_FACE, face_tag::J_FACE, face_tag::J_FACE, face_tag::J_FACE,
            face_tag::K_FACE, face_tag::K_FACE, face_tag::K_FACE, face_tag::K_FACE,
        };

        BOOST_CHECK_EQUAL_COLLECTIONS(out.face_tag,
                                      out.face_tag + out.number_of_faces,
                                      expect.begin(), expect.end());
    }

    // Vertex coordinates
    {
        const auto expect = std::array {
            // Back left pillar
            std::array { 0.0, 0.0, 1.0 }, // 0
            std::array { 0.0, 0.0, 2.0 }, // 1

            // Back centre pillar
            std::array { 1.0, 0.0, 1.0 }, // 2
            std::array { 1.0, 0.0, 1.5 }, // 3
            std::array { 1.0, 0.0, 2.0 }, // 4

            // Back right pillar
            std::array { 2.0, 0.0, 1.5 }, // 5
            std::array { 2.0, 0.0, 2.0 }, // 6

            // -------------------------------

            // Front left pillar
            std::array { 0.0, 1.0, 1.0 }, // 7
            std::array { 0.0, 1.0, 2.0 }, // 8

            // Front centre pillar
            std::array { 1.0, 1.0, 0.5 }, //  9
            std::array { 1.0, 1.0, 1.0 }, // 10
            std::array { 1.0, 1.0, 2.0 }, // 11

            // Front right pillar
            std::array { 2.0, 1.0, 1.0 }, // 12
            std::array { 2.0, 1.0, 2.0 }, // 13

            // --------------------------------

            // Fault node
            std::array { 1.0, 0.5, 1.0 }, // 14
        };

        BOOST_REQUIRE_EQUAL(expect.size(), static_cast<std::size_t>(out.number_of_nodes));

        auto vertexId = 0*expect.size();
        const auto* coord = out.node_coordinates;
        for (const auto& vertex : expect) {
            auto coordCompId = 0*vertex.size();
            for (const auto& component : vertex) {
                BOOST_TEST_MESSAGE("Vertex " << vertexId <<
                                   ", component " << coordCompId);
                BOOST_CHECK_CLOSE(*coord++, component, 1.0e-8);

                ++coordCompId;
            }

            ++vertexId;
        }
    }

    // Face vertices (offsets)
    {
        const auto expect = std::array {
            // I faces (0..4)
            0, 4, 7, 10, 15,

            // J faces (5..9)
            19, 23, 27, 31,

            // K faces (10..13)
            35, 40, 44, 49,

            // End
            53,
        };

        BOOST_CHECK_EQUAL_COLLECTIONS(out.face_node_ptr, out.face_node_ptr + out.number_of_faces + 1,
                                      expect.begin(), expect.end());
    }

    // Face vertices (vertex IDs)
    {
        const auto expect = std::array {
            0, 7, 8, 1,         // 0, I-min (left)
            14, 9, 10,          // 1, I-min (right)
            2, 14, 3,           // 2, I-max (1, left)
            3, 14, 10, 11, 4,   // 3, I-max (2, left)
            5, 12, 13, 6,       // 4, I-max (right)

            // --------------------------------------

            2, 0, 1, 4,         // 5, J-min (left)
            5, 3, 4, 6,         // 6, J-max (left)
            10, 7, 8, 11,       // 7, J-min (right)
            12, 9, 11, 13,      // 8, J-max (right)

            // --------------------------------------

            0, 2, 14, 10, 7,    //  9, K-min (left)
            1, 4, 11, 8,        // 10, K-max (left)
            3, 5, 12, 9, 14,    // 11, K-min (right)
            4, 6, 13, 11,       // 12, K-max (right)
        };

        BOOST_REQUIRE_EQUAL(expect.size(), static_cast<std::size_t>(out.face_node_ptr[out.number_of_faces]));

        BOOST_CHECK_EQUAL_COLLECTIONS(out.face_nodes, out.face_nodes + out.face_node_ptr[out.number_of_faces],
                                      expect.begin(), expect.end());
    }
}

BOOST_AUTO_TEST_SUITE_END()     // Edge_Conformal_Processing
