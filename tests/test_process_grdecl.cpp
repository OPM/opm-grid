#include <config.h>

#define NVERBOSE // to suppress our messages when throwing

#define BOOST_TEST_MODULE Process_Grdecl
#include <boost/test/unit_test.hpp>

#include <opm/grid/cpgpreprocess/preprocess.h>

#include <array>
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

        TestGrid& process()
        {
            auto input = grdecl{};

            std::copy(this->dims_.begin(), this->dims_.end(), input.dims);

            input.coord  = this->coord_.data();
            input.zcorn  = this->zcorn_.data();
            input.actnum = this->actnum_.data();

            this->g_.emplace(processed_grid{});

            this->status_ =
                process_grdecl(&input,
                               this->ztol_,
                               nullptr,
                               &*this->g_,
                               this->pinch_active_);

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
} // Anonymous namespace

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
