/*
  Copyright 2017 SINTEF ICT, Applied Mathematics.
  Copyright 2017 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#include <config.h>

#define NVERBOSE  // Suppress own messages when throw()ing

#define BOOST_TEST_MODULE TEST_RepairZCORN

#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif

/* --- our own headers --- */

#include <opm/grid/RepairZCORN.hpp>

#include <algorithm>
#include <cstddef>
#include <vector>

namespace {
    template <class Coll1, class Coll2>
    void check_is_close(const Coll1& c1, const Coll2& c2)
    {
        BOOST_REQUIRE_EQUAL(c1.size(), c2.size());

        if (! c1.empty()) {
            auto i1 = c1.begin(), e1 = c1.end();
            auto i2 = c2.begin();

            for (; i1 != e1; ++i1, ++i2) {
                BOOST_CHECK_CLOSE(*i1, *i2, 1.0e-10);
            }
        }
    }
} // Namespace anonymous

BOOST_AUTO_TEST_SUITE (Repair_AllActive)

BOOST_AUTO_TEST_CASE (NoChange)
{
    const auto cartDims = std::vector<int>{ 1, 1, 2 };
    const auto actnum   = std::vector<int>{};  // empty => all active

    auto zcorn = std::vector<double> {
        0.0, 0.0,
        0.0, 0.0,
        1.0, 1.0,
        1.0, 1.0,

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    const auto expect = zcorn;

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), false);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{0});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{0});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_CASE (TopBelowBottom)
{
    const auto cartDims = std::vector<int>{ 1, 1, 2 };
    const auto actnum   = std::vector<int>{};  // empty => all active

    auto zcorn = std::vector<double> {
        0.0, 0.0,
        0.0, 1.5,
        1.0, 1.0,
        1.0, 1.0,

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    const auto expect = std::vector<double> {
        0.0, 0.0,
        0.0, 1.0,               // Top <- Bottom at c(1,1)
        1.0, 1.0,
        1.0, 1.0,

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), false);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{1});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{0});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_CASE (BottomBelowLowerTop)
{
    const auto cartDims = std::vector<int>{ 1, 1, 2 };
    const auto actnum   = std::vector<int>{};  // empty => all active

    auto zcorn = std::vector<double> {
        0.0, 0.0,
        0.0, 1.0,
        1.0, 1.0,
        1.0, 1.5,               // Bottom below lower top at c(1,1)

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    const auto expect = std::vector<double> {
        0.0, 0.0,
        0.0, 1.0,
        1.0, 1.0,
        1.0, 1.5,

        1.0, 1.0,
        1.0, 1.5,               // Lower Top <- Upper Bottom at c(1,1)
        2.0, 2.0,
        2.0, 2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), false);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{0});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{1});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_CASE (TBB_And_BBLT)
{
    const auto cartDims = std::vector<int>{ 1, 1, 2 };
    const auto actnum   = std::vector<int>{};  // empty => all active

    auto zcorn = std::vector<double> {
        0.0, 0.0,
        0.0, 1.6,
        1.0, 1.0,
        1.0, 1.5,               // Bottom below lower top at c(1,1)

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    const auto expect = std::vector<double> {
        0.0, 0.0,
        0.0, 1.5,
        1.0, 1.0,
        1.0, 1.5,

        1.0, 1.0,
        1.0, 1.5,               // Lower Top <- Upper Bottom at c(1,1)
        2.0, 2.0,
        2.0, 2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), false);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{1});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{1});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_CASE (Elevation)
{
    const auto cartDims = std::vector<int>{ 1, 1, 2 };
    const auto actnum   = std::vector<int>{};  // empty => all active

    auto zcorn = std::vector<double> {
         0.0,  0.0,
         0.0,  0.0,
        -1.0, -1.0,
        -1.0, -1.0,

        -1.0, -1.0,
        -1.0, -1.0,
        -2.0, -2.0,
        -2.0, -2.0,
    };

    const auto expect = std::vector<double> {
        0.0, 0.0,
        0.0, 0.0,
        1.0, 1.0,
        1.0, 1.0,

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), true);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{0});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{0});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_CASE (TBB_And_Elevation)
{
    const auto cartDims = std::vector<int>{ 1, 1, 2 };
    const auto actnum   = std::vector<int>{};  // empty => all active

    auto zcorn = std::vector<double> {
         0.0,  0.0,
         0.0, -1.5,
        -1.0, -1.0,
        -1.0, -1.0,

        -1.0, -1.0,
        -1.0, -1.0,
        -2.0, -2.0,
        -2.0, -2.0,
    };

    const auto expect = std::vector<double> {
        0.0, 0.0,
        0.0, 1.0,               // Top <- Bottom at c(1,1)
        1.0, 1.0,
        1.0, 1.0,

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), true);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{1});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{0});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_CASE (BBLT_And_Elevation)
{
    const auto cartDims = std::vector<int>{ 1, 1, 2 };
    const auto actnum   = std::vector<int>{};  // empty => all active

    auto zcorn = std::vector<double> {
         0.0,  0.0,
         0.0, -1.0,
        -1.0, -1.0,
        -1.0, -1.5,             // Bottom below lower top at c(1,1)

        -1.0, -1.0,
        -1.0, -1.0,
        -2.0, -2.0,
        -2.0, -2.0,
    };

    const auto expect = std::vector<double> {
        0.0, 0.0,
        0.0, 1.0,
        1.0, 1.0,
        1.0, 1.5,

        1.0, 1.0,
        1.0, 1.5,               // Lower Top <- Upper Bottom at c(1,1)
        2.0, 2.0,
        2.0, 2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), true);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{0});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{1});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_CASE (TBB_And_BBLT_And_Elevation)
{
    const auto cartDims = std::vector<int>{ 1, 1, 2 };
    const auto actnum   = std::vector<int>{};  // empty => all active

    auto zcorn = std::vector<double> {
         0.0,  0.0,
         0.0, -1.6,             // Top below bottom at c(1,1)
        -1.0, -1.0,
        -1.0, -1.5,             // Bottom below lower top at c(1,1)

        -1.0, -1.0,
        -1.0, -1.0,
        -2.0, -2.0,
        -2.0, -2.0,
    };

    const auto expect = std::vector<double> {
        0.0, 0.0,
        0.0, 1.5,               // Top <- Bottom at c(1,1)
        1.0, 1.0,
        1.0, 1.5,

        1.0, 1.0,
        1.0, 1.5,               // Lower Top <- Upper Bottom at c(1,1)
        2.0, 2.0,
        2.0, 2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), true);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{1});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{1});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_SUITE_END()

// ======================================================================

BOOST_AUTO_TEST_SUITE (Repair_With_Inactive)

BOOST_AUTO_TEST_CASE (NoChange)
{
    const auto cartDims = std::vector<int>{ 1, 1, 3 };
    const auto actnum   = std::vector<int>{ 1, 0, 1 };

    auto zcorn = std::vector<double> {
          0.0,  0.0,
          0.0,  0.0,
          1.0,  1.0,
          1.0,  1.0,

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

          1.0,  1.0,
          1.0,  1.0,
          2.0,  2.0,
          2.0,  2.0,
    };

    const auto expect = zcorn;

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), false);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{0});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{0});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_CASE (TopBelowBottom)
{
    const auto cartDims = std::vector<int>{ 1, 1, 3 };
    const auto actnum   = std::vector<int>{ 1, 0, 1 };

    auto zcorn = std::vector<double> {
        0.0, 0.0,
        0.0, 1.5,
        1.0, 1.0,
        1.0, 1.0,

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    const auto expect = std::vector<double> {
        0.0, 0.0,
        0.0, 1.0,               // Top <- Bottom at c(1,1)
        1.0, 1.0,
        1.0, 1.0,

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), false);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{1});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{0});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_CASE (BottomBelowLowerTop)
{
    const auto cartDims = std::vector<int>{ 1, 1, 3 };
    const auto actnum   = std::vector<int>{ 1, 0, 1 };

    auto zcorn = std::vector<double> {
        0.0, 0.0,
        0.0, 1.0,
        1.0, 1.0,
        1.0, 1.5,               // Bottom below lower top at c(1,1)

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    const auto expect = std::vector<double> {
        0.0, 0.0,
        0.0, 1.0,
        1.0, 1.0,
        1.0, 1.5,

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        1.0, 1.0,
        1.0, 1.5,               // Lower Top <- Upper Bottom at c(1,1)
        2.0, 2.0,
        2.0, 2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), false);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{0});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{1});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_CASE (TBB_And_BBLT)
{
    const auto cartDims = std::vector<int>{ 1, 1, 3 };
    const auto actnum   = std::vector<int>{ 1, 0, 1 };

    auto zcorn = std::vector<double> {
        0.0, 0.0,
        0.0, 1.6,
        1.0, 1.0,
        1.0, 1.5,               // Bottom below lower top at c(1,1)

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    const auto expect = std::vector<double> {
        0.0, 0.0,
        0.0, 1.5,
        1.0, 1.0,
        1.0, 1.5,

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        1.0, 1.0,
        1.0, 1.5,               // Lower Top <- Upper Bottom at c(1,1)
        2.0, 2.0,
        2.0, 2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), false);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{1});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{1});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_CASE (Elevation)
{
    const auto cartDims = std::vector<int>{ 1, 1, 3 };
    const auto actnum   = std::vector<int>{ 1, 0, 1 };

    auto zcorn = std::vector<double> {
         0.0,  0.0,
         0.0,  0.0,
        -1.0, -1.0,
        -1.0, -1.0,

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        -1.0, -1.0,
        -1.0, -1.0,
        -2.0, -2.0,
        -2.0, -2.0,
    };

    const auto expect = std::vector<double> {
        0.0, 0.0,
        0.0, 0.0,
        1.0, 1.0,
        1.0, 1.0,

          0.0,  1.0,
       - 10.0,  0.0,
          2.0,  0.0,
          0.0,  0.0,

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), true);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{0});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{0});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_CASE (TBB_And_Elevation)
{
    const auto cartDims = std::vector<int>{ 1, 1, 3 };
    const auto actnum   = std::vector<int>{ 1, 0, 1 };

    auto zcorn = std::vector<double> {
         0.0,  0.0,
         0.0, -1.5,
        -1.0, -1.0,
        -1.0, -1.0,

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        -1.0, -1.0,
        -1.0, -1.0,
        -2.0, -2.0,
        -2.0, -2.0,
    };

    const auto expect = std::vector<double> {
        0.0, 0.0,
        0.0, 1.0,               // Top <- Bottom at c(1,1)
        1.0, 1.0,
        1.0, 1.0,

          0.0,  1.0,
       - 10.0,  0.0,
          2.0,  0.0,
          0.0,  0.0,

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), true);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{1});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{0});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_CASE (BBLT_And_Elevation)
{
    const auto cartDims = std::vector<int>{ 1, 1, 3 };
    const auto actnum   = std::vector<int>{ 1, 0, 1 };

    auto zcorn = std::vector<double> {
         0.0,  0.0,
         0.0, -1.0,
        -1.0, -1.0,
        -1.0, -1.5,             // Bottom below lower top at c(1,1)

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        -1.0, -1.0,
        -1.0, -1.0,
        -2.0, -2.0,
        -2.0, -2.0,
    };

    const auto expect = std::vector<double> {
        0.0, 0.0,
        0.0, 1.0,
        1.0, 1.0,
        1.0, 1.5,

          0.0,  1.0,
       - 10.0,  0.0,
          2.0,  0.0,
          0.0,  0.0,

        1.0, 1.0,
        1.0, 1.5,               // Lower Top <- Upper Bottom at c(1,1)
        2.0, 2.0,
        2.0, 2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), true);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{0});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{1});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_CASE (TBB_And_BBLT_And_Elevation)
{
    const auto cartDims = std::vector<int>{ 1, 1, 3 };
    const auto actnum   = std::vector<int>{ 1, 0, 1 };

    auto zcorn = std::vector<double> {
         0.0,  0.0,
         0.0, -1.6,             // Top below bottom at c(1,1)
        -1.0, -1.0,
        -1.0, -1.5,             // Bottom below lower top at c(1,1)

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        -1.0, -1.0,
        -1.0, -1.0,
        -2.0, -2.0,
        -2.0, -2.0,
    };

    const auto expect = std::vector<double> {
        0.0, 0.0,
        0.0, 1.5,               // Top <- Bottom at c(1,1)
        1.0, 1.0,
        1.0, 1.5,

          0.0,  1.0,
       - 10.0,  0.0,
          2.0,  0.0,
          0.0,  0.0,

        1.0, 1.0,
        1.0, 1.5,               // Lower Top <- Upper Bottom at c(1,1)
        2.0, 2.0,
        2.0, 2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), true);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{1});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{1});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_SUITE_END()

// ======================================================================

BOOST_AUTO_TEST_SUITE (Repair_With_NoBottomNeigh)

BOOST_AUTO_TEST_CASE (NoChange)
{
    const auto cartDims = std::vector<int>{ 1, 1, 3 };
    const auto actnum   = std::vector<int>{ 1, 0, 0 };

    auto zcorn = std::vector<double> {
          0.0,  0.0,
          0.0,  0.0,
          1.0,  1.0,
          1.0,  1.0,

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

          1.0,  1.0,
          1.0,  1.0,
          2.0,  2.0,
          2.0,  2.0,
    };

    const auto expect = zcorn;

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), false);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{0});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{0});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_CASE (TopBelowBottom)
{
    const auto cartDims = std::vector<int>{ 1, 1, 3 };
    const auto actnum   = std::vector<int>{ 1, 0, 0 };

    auto zcorn = std::vector<double> {
        0.0, 0.0,
        0.0, 1.5,
        1.0, 1.0,
        1.0, 1.0,

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    const auto expect = std::vector<double> {
        0.0, 0.0,
        0.0, 1.0,               // Top <- Bottom at c(1,1)
        1.0, 1.0,
        1.0, 1.0,

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), false);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{1});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{0});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_CASE (BottomBelowLowerTop)
{
    const auto cartDims = std::vector<int>{ 1, 1, 3 };
    const auto actnum   = std::vector<int>{ 1, 0, 0 };

    auto zcorn = std::vector<double> {
        0.0, 0.0,
        0.0, 1.0,
        1.0, 1.0,
        1.0, 1.5,               // Bottom below lower top at c(1,1)

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    const auto expect = std::vector<double> {
        0.0, 0.0,
        0.0, 1.0,
        1.0, 1.0,
        1.0, 1.5,

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), false);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{0});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{0});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_CASE (TBB_And_BBLT)
{
    const auto cartDims = std::vector<int>{ 1, 1, 3 };
    const auto actnum   = std::vector<int>{ 1, 0, 0 };

    auto zcorn = std::vector<double> {
        0.0, 0.0,
        0.0, 1.6,
        1.0, 1.0,
        1.0, 1.5,               // Bottom below lower top at c(1,1)

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    const auto expect = std::vector<double> {
        0.0, 0.0,
        0.0, 1.5,
        1.0, 1.0,
        1.0, 1.5,

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), false);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{1});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{0});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_CASE (Elevation)
{
    const auto cartDims = std::vector<int>{ 1, 1, 3 };
    const auto actnum   = std::vector<int>{ 1, 0, 0 };

    auto zcorn = std::vector<double> {
         0.0,  0.0,
         0.0,  0.0,
        -1.0, -1.0,
        -1.0, -1.0,

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        -1.0, -1.0,
        -1.0, -1.0,
        -2.0, -2.0,
        -2.0, -2.0,
    };

    const auto expect = std::vector<double> {
        0.0, 0.0,
        0.0, 0.0,
        1.0, 1.0,
        1.0, 1.0,

          0.0,  1.0,
       - 10.0,  0.0,
          2.0,  0.0,
          0.0,  0.0,

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), true);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{0});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{0});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

// This test (Repair_With_NoBottomNeigh/TBB_And_Elevation) verifies the
// current behaviour of class RepairZCORN.  The sanitised ZCORN array
// results are arguably wrong.
BOOST_AUTO_TEST_CASE (TBB_And_Elevation)
{
    const auto cartDims = std::vector<int>{ 1, 1, 3 };
    const auto actnum   = std::vector<int>{ 1, 0, 0 };

    auto zcorn = std::vector<double> {
         0.0,  0.0,
         0.0, -1.5,
        -1.0, -1.0,
        -1.0, -1.0,

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        -1.0, -1.0,
        -1.0, -1.0,
        -2.0, -2.0,
        -2.0, -2.0,
    };

    const auto expect = std::vector<double> {
        -1.0, -1.0,
        -1.0, -1.5,
        -1.0, -1.0,
        -1.0, -1.0,

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        -1.0, -1.0,
        -1.0, -1.0,
        -2.0, -2.0,
        -2.0, -2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    // Only active input cell is twisted and doesn't allow determining sign
    // of ZCORN delta.  This is a deficiency of the current implementation.
    //
    // Hopefully this case does not occur too often in real input decks.
    BOOST_CHECK_EQUAL(repair.switchedToDepth(), false);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{3});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{0});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_CASE (BBLT_And_Elevation)
{
    const auto cartDims = std::vector<int>{ 1, 1, 3 };
    const auto actnum   = std::vector<int>{ 1, 0, 0 };

    auto zcorn = std::vector<double> {
         0.0,  0.0,
         0.0, -1.0,
        -1.0, -1.0,
        -1.0, -1.5,

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        -1.0, -1.0,
        -1.0, -1.0,
        -2.0, -2.0,
        -2.0, -2.0,
    };

    const auto expect = std::vector<double> {
        0.0, 0.0,
        0.0, 1.0,
        1.0, 1.0,
        1.0, 1.5,

          0.0,  1.0,
       - 10.0,  0.0,
          2.0,  0.0,
          0.0,  0.0,

        1.0, 1.0,
        1.0, 1.0,
        2.0, 2.0,
        2.0, 2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    BOOST_CHECK_EQUAL(repair.switchedToDepth(), true);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{0});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{0});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

// This test (Repair_With_NoBottomNeigh/TBB_And_BBLT_And_Elevation) verifies
// the current behaviour of class RepairZCORN.  The sanitised ZCORN array
// results are arguably wrong.
BOOST_AUTO_TEST_CASE (TBB_And_BBLT_And_Elevation)
{
    const auto cartDims = std::vector<int>{ 1, 1, 3 };
    const auto actnum   = std::vector<int>{ 1, 0, 0 };

    auto zcorn = std::vector<double> {
         0.0,  0.0,
         0.0, -1.6,             // Top below bottom at c(1,1)
        -1.0, -1.0,
        -1.0, -1.5,             // Bottom below lower top at c(1,1)

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        -1.0, -1.0,
        -1.0, -1.0,
        -2.0, -2.0,
        -2.0, -2.0,
    };

    const auto expect = std::vector<double> {
        -1.0, -1.0,
        -1.0, -1.6,
        -1.0, -1.0,
        -1.0, -1.5,

          0.0, -1.0,
         10.0,  0.0,
       -  2.0,  0.0,
          0.0,  0.0,

        -1.0, -1.0,
        -1.0, -1.0,
        -2.0, -2.0,
        -2.0, -2.0,
    };

    auto repair = ::Opm::UgGridHelpers::RepairZCORN{
        std::move(zcorn), actnum, cartDims
    };

    // Only active input cell is twisted and doesn't allow determining sign
    // of ZCORN delta.  This is a deficiency of the current implementation.
    //
    // Hopefully this case does not occur too often in real input decks.
    BOOST_CHECK_EQUAL(repair.switchedToDepth(), false);

    {
        const auto& tbb = repair.statTopBelowBottom();

        BOOST_CHECK_EQUAL(tbb.cells  , std::size_t{1});
        BOOST_CHECK_EQUAL(tbb.corners, std::size_t{3});
    }

    {
        const auto& bblt = repair.statBottomBelowLowerTop();

        BOOST_CHECK_EQUAL(bblt.cells  , std::size_t{0});
        BOOST_CHECK_EQUAL(bblt.corners, std::size_t{0});
    }

    zcorn = repair.destructivelyGrabSanitizedValues();

    check_is_close(zcorn, expect);
}

BOOST_AUTO_TEST_SUITE_END()
