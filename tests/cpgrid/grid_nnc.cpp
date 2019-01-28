/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#include <config.h>

#define NVERBOSE

#define BOOST_TEST_MODULE CpGridWithNNC

#include <boost/test/unit_test.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/grid/CpGrid.hpp>

struct Fixture
{
    Fixture()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(m_argc, m_argv);
    }
};

BOOST_AUTO_TEST_SUITE(ConstructingWithNNC)

BOOST_FIXTURE_TEST_CASE(Construct, Fixture)
{
#if HAVE_ECL_INPUT
    Opm::Parser parser;
    Opm::EclipseState es(parser.parseFile("CORNERPOINT_ACTNUM.DATA"));
    Dune::CpGrid grid;
    grid.processEclipseFormat(es.getInputGrid(), false, false, false, {}, es.getInputNNC());
#endif
}

BOOST_AUTO_TEST_SUITE_END()
