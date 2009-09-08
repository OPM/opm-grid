# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

check_PROGRAMS =
noinst_PROGRAMS = mimetic_ipeval_test \
                  mimetic_solver_test \
                  mimetic_periodic_test \
		  parsolver_test

LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) \
        $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)   \
        $(BOOST_SYSTEM_LIB)                             \
        ../../../grid/cpgrid/libcpgrid.la               \
        ../../../grid/preprocess/libpreprocess.la       \
        ../../../grid/common/libtemporaryhack.la        \
        ../../../common/param/libparam.la -llapack -lblas

mimetic_ipeval_test_SOURCES  = mimetic_ipeval_test.cpp
mimetic_ipeval_test_CXXFLAGS = $(DUNE_CPPFLAGS) $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)

mimetic_solver_test_SOURCES  = mimetic_solver_test.cpp
mimetic_solver_test_CXXFLAGS = $(DUNE_CPPFLAGS) $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)

mimetic_periodic_test_SOURCES  = mimetic_periodic_test.cpp
mimetic_periodic_test_CXXFLAGS = $(DUNE_CPPFLAGS) $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)

parsolver_test_SOURCES  = parsolver_test.cpp
parsolver_test_CXXFLAGS = $(DUNE_CPPFLAGS) $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
