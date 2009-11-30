# $Date$
# $Revision$

check_PROGRAMS  = entityrep_test orientedentitytable_test entity_test
#noinst_PROGRAMS = buildcpgrid_test

AM_CXXFLAGS = $(DUNE_CPPFLAGS) $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
AM_LDFLAGS  = $(DUNE_LDFLAGS) $(MPI_LDFLAGS) $(BOOST_LDFLAGS)

entity_test_SOURCES  = entity_test.cpp
entity_test_LDADD    = $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

entityrep_test_SOURCES  = entityrep_test.cpp
entityrep_test_LDADD    = $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

orientedentitytable_test_SOURCES = orientedentitytable_test.cpp
orientedentitytable_test_LDADD   = $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

#buildcpgrid_test_SOURCES  = BuildCpGrid.hpp buildcpgrid_test.cpp
#buildcpgrid_test_LDADD    = $(DUNE_LIBS) $(BOOST_SYSTEM_LIB) \
#                            ../../lib/libdunecornerpoint.la

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
