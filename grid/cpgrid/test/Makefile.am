# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

check_PROGRAMS = entityrep_test orientedentitytable_test entity_test
noinst_PROGRAMS = buildcpgrid_test

entity_test_SOURCES = entity_test.cpp
entity_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
entity_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) \
                    $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

entityrep_test_SOURCES = entityrep_test.cpp
entityrep_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
entityrep_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) \
                    $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

orientedentitytable_test_SOURCES = orientedentitytable_test.cpp
orientedentitytable_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
orientedentitytable_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) \
                    $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)


buildcpgrid_test_SOURCES = BuildCpGrid.hpp buildcpgrid_test.cpp
buildcpgrid_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
buildcpgrid_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS) \
                  $(BOOST_LDFLAGS) $(BOOST_SYSTEM_LIB) \
                  ../libcpgrid.la \
                  ../../common/libtemporaryhack.la \
                  ../../preprocess/libpreprocess.la \
                  ../../../common/param/libparam.la

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules

