# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

check_PROGRAMS = entity_test
noinst_PROGRAMS = BuildCpGrid

entity_test_SOURCES = entity_test.cpp
entity_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
entity_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) \
                    $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

BuildCpGrid_SOURCES = BuildCpGrid.cpp
BuildCpGrid_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
BuildCpGrid_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS) \
                  ../libcpgrid.la \
                  ../../common/libtemporaryhack.la \
                  ../../preprocess/libpreprocess.la \
                  ../../../common/param/libparam.la

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules

