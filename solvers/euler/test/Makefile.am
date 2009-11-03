# $Date$
# $Revision: duneproject 5489 2009-03-25 11:19:24Z sander $

noinst_PROGRAMS = euler_upstream_test

AM_CPPFLAGS += $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
AM_LDFLAGS  += $(MPI_LDFLAGS) $(BOOST_LDFLAGS)

LDADD = $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
        $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)   \
        ../../../lib/libdunecornerpoint.la            \
        ../../../grid/preprocess/libpreprocess.la     \
        -llapack -lblas

euler_upstream_test_SOURCES = euler_upstream_test.cpp

#TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
