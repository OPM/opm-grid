# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

check_PROGRAMS =
noinst_PROGRAMS = multiscale_test

LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) \
        $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)   \
        $(BOOST_SYSTEM_LIB)                             \
        ../../../grid/cpgrid/libcpgrid.la               \
        ../../../grid/preprocess/libpreprocess.la       \
        ../../../grid/common/libtemporaryhack.la        \
        ../../../common/param/libparam.la -llapack -lblas

multiscale_test_SOURCES  = multiscale_test.cpp
multiscale_test_CXXFLAGS = $(DUNE_CPPFLAGS) $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
