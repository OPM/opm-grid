# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

check_PROGRAMS =
noinst_PROGRAMS = mimetic_ipeval_test

mimetic_ipeval_test_SOURCES  = mimetic_ipeval_test.cpp
mimetic_ipeval_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
mimetic_ipeval_test_LDADD    = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) \
                               $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)   \
                               $(BOOST_SYSTEM_LIB)                             \
                               ../../../grid/cpgrid/libcpgrid.la               \
                               ../../../grid/preprocess/libpreprocess.la       \
                               ../../../common/param/libparam.la -llapack -lblas

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules

