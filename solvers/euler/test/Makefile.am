# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

noinst_PROGRAMS = euler_upstream_test

euler_upstream_test_SOURCES = euler_upstream_test.cpp
euler_upstream_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
euler_upstream_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) \
                    $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
                    $(BOOST_SYSTEM_LIB) \
                    ../../../grid/cpgrid/libcpgrid.la \
                    ../../../grid/preprocess/libpreprocess.la \
                    ../../../grid/common/libtemporaryhack.la \
                    ../../../common/param/libparam.la -llapack -lblas

#TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
