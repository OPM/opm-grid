# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

check_PROGRAMS = partition_test
noinst_PROGRAMS =

LDADD =  $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) \
         $(DUNE_LIBS)  $(BOOST_SYSTEM_LIB) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
         ../../cpgrid/libcpgrid.la \
         ../../preprocess/libpreprocess.la \
         ../../../common/param/libparam.la \
         ../../common/libtemporaryhack.la \
         ../../../common/libcommon_hack.la


partition_test_SOURCES = partition_test.cpp
partition_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
