# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

check_PROGRAMS = gie_test

gie_test_SOURCES = gie_test.cpp
gie_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
gie_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) \
                    $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules

