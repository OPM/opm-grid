# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

check_PROGRAMS = sparsetable_test sparsevector_test

sparsetable_test_SOURCES = sparsetable_test.cpp
sparsetable_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
sparsetable_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) \
                         $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

sparsevector_test_SOURCES = sparsevector_test.cpp
sparsevector_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
sparsevector_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) \
                          $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules

