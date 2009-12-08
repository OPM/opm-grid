# $Date$
# $Revision$

check_PROGRAMS = sparsetable_test sparsevector_test monotcubicinterpolator_test
noinst_PROGRAMS = unit_test

sparsetable_test_SOURCES = sparsetable_test.cpp
sparsetable_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
sparsetable_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) \
                         $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

sparsevector_test_SOURCES = sparsevector_test.cpp
sparsevector_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
sparsevector_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) \
                          $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

monotcubicinterpolator_test_SOURCES = monotcubicinterpolator_test.cpp
monotcubicinterpolator_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
monotcubicinterpolator_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS) ../libcommon.la

unit_test_SOURCES = unit_test.cpp
unit_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
unit_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS)

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules

