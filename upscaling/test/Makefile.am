# $Date$
# $Revision$

check_PROGRAMS = 
noinst_PROGRAMS = upscaling_test steadystate_test

AM_CPPFLAGS += $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
AM_LDFLAGS  += $(MPI_LDFLAGS) $(BOOST_LDFLAGS)
LDADD = $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
        $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)   \
        ../../lib/libdunecornerpoint.la               \
        $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) 


upscaling_test_SOURCES = upscaling_test.cpp

steadystate_test_SOURCES = steadystate_test.cpp

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
