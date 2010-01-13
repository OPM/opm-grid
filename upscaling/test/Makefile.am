# $Date$
# $Revision$

check_PROGRAMS = 
noinst_PROGRAMS = upscaling_test steadystate_test #upscale_perm upscale_relperm

AM_CPPFLAGS += $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS) $(SUPERLU_CPPFLAGS)
AM_LDFLAGS  += $(MPI_LDFLAGS) $(BOOST_LDFLAGS) $(SUPERLU_LDFLAGS)
LDADD = $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
        $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)   \
        ../../lib/libdunecornerpoint.la               \
        $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SUPERLU_LIBS)


upscaling_test_SOURCES = upscaling_test.cpp

steadystate_test_SOURCES = steadystate_test.cpp

#upscale_perm_SOURCES = upscale_perm.C

#upscale_relperm_SOURCES = upscale_relperm.C

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
