# $Date$
# $Revision$

check_PROGRAMS =
noinst_PROGRAMS = multiscale_test

AM_CPPFLAGS += $(DUNEMPICPPFLAGS) $(BOOST_CPPFLAGS) $(SUPERLU_CPPFLAGS)
AM_LDFLAGS  += $(DUNEMPILDFLAGS) $(BOOST_LDFLAGS) $(SUPERLU_LDFLAGS)
LDADD = $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
        $(BOOST_SYSTEM_LIB) $(DUNEMPILIBS)            \
        ../../../lib/libdunecornerpoint.la            \
        $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SUPERLU_LIBS)

multiscale_test_SOURCES  = multiscale_test.cpp

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
