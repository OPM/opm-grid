# $Date$
# $Revision$

check_PROGRAMS =
noinst_PROGRAMS = mimetic_ipeval_test \
                  mimetic_solver_test \
                  mimetic_aniso_solver_test \
                  mimetic_periodic_test \
		  parsolver_test \
		  istl_test \
		  spe10_test

AM_CPPFLAGS += $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS) $(SUPERLU_CPPFLAGS)
AM_LDFLAGS  += $(MPI_LDFLAGS) $(BOOST_LDFLAGS) $(SUPERLU_LDFLAGS)
LDADD = $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)   \
        $(BOOST_SYSTEM_LIB)                             \
        ../../../lib/libdunecornerpoint.la              \
        $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SUPERLU_LIBS)

# Note: Need to fill _SOURCES because default is "<base>.c"
mimetic_ipeval_test_SOURCES = mimetic_ipeval_test.cpp

mimetic_solver_test_SOURCES = mimetic_solver_test.cpp

mimetic_aniso_solver_test_SOURCES = mimetic_aniso_solver_test.cpp

mimetic_periodic_test_SOURCES = mimetic_periodic_test.cpp

parsolver_test_SOURCES = parsolver_test.cpp

istl_test_SOURCES = istl_test.cpp

spe10_test_SOURCES = spe10_test.cpp

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
