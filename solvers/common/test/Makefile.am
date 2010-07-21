# $Date$
# $Revision$

check_PROGRAMS = boundaryconditions_test nonuniformtablelinear_test
noinst_PROGRAMS = gie_test simulator_test periodic_test aniso_simulator_test implicitcap_test aniso_implicitcap_test rockjfunc_test

# Common flags
AM_CPPFLAGS += $(DUNEMPICPPFLAGS) $(BOOST_CPPFLAGS) $(SUPERLU_CPPFLAGS)
AM_LDFLAGS  += $(DUNEMPILDFLAGS) $(BOOST_LDFLAGS) $(SUPERLU_LDFLAGS)
LDADD = $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)   \
	$(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)     \
	$(SUPERLU_LIBS) $(DUNEMPILIBS) \
	../../../lib/libdunecornerpoint.la

# Additional quirks for individual programs
boundaryconditions_test_SOURCES = boundaryconditions_test.cpp

nonuniformtablelinear_test_SOURCES = nonuniformtablelinear_test.cpp

periodic_test_SOURCES = periodic_test.cpp

gie_test_SOURCES = gie_test.cpp

rockjfunc_test_SOURCES = rockjfunc_test.cpp
rockjfunc_test_LDADD   = $(LDADD) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SUPERLU_LIBS)

simulator_test_SOURCES = simulator_test.cpp
simulator_test_LDADD   = $(LDADD) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SUPERLU_LIBS) # -ltbb

aniso_simulator_test_SOURCES = aniso_simulator_test.cpp
aniso_simulator_test_LDADD   = $(LDADD) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SUPERLU_LIBS)

implicitcap_test_SOURCES = implicitcap_test.cpp
implicitcap_test_LDADD   = $(LDADD) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SUPERLU_LIBS)

aniso_implicitcap_test_SOURCES = aniso_implicitcap_test.cpp
aniso_implicitcap_test_LDADD   = $(LDADD) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SUPERLU_LIBS)

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
