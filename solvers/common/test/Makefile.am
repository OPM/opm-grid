# $Date$
# $Revision$

check_PROGRAMS = boundaryconditions_test nonuniformtablelinear_test
noinst_PROGRAMS = gie_test simulator_test periodic_test aniso_simulator_test

# Common flags
AM_CPPFLAGS += $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
AM_LDFLAGS  += $(MPI_LDFLAGS) $(BOOST_LDFLAGS)
LDADD = $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)   \
	$(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)     \
	../../../lib/libdunecornerpoint.la

# Additional quirks for individual programs
boundaryconditions_test_SOURCES = boundaryconditions_test.cpp

nonuniformtablelinear_test_SOURCES = nonuniformtablelinear_test.cpp

periodic_test_SOURCES = periodic_test.cpp

gie_test_SOURCES = gie_test.cpp

simulator_test_SOURCES = simulator_test.cpp
simulator_test_LDADD   = $(LDADD) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) # -ltbb
#simulator_test_LDADD   = $(LDADD) -llapack -lblas # -ltbb

aniso_simulator_test_SOURCES = aniso_simulator_test.cpp
aniso_simulator_test_LDADD   = $(LDADD) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) 
#aniso_simulator_test_LDADD   = $(LDADD) -llapack -lblas

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
