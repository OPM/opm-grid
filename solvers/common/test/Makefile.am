check_PROGRAMS = boundaryconditions_test nonuniformtablelinear_test
noinst_PROGRAMS = gie_test simulator_test periodic_test aniso_simulator_test

# Common flags
AM_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) \
	$(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)   \
	$(BOOST_SYSTEM_LIB)                             \
	../../../grid/cpgrid/libcpgrid.la               \
	../../../grid/preprocess/libpreprocess.la       \
	../../../grid/common/libtemporaryhack.la        \
	../../../common/param/libparam.la

# Additional quirks for individual programs
boundaryconditions_test_SOURCES = boundaryconditions_test.cpp

nonuniformtablelinear_test_SOURCES = nonuniformtablelinear_test.cpp

periodic_test_SOURCES = periodic_test.cpp

gie_test_SOURCES                = gie_test.cpp

simulator_test_SOURCES          = simulator_test.cpp
simulator_test_LDADD            = $(LDADD) -llapack -lblas

aniso_simulator_test_SOURCES          = aniso_simulator_test.cpp
aniso_simulator_test_LDADD            = $(LDADD) -llapack -lblas

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
