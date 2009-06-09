# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

check_PROGRAMS = cpgrid_test

noinst_PROGRAMS = grid_read_test finitevolume_test

cpgrid_test_SOURCES = cpgrid_test.cpp
cpgrid_test_CXXFLAGS = $(MPI_CPPFLAGS)
cpgrid_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS)

finitevolume_test_SOURCES = finitevolume.cc evolve.hh initialize.hh transportproblem2.hh unitcube.hh vtkout.hh
finitevolume_test_CXXFLAGS = $(MPI_CPPFLAGS)
finitevolume_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS) ../cpgrid/libcpgrid.la ../../common/param/libparam.la -lboost_system

grid_read_test_SOURCES = grid_read_test.cpp
grid_read_test_CXXFLAGS = $(MPI_CPPFLAGS)
grid_read_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS) ../cpgrid/libcpgrid.la ../../common/param/libparam.la -lboost_system

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules

