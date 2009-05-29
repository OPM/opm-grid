# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

noinst_PROGRAMS = cpgrid_test

cpgrid_test_SOURCES = cpgrid_test.cpp
cpgrid_test_CXXFLAGS = $(MPI_CPPFLAGS)
cpgrid_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS)

include $(top_srcdir)/am/global-rules

