# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

check_PROGRAMS = entity_test

entity_test_SOURCES = entity_test.cpp
entity_test_CXXFLAGS = $(MPI_CPPFLAGS)
entity_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS) -lboost_unit_test_framework

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules

