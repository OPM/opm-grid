# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

check_PROGRAMS = sparsetable_test

sparsetable_test_SOURCES = sparsetable_test.cpp
sparsetable_test_CXXFLAGS = $(MPI_CPPFLAGS)
sparsetable_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS) -lboost_unit_test_framework

TEST = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules

