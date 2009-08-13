# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

check_PROGRAMS = cpgrid_test mapper_test

noinst_PROGRAMS = grid_read_test finitevolume_test old_new_cpgrid_test #dumux_test

LDADD =  $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) \
         $(DUNE_LIBS)  $(BOOST_SYSTEM_LIB) \
         ../cpgrid/libcpgrid.la \
	 ../preprocess/libpreprocess.la \
         ../../common/param/libparam.la \
         ../common/libtemporaryhack.la \
         ../../common/libcommon_hack.la

#CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)

cpgrid_test_SOURCES = cpgrid_test.cpp
cpgrid_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)

mapper_test_SOURCES = mapper_test.cpp
mapper_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)

finitevolume_test_SOURCES = finitevolume.cc evolve.hh initialize.hh \
                            transportproblem2.hh unitcube.hh vtkout.hh
finitevolume_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)

grid_read_test_SOURCES = grid_read_test.cpp
grid_read_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)

old_new_cpgrid_test_SOURCES = old_new_cpgrid_test.cpp
old_new_cpgrid_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)

#dumux_test_SOURCES = dumux_test.cpp
#dumux_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS) -I/Users/atgeirr/openrs/dune-mux-1.0-rc1
#dumux_test_LDADD = $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) \
#                       $(DUNE_LIBS) ../cpgrid/libcpgrid.la \
#                       ../../common/param/libparam.la $(BOOST_SYSTEM_LIB)

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules

