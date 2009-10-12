# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

check_PROGRAMS = cpgrid_test mapper_test

noinst_PROGRAMS = finitevolume_test make_vtk_test partition_test max_zdist_test

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

max_zdist_test_SOURCES = max_zdist_test.cpp

finitevolume_test_SOURCES = finitevolume.cc evolve.hh initialize.hh \
                            transportproblem2.hh unitcube.hh vtkout.hh
finitevolume_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)

make_vtk_test_SOURCES = make_vtk_test.cpp
make_vtk_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
make_vtk_test_LDADD = $(LDADD) -llapack -lblas


partition_test_SOURCES = partition_test.cpp
partition_test_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
partition_test_LDADD = $(LDADD) -llapack -lblas


TESTS = $(check_PROGRAMS) partition_test

include $(top_srcdir)/am/global-rules
