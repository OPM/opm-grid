# $Date$
# $Revision$

check_PROGRAMS = cpgrid_test mapper_test partition_test

#noinst_PROGRAMS = finitevolume_test make_vtk_test partition_test max_zdist_test

AM_CPPFLAGS += $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
AM_LDFLAGS  += $(MPI_LDFLAGS) $(BOOST_LDFLAGS)
LDADD =  $(DUNE_LIBS) $(BOOST_SYSTEM_LIB) \
         ../../lib/libdunecornerpoint.la

cpgrid_test_SOURCES = cpgrid_test.cpp

mapper_test_SOURCES = mapper_test.cpp

#max_zdist_test_SOURCES = max_zdist_test.cpp

#finitevolume_test_SOURCES = finitevolume.cc evolve.hh initialize.hh \
#                            transportproblem2.hh unitcube.hh vtkout.hh

#make_vtk_test_SOURCES = make_vtk_test.cpp
#make_vtk_test_LDADD = $(LDADD) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
#make_vtk_test_LDADD = $(LDADD) -llapack -lblas


partition_test_SOURCES = partition_test.cpp
partition_test_LDADD = $(LDADD) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
#partition_test_LDADD = $(LDADD) -llapack -lblas


TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
