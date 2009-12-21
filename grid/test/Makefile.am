# $Date$
# $Revision$

check_PROGRAMS = cpgrid_test mapper_test partition_test

#noinst_PROGRAMS = finitevolume_test make_vtk_test max_zdist_test grdecl_to_legacy_test
noinst_PROGRAMS = grdecl_to_legacy_test make_vtk_test

AM_CPPFLAGS += $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
AM_LDFLAGS  += $(MPI_LDFLAGS) $(BOOST_LDFLAGS)
# Here we depend on libraries in other parts of
# the tree that should be built before this one,
# but if you do not run make from the top you may
# run into trouble. Recursive make considered harmful...
LDADD =  $(DUNE_LIBS) $(BOOST_SYSTEM_LIB) \
         ../cpgrid/libcpgrid.la ../preprocess/libpreprocess.la \
         ../common/libgrid.la  ../../common/libcommon.la \
         ../../common/param/libparam.la

cpgrid_test_SOURCES = cpgrid_test.cpp

mapper_test_SOURCES = mapper_test.cpp

#max_zdist_test_SOURCES = max_zdist_test.cpp

#finitevolume_test_SOURCES = finitevolume.cc evolve.hh initialize.hh \
#                            transportproblem2.hh unitcube.hh vtkout.hh

make_vtk_test_SOURCES = make_vtk_test.cpp
make_vtk_test_LDADD = $(LDADD) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)


partition_test_SOURCES = partition_test.cpp
partition_test_LDADD = $(LDADD) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

grdecl_to_legacy_test_SOURCES = grdecl_to_legacy_test.cpp
grdecl_to_legacy_test_LDADD = $(LDADD) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
