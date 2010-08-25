# $Date$
# $Revision$

check_PROGRAMS = cpgrid_test mapper_test partition_test

#noinst_PROGRAMS = finitevolume_test make_vtk_test max_zdist_test grdecl_to_legacy_test
noinst_PROGRAMS = grdecl2vtu check_grid_normals

AM_CPPFLAGS += $(DUNEMPICPPFLAGS) $(BOOST_CPPFLAGS)
AM_LDFLAGS  += $(DUNEMPILDFLAGS) $(BOOST_LDFLAGS)
# Here we depend on libraries in other parts of
# the tree that should be built before this one,
# but if you do not run make from the top you may
# run into trouble. Recursive make considered harmful...
LDADD =  $(DUNE_LIBS) $(DUNEMPILIBS) $(BOOST_SYSTEM_LIB) \
         ../cpgrid/libcpgrid.la ../preprocess/libpreprocess.la \
         ../common/libgrid.la  ../../common/libcommon.la \
         ../../common/param/libparam.la

cpgrid_test_SOURCES = cpgrid_test.cpp

mapper_test_SOURCES = mapper_test.cpp

grdecl2vtu_SOURCES = grdecl2vtu.cpp

check_grid_normals_SOURCES = check_grid_normals.cpp

partition_test_SOURCES = partition_test.cpp

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
