# $Date$
# $Revision$

# The dune-cornerpoint library
lib_LTLIBRARIES = libdunecornerpoint.la

libdunecornerpoint_la_SOURCES =

# Dummy C++ source to cause C++ linking.
nodist_EXTRA_libdunecornerpoint_la_SOURCES = dummy.cpp

# Collect constituent parts into cpg support library.
libdunecornerpoint_la_LIBADD = \
   ../dune/grid/cpgrid/libcpgrid.la \
   ../dune/grid/preprocess/libcpgpreprocess_noinst.la \
   ../dune/grid/common/libgrid.la \
   ../dune/common/param/libparam.la \
   ../dune/common/libcommon.la

include $(top_srcdir)/am/global-rules
