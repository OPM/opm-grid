# $Date$
# $Revision$

# The dune-cornerpoint library
lib_LTLIBRARIES = libdunecornerpoint.la

libdunecornerpoint_la_SOURCES =

# Dummy C++ source to cause C++ linking.
nodist_EXTRA_libdunecornerpoint_la_SOURCES = dummy.cpp

# Collect constituent parts into cpg support library.
libdunecornerpoint_la_LIBADD = \
   ../grid/cpgrid/libcpgrid.la \
   ../grid/preprocess/libpreprocess.la \
   ../grid/common/libgrid.la \
   ../common/param/libparam.la \
   ../common/libcommon.la

include $(top_srcdir)/am/global-rules
