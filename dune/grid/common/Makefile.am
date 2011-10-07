# $Date$
# $Revision$

SUBDIRS = . test

commondir = $(includedir)/dune/grid/common
common_HEADERS = GeometryHelpers.hpp GridPartitioning.hpp Volumes.hpp

noinst_LTLIBRARIES = libgrid.la

libgrid_la_SOURCES = GeometryHelpers.cpp GridPartitioning.cpp

libgrid_la_CPPFLAGS = $(DUNE_CPPFLAGS) $(BOOST_CPPFLAGS)
libgrid_la_LDFLAGS  = $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) $(BOOST_DATE_TIME_LIB)


include $(top_srcdir)/am/global-rules
