# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

SUBDIRS = test

cpgriddir = $(includedir)/dune/grid/cpgrid

cpgrid_HEADERS = Entity.hpp Geometry.hpp GridView.hpp Iterators.hpp Indexsets.hpp DefaultGeometryPolicy.hpp

lib_LTLIBRARIES = libcpgrid.la

libcpgrid_la_SOURCES = CpGrid.cpp
libcpgrid_la_CXXFLAGS = $(MPI_CPPFLAGS)

include $(top_srcdir)/am/global-rules

