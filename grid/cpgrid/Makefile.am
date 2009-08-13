# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

SUBDIRS = test

cpgriddir = $(includedir)/dune/grid/cpgrid

cpgrid_HEADERS = Entity.hpp EntityRep.hpp Geometry.hpp \
                 Intersection.hpp Iterators.hpp Indexsets.hpp \
                 DefaultGeometryPolicy.hpp

lib_LTLIBRARIES = libcpgrid.la

libcpgrid_la_SOURCES = CpGrid.cpp readSintefLegacyFormat.cpp \
                       readEclipseFormat.cpp

libcpgrid_la_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
libcpgrid_la_LIBADD = ../../common/libcommon_hack.la

include $(top_srcdir)/am/global-rules

