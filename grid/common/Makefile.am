# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

SUBDIRS = . test

commondir = $(includedir)/dune/grid/common
common_HEADERS = EclipseGridInspector.hpp EclipseGridParser.hpp GeometryHelpers.hpp Volumes.hpp


lib_LTLIBRARIES = libtemporaryhack.la

libtemporaryhack_la_SOURCES = GeometryHelpers.hpp GeometryHelpers.cpp \
                              EclipseGridInspector.hpp EclipseGridInspector.cpp \
                              EclipseGridParser.hpp EclipseGridParser.cpp \
                              GridPartitioning.hpp GridPartitioning.cpp

libtemporaryhack_la_CPPFLAGS = $(DUNE_CPPFLAGS) $(BOOST_CPPFLAGS)
libtemporaryhack_la_LIBADD   = $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) $(BOOST_DATE_TIME_LIB)



include $(top_srcdir)/am/global-rules

