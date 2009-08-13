# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

SUBDIRS = . param

common_HEADERS = Average.hpp SparseTable.hpp StopWatch.hpp \
                 ErrorMacros.hpp SparseVector.hpp Units.hpp

commondir = $(includedir)/dune/common


lib_LTLIBRARIES = libcommon_hack.la

libcommon_hack_la_SOURCES = StopWatch.cpp
libcommon_hack_la_CPPFLAGS = $(DUNE_CPPFLAGS) $(BOOST_CPPFLAGS)
libcommon_hack_la_LIBADD   = $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) $(BOOST_DATE_TIME_LIB)



include $(top_srcdir)/am/global-rules

