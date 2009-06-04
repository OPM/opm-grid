# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

SUBDIRS = . test

paramdir = $(includedir)/dune/common/param

lib_LTLIBRARIES = libparam.la

libparam_la_SOURCES = Parameter.hpp ParameterGroup.hpp ParameterGroup_impl.hpp \
                      ParameterMapItem.hpp ParameterRequirement.hpp \
                      ParameterStrings.hpp ParameterTools.hpp ParameterXML.hpp \
                      Parameter.cpp ParameterGroup.cpp ParameterTools.cpp ParameterXML.cpp
libparam_la_CPPFLAGS = $(libxml2_CFLAGS)
libparam_la_LIBADD   = $(libxml2_LIBS)

include $(top_srcdir)/am/global-rules

