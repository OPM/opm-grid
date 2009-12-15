# $Date$
# $Revision$

SUBDIRS = . test

paramdir = $(includedir)/dune/common/param

param_HEADERS = Parameter.hpp ParameterGroup.hpp ParameterGroup_impl.hpp \
                ParameterMapItem.hpp ParameterRequirement.hpp \
                ParameterStrings.hpp ParameterTools.hpp ParameterXML.hpp

noinst_LTLIBRARIES = libparam.la

libparam_la_SOURCES = Parameter.cpp ParameterGroup.cpp ParameterTools.cpp \
                      ParameterXML.cpp

libparam_la_CPPFLAGS = $(libxml2_CFLAGS) $(BOOST_CPPFLAGS)
libparam_la_LDFLAGS  = $(libxml2_LDFLAGS) $(BOOST_LDFLAGS)
libparam_la_LIBADD   = $(libxml2_LIBS) $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)

include $(top_srcdir)/am/global-rules
