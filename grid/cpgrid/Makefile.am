# $Date$
# $Revision$

SUBDIRS = . test

cpgriddir = $(includedir)/dune/grid/cpgrid

cpgrid_HEADERS = DefaultGeometryPolicy.hpp Entity.hpp EntityRep.hpp \
                 Geometry.hpp Intersection.hpp Iterators.hpp \
                 Indexsets.hpp OrientedEntityTable.hpp

noinst_LTLIBRARIES = libcpgrid.la

libcpgrid_la_SOURCES = CpGrid.cpp readSintefLegacyFormat.cpp \
                       readEclipseFormat.cpp

libcpgrid_la_CXXFLAGS = $(MPI_CPPFLAGS) $(BOOST_CPPFLAGS)
libcpgrid_la_LDFLAGS = ../../common/libcommon.la

include $(top_srcdir)/am/global-rules
