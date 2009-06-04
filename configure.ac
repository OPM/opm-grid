# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE([foreign])

AC_PROG_LIBTOOL
PKG_PROG_PKG_CONFIG([0.21])
PKG_CHECK_MODULES([libxml2], [libxml-2.0 >= 2.0])

AC_CONFIG_SRCDIR([grid/CpGrid.hpp])
AM_CONFIG_HEADER([config.h])

# we need no more than the standard DE-stuff
# this module depends on dune-common dune-grid
# this implies checking for [dune-common], [dune-grid]
DUNE_CHECK_ALL

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

# Boost support
AC_SUBST(BOOST_LIBDIR)
AC_SUBST(BOOST_UNIT_TEST)
AC_SUBST(BOOST_ROOT)

AC_CONFIG_FILES([
  Makefile
  common/Makefile
  common/param/Makefile
  grid/Makefile
  grid/test/Makefile
  grid/common/Makefile
  grid/common/test/Makefile
  grid/cpgrid/Makefile
  grid/cpgrid/test/Makefile
  m4/Makefile
  dune-cornerpoint.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
