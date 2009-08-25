# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.50])

DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE([foreign])

#LT_INIT
#LT_LANG([C++])
#LT_OUTPUT

#AC_PROG_SED
PKG_PROG_PKG_CONFIG([0.21])
PKG_CHECK_MODULES([libxml2], [libxml-2.0 >= 2.0])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([grid/CpGrid.hpp])
AM_CONFIG_HEADER([config.h])

# we need no more than the standard DE-stuff
# this module depends on dune-common dune-grid
# this implies checking for [dune-common], [dune-grid]
DUNE_CHECK_ALL

# Boost support.
# NOTE: AX_BOOST_* macros for compiled Boost libraries such as
#       Boost.Filesystem *MUST* be called before LIBS="$DUNE_LIBS"
#       below.  This is due to the macros' explicit use of the $LIBS
#       variable.
AX_BOOST_BASE([[1.37]])
AX_BOOST_DATE_TIME
AX_BOOST_FILESYSTEM
AX_BOOST_SYSTEM
AX_BOOST_UNIT_TEST_FRAMEWORK

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

AC_CONFIG_FILES([
  dune-cornerpoint.pc
  Makefile
  m4/Makefile
  common/Makefile
  common/param/Makefile
  common/test/Makefile
  grid/Makefile
  grid/test/Makefile
  grid/common/Makefile
  grid/common/test/Makefile
  grid/preprocess/Makefile
  grid/cpgrid/Makefile
  grid/cpgrid/test/Makefile
  solvers/Makefile
  solvers/common/Makefile
  solvers/common/test/Makefile
  solvers/euler/Makefile
  solvers/euler/test/Makefile
  solvers/mimetic/Makefile
  solvers/mimetic/test/Makefile
  upscaling/Makefile
  upscaling/test/Makefile
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
