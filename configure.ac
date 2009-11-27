# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.50])

DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE([foreign])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([grid/CpGrid.hpp])
AM_CONFIG_HEADER([config.h])

# Make dune -> . symlink
DUNE_SYMLINK

# Check all dune-module stuff
DUNE_CHECK_ALL

# Set up flags for automated build system
DUNE_AUTOBUILD_FLAGS

# implicitly set the Dune-flags everywhere
AC_SUBST([AM_CPPFLAGS], '$(DUNE_CPPFLAGS) -I$(top_srcdir)')
AC_SUBST([AM_LDFLAGS], '$(DUNE_LDFLAGS) $(DUNE_LIBS)')

AC_CONFIG_FILES([
  Makefile
  lib/Makefile
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
  solvers/multiscale/Makefile
  solvers/multiscale/test/Makefile
  upscaling/Makefile
  upscaling/test/Makefile
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  dune-cornerpoint.pc
])
AC_OUTPUT

DUNE_SUMMARY_ALL
