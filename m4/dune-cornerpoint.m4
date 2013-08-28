dnl -*- autoconf -*-

dnl locate dune-cornerpoint library itself; this macro is called by every
dnl module that depends on dune-cornerpoint.
AC_DEFUN([DUNE_CORNERPOINT_CHECK_MODULE],
[
 OPM_CHECK_PKG_MODULE([dune-cornerpoint],[1.0],[DUNE module supporting grids in a corner-point format])
])

dnl find all prerequisites of dune-cornerpoint; nothing to do here since
dnl this is done by the CMake module and then stored in the -config file.
AC_DEFUN([DUNE_CORNERPOINT_CHECKS],[])
