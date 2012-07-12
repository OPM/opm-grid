# $Date$
# $Revision$

# Additional checks needed to build the module
AC_DEFUN([DUNE_CORNERPOINT_CHECKS],
[
dnl     Boost support.

        OPM_BOOST_BASE
        AX_BOOST_SYSTEM
        AX_BOOST_DATE_TIME
        AX_BOOST_FILESYSTEM
        AX_BOOST_UNIT_TEST_FRAMEWORK

        # Add Boost support to module dependencies
        DUNE_ADD_MODULE_DEPS([DUNE_CORNERPOINT],dnl
                             [DUNE_CORNERPOINT],dnl
          [$OPM_BOOST_CPPFLAGS],dnl
          [$OPM_BOOST_LDFLAGS],dnl
          [[$BOOST_DATE_TIME_LIB]dnl
           [$BOOST_FILESYSTEM_LIB]dnl
           [$BOOST_SYSTEM_LIB]])dnl

        DUNE_DEFINE_GRIDTYPE([CPGRID],[(GRIDDIM == 3) && (WORLDDIM == 3)],dnl
                             [Dune::CpGrid], [dune/grid/CpGrid.hpp],dnl
                             [dune/grid/cpgrid/dgfparser.hh])
])



# Additional checks needed to find the module
AC_DEFUN([DUNE_CORNERPOINT_CHECK_MODULE],
[
        DUNE_CHECK_MODULES([dune-cornerpoint],
                           [grid/CpGrid.hpp],
                           [Dune::CpGrid g;])
])
