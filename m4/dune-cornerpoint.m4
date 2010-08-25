# $Date$
# $Revision$

# Additional checks needed to build the module
AC_DEFUN([DUNE_CORNERPOINT_CHECKS],
[
        # LibXML2 support.
        AM_PATH_XML2([2.0.0])

        # Boost support.
        AX_BOOST_BASE([1.37])
        AX_BOOST_DATE_TIME
        AX_BOOST_FILESYSTEM
        AX_BOOST_SYSTEM
        AX_BOOST_UNIT_TEST_FRAMEWORK

        # Add Boost support to module dependencies
        DUNE_ADD_MODULE_DEPS([DUNE_CORNERPOINT],dnl
                             [DUNE_CORNERPOINT],dnl
          [$BOOST_CPPFLAGS],dnl
          [$BOOST_LDFLAGS],dnl
          [[$BOOST_DATE_TIME_LIB]dnl
           [$BOOST_FILESYSTEM_LIB]dnl
           [$BOOST_SYSTEM_LIB]])dnl
])

# Additional checks needed to find the module
AC_DEFUN([DUNE_CORNERPOINT_CHECK_MODULE],
[
        DUNE_CHECK_MODULES([dune-cornerpoint],
                           [common/StopWatch.hpp],
                           [Dune::time::StopWatch s;])
])
