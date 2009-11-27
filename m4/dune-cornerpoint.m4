# $Date$
# $Revision$

# Additional checks needed to build the module
AC_DEFUN([DUNE_CORNERPOINT_CHECKS],
[
        # LibXML2 support.
        PKG_CHECK_MODULES([libxml2], [libxml-2.0 >= 2.0])

        # BLAS and LAPACK support.
        #
        # NOTE: ACX_LAPACK internally AC_REQUIRE's ACX_BLAS which,
        # subsequently, AC_REQUIRE's AC_F77_LIBRARY_LDFLAGS which sets
        # the $FLIBS macro.
        AC_REQUIRE([AC_F77_WRAPPERS])
        AC_REQUIRE([ACX_LAPACK])

        # Boost support.
        AX_BOOST_BASE([1.37])
        AX_BOOST_DATE_TIME
        AX_BOOST_FILESYSTEM
        AX_BOOST_SYSTEM
        AX_BOOST_UNIT_TEST_FRAMEWORK

        # Additional summary entries.
        DUNE_ADD_SUMMARY_ENTRY([BLAS], [$acx_blas_ok])
        DUNE_ADD_SUMMARY_ENTRY([LAPACK], [$acx_lapack_ok])
])

# Additional checks needed to find the module
AC_DEFUN([DUNE_CORNERPOINT_CHECK_MODULE],
[
        DUNE_CHECK_MODULES([dune-cornerpoint],
                           [grid/CpGrid.hpp])
])
