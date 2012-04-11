# $Date$
# $Revision$

# Additional checks needed to build the module
AC_DEFUN([DUNE_CORNERPOINT_CHECKS],
[
dnl        # Boost support.
dnl
dnl     ISTL already configures Boost during checking for Boost.Fusion
dnl     support.  Admittedly, this test uses an unadorned AX_BOOST_BASE
dnl     call so we're only guaranteed to have Boost >= 1.20.0, but in
dnl     practice any recent Linux distro will have Boost >= 1.37.  Bank
dnl     on that...
dnl
dnl     AX_BOOST_BASE([1.37])
dnl
        AX_BOOST_SYSTEM
        AX_BOOST_DATE_TIME
        AX_BOOST_FILESYSTEM
        AX_BOOST_UNIT_TEST_FRAMEWORK

        dnl Check for OPM-Core

        dnl Cater to the possibility that the OPM-Core library is
        dnl installed and functional, but might depend on UMFPACK and
        dnl does not transitively convey that dependency.
        dnl
        AC_SEARCH_LIBS([amd_free],             [amd])
        AC_SEARCH_LIBS([camd_free],            [camd])
        AC_SEARCH_LIBS([colamd_set_defaults],  [colamd])
        AC_SEARCH_LIBS([ccolamd_set_defaults], [ccolamd])
        AC_SEARCH_LIBS([cholmod_l_start],      [cholmod])
        AC_SEARCH_LIBS([umfpack_dl_solve],     [umfpack],dnl
                       ,,[[$BLAS_LIBS] [$LIBS] [$FLIBS]])

        AC_CHECK_HEADERS([opm/core/grid/cornerpoint_grid.h],dnl
                         [opmcore_header=yes], [opmcore_header=no])

        AC_SEARCH_LIBS([create_grid_cornerpoint], [opmcore],dnl
                       [opmcore_lib=yes], [opmcore_lib=no],dnl
                       [[$BOOST_LDFLAGS]dnl
                        [$BOOST_SYSTEM_LIB]dnl
                        [$LAPACK_LIBS] [$BLAS_LIBS] [$LIBS] [$FLIBS]])

        AS_IF([test "$opmcore_header" != "yes" -o "$opmcore_lib" != "yes"],dnl
              [AC_MSG_ERROR([No suitable OPM-Core library found!])],dnl
              [:])

        # Add Boost support to module dependencies
        DUNE_ADD_MODULE_DEPS([DUNE_CORNERPOINT],dnl
                             [DUNE_CORNERPOINT],dnl
          [$BOOST_CPPFLAGS],dnl
          [$BOOST_LDFLAGS],dnl
          [[$BOOST_DATE_TIME_LIB]dnl
           [$BOOST_FILESYSTEM_LIB]dnl
           [$BOOST_SYSTEM_LIB]])dnl

        DUNE_DEFINE_GRIDTYPE([CPGRID],[(GRIDDIM == 3) && (WORLDDIM == 3)],dnl
                             [Dune::CpGrid], [dune/grid/CpGrid.hpp],dnl
                             [dune/grid/cpgrid/dgfparser.hh])

        DUNE_ADD_SUMMARY_ENTRY([opm-core], [$opmcore_lib])
])



# Additional checks needed to find the module
AC_DEFUN([DUNE_CORNERPOINT_CHECK_MODULE],
[
        DUNE_CHECK_MODULES([dune-cornerpoint],
                           [grid/CpGrid.hpp],
                           [Dune::CpGrid g;])
])
