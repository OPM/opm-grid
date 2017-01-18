# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

###########################################################################
#                                                                         #
# Note: The bulk of the build system is located in the cmake/ directory.  #
#       This file only contains the specializations for this particular   #
#       project. Most likely you are interested in editing one of these   #
#       files instead:                                                    #
#                                                                         #
#       dune.module                              Name and version number  #
#       CMakeLists_files.cmake                   Path of source files     #
#       cmake/Modules/${project}-prereqs.cmake   Dependencies             #
#                                                                         #
###########################################################################

cmake_minimum_required (VERSION 2.8)
# project information is in dune.module. Read this file and set variables.
# we cannot generate dune.module since it is read by dunecontrol before
# the build starts, so it makes sense to keep the data there then.

# if true the OPM build system is used, otherwise the DUNE build system
set(USE_OPM_BUILDSYSTEM TRUE)

set( OPM_COMMON_ROOT "" CACHE PATH "Root directory containing OPM related cmake modules")
option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)

if(NOT OPM_COMMON_ROOT)
   find_package(opm-common QUIET)
endif()

if (opm-common_FOUND)
   include(OpmInit)
else()
   unset(opm-common_FOUND)

   if (NOT OPM_COMMON_ROOT AND SIBLING_SEARCH)
      set(OPM_COMMON_ROOT ${PROJECT_SOURCE_DIR}/../opm-common)
   endif()
   if (OPM_COMMON_ROOT)
      list( APPEND CMAKE_MODULE_PATH "${OPM_COMMON_ROOT}/cmake/Modules")
      include (OpmInit OPTIONAL RESULT_VARIABLE OPM_INIT)
      set( OPM_MACROS_ROOT ${OPM_COMMON_ROOT} )
   endif()

   if (NOT OPM_INIT)
      message( "" )
      message( " /---------------------------------------------------------------------------------\\")
      message( " |  Could not locate the opm build macros. The opm build macros                    |")
      message( " |  are in a separate repository - instructions to proceed:                        |")
      message( " |                                                                                 |")
      message( " |    1. Clone the repository: git clone git@github.com:OPM/opm-common.git         |")
      message( " |                                                                                 |")
      message( " |    2. Run cmake in the current project with -DOPM_COMMON_ROOT=<path>/opm-common |")
      message( " |                                                                                 |")
      message( " |    3. If only using DUNE not further action needed!                             |")
      message( " |                                                                                 |")
      message( " \\---------------------------------------------------------------------------------/")
      message( "" )
      set(USE_OPM_BUILDSYSTEM FALSE)
   endif()

endif()

# true if opm-common was found
if(USE_OPM_BUILDSYSTEM)

# not the same location as most of the other projects? this hook overrides
macro (dir_hook)
endmacro (dir_hook)

# project information is in dune.module. Read this file and set variables.
# we cannot generate dune.module since it is read by dunecontrol before
# the build starts, so it makes sense to keep the data there then.
include (OpmInit)

# list of prerequisites for this particular project; this is in a
# separate file (in cmake/Modules sub-directory) because it is shared
# with the find module
include (${project}-prereqs)

# read the list of components from this file (in the project directory);
# it should set various lists with the names of the files to include
include (CMakeLists_files.cmake)

macro (config_hook)
  opm_need_version_of ("dune-common")
  opm_need_version_of ("dune-grid")
  find_file(_HAVE_DUNE_GRID_CHECKS checkpartition.cc HINTS ${dune-grid_INCLUDE_DIR}
    PATH_SUFFIXES dune/grid/test)
  if(_HAVE_DUNE_GRID_CHECKS)
    set(HAVE_DUNE_GRID_CHECKS 1)
  else(_HAVE_DUNE_GRID_CHECK)
    set(HAVE_DUNE_GRID_CHECKS 0)
  endif(_HAVE_DUNE_GRID_CHECKS)
  list (APPEND ${project}_CONFIG_IMPL_VARS
    HAVE_DUNE_GRID_CHECKS
    )
endmacro (config_hook)

macro (prereqs_hook)
endmacro (prereqs_hook)

macro (sources_hook)
endmacro (sources_hook)

macro (fortran_hook)
endmacro (fortran_hook)

macro (files_hook)
endmacro (files_hook)

macro (tests_hook)
endmacro (tests_hook)

# all setup common to the OPM library modules is done here
include (OpmLibMain)

# end if OPM build system part
else()
###########################################################
#
# DUNE build system
# if opm-common has not been found, use dune-build system
#
###########################################################
  project("opm-grid" C CXX)

  # guess dune-common build dir
  if(NOT (dune-common_DIR OR dune-common_ROOT OR
        "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
        string(REPLACE ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
        ${PROJECT_BINARY_DIR})
  endif()

  #find dune-common and set the module path
  find_package(dune-common REQUIRED)

  list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
    "${PROJECT_SOURCE_DIR}/cmake/modules")

  #include the dune macros
  include(DuneMacros)

  # start a dune project with information from dune.module
  dune_project()
  dune_enable_all_packages()

  # get the prerequisite Boost libraries
  #find_package(Boost 1.44.0 COMPONENTS filesystem date_time system unit_test_framework regex ${OPM_PARSER_QUIET})
  set(Boost_USE_STATIC_LIBS OFF)
  set(BOOST_ALL_DYN_LINK    ON)
  find_package(Boost 1.44.0 COMPONENTS filesystem date_time system unit_test_framework regex)

  if( Boost_FOUND )
    message( "Boost ${Boost_LIBRARIES} ${Boost_INCLUDE_DIRS}")
    # register Boost to all pkg flags
    dune_register_package_flags(LIBRARIES "${Boost_LIBRARIES}"
     INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
    add_definitions( "-DHAVE_DYNAMIC_BOOST_TEST=1" )
  endif()

  # include source file lists
  include(CMakeLists_files.cmake)

  # libopmgrid
  dune_add_library(opmgrid ${MAIN_SOURCE_FILES} ADD_LIBS ${DUNE_LIBS} )

  # register opmgrid to all pkg flags
  dune_register_package_flags(LIBRARIES "${OPM_GRID_LIBRARY}"
                              INCLUDE_DIRS "${OPM_GRID_INCLUDEDIRS}")

  # set output dir of binaries to bin
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

  if( DUNE_GRID_VERSION VERSION_LESS "2.5" )
    message( "libraries = ${LIBRARIES} ${ALL_PKG_LIBS}" )
    add_executable( grid_test tests/grid_test.cc )
    target_link_libraries( grid_tes opmgrid )
    add_dune_all_flags( grid_test )
  else()
    foreach( TESTSRCFILE ${TEST_SOURCE_FILES} ${EXAMPLE_SOURCE_FILES})
      set(TESTPROG "${TESTSRCFILE}")
      string(REGEX REPLACE "examples/" "" TESTPROG "${TESTPROG}" )
      string(REGEX REPLACE "finitevolume/" "" TESTPROG "${TESTPROG}" )
      string(REGEX REPLACE "tests/" "" TESTPROG "${TESTPROG}" )
      string(REGEX REPLACE "cpgrid/" "" TESTPROG "${TESTPROG}" )
      string(REGEX REPLACE ".cc" "" TESTPROG "${TESTPROG}" )
      string(REGEX REPLACE ".cpp" "" TESTPROG "${TESTPROG}" )
      # message( "add test ${TESTPROG} with file ${TESTSRCFILE} ")
      dune_add_test( NAME ${TESTPROG} SOURCES ${TESTSRCFILE} LINK_LIBRARIES opmgrid)
    endforeach()
  endif()

  # finalize the dune project, e.g., generate config.h etc.
  finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
endif()
