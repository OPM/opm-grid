# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# set up project and specify the minimum cmake version
project("opm-grid" C CXX)

cmake_minimum_required(VERSION 3.0)

# Sibling build
option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
if(SIBLING_SEARCH)
  # guess the sibling dir
  get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
  get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
  get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
  get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
  foreach(MOD dune-common dune-geometry dune-grid dune-istl opm-common libecl opm-parser)
    # do not overwrite the location of a module if it was explicitly
    # specified by the user
    if(${MOD}_DIR)
      continue()
    endif()

    # Try various possible locations for the build directory of the dependency
    foreach(BUILD_DIR "${_leaf_dir_name}" "build-cmake" "build" "../${MOD}-build" "../build-${MOD}" "../build" ".")
      if(EXISTS "${_modules_dir}/${MOD}/${BUILD_DIR}/${MOD}-config.cmake")
        set(${MOD}_DIR ${_modules_dir}/${MOD}/${BUILD_DIR})
        break()
      endif()
    endforeach()
  endforeach()
endif()

if(dune-common_DIR AND NOT IS_DIRECTORY ${dune-common_DIR})
  message(WARNING "Value ${dune-common_DIR} passed to variable"
    " dune-common_DIR is not a directory")
endif()
if(opm-common_DIR AND NOT IS_DIRECTORY ${opm-common_DIR})
  message(WARNING "Value ${opm-common_DIR} passed to variable"
    " opm-common_DIR is not a directory")
endif()

# Set CMP0053 (how to handle escape sequences in strings) to the new
# behavior to avoid a pretty annoying cmake warning if a library is
# defined in the toplevel CMakeLists.txt. This should probably be
# considered to be a bug in the dune build system. Note that the old
# behaviour will most likely also work fine, but the result of setting
# this policy to NEW is most likely what's intended.
if (POLICY CMP0053)
  cmake_policy(SET CMP0053 NEW)
endif()

# find the build system (i.e., dune-common) and set cmake's module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
  "${PROJECT_SOURCE_DIR}/cmake/modules")

# include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

# include the OPM cmake macros
include(OpmMacros)

# find the packages needed to compile the module
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

# recursively mark all header files beneath the "opm" and "dune"
# directories for installation.
opm_recusive_export_all_headers("opm" "dune")

# the cmake modules get a special treatment
opm_export_cmake_modules()

# we want all features detected by the build system to be enabled,
# thank you!
dune_enable_all_packages()

# add libgrid.a and specify that it encompasses all source files which
# are located somewhere beneath the "opm" and "dune" directories. Note
# that if the tests and the source files of the library were
# intermingled, this would be not so easy, i.e. it would require a
# separate CMakeLists.txt file in each subdirectory and a call to
# add_subdirectory() for each of them.
opm_recursive_add_library("opmgrid" "opm" "dune")

# add all unit tests
#opm_add_test(grid_test LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_test(test_geom2d LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
#opm_add_test(test_repairzcorn LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_test(p2pcommunicator_test LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_test(test_gridutilities LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_test(test_sparsetable LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_test(test_cartgrid LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_test(test_minpvprocessor LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_test(test_ug LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_test(test_column_extract LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_test(test_quadratures LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_test(distribution_test
  SOURCES "tests/cpgrid/distribution_test.cpp"
  LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_test(geometry_test
  SOURCES "tests/cpgrid/geometry_test.cpp"
  LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_test(entityrep_test
  SOURCES "tests/cpgrid/entityrep_test.cpp"
  LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_test(orientedentitytable_test
  SOURCES "tests/cpgrid/orientedentitytable_test.cpp"
  LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_test(entity_test
  SOURCES "tests/cpgrid/entity_test.cpp"
  LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_test(partition_iterator_test
  SOURCES "tests/cpgrid/partition_iterator_test.cpp"
  LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_test(facetag_test
  SOURCES "tests/cpgrid/facetag_test.cpp"
  LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_test(zoltan_test
  SOURCES "tests/cpgrid/zoltan_test.cpp"
  CONDITION ZOLTAN_FOUND
  LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")

opm_recusive_copy_testdata("tests/*.data" "tests/*.DATA")

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
