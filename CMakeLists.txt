# -*- mode: cmake; tab-width: 2; indent-tabs-mode: nil; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 expandtab softtabstop=2 nowrap:

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

cmake_minimum_required (VERSION 3.23)

# Mandatory call to project
project(opm-grid C CXX Fortran)

# project information is in dune.module. Read this file and set variables.
# we cannot generate dune.module since it is read by dunecontrol before
# the build starts, so it makes sense to keep the data there then.

option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
option(REQUIRE_ZOLTAN "Require Zoltan to be found (needed for productive run" ON)
option(USE_OPM_COMMON "Use data structures from opm-common?" ON)

macro(opm-grid_prereqs_hook)
  if(USE_OPM_COMMON)
    target_link_libraries(opmgrid PUBLIC opmcommon)
  endif()

  if(NOT ZOLTAN_FOUND AND MPI_C_FOUND AND REQUIRE_ZOLTAN)
    message(SEND_ERROR
      "opm-grid with MPI support requires the package ZOLTAN. "
      "Please install it (e.g. from http://www.cs.sandia.gov/zoltan/.)"
    )
  endif()

  if(MPI_FOUND)
    target_link_libraries(opmgrid PUBLIC MPI::MPI_C MPI::MPI_CXX)
  endif()

  target_link_libraries(opmgrid PUBLIC dunecommon dunegrid)

  if(TARGET Dune::ISTL)
    target_link_libraries(opmgrid PUBLIC Dune::ISTL)
  endif()

  if(TARGET duneuggrid)
    target_link_libraries(opmgrid PUBLIC duneuggrid)
  endif()

  if(TARGET Zoltan::Zoltan)
    target_link_libraries(opmgrid PUBLIC Zoltan::Zoltan)
  endif()

  if(TARGET METIS::METIS)
    target_link_libraries(opmgrid PUBLIC METIS::METIS)
  endif()
endmacro()

macro(opm-grid_language_hook)
  include(FortranCInterface)
  fortrancinterface_header(${PROJECT_BINARY_DIR}/FCMacros.h MACRO_NAMESPACE FC_)
endmacro()

macro(opm-grid_config_hook)
  opm_need_version_of ("dune-grid")
endmacro()

macro(opm-grid_targets_hook)
  target_include_directories(opmgrid PRIVATE ${PROJECT_BINARY_DIR})
endmacro()

macro(opm-grid_tests_hook)
  # Additional parallel tests
  if(MPI_FOUND)
    set(tests_4proc
      addLgrsOnDistributedGrid_test
      autoRefine_test
      distribution_test
      level_and_grid_cartesianIndexMappers_test
      logicalCartesianSize_and_refinement_test
      test_communication_utils
      id_entity_entityrep_test
    )
    if(USE_OPM_COMMON)
      list(APPEND tests_4proc
        communicate_distributed_grid_with_lgrs_test
        distribute_level_zero_from_grid_with_lgrs_test
        distribute_level_zero_from_grid_with_lgrs_and_wells_test
        grid_global_id_set_test
        lgr_cell_id_sync_test
        mapLevelIndicesToCartesianOutputOrder_test
      )
    endif()

    list(APPEND test_4proc lgr_with_inactive_parent_cells_test)
    foreach(procs 3 4)
      add_test(
        NAME
          graphofgrid_parallel${procs}
        COMMAND
          ${MPIEXEC_EXECUTABLE}
          ${MPIEXEC_NUMPROC_FLAG}
          ${procs}
          $<TARGET_FILE:test_graphofgrid_parallel>
      )
      set_tests_properties(graphofgrid_parallel${procs} PROPERTIES PROCESSORS ${procs})
    endforeach()

    foreach(test4 ${tests_4proc})
      add_test(
        NAME
          ${test4}_parallel
        COMMAND
          ${MPIEXEC_EXECUTABLE}
          ${MPIEXEC_NUMPROC_FLAG}
          4
          $<TARGET_FILE:${test4}>
      )
      set_tests_properties(${test4}_parallel PROPERTIES PROCESSORS 4)
    endforeach()

    if(HAVE_OPM_TESTS AND USE_OPM_COMMON)
      opm_add_test(cpgrid_aquifer_test
        ONLY_COMPILE
        SOURCES
          tests/cpgrid/cpgrid_aquifer_test.cpp
        LIBRARIES
          opmcommon
          opmgrid
          Boost::unit_test_framework
      )
      add_test(
        NAME
          cpgrid_aquifer_parallel
        COMMAND
          ${MPIEXEC_EXECUTABLE}
          ${MPIEXEC_NUMPROC_FLAG}
          4
          $<TARGET_FILE:cpgrid_aquifer_test>
          -- ${OPM_TESTS_ROOT}/aquifer-num/3D_2AQU_NUM.DATA
      )
    endif()
  endif()
endmacro()

macro(opm-grid_install_hook)
  install(
    DIRECTORY
      doc/man1
    DESTINATION
      ${CMAKE_INSTALL_MANDIR}
    FILES_MATCHING PATTERN
      "*.1"
  )
endmacro()

# NOTE: Inline sibling search: duplicated from opm-common/cmake/Modules/OpmSiblingSearch.cmake
#   because we need to locate opm-common before find_package(opm-common) makes that module available.
# NOTE: If you modify this logic, you must also modify the logic in opm-simulators/CMakeLists.txt and in
#   opm-common/cmake/Modules/OpmSiblingSearch.cmake.
# NOTE: See description of logic in opm-common/cmake/Modules/OpmSiblingSearch.cmake for more details.
if(SIBLING_SEARCH AND NOT opm-common_DIR)
  # Strategy 1: Compute relative path from source to build dir, apply to sibling.
  # Handles both flat (repo/build) and nested (repo/builds/debug) layouts.
  get_filename_component(_source_parent "${PROJECT_SOURCE_DIR}" DIRECTORY)
  file(RELATIVE_PATH _build_rel "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
  set(_sibling_build "${_source_parent}/opm-common/${_build_rel}")
  if(NOT _build_rel MATCHES "^\\.\\." AND IS_DIRECTORY "${_sibling_build}")
    set(opm-common_DIR "${_sibling_build}")
  else()
    # Fallback strategies for out-of-source or non-standard layouts
    get_filename_component(_leaf_dir_name "${PROJECT_BINARY_DIR}" NAME)
    get_filename_component(_parent_full_dir "${PROJECT_BINARY_DIR}" DIRECTORY)
    string(REPLACE "${PROJECT_NAME}" opm-common _opm_common_leaf "${_leaf_dir_name}")
    if(NOT _leaf_dir_name STREQUAL _opm_common_leaf
        AND IS_DIRECTORY "${_parent_full_dir}/${_opm_common_leaf}")
      set(opm-common_DIR "${_parent_full_dir}/${_opm_common_leaf}")
    elseif(IS_DIRECTORY "${_parent_full_dir}/opm-common" AND
        EXISTS "${_parent_full_dir}/opm-common/CMakeCache.txt")
      set(opm-common_DIR "${_parent_full_dir}/opm-common")
    endif()
  endif()
endif()
if(opm-common_DIR AND NOT IS_DIRECTORY ${opm-common_DIR})
  message(WARNING "Value ${opm-common_DIR} passed to variable"
    " opm-common_DIR is not a directory")
endif()

find_package(opm-common REQUIRED)

# project information is in dune.module. Read this file and set variables.
# we cannot generate dune.module since it is read by dunecontrol before
# the build starts, so it makes sense to keep the data there then.
include (OpmInit)

# Search opm-tests to use some of the data sets for testing
include(Findopm-tests)

# all setup common to the OPM library modules is done here
include (OpmLibMain)

