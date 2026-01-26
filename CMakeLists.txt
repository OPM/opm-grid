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

cmake_minimum_required (VERSION 3.23)

# Mandatory call to project
project(opm-grid C CXX)

# project information is in dune.module. Read this file and set variables.
# we cannot generate dune.module since it is read by dunecontrol before
# the build starts, so it makes sense to keep the data there then.

option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
option(REQUIRE_ZOLTAN "Require Zoltan to be found (needed for productive run" ON)

if(SIBLING_SEARCH AND NOT opm-common_DIR)
  # guess the sibling dir
  get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
  get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
  get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
  #Try if <module-name>/<build-dir> is used
  get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
  if(IS_DIRECTORY ${_modules_dir}/opm-common/${_leaf_dir_name})
    set(opm-common_DIR ${_modules_dir}/opm-common/${_leaf_dir_name})
  else()
    string(REPLACE ${PROJECT_NAME} opm-common _opm_common_leaf ${_leaf_dir_name})
    if(NOT _leaf_dir_name STREQUAL _opm_common_leaf
        AND IS_DIRECTORY ${_parent_full_dir}/${_opm_common_leaf})
      # We are using build directories named <prefix><module-name><postfix>
      set(opm-common_DIR ${_parent_full_dir}/${_opm_common_leaf})
    elseif(IS_DIRECTORY ${_parent_full_dir}/opm-common)
      # All modules are in a common build dir
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
OpmSetPolicies()
# not the same location as most of the other projects? this hook overrides
macro (dir_hook)
endmacro (dir_hook)


# list of prerequisites for this particular project; this is in a
# separate file (in cmake/Modules sub-directory) because it is shared
# with the find module
include (${project}-prereqs)

# read the list of components from this file (in the project directory);
# it should set various lists with the names of the files to include
include (CMakeLists_files.cmake)

# Search opm-tests to use some of the data sets for testing
include(Findopm-tests)

macro (config_hook)
	opm_need_version_of ("dune-common")
	opm_need_version_of ("dune-geometry")
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
	if(NOT ZOLTAN_FOUND AND MPI_C_FOUND AND REQUIRE_ZOLTAN)
		message(SEND_ERROR "opm-grid with MPI support requires the package ZOLTAN."
			"Please install it (e.g. from http://www.cs.sandia.gov/zoltan/.)")
	endif(NOT ZOLTAN_FOUND AND MPI_C_FOUND AND REQUIRE_ZOLTAN)
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

# Additional parallel tests
if(MPI_FOUND)
  if(NOT MPIEXEC_EXECUTABLE)
    set(MPIEXEC_EXECUTABLE ${MPIEXEC})
  endif()

  set(tests_4proc
    addLgrsOnDistributedGrid_test
    autoRefine_test
    communicate_distributed_grid_with_lgrs_test
    distribute_level_zero_from_grid_with_lgrs_and_wells_test
    distribute_level_zero_from_grid_with_lgrs_test
    distribution_test
    grid_global_id_set_test
    level_and_grid_cartesianIndexMappers_test
    lgr_cell_id_sync_test
    logicalCartesianSize_and_refinement_test
    mapLevelIndicesToCartesianOutputOrder_test
    test_communication_utils
    id_entity_entityrep_test
  )

  if(Boost_VERSION_STRING VERSION_GREATER 1.53)
    list(APPEND test_4proc lgr_with_inactive_parent_cells_test)
     foreach(procs 3 4)
       add_test(NAME graphofgrid_parallel${procs} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${procs} $<TARGET_FILE:test_graphofgrid_parallel>)
       set_tests_properties(graphofgrid_parallel${procs} PROPERTIES PROCESSORS ${procs})
     endforeach()
  endif()

  foreach(test4 ${tests_4proc})
    add_test(NAME ${test4}_parallel COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:${test4}>)
    set_tests_properties(${test4}_parallel PROPERTIES PROCESSORS 4)
  endforeach()
endif()

if(MPI_FOUND AND HAVE_OPM_TESTS AND HAVE_ECL_INPUT)
  # opm_add_test would need a driver script to run the
  # parallel test (tests/run-parallel-unitTest.sh which
  # is in opm-simulators and needs to be moved) until
  # then we use opm_add_test and an additional call to add_test
  # to run the test in parall. ONLY_COMPILE does not work because
  # somehow -DBOOST_TEST_DYN_LINK is skipped and there is no main.
  opm_add_test(cpgrid_aquifer_test
    SOURCES tests/cpgrid/cpgrid_aquifer_test.cpp
    DEPENDS opmgrid
    #ONLY_COMPILE
    PROCESSORS 4
    CONDITION MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
    DRIVER_ARGS -np 4
    TEST_ARGS -- ${OPM_TESTS_ROOT}/aquifer-num/3D_2AQU_NUM.DATA
    LIBRARIES opmgrid ${Boost_LIBRARIES})
  add_test(NAME cpgrid_aquifer_parallel COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:cpgrid_aquifer_test> -- ${OPM_TESTS_ROOT}/aquifer-num/3D_2AQU_NUM.DATA)
endif()

install(DIRECTORY doc/man1 DESTINATION ${CMAKE_INSTALL_MANDIR}
  FILES_MATCHING PATTERN "*.1")
