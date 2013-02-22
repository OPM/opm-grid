# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# key information about the library
set (project "opm-core")
set (${project}_NAME "${project}")
set (${project}_DESCRIPTION "Open Porous Media Initiative Core Library")
set (${project}_DIR "opm")
set (${project}_VERSION_MAJOR 1)
set (${project}_VERSION_MINOR 0)
set (doxy_dir "Documentation")

# defines that must be present in config.h for our headers
set (${project}_CONFIG_VAR
	HAVE_DYNAMIC_BOOST_TEST
	)

# dependencies
set (${project}_DEPS
	# compile with C99 support if available
	"C99"
	# compile with C++0x/11 support if available
	"CXX11Features"
	# various runtime library enhancements
	"Boost 1.39.0
		COMPONENTS date_time filesystem system unit_test_framework REQUIRED"
	# DUNE dependency
	"dune-common"
	)

# C++ project
cmake_minimum_required (VERSION 2.8)
project (${${project}_NAME})
enable_language (C)
enable_language (CXX)

# additional search modules
set (${project}_MODULE_DIR "${PROJECT_SOURCE_DIR}/cmake/Modules")
list (APPEND CMAKE_MODULE_PATH ${${project}_MODULE_DIR})

# print system information to better pinpoint issues from log alone
include (UseSystemInfo)
system_info ()

# very early try to print repo id (to pinpoint version if something goes wrong)
include (UseVCSInfo)
vcs_info ()

# include special
if (CMAKE_VERSION VERSION_LESS "2.8.7")
	message (STATUS "Enabling backward compatibility modules for CMake ${CMAKE_VERSION}")
	list (APPEND CMAKE_MODULE_PATH "${${project}_MODULE_DIR}/compat-2.8.7")
endif (CMAKE_VERSION VERSION_LESS "2.8.7")

# default settings: build static debug library
include (OpmDefaults)
opm_defaults (${project})
message (STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# use tricks to do faster builds
include (UseFastBuilds)

# precompiled headers
include (UsePrecompHeaders)

# macro to set standard variables (INCLUDE_DIRS, LIBRARIES etc.)
include (OpmFind)
find_and_append_package_list_to (${project} ${${project}_DEPS})

# remove the dependency on the testing framework from the main library;
# it is not possible to query for Boost twice with different components.
list (REMOVE_ITEM ${project}_LIBRARIES ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# don't import more libraries than we need to
include (UseOnlyNeeded)

# put debug information into every executable
include (UseDebugSymbols)

# optimize full if we're not doing a debug build
include (UseOptimization)

# turn on all warnings
include (UseWarnings)

# detect if Boost is in a shared library
include (UseDynamicBoost)

# needed for Debian installation scheme
include (UseMultiArch)

# this module contains code to figure out which files is where
include (OpmFiles)
opm_auto_dirs ()

# put libraries in lib/
opm_out_dirs ()

# identify the compilation units in the library
opm_sources (${project})

# enumerate all testing programs in test/ directory
opm_find_tests ()

# tutorial programs are found in the tutorials/ directory
opm_find_tutorials ()

# example programs are found in the examples/ directory
opm_find_examples ()

# create configuration header which describes available features
# necessary to compile this library. singular version is the names that
# is required by this project alone, plural version transitively
# includes the necessary defines by the dependencies
include (ConfigVars)
list (APPEND ${project}_CONFIG_VARS ${${project}_CONFIG_VAR})
set (CONFIG_H "${PROJECT_BINARY_DIR}/config.h")
configure_vars (
	FILE  CXX  ${CONFIG_H}
	WRITE ${${project}_CONFIG_VARS}
	)

# compile main library; pull in all required includes and libraries
include (OpmCompile)
opm_compile (${project})

# installation target: copy the library together with debug and
# configuration files to system directories
include (OpmInstall)
opm_install (${project})
message (STATUS "This build defaults to installing in ${CMAKE_INSTALL_PREFIX}")

# installation of CMake modules to help user programs locate the library
include (OpmProject)
opm_cmake_config (${project})

# routines to build satellites such as tests, tutorials and samples
include (OpmSatellites)

# tutorial programs are found in the tutorials/ directory
opm_compile_satellites (${project} tutorial "" "")
opm_compile_satellites (${project} examples "" "")

# infrastructure for testing
enable_testing ()
include (CTest)

# make datafiles necessary for tests available in output directory
opm_data (tests datafiles "${tests_DIR}" "*.xml")
opm_compile_satellites (${project} tests "" "${tests_REGEXP}")

# use this target to run all tests
add_custom_target (check
	COMMAND ${CMAKE_CTEST_COMMAND}
	DEPENDS tests
	COMMENT "Checking if library is functional"
	VERBATIM
	)

# generate documentation from source code with Doxygen;
# setup install target for this documentation
include (OpmDoc)
opm_doc (${project} ${doxy_dir})

# provide compatibility with using this build in dunecontrol
include (DuneCompat)
include (LibtoolArchives)
configure_la (${project} ${${project}_TARGET} ${project}_LIBTOOL_ARCHIVE)

### clean in-source builds ###
include (OpmDistClean)
opm_dist_clean (${project})

# smart wrapper that auto-parallelizes builds
file (COPY
	GNUmakefile
	DESTINATION ${PROJECT_BINARY_DIR}
	)
