# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# This file sets up six lists:
#	MAIN_SOURCE_FILES     List of compilation units which will be included in
#	                      the library. If it isn't on this list, it won't be
#	                      part of the library. Please try to keep it sorted to
#	                      maintain sanity.
#
#	TEST_SOURCE_FILES     List of programs that will be run as unit tests.
#
#	TEST_DATA_FILES       Files from the source three that should be made
#	                      available in the corresponding location in the build
#	                      tree in order to run tests there.
#
#	EXAMPLE_SOURCE_FILES  Other programs that will be compiled as part of the
#	                      build, but which is not part of the library nor is
#	                      run as tests.
#
#	PUBLIC_HEADER_FILES   List of public header files that should be
#	                      distributed together with the library. The source
#	                      files can of course include other files than these;
#	                      you should only add to this list if the *user* of
#	                      the library needs it.
#
# NOT_UNIT_TEST_SOURCE_FILES List of programms that will be build during 
#                            make non-unit-tests

# originally generated with the command:
# find dune -name '*.c*' -printf '\t%p\n' | sort
list (APPEND MAIN_SOURCE_FILES
	dune/grid/common/GeometryHelpers.cpp
	dune/grid/common/GridPartitioning.cpp
	dune/grid/cpgrid/CpGrid.cpp
	dune/grid/cpgrid/readEclipseFormat.cpp
	dune/grid/cpgrid/readSintefLegacyFormat.cpp
	dune/grid/cpgrid/writeSintefLegacyFormat.cpp
	)

# originally generated with the command:
# find tests/not-unit/ -name \*.cpp -o \*.cc
list (APPEND NOT_UNIT_TEST_SOURCE_FILES
  tests/not-unit/finitevolume.cc
  tests/not-unit/partition_test.cpp
#  tests/not-unit/dumux_test.cpp
  tests/not-unit/mapper_test.cpp
  tests/not-unit/cpgrid/buildcpgrid_test.cpp
#  tests/not-unit/cpgrid_test.cpp
  tests/not-unit/max_zdist_test.cpp
  tests/not-unit/check_grid_normals.cpp
  )

# originally generated with the command:
# find tests -name '*.cpp' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_SOURCE_FILES
	tests/cpgrid/entityrep_test.cpp
	tests/cpgrid/entity_test.cpp
	tests/cpgrid/geometry_test.cpp
	tests/cpgrid/orientedentitytable_test.cpp
	)

# originally generated with the command:
# find tests -name '*.xml' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_DATA_FILES
	)

# originally generated with the command:
# find tutorials examples -name '*.c*' -printf '\t%p\n' | sort
list (APPEND EXAMPLE_SOURCE_FILES
	examples/grdecl2vtu.cpp
	)

# programs listed here will not only be compiled, but also marked for
# installation
list (APPEND PROGRAM_SOURCE_FILES
	examples/grdecl2vtu.cpp
	)

# originally generated with the command:
# find dune -name '*.h*' -a ! -name '*-pch.hpp' -printf '\t%p\n' | sort
list (APPEND PUBLIC_HEADER_FILES
	dune/grid/common/GeometryHelpers.hpp
	dune/grid/common/GridAdapter.hpp
	dune/grid/common/GridPartitioning.hpp
	dune/grid/common/Volumes.hpp
	dune/grid/cpgrid/DefaultGeometryPolicy.hpp
	dune/grid/cpgrid/dgfparser.hh
	dune/grid/cpgrid/Entity.hpp
	dune/grid/cpgrid/EntityRep.hpp
	dune/grid/cpgrid/Geometry.hpp
	dune/grid/CpGrid.hpp
	dune/grid/cpgrid/Indexsets.hpp
	dune/grid/cpgrid/Intersection.hpp
	dune/grid/cpgrid/Iterators.hpp
	dune/grid/cpgrid/OrientedEntityTable.hpp
	dune/grid/cpgrid/PersistentContainer.hpp
	)
