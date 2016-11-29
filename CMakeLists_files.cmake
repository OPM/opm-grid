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
# ATTIC_FILES           Unmaintained files. This for the projects developers
#                       only. Don't expect these files to build.

# originally generated with the command:
# find dune -name '*.c*' -printf '\t%p\n' | sort
list (APPEND MAIN_SOURCE_FILES
	dune/grid/cpgrid/Intersection.cpp
	dune/grid/cpgrid/CpGridData.cpp
	dune/grid/cpgrid/CpGrid.cpp
	dune/grid/cpgrid/GridHelpers.cpp
	dune/grid/cpgrid/PartitionTypeIndicator.cpp
	dune/grid/cpgrid/processEclipseFormat.cpp
	dune/grid/cpgrid/readSintefLegacyFormat.cpp
	dune/grid/cpgrid/writeSintefLegacyFormat.cpp
	dune/grid/common/GeometryHelpers.cpp
	dune/grid/common/GridPartitioning.cpp
	dune/grid/common/WellConnections.cpp
	dune/grid/common/ZoltanGraphFunctions.cpp
	dune/grid/common/ZoltanPartition.cpp
        opm/core/grid/GridHelpers.cpp
        opm/core/grid/GridManager.cpp
        opm/core/grid/GridUtilities.cpp
        opm/core/grid/cart_grid.c
        opm/core/grid/cornerpoint_grid.c
        opm/core/grid/cpgpreprocess/facetopology.c
        opm/core/grid/cpgpreprocess/geometry.c
        opm/core/grid/cpgpreprocess/preprocess.c
        opm/core/grid/cpgpreprocess/uniquepoints.c
        opm/core/grid/grid.c
        opm/core/grid/grid_equal.cpp
        opm/core/utility/StopWatch.cpp
	)

# originally generated with the command:
# find tests/not-unit/ -name \*.cpp -o \*.cc
list (APPEND ATTIC_FILES
	attic/partition_test.cpp
#	attic/dumux_test.cpp
	attic/mapper_test.cpp
	attic/buildcpgrid_test.cpp
#	attic/cpgrid_test.cpp
	attic/check_grid_normals.cpp
  )

# originally generated with the command:
# find tests -name '*.cpp' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_SOURCE_FILES
  tests/cpgrid/distribution_test.cpp
	tests/cpgrid/entityrep_test.cpp
	tests/cpgrid/entity_test.cpp
	tests/cpgrid/facetag_test.cpp
	tests/cpgrid/geometry_test.cpp
	tests/cpgrid/orientedentitytable_test.cpp
	tests/cpgrid/partition_iterator_test.cpp
	tests/cpgrid/zoltan_test.cpp
	tests/grid_test.cc
	tests/p2pcommunicator_test.cc
        tests/test_sparsetable.cpp
	)

# originally generated with the command:
# find tests -name '*.xml' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_DATA_FILES
	)

# originally generated with the command:
# find tutorials examples -name '*.c*' -printf '\t%p\n' | sort
list (APPEND EXAMPLE_SOURCE_FILES
	examples/grdecl2vtu.cpp
	examples/finitevolume/finitevolume.cc
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
	dune/grid/common/p2pcommunicator.hh
	dune/grid/common/p2pcommunicator_impl.hh
	dune/grid/cpgrid/CartesianIndexMapper.hpp
	dune/grid/cpgrid/CpGridData.hpp
	dune/grid/cpgrid/DefaultGeometryPolicy.hpp
	dune/grid/cpgrid/dgfparser.hh
	dune/grid/cpgrid/Entity2IndexDataHandle.hpp
	dune/grid/cpgrid/Entity.hpp
	dune/grid/cpgrid/EntityRep.hpp
	dune/grid/cpgrid/Geometry.hpp
	dune/grid/cpgrid/GlobalIdMapping.hpp
	dune/grid/cpgrid/GridHelpers.hpp
	dune/grid/CpGrid.hpp
	dune/grid/cpgrid/Indexsets.hpp
	dune/grid/cpgrid/Intersection.hpp
	dune/grid/cpgrid/Iterators.hpp
	dune/grid/cpgrid/OrientedEntityTable.hpp
	dune/grid/cpgrid/PartitionIteratorRule.hpp
	dune/grid/cpgrid/PartitionTypeIndicator.hpp
	dune/grid/cpgrid/PersistentContainer.hpp
	dune/grid/common/CartesianIndexMapper.hpp
	dune/grid/common/WellConnections.hpp
	dune/grid/common/ZoltanGraphFunctions.hpp
	dune/grid/common/ZoltanPartition.hpp
	dune/grid/polyhedralgrid/capabilities.hh
	dune/grid/polyhedralgrid/cartesianindexmapper.hh
	dune/grid/polyhedralgrid/declaration.hh
	dune/grid/polyhedralgrid/dgfparser.hh
	dune/grid/polyhedralgrid/entity.hh
	dune/grid/polyhedralgrid/entitypointer.hh
	dune/grid/polyhedralgrid/entityseed.hh
	dune/grid/polyhedralgrid/geometry.hh
	dune/grid/polyhedralgrid/gridhelpers.hh
	dune/grid/polyhedralgrid/grid.hh
	dune/grid/polyhedralgrid/gridview.hh
	dune/grid/polyhedralgrid.hh
	dune/grid/polyhedralgrid/idset.hh
	dune/grid/polyhedralgrid/indexset.hh
	dune/grid/polyhedralgrid/intersection.hh
	dune/grid/polyhedralgrid/intersectioniterator.hh
	dune/grid/polyhedralgrid/iterator.hh
	dune/grid/polyhedralgrid/persistentcontainer.hh
        opm/core/grid.h
        opm/core/grid/CellQuadrature.hpp
        opm/core/grid/ColumnExtract.hpp
        opm/core/grid/FaceQuadrature.hpp
        opm/core/grid/GridHelpers.hpp
        opm/core/grid/GridManager.hpp
        opm/core/grid/GridUtilities.hpp
        opm/core/grid/MinpvProcessor.hpp
        opm/core/grid/PinchProcessor.hpp
        opm/core/grid/cart_grid.h
        opm/core/grid/cornerpoint_grid.h
        opm/core/grid/cpgpreprocess/facetopology.h
        opm/core/grid/cpgpreprocess/geometry.h
        opm/core/grid/cpgpreprocess/preprocess.h
        opm/core/grid/cpgpreprocess/uniquepoints.h
        opm/core/utility/SparseTable.hpp
        opm/core/utility/StopWatch.hpp
	)
