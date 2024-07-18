# defines that must be present in config.h for our headers
set (opm-grid_CONFIG_VAR
  DUNE_GRID_VERSION_MAJOR
  DUNE_GRID_VERSION_MINOR
  DUNE_GRID_VERSION_REVISION
  DUNE_GEOMETRY_VERSION_MAJOR
  DUNE_GEOMETRY_VERSION_MINOR
  DUNE_GEOMETRY_VERSION_REVISION
  DUNE_COMMON_VERSION_MAJOR
  DUNE_COMMON_VERSION_MINOR
  DUNE_COMMON_VERSION_REVISION
  HAVE_DUNE_ISTL
  HAVE_METIS
  HAVE_MPI
  HAVE_PTSCOTCH
  IS_SCOTCH_METIS_HEADER
  HAVE_ZOLTAN
  HAVE_OPM_COMMON
  HAVE_ECL_INPUT
  )

# dependencies
set (opm-grid_DEPS
  # compile with C99 support if available
  "C99"
  # various runtime library enhancements
  "MPI"
  "dune-common REQUIRED"
  "dune-grid REQUIRED"
  "dune-istl"
  "opm-common REQUIRED"
  "ZOLTAN"
  "PTScotch"
  "Scotch"
  "METIS"
  )

find_package_deps(opm-grid)
