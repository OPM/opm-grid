# defines that must be present in config.h for our headers
set (opm-grid_CONFIG_VAR
  HAVE_DUNE_ISTL
  HAVE_METIS
  HAVE_MPI
  IS_SCOTCH_METIS_HEADER
  HAVE_ZOLTAN
)

# dependencies
set (opm-grid_DEPS
  # various runtime library enhancements
  "MPI"
  "dune-common REQUIRED"
  "dune-grid REQUIRED"
  "dune-istl"
  "ZOLTAN"
  "PTScotch"
  "Scotch"
  "METIS"
  )

find_package_deps(opm-grid)
