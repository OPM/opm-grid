# These packages are always required
find_package(dune-common REQUIRED)
find_package(dune-grid REQUIRED)

# If the target is created, it means we are used from the config file.
# Use compile definitions to decide which packages are required.
if(TARGET opmgrid)
  get_property(opm-grid_COMPILE_DEFINITIONS TARGET opmgrid PROPERTY INTERFACE_COMPILE_DEFINITIONS)
  get_property(opm-grid_LIBS TARGET opmgrid PROPERTY INTERFACE_LINK_LIBRARIES)
  if(opm-grid_COMPILE_DEFINITIONS MATCHES HAVE_MPI)
    find_package(MPI REQUIRED)
  endif()
  if(opm-grid_COMPILE_DEFINITIONS MATCHES HAVE_DUNE_ISTL)
    find_package(dune-istl REQUIRED)
  endif()
  if(opm-grid_COMPILE_DEFINITIONS MATCHES HAVE_ZOLTAN)
    find_package(ZOLTAN REQUIRED)
  endif()
  if(opm-grid_COMPILE_DEFINITIONS MATCHES HAVE_METIS)
    find_package(METIS REQUIRED)
  endif()
  if(opm-grid_LIBS MATCHES duneuggrid)
    find_package(dune-uggrid REQUIRED)
  endif()
else()
  # Used in the opm-grid build system. These are optional
  # so no required flag here.
  find_package(dune-uggrid)
  find_package(MPI)
  find_package(dune-istl)
  find_package(ZOLTAN)
  find_package(METIS)
endif()
