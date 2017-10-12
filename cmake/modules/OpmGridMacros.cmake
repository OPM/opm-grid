# .. cmake_module::
#
# This module's content is executed whenever a Dune module requires or
# suggests opm-grid!
#

find_package(ZOLTAN)
if(ZOLTAN_FOUND)
  dune_register_package_flags(
    LIBRARIES "${ZOLTAN_LIBRARIES}"
    INCLUDE_DIRS "${ZOLTAN_INCLUDE_DIRS}")
endif()

find_package(Boost
  COMPONENTS filesystem system date_time
  REQUIRED)
