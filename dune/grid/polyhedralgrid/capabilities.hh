// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_CAPABILITIES_HH
#define DUNE_POLYHEDRALGRID_CAPABILITIES_HH

//- dune-grid includes
#include <dune/grid/common/capabilities.hh>

//- dune-metagrid includes
#include <dune/grid/polyhedralgrid/declaration.hh>

namespace Dune
{

  namespace Capabilities
  {

    // Capabilities from dune-grid
    // ---------------------------

    template< int dim, int dimworld >
    struct hasSingleGeometryType< PolyhedralGrid< dim, dimworld > >
    {
      static const bool v = false;
      static const unsigned int topologyId = ~0u;
    };


    template< int dim, int dimworld >
    struct isCartesian< PolyhedralGrid< dim, dimworld > >
    {
      static const bool v = false;
    };


    template< int dim, int dimworld, int codim >
    struct hasEntity< PolyhedralGrid< dim, dimworld >, codim >
    {
      static const bool v = (codim == 0 || codim == 1 || codim == dim);
    };


#if ! DUNE_VERSION_NEWER(DUNE_GRID, 2, 5)
    template< int dim, int dimworld >
    struct isParallel< PolyhedralGrid< dim, dimworld > >
    {
        static const bool v = false;
    };
#endif


    template< int dim, int dimworld, int codim >
    struct canCommunicate< PolyhedralGrid< dim, dimworld >, codim >
    {
        static const bool v = false;
    };


    template< int dim, int dimworld >
    struct hasBackupRestoreFacilities< PolyhedralGrid< dim, dimworld > >
    {
      static const bool v = false;
    };

    template< int dim, int dimworld >
    struct isLevelwiseConforming< PolyhedralGrid< dim, dimworld > >
    {
      static const bool v = false;
    };

    template< int dim, int dimworld >
    struct isLeafwiseConforming< PolyhedralGrid< dim, dimworld > >
    {
      static const bool v = false;
    };

    template< int dim, int dimworld >
    struct threadSafe< PolyhedralGrid< dim, dimworld > >
    {
      static const bool v = false;
    };

    template< int dim, int dimworld >
    struct viewThreadSafe< PolyhedralGrid< dim, dimworld > >
    {
      static const bool v = false;
    };

  } // namespace Capabilities

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_CAPABILITIES_HH
