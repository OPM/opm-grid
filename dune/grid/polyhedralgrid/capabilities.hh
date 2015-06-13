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
      static const bool v = hasSingleGeometryType< HostGrid >::v;
      static const unsigned int topologyId = hasSingleGeometryType< HostGrid >::topologyId;
    };


    template< int dim, int dimworld >
    struct isCartesian< PolyhedralGrid< dim, dimworld > >
    {
      static const bool v = isCartesian< HostGrid >::v;
    };


    template< int dim, int dimworld, int codim >
    struct hasEntity< PolyhedralGrid< dim, dimworld >, codim >
    {
      static const bool v = hasEntity< HostGrid, codim >::v;
    };


    template< int dim, int dimworld >
    struct isParallel< PolyhedralGrid< dim, dimworld > >
    {
      static const bool v = isParallel< HostGrid >::v;
    };


    template< int dim, int dimworld, int codim >
    struct canCommunicate< PolyhedralGrid< dim, dimworld >, codim >
    {
      static const bool v = canCommunicate< HostGrid, codim >::v;
    };


    template< int dim, int dimworld >
    struct hasBackupRestoreFacilities< PolyhedralGrid< dim, dimworld > >
    {
      static const bool v = hasBackupRestoreFacilities< HostGrid >::v;
    };

    template< int dim, int dimworld >
    struct isLevelwiseConforming< PolyhedralGrid< dim, dimworld > >
    {
      static const bool v = isLevelwiseConforming< HostGrid >::v;
    };

    template< int dim, int dimworld >
    struct isLeafwiseConforming< PolyhedralGrid< dim, dimworld > >
    {
      static const bool v = isLeafwiseConforming< HostGrid >::v;
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
