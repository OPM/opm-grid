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

    template< class HostGrid >
    struct hasSingleGeometryType< PolyhedralGrid< HostGrid > >
    {
      static const bool v = hasSingleGeometryType< HostGrid >::v;
      static const unsigned int topologyId = hasSingleGeometryType< HostGrid >::topologyId;
    };


    template< class HostGrid >
    struct isCartesian< PolyhedralGrid< HostGrid > >
    {
      static const bool v = isCartesian< HostGrid >::v;
    };


    template< class HostGrid, int codim >
    struct hasEntity< PolyhedralGrid< HostGrid >, codim >
    {
      static const bool v = hasEntity< HostGrid, codim >::v;
    };


    template< class HostGrid >
    struct isParallel< PolyhedralGrid< HostGrid > >
    {
      static const bool v = isParallel< HostGrid >::v;
    };


    template< class HostGrid, int codim >
    struct canCommunicate< PolyhedralGrid< HostGrid >, codim >
    {
      static const bool v = canCommunicate< HostGrid, codim >::v;
    };


    template< class HostGrid >
    struct hasBackupRestoreFacilities< PolyhedralGrid< HostGrid > >
    {
      static const bool v = hasBackupRestoreFacilities< HostGrid >::v;
    };

    template< class HostGrid >
    struct isLevelwiseConforming< PolyhedralGrid< HostGrid > >
    {
      static const bool v = isLevelwiseConforming< HostGrid >::v;
    };

    template< class HostGrid >
    struct isLeafwiseConforming< PolyhedralGrid< HostGrid > >
    {
      static const bool v = isLeafwiseConforming< HostGrid >::v;
    };

    template< class HostGrid >
    struct threadSafe< PolyhedralGrid< HostGrid > >
    {
      static const bool v = false;
    };

    template< class HostGrid >
    struct viewThreadSafe< PolyhedralGrid< HostGrid > >
    {
      static const bool v = false;
    };

  } // namespace Capabilities

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_CAPABILITIES_HH
