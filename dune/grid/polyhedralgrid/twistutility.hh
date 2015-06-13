#ifndef DUNE_POLYHEDRALGRID_TWISTUTILITY_HH
#define DUNE_POLYHEDRALGRID_TWISTUTILITY_HH

//- C++ includes
#include <cassert>

//- dune-metagrid includes
#include <dune/grid/polyhedralgrid/declaration.hh>

#if HAVE_DUNE_FEM
//- dune-fem includes
#include <dune/fem/quadrature/caching/twistutility.hh>

namespace Dune
{

  namespace Fem
  {

    // Specialization for PolyhedralGrid
    // -------------------------

    template< class HostGrid >
    struct TwistUtility< PolyhedralGrid< HostGrid > >
    {
      typedef PolyhedralGrid< HostGrid > GridType;

      typedef typename GridType::Traits::LeafIntersectionIterator  LeafIntersectionIterator;
      typedef typename LeafIntersectionIterator::Intersection  LeafIntersection;
      typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
      typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

      static const int dimension = GridType::dimension;

      typedef TwistUtility< HostGrid > HostTwistUtilityType;

    public:
      template< class Intersection >
      static int twistInSelf ( const GridType &grid, const Intersection &it )
      {
        return HostTwistUtilityType::twistInSelf( grid.hostGrid(),
            GridType::getRealImplementation( it ).hostIntersection() );
      }

      template< class Intersection >
      static int twistInNeighbor ( const GridType &grid, const Intersection &it )
      {
        return HostTwistUtilityType::twistInNeighbor( grid.hostGrid(),
            GridType::getRealImplementation( it ).hostIntersection() );
      }

      template< class Intersection >
      static inline GeometryType
      elementGeometry ( const Intersection &intersection, const bool inside)
      {
        return HostTwistUtilityType::elementGeometry(
            GridType::getRealImplementation( intersection ).hostIntersection(),
            inside );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_DUNE_FEM

#endif // #ifndef DUNE_POLYHEDRALGRID_TWISTUTILITY_HH
