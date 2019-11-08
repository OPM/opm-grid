// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_GRIDVIEW_HH
#define DUNE_POLYHEDRALGRID_GRIDVIEW_HH

//- dune-common includes
#include <dune/common/typetraits.hh>

//- dune-grid includes
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridview.hh>

//- polyhedralgrid includes
#include <opm/grid/polyhedralgrid/indexset.hh>
#include <opm/grid/polyhedralgrid/intersection.hh>
#include <opm/grid/polyhedralgrid/intersectioniterator.hh>
#include <opm/grid/polyhedralgrid/iterator.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< int dim, int dimworld, typename coord_t, PartitionIteratorType defaultpitype >
  class PolyhedralGridView;

  template< int dim, int dimworld, typename coord_t, PartitionIteratorType ptype >
  struct PolyhedralGridViewTraits;


  // PolyhedralGridView
  // ------------------

  template< int dim, int dimworld, typename coord_t, PartitionIteratorType defaultpitype >
  class PolyhedralGridView
  {
    typedef PolyhedralGridView< dim, dimworld, coord_t, defaultpitype > This;

  public:
    typedef PolyhedralGridViewTraits< dim, dimworld, coord_t, defaultpitype > Traits;

    typedef typename Traits::Grid Grid;
    typedef typename Traits::IndexSet IndexSet;
    typedef typename Traits::Intersection Intersection;
    typedef typename Traits::IntersectionIterator IntersectionIterator;
    typedef typename Traits::CollectiveCommunication CollectiveCommunication;

    template< int codim >
    struct Codim
    : public Traits::template Codim< codim >
    {};

    static const bool conforming = Traits :: conforming;
    static const PartitionIteratorType pitype = Traits :: pitype;

    PolyhedralGridView ( const Grid &grid, const int level = 0 )
    : grid_( &grid )
    {
    }

    const Grid &grid () const
    {
      assert( grid_ );
      return *grid_;
    }

    const IndexSet &indexSet () const
    {
        return grid().leafIndexSet();
    }

    int size ( int codim ) const
    {
      return grid().size( codim );
    }

    int size ( const GeometryType &type ) const
    {
      return grid().size( type );
    }

    template< int codim >
    typename Codim< codim >::Iterator begin () const
    {
      return begin< codim, defaultpitype >();
    }

    template< int codim, PartitionIteratorType pit >
    typename Codim< codim >::template Partition< pit >::Iterator begin () const
    {
      typedef typename Traits::template Codim< codim >::template Partition< pit >::IteratorImpl Impl;
      return Impl( grid().extraData(), true );
    }

    template< int codim >
    typename Codim< codim >::Iterator end () const
    {
      return end< codim, defaultpitype >();
    }

    template< int codim, PartitionIteratorType pit >
    typename Codim< codim >::template Partition< pit >::Iterator end () const
    {
      typedef typename Traits::template Codim< codim >::template Partition< pit >::IteratorImpl Impl;
      return Impl( grid().extraData(), false );
    }

    IntersectionIterator ibegin ( const typename Codim< 0 >::Entity &entity ) const
    {
      typedef typename Traits::IntersectionIteratorImpl IntersectionIteratorImpl;
      return IntersectionIteratorImpl( grid().extraData(), entity.seed(), true);
    }

    IntersectionIterator iend ( const typename Codim< 0 >::Entity &entity ) const
    {
      typedef typename Traits::IntersectionIteratorImpl IntersectionIteratorImpl;
      return IntersectionIteratorImpl( grid().extraData(), entity.seed(), false);
    }

    const CollectiveCommunication &comm () const
    {
      return grid().comm();
    }

    int overlapSize ( int codim ) const
    {
      return grid().overlapSize( codim );
    }

    int ghostSize ( int codim ) const
    {
      return grid().ghostSize( codim );
    }

    template< class DataHandle, class Data >
    void communicate ( CommDataHandleIF< DataHandle, Data > &dataHandle,
                       InterfaceType interface,
                       CommunicationDirection direction ) const
    {
    }

  protected:
    const Grid *grid_;
  };

  // PolyhedralGridViewTraits
  // ------------------------

  template< int dim, int dimworld, typename coord_t, PartitionIteratorType ptype >
  struct PolyhedralGridViewTraits
  {
    typedef PolyhedralGrid< dim, dimworld, coord_t > Grid;
    static const PartitionIteratorType pitype = ptype;

    typedef PolyhedralGridView< dim, dimworld, coord_t, pitype > GridViewImp;
    typedef PolyhedralGridIndexSet< Grid::dimension, Grid::dimensionworld, coord_t > IndexSet;

    typedef PolyhedralGridIntersection< const Grid > IntersectionImpl;
    typedef PolyhedralGridIntersectionIterator< const Grid > IntersectionIteratorImpl;

#if DUNE_VERSION_NEWER(DUNE_GRID,2,3)
    typedef Dune::Intersection< const Grid, IntersectionImpl > Intersection;
    typedef Dune::IntersectionIterator< const Grid, IntersectionIteratorImpl, IntersectionImpl > IntersectionIterator;
#else
    typedef Dune::Intersection< const Grid, PolyhedralGridIntersection > Intersection;
    typedef Dune::IntersectionIterator< const Grid, PolyhedralGridIntersectionIterator, PolyhedralGridIntersection > IntersectionIterator;
#endif

    typedef typename Grid::CollectiveCommunication CollectiveCommunication;

    template< int codim >
    struct Codim
    {
      typedef typename Grid::Traits::template Codim< codim >::Entity Entity;
      typedef typename Grid::Traits::template Codim< codim >::EntityPointer EntityPointer;

      typedef typename Grid::template Codim< codim >::Geometry Geometry;
      typedef typename Grid::template Codim< codim >::LocalGeometry LocalGeometry;

      template< PartitionIteratorType pit >
      struct Partition
      {
        typedef PolyhedralGridIterator< codim, const Grid, pit > IteratorImpl;
        typedef Dune::EntityIterator< codim, const Grid, IteratorImpl > Iterator;
      };

      typedef typename Partition< pitype >::Iterator Iterator;
    };

    static const bool conforming = false;
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_GRIDVIEW_HH
