// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_GRIDVIEW_HH
#define DUNE_POLYHEDRALGRID_GRIDVIEW_HH

//- dune-common includes
#include <dune/common/typetraits.hh>

//- dune-grid includes
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridview.hh>

//- dune-metagrid includes
#include <dune/grid/polyhedralgrid/indexset.hh>
#include <dune/grid/polyhedralgrid/intersection.hh>
#include <dune/grid/polyhedralgrid/intersectioniterator.hh>
#include <dune/grid/polyhedralgrid/iterator.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< int dim, int dimworld, PartitionIteratorType defaultpitype >
  class PolyhedralGridView;

  template< int dim, int dimworld, PartitionIteratorType ptype >
  struct PolyhedralGridViewTraits;


  // PolyhedralGridView
  // ----------

  template< int dim, int dimworld, PartitionIteratorType defaultpitype >
  class PolyhedralGridView
  {
    typedef PolyhedralGridView< dim, dimworld, defaultpitype > This;

  public:
    typedef PolyhedralGridViewTraits< dim, dimworld, defaultpitype > Traits;

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
    static const PartitionIteratorType pitype = Traits :: defaultpitype;

    PolyhedralGridView ( const Grid &grid, int level )
    : grid_( &grid )
    {}

    const Grid &grid () const
    {
      assert( grid_ );
      return *grid_;
    }

    const IndexSet &indexSet () const
    {
        return grid().indexSet();
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
      return Impl( grid().extraData() );
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
      return Impl( grid().extraData() );
    }

    IntersectionIterator ibegin ( const typename Codim< 0 >::Entity &entity ) const
    {
      typedef typename Traits::IntersectionIteratorImpl IntersectionIteratorImpl;
      return IntersectionIteratorImpl( grid().extraData(), entity, true);
    }

    IntersectionIterator iend ( const typename Codim< 0 >::Entity &entity ) const
    {
      typedef typename Traits::IntersectionIteratorImpl IntersectionIteratorImpl;
      return IntersectionIteratorImpl( grid().extraData(), entity, false);
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
#warning TODO
#if 0
      typedef CommDataHandleIF< DataHandle, Data > DataHandleIF;
      typedef PolyhedralGridDataHandle< DataHandleIF, Grid > WrappedDataHandle;

      WrappedDataHandle wrappedDataHandle( grid().extraData(), dataHandle );
      hostGridView().communicate( wrappedDataHandle, interface, direction );
#endif
    }

  protected:
    const Grid *grid_;
  };

  // PolyhedralGridViewTraits
  // ----------------

  template< int dim, int dimworld, PartitionIteratorType ptype >
  struct PolyhedralGridViewTraits
  {
    typedef PolyhedralGrid< dim, dimworld > Grid;
    static const PartitionIteratorType pitype = ptype;
    //friend class PolyhedralGridView< pitype >;

    typedef PolyhedralGridView< dim, dimworld, pitype > GridViewImp;
    typedef PolyhedralGridIndexSet< Grid::dimension, Grid::dimensionworld > IndexSet;

    typedef PolyhedralGridIntersection< const Grid > IntersectionImpl;
    typedef Dune::Intersection< const Grid, IntersectionImpl > Intersection;

    typedef PolyhedralGridIntersectionIterator< const Grid > IntersectionIteratorImpl;
    typedef Dune::IntersectionIterator< const Grid, IntersectionIteratorImpl, IntersectionImpl > IntersectionIterator;

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
        typedef PolyhedralGridIterator< const Grid, codim > IteratorImpl;
        typedef Dune::EntityIterator< codim, const Grid, IteratorImpl > Iterator;
      };

      typedef typename Partition< pitype >::Iterator Iterator;
    };

    static const bool conforming = false;
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_GRIDVIEW_HH
