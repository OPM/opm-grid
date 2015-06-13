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

  template< class Grid, PartitionIteratorType defaultpitype >
  class PolyhedralGridView;

  template< class Grid, PartitionIteratorType ptype >
  struct PolyhedralGridViewTraits;


  // PolyhedralGridView
  // ----------

  template< class Grid, PartitionIteratorType defaultpitype >
  class PolyhedralGridView
  {
    typedef PolyhedralGridView< Grid, defaultpitype > This;

  public:
    typedef PolyhedralGridViewTraits< Grid, defaultpitype > Traits;

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
        return grid_->indexSet();
    }

    int size ( int codim ) const
    {
      return grid_->size( codim );
    }

    int size ( const GeometryType &type ) const
    {
      return grid_->size( type );
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
      return IntersectionIteratorImpl( grid().extraData(), entity, true) );
    }

    IntersectionIterator iend ( const typename Codim< 0 >::Entity &entity ) const
    {
      typedef typename Traits::IntersectionIteratorImpl IntersectionIteratorImpl;
      return IntersectionIteratorImpl( grid().extraData(), entity, false) );
    }

    const CollectiveCommunication &comm () const
    {
      return grid_->.comm();
    }

    int overlapSize ( int codim ) const
    {
      return grid_->.overlapSize( codim );
    }

    int ghostSize ( int codim ) const
    {
      return grid_->.ghostSize( codim );
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

  template< class Grid, PartitionIteratorType ptype >
  struct PolyhedralGridViewTraits
  {
    static const PartitionIteratorType pitype = ptype;
    friend class PolyhedralGridView< pitype >;

    typedef PolyhedralGridView< Grid, pitype > GridViewImp;
    typedef PolyhedralGrid< dim, dimworld > Grid;
    typedef PolyhedralGridIndexSet< Grid::dimension, Grid::dimensionworld > IndexSet;

    typedef PolyhedralGridIntersection< const Grid, Grid::Intersection > IntersectionImpl;
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
        typedef PolyhedralGridIterator< const Grid > IteratorImpl;
        typedef Dune::EntityIterator< codim, const Grid, IteratorImpl > Iterator;
      };

      typedef typename Partition< pitype >::Iterator Iterator;
    };

    static const bool conforming = false;
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_GRIDVIEW_HH
