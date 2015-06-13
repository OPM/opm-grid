#ifndef DUNE_POLYHEDRALGRID_GRIDVIEW_HH
#define DUNE_POLYHEDRALGRID_GRIDVIEW_HH

//- dune-common includes
#include <dune/common/typetraits.hh>

//- dune-grid includes
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridview.hh>

//- dune-metagrid includes
#include <dune/grid/polyhedralgrid/datahandle.hh>
#include <dune/grid/polyhedralgrid/indexset.hh>
#include <dune/grid/polyhedralgrid/intersection.hh>
#include <dune/grid/polyhedralgrid/intersectioniterator.hh>
#include <dune/grid/polyhedralgrid/iterator.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class HostGridView, PartitionIteratorType pitype >
  class PolyhedralGridView;


  // PolyhedralGridView
  // ----------

  template< class GVTraits >
  class PolyhedralGridViewBasic
  {
    typedef PolyhedralGridViewBasic< GVTraits > This;

  public:
    typedef GVTraits Traits;

    typedef typename Traits::HostGridView HostGridView;

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

    PolyhedralGridViewBasic ( const Grid &grid, const HostGridView &hostGridView )
    : grid_( &grid ),
      hostGridView_( hostGridView )
    {}

    const Grid &grid () const
    {
      assert( grid_ );
      return *grid_;
    }

    const IndexSet &indexSet () const
    {
      if( !indexSet_ )
        indexSet_ = IndexSet( hostGridView().indexSet() );
      return indexSet_;
    }

    int size ( int codim ) const
    {
      return hostGridView().size( codim );
    }

    int size ( const GeometryType &type ) const
    {
      return hostGridView().size( type );
    }

    template< int codim >
    typename Codim< codim >::Iterator begin () const
    {
      return begin< codim, pitype >();
    }

    template< int codim, PartitionIteratorType pit >
    typename Codim< codim >::template Partition< pit >::Iterator begin () const
    {
      typedef typename Traits::template Codim< codim >::template Partition< pit >::IteratorImpl Impl;
      return Impl( grid().extraData(), hostGridView().template begin< codim, pit >() );
    }

    template< int codim >
    typename Codim< codim >::Iterator end () const
    {
      return end< codim, pitype >();
    }

    template< int codim, PartitionIteratorType pit >
    typename Codim< codim >::template Partition< pit >::Iterator end () const
    {
      typedef typename Traits::template Codim< codim >::template Partition< pit >::IteratorImpl Impl;
      return Impl( grid().extraData(), hostGridView().template end< codim, pit >() );
    }

    IntersectionIterator ibegin ( const typename Codim< 0 >::Entity &entity ) const
    {
      typedef typename Traits::IntersectionIteratorImpl IntersectionIteratorImpl;
      return IntersectionIteratorImpl( grid().extraData(), hostGridView().ibegin( Grid::getRealImplementation( entity ).hostEntity() ) );
    }

    IntersectionIterator iend ( const typename Codim< 0 >::Entity &entity ) const
    {
      typedef typename Traits::IntersectionIteratorImpl IntersectionIteratorImpl;
      return IntersectionIteratorImpl( grid().extraData(), hostGridView().iend( Grid::getRealImplementation( entity ).hostEntity() ) );
    }

    const CollectiveCommunication &comm () const
    {
      return hostGridView().comm();
    }

    int overlapSize ( int codim ) const
    {
      return hostGridView().overlapSize( codim );
    }

    int ghostSize ( int codim ) const
    {
      return hostGridView().ghostSize( codim );
    }

    template< class DataHandle, class Data >
    void communicate ( CommDataHandleIF< DataHandle, Data > &dataHandle,
                       InterfaceType interface,
                       CommunicationDirection direction ) const
    {
      typedef CommDataHandleIF< DataHandle, Data > DataHandleIF;
      typedef PolyhedralGridDataHandle< DataHandleIF, Grid > WrappedDataHandle;

      WrappedDataHandle wrappedDataHandle( grid().extraData(), dataHandle );
      hostGridView().communicate( wrappedDataHandle, interface, direction );
    }

    const HostGridView &hostGridView () const { return hostGridView_; }

  protected:
    const Grid *grid_;
    HostGridView hostGridView_;
    mutable IndexSet indexSet_;
  };

  // PolyhedralGridViewTraits
  // ----------------

  template< class HGV, PartitionIteratorType ptype >
  struct PolyhedralGridViewTraits
  {
    static const PartitionIteratorType pitype = ptype;

    friend class PolyhedralGridView< HGV, pitype >;

    typedef HGV HostGridView;

    typedef typename HostGridView::Grid HostGrid;

    typedef PolyhedralGridView< HostGridView, pitype > GridViewImp;

    typedef PolyhedralGrid< HostGrid > Grid;

    typedef PolyhedralGridIndexSet< const Grid, typename HostGridView::IndexSet > IndexSet;

    typedef PolyhedralGridIntersection< const Grid, typename HostGridView::Intersection > IntersectionImpl;
    typedef Dune::Intersection< const Grid, IntersectionImpl > Intersection;

    typedef PolyhedralGridIntersectionIterator< const Grid, typename HostGridView::IntersectionIterator > IntersectionIteratorImpl;
    typedef Dune::IntersectionIterator< const Grid, IntersectionIteratorImpl, IntersectionImpl > IntersectionIterator;

    typedef typename HostGridView::CollectiveCommunication CollectiveCommunication;

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
        typedef typename HostGridView::template Codim< codim >::template Partition< pit > HostPartition;

        typedef PolyhedralGridIterator< const Grid, typename HostPartition::Iterator > IteratorImpl;
        typedef Dune::EntityIterator< codim, const Grid, IteratorImpl > Iterator;
      };

      typedef typename Partition< pitype >::Iterator Iterator;
    };

    static const bool conforming = HostGridView::conforming;
  };


  template< class HGV, PartitionIteratorType pitype >
  class PolyhedralGridView : public PolyhedralGridViewBasic< PolyhedralGridViewTraits< HGV, pitype > >
  {
    typedef PolyhedralGridView< HGV, pitype > This;
    typedef PolyhedralGridViewBasic< PolyhedralGridViewTraits< HGV, pitype > > Base;

  public:
    typedef typename Base::HostGridView HostGridView;
    typedef typename Base::Grid Grid;

    PolyhedralGridView ( const Grid &grid, const HostGridView &hostGridView )
    : Base( grid, hostGridView )
    {}
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_GRIDVIEW_HH
