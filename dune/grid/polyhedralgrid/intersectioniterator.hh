#ifndef DUNE_POLYHEDRALGRID_INTERSECTIONITERATOR_HH
#define DUNE_POLYHEDRALGRID_INTERSECTIONITERATOR_HH

#include <dune/grid/polyhedralgrid/entitypointer.hh>
#include <dune/grid/polyhedralgrid/intersection.hh>

namespace Dune
{

  // PolyhedralGridIntersectionIterator
  // --------------------------

  template< class Grid, class HostIntersectionIterator >
  class PolyhedralGridIntersectionIterator
  {
  protected:
    typedef PolyhedralGridIntersectionIterator< Grid, HostIntersectionIterator > This;

    typedef typename remove_const< Grid >::type::Traits Traits;

    static const bool isLeafIntersection =
      is_same< HostIntersectionIterator,
               typename Grid::HostGrid::Traits::LeafIntersectionIterator > :: value ;
  public:
    typedef typename conditional< isLeafIntersection,
                                  typename Traits :: LeafIntersection,
                                  typename Traits :: LevelIntersection > :: type  Intersection ;
    typedef typename Intersection :: Implementation IntersectionImpl ;

    typedef typename Traits :: ExtraDataType ExtraData;

    typedef typename Grid::template Codim< 0 >::EntityPointer EntityPointer;

    PolyhedralGridIntersectionIterator ( ExtraData data, const HostIntersectionIterator &hostIterator )
    : intersection_( IntersectionImpl( data ) ),
      hostIterator_( hostIterator )
    {}

    PolyhedralGridIntersectionIterator ( const This &other )
    : intersection_( IntersectionImpl( other.data() ) ),
      hostIterator_( other.hostIterator_ )
    {}

    const This &operator= ( const This &other )
    {
      intersectionImpl() = IntersectionImpl( other.data() );
      hostIterator_ = other.hostIterator_;
      return *this;
    }

    bool equals ( const This &other ) const
    {
      return (hostIterator_ == other.hostIterator_);
    }

    void increment ()
    {
      ++hostIterator_;
      intersectionImpl() = IntersectionImpl( data() );
    }

    const Intersection &dereference () const
    {
      if( !intersectionImpl() )
        intersectionImpl() = IntersectionImpl( data(), *hostIterator_ );
      return intersection_;
    }

    ExtraData data() const { return intersectionImpl().data(); }

  protected:
    IntersectionImpl &intersectionImpl () const
    {
      return Grid::getRealImplementation( intersection_ );
    }

    mutable Intersection intersection_;
    HostIntersectionIterator hostIterator_;
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_INTERSECTIONITERATOR_HH
