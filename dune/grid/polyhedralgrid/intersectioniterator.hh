#ifndef DUNE_POLYHEDRALGRID_INTERSECTIONITERATOR_HH
#define DUNE_POLYHEDRALGRID_INTERSECTIONITERATOR_HH

#include <dune/grid/polyhedralgrid/entitypointer.hh>
#include <dune/grid/polyhedralgrid/intersection.hh>

namespace Dune
{

  // PolyhedralGridIntersectionIterator
  // --------------------------

  template< class Grid >
  class PolyhedralGridIntersectionIterator
  {
  protected:
    typedef PolyhedralGridIntersectionIterator< Grid > This;

    typedef typename remove_const< Grid >::type::Traits Traits;
    static const bool isLeafIntersection = true;

  public:
    typedef typename conditional< isLeafIntersection,
                                  typename Traits :: LeafIntersection,
                                  typename Traits :: LevelIntersection > :: type  Intersection ;
    typedef typename Intersection :: Implementation IntersectionImpl ;

    typedef typename Traits :: ExtraDataType ExtraData;

    typedef typename Grid::template Codim< 0 >::EntityPointer EntityPointer;

    PolyhedralGridIntersectionIterator ( ExtraData data )
    : intersection_( IntersectionImpl( data ) )
    {}

    PolyhedralGridIntersectionIterator ( const This &other ) = default;
    const This &operator= ( const This &other ) = default;

    bool equals ( const This &other ) const
    {
        return seed_ == other_.seed_ && intersection_.faceIdx_ == other.intersection_.faceIdx_;
    }

    void increment ()
    {
      ++intersection_.faceIdx_;
    }

    const Intersection &dereference () const
    {
      return intersection_;
    }

    ExtraData data() const { return intersectionImpl().data(); }

  protected:
    IntersectionImpl &intersectionImpl () const
    {
      return Grid::getRealImplementation( intersection_ );
    }

    mutable Intersection intersection_;
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_INTERSECTIONITERATOR_HH
