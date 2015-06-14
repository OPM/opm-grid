// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_INTERSECTIONITERATOR_HH
#define DUNE_POLYHEDRALGRID_INTERSECTIONITERATOR_HH

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

    typedef typename Grid::Traits Traits;
    typedef typename Traits::template Codim<0>::Entity Element;
    static const bool isLeafIntersection = true;

  public:
    typedef typename conditional< isLeafIntersection,
                                  typename Traits :: LeafIntersection,
                                  typename Traits :: LevelIntersection > :: type  Intersection ;
    typedef typename Intersection :: Implementation IntersectionImpl ;

    typedef typename Traits :: ExtraData ExtraData;

    typedef typename Grid::template Codim< 0 >::EntityPointer EntityPointer;

    PolyhedralGridIntersectionIterator ( ExtraData data, const Element& elem, bool isBegin )
      : intersection_( IntersectionImpl( data, elem.seed(), isBegin?0:-1 ) )
    {}

    PolyhedralGridIntersectionIterator ( const This &other ) = default;
    This &operator= ( const This &other ) = default;

    bool equals ( const This &other ) const
    {
        return intersection_ == other.intersection_;
    }

    void increment ()
    {
      ++(data()->getRealImplementation(intersection_).intersectionIdx_);
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
