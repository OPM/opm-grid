// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_ITERATOR_HH
#define DUNE_POLYHEDRALGRID_ITERATOR_HH

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/entityiterator.hh>

#include <dune/grid/polyhedralgrid/entitypointer.hh>

namespace Dune
{

  // PolyhedralGridIterator
  // --------------

  template< class Grid, int codim >
  class PolyhedralGridIterator
  : public PolyhedralGridEntityPointer< Grid, codim >
  {
    typedef PolyhedralGridIterator< Grid, codim > This;
    typedef PolyhedralGridEntityPointer< Grid, codim > Base;

  protected:
    typedef typename Base::ExtraData ExtraData;
    using Base :: entityImpl;

  public:
    typedef typename Grid::Traits::template Codim<codim>::EntitySeed EntitySeed;
    typedef typename Grid::Traits::template Codim<codim>::EntityImpl EntityImpl;

    PolyhedralGridIterator ( ExtraData data, const bool beginIterator )
    : Base( data )
    {
      if( beginIterator )
        entityImpl() = EntityImpl( data, EntitySeed( 0 ) );
      else
        entityImpl() = EntityImpl( data );
    }

    /** \brief increment */
    void increment ()
    {
      int index = entityImpl().seed().index();
      ++index;

      if( index >= entityImpl().data()->size( 0 ) )
        entityImpl() = EntityImpl( entityImpl().data() );
      else
        entityImpl() = EntityImpl( entityImpl().data(), EntitySeed( index ) );
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_ITERATOR_HH
