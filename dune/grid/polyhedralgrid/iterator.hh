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

  template< class Grid >
  class PolyhedralGridIterator
  : public PolyhedralGridEntityPointer< Grid >
  {
    typedef PolyhedralGridIterator< Grid > This;
    typedef PolyhedralGridEntityPointer< Grid > Base;

  protected:
    typedef typename Base::ExtraData ExtraData;

  public:
    PolyhedralGridIterator ( ExtraData data, const bool beginIterator )
    : Base( data )
    {
      if( beginIterator )
        entity_ = Entity( data, EntitySeed( 0 ) );
      else
        entity_ = Entity( data );
    }

    /** \brief increment */
    void increment ()
    {
      const int index = entity_.seed().index() + 1 ;
      if( index >= data->size( 0 ) )
        entity_ = Entity( entity_.data() );
      else
        entity_ = Entity( entity_.data(), EntitySeed( index ) );
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_ITERATOR_HH
