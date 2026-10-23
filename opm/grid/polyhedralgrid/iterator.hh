// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_ITERATOR_HH
#define DUNE_POLYHEDRALGRID_ITERATOR_HH

#include <algorithm>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/entityiterator.hh>

#include <opm/grid/polyhedralgrid/entitypointer.hh>

namespace Dune
{

  // PolyhedralGridIterator
  // --------------

  template< int codim, class Grid, PartitionIteratorType pitype >
  class PolyhedralGridIterator
  : public PolyhedralGridEntityPointer< codim, Grid >
  {
    typedef PolyhedralGridIterator< codim, Grid, pitype > This;
    typedef PolyhedralGridEntityPointer< codim, Grid > Base;

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
        setIndex( data, firstValid( data, 0 ) );
    }

    /** \brief increment */
    void increment ()
    {
      const auto data = entityImpl().data();
      setIndex( data, firstValid( data, entityImpl().seed().index() + 1 ) );
    }

  protected:
    // Does an entity of the given partition type belong to the partition of the iterator?
    static constexpr bool contains ( const PartitionType type )
    {
      switch( pitype )
      {
      case Interior_Partition:       return type == InteriorEntity;
      case InteriorBorder_Partition: return type == InteriorEntity || type == BorderEntity;
      case Overlap_Partition:        return type != FrontEntity && type != GhostEntity;
      case OverlapFront_Partition:   return type != GhostEntity;
      case Ghost_Partition:          return type == GhostEntity;
      default:                       return true;
      }
    }

    // The first index not smaller than start of an entity in the partition, or the
    // number of entities if there is none.
    static int firstValid ( ExtraData data, int index )
    {
      const int size = data->size( codim );
      if constexpr( pitype != All_Partition )
      {
        while( index < size && ! contains( data->template entityPartitionType< codim >( index ) ) )
          ++index;
      }
      return std::min( index, size );
    }

    void setIndex ( ExtraData data, const int index )
    {
      if( index >= data->size( codim ) )
        entityImpl() = EntityImpl( data );
      else
        entityImpl() = EntityImpl( data, EntitySeed( index ) );
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_ITERATOR_HH
