// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_IDSET_HH
#define DUNE_POLYHEDRALGRID_IDSET_HH

#include <dune/grid/common/indexidset.hh>

namespace Dune
{
  template <int dim, int dimworld, typename coord_t> class PolyhedralGrid;

  // PolyhedralGridIdSet
  // -------------------

  template< int dim, int dimworld, typename coord_t >
  class PolyhedralGridIdSet
      : public IdSet< PolyhedralGrid< dim, dimworld, coord_t >, PolyhedralGridIdSet< dim, dimworld, coord_t >, std::size_t /*IdType=size_t*/ >
  {
  public:
    typedef PolyhedralGrid<  dim, dimworld, coord_t > Grid;
    typedef typename std::remove_const< Grid >::type::Traits Traits;
    typedef std::size_t IdType;

    typedef PolyhedralGridIdSet< dim, dimworld, coord_t > This;
    typedef IdSet< Grid, This, IdType > Base;

    explicit PolyhedralGridIdSet (const Grid& grid)
        : grid_( grid )
    {
      update();
    }

    //! Recompute the cached data, to be called when the grid view changes.
    void update()
    {
      globalCellPtr_ = grid_.globalCellPtr();
      globalFacePtr_ = grid_.globalFacePtr();
      globalNodePtr_ = grid_.globalNodePtr();
      // The ids of cells, faces and nodes are in separate ranges, which are the same
      // in all views of the grid: [0, cells), [cells, cells+faces), then the nodes.
      const auto ranges = grid_.idRanges();
      codimOffset_[ 0 ] = 0;
      for( int i=1; i<dim; ++i )
      {
        codimOffset_[ i ] = ranges[ 0 ];
      }
      codimOffset_[ dim ] = ranges[ 0 ] + ranges[ 1 ];
    }

    //! id meethod for entity and specific codim
    template< int codim >
    IdType id ( const typename Traits::template Codim< codim >::Entity &entity ) const
    {
      const int index = entity.seed().index();
      const int* global = codim == 0 ? globalCellPtr_
                        : ( codim == 1 ? globalFacePtr_
                        : ( codim == dim ? globalNodePtr_ : nullptr ) );
      return codimOffset_[ codim ] + IdType( global ? global[ index ] : index );
    }

    //! id method of all entities
    template< class Entity >
    IdType id ( const Entity &entity ) const
    {
      return id< Entity::codimension >( entity );
    }

    //! id method of all entities
    template< class IntersectionImpl >
    IdType id ( const Dune::Intersection< const Grid, IntersectionImpl >& intersection ) const
    {
      return intersection.impl().id();
    }

    //! subId method for entities
    template< class Entity >
    IdType subId ( const Entity &entity, int i, unsigned int codim ) const
    {
      if( codim == 0 )
        return id( entity );
      else if ( codim == 1 )
        return id( entity.template subEntity< 1 >( i ) );
      else if ( codim == dim )
      {
        return id( entity.template subEntity< dim >( i ) );
      }
      else
      {
        DUNE_THROW(NotImplemented,"codimension not available");
        return IdType( -1 );
      }
    }

  protected:
    const Grid& grid_;
    const int* globalCellPtr_ = nullptr;
    const int* globalFacePtr_ = nullptr;
    const int* globalNodePtr_ = nullptr;
    IdType codimOffset_[ dim+1 ];
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_IDSET_HH
