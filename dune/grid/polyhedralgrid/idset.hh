// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_IDSET_HH
#define DUNE_POLYHEDRALGRID_IDSET_HH

#include <dune/common/nullptr.hh>
#include <dune/grid/common/indexidset.hh>

namespace Dune
{

  // PolyhedralGridIdSet
  // -----------

  template< int dim, int dimworld >
  class PolyhedralGridIdSet
      : public IdSet< PolyhedralGrid< dim, dimworld >, PolyhedralGridIdSet< dim, dimworld >, /*IdType=*/int >
  {
  public:
    typedef PolyhedralGrid<  dim, dimworld > Grid;
    typedef typename remove_const< Grid >::type::Traits Traits;
    typedef typename Traits::Index  IdType;

    typedef PolyhedralGridIdSet< dim, dimworld > This;
    typedef IdSet< Grid, PolyhedralGridIdSet< dim, dimworld >, IdType > Base;

    PolyhedralGridIdSet (const Grid& grid)
        : grid_(grid)
    {}

    PolyhedralGridIdSet ( const This &other ) = default;

    const This &operator= ( const This &other ) = default;

    //! id meethod for entity and specific codim
    template< int codim >
    IdType id ( const typename Traits::template Codim< codim >::Entity &entity ) const
    {
        if (codim == 0)
            return grid_.cartesianElementIndexIndex(entity.seed());
        else
            return entity.seed();
    }

    //! id method of all entities
    template< class Entity >
    IdType id ( const Entity &entity ) const
    {
      return id< Entity::codimension >( entity );
    }

    //! subId method for entities
    template< class Entity >
    IdType subId ( const Entity &entity, int i, unsigned int codim ) const
    {
#warning TODO
        return 0;
    }

  protected:
    const Grid& grid_;
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_IDSET_HH
