// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_INDEXSET_HH
#define DUNE_POLYHEDRALGRID_INDEXSET_HH

#include <vector>

#include <dune/common/nullptr.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/indexidset.hh>

#include <dune/grid/polyhedralgrid/declaration.hh>

namespace Dune
{

  // PolyhedralGridIndexSet
  // --------------

  template< int dim, int dimworld >
  class PolyhedralGridIndexSet
      : public IndexSet< PolyhedralGrid< dim, dimworld >, PolyhedralGridIndexSet< dim, dimworld >, int >
  {
    typedef PolyhedralGrid<dim, dimworld> GridType;

  protected:
    typedef PolyhedralGridIndexSet< dim, dimworld > This;
      typedef IndexSet< GridType, This, int > Base;

    typedef typename remove_const< GridType >::type::Traits Traits;

  public:
    static const int dimension = Traits::dimension;

    typedef typename Base::IndexType IndexType;

    PolyhedralGridIndexSet ( const GridType& grid )
        : grid_(&grid)
    {
        GeometryType t;
        t.makeCube(/*dim=*/3);
        geomTypes_[/*codim=*/0].push_back(t);

        t.makeCube(/*dim=*/0);
        geomTypes_[/*codim=*/3].push_back(t);
    }

    //PolyhedralGridIndexSet( const This &other ) = default;
    //const This &operator= ( const This &other ) = default;

    template< class Entity >
    IndexType index ( const Entity &entity ) const
    {
      return index< Entity::codimension >( entity );
    }

    template< int cd >
    IndexType index ( const typename Traits::template Codim< cd >::Entity &entity ) const
    {
#warning TODO
      return 0;
      //return grid_->getRealImplementation(entity).template index< cd >( );
    }

    template< class Entity >
    IndexType subIndex ( const Entity &entity, int i, unsigned int codim ) const
    {
      return subIndex< Entity::codimension >( entity, i, codim );
    }

    IndexType size ( GeometryType type ) const
    {
      return grid().size( type );
    }

    int size ( int codim ) const
    {
      return grid().size( codim );
    }

    template< class Entity >
    bool contains ( const Entity &entity ) const
    {
        return index(entity) >= 0 && index(entity) < size(Entity::codimension);
    }

    const std::vector< GeometryType > &geomTypes ( int codim ) const
    {
        static std::vector<GeometryType> gt{{GeometryType::cube}};
        return gt;
    }

    const GridType& grid() const { assert( grid_ ); return *grid_; }

  protected:
    const GridType *grid_;
    std::vector<GeometryType> geomTypes_[4];
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_INDEXSET_HH
