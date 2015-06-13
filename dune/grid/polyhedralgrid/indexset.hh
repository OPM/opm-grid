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

  template< class Grid, class HostIndexSet >
  class PolyhedralGridIndexSet
  : public IndexSet< Grid, PolyhedralGridIndexSet< Grid, HostIndexSet >, typename HostIndexSet::IndexType >
  {
  protected:
    typedef PolyhedralGridIndexSet< Grid, HostIndexSet > This;
    typedef IndexSet< Grid, PolyhedralGridIndexSet< Grid, HostIndexSet >, typename HostIndexSet::IndexType > Base;

    typedef typename remove_const< Grid >::type::Traits Traits;

    typedef typename Traits::HostGrid HostGrid;

  public:
    static const int dimension = Traits::dimension;

    typedef typename Base::IndexType IndexType;

    PolyhedralGridIndexSet ()
    : hostIndexSet_( nullptr )
    {}

    explicit PolyhedralGridIndexSet ( const HostIndexSet &hostIndexSet )
    : hostIndexSet_( &hostIndexSet )
    {}

    PolyhedralGridIndexSet ( const This &other )
    : hostIndexSet_( other.hostIndexSet_ )
    {}

    const This &operator= ( const This &other )
    {
      hostIndexSet_ = other.hostIndexSet_;
      return *this;
    }

    template< class Entity >
    IndexType index ( const Entity &entity ) const
    {
      return index< Entity::codimension >( entity );
    }

    template< int cd >
    IndexType index ( const typename Traits::template Codim< cd >::Entity &entity ) const
    {
      return index( Grid::template getHostEntity< cd >( entity ) );
    }

    template< int cd >
    IndexType index ( const typename Traits::HostGrid::template Codim< cd >::Entity &entity ) const
    {
      return hostIndexSet().index( entity );
    }

    template< class Entity >
    IndexType subIndex ( const Entity &entity, int i, unsigned int codim ) const
    {
      return subIndex< Entity::codimension >( entity, i, codim );
    }

    template< int cd >
    IndexType subIndex ( const typename Traits::template Codim< cd >::Entity &entity, int i, unsigned int codim ) const
    {
      return subIndex( Grid::template getHostEntity< cd >( entity ), i, codim );
    }

    template< int cd >
    IndexType subIndex ( const typename Traits::HostGrid::template Codim< cd >::Entity &entity, int i, unsigned int codim ) const
    {
      return hostIndexSet().subIndex( entity, i, codim );
    }

    IndexType size ( GeometryType type ) const
    {
      return hostIndexSet().size( type );
    }

    int size ( int codim ) const
    {
      return hostIndexSet().size( codim );
    }

    template< class Entity >
    bool contains ( const Entity &entity ) const
    {
      static const int cc = Entity::codimension;
      return hostIndexSet().contains( Grid::template getHostEntity< cc >( entity ) );
    }

    const std::vector< GeometryType > &geomTypes ( int codim ) const
    {
      return hostIndexSet().geomTypes( codim );
    }

    operator bool () const { return bool( hostIndexSet_ ); }

  protected:
    const HostIndexSet &hostIndexSet () const
    {
      assert( hostIndexSet_ );
      return *hostIndexSet_;
    }

    const HostIndexSet *hostIndexSet_;
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_INDEXSET_HH
