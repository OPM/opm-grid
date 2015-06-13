#ifndef DUNE_POLYHEDRALGRID_ENTITYSEED_HH
#define DUNE_POLYHEDRALGRID_ENTITYSEED_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/common/entityseed.hh>

namespace Dune
{

  template< int codim, class Grd >
  class PolyhedralGridEntitySeed
  {
    typedef typename remove_const< Grd >::type::Traits Traits;

  public:
    static const int codimension = codim;
    static const int dimension = Traits::dimension;
    static const int mydimension = dimension - codimension;
    static const int dimensionworld = Traits::dimensionworld;


    typedef typename Traits::Grid Grid;
    typedef typename Traits::template Codim< codim >::Entity Entity;
    typedef typename Traits :: Index Index ;

    typedef typename Traits::HostGrid HostGrid;
    typedef typename HostGrid::template Codim< codim >::EntitySeed HostEntitySeed;

    static const Index defaulIndex = -1;

    explicit PolyhedralGridEntitySeed ( const Index& index )
      : index_( index )
    {}

    PolyhedralGridEntitySeed ()
      : index_( defaultIndex )
    {}

    bool isValid() const { return index_ != defaultIndex; }

  protected:
    Index index_;
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_ENTITYSEED_HH
