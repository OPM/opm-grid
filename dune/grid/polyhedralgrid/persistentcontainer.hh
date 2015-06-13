#ifndef DUNE_POLYHEDRALGRID_PERSISTENTCONTAINER_HH
#define DUNE_POLYHEDRALGRID_PERSISTENTCONTAINER_HH

#include <dune/grid/utility/persistentcontainerwrapper.hh>

#include <dune/grid/polyhedralgrid/declaration.hh>

namespace Dune
{

  // PersistentContainer for PolyhedralGrid
  // ------------------------------

  template< int dim, int dimworld, class T >
  class PersistentContainer< PolyhedralGrid< dim, dimworld >, T >
  : public PersistentContainerWrapper< PolyhedralGrid< dim, dimworld >, T >
  {
    typedef PersistentContainerWrapper< PolyhedralGrid< dim, dimworld >, T > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
    : Base( grid, codim, value )
    {}
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_PERSISTENTCONTAINER_HH
