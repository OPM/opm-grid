#ifndef DUNE_CPGRID_PERSISTENTCONTAINER_HH
#define DUNE_CPGRID_PERSISTENTCONTAINER_HH

#include <dune/common/version.hh>

#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/grid/CpGrid.hpp>

#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 3)
#include <dune/grid/utility/persistentcontainervector.hh>
#endif

namespace Dune
{
  // PersistentContainer for CpGrid
  // -------------------------------
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 3)
  template< class Data >
  class PersistentContainer< CpGrid, Data >
  : public PersistentContainerVector< CpGrid,
                                      typename CpGrid::Traits::LeafIndexSet,
                                      std::vector<Data> >
#else
  template< class Data, class Allocator >
  class PersistentContainer< CpGrid, Data, Allocator >
  : public PersistentContainerVector< CpGrid, 
                                      typename CpGrid::Traits::LeafIndexSet,
                                      std::vector<Data,Allocator> >
#endif
  {
  public:
    typedef CpGrid  GridType;
  private:
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 3)
    typedef PersistentContainerVector< GridType, typename GridType::Traits::LeafIndexSet, std::vector<Data> > BaseType;
#else
    typedef PersistentContainerVector< GridType, typename GridType::Traits::LeafIndexSet, std::vector<Data,Allocator> > BaseType;
#endif

  public:
    //! Constructor filling the container with values using the default constructor 
    //! Depending on the implementation this could be achieved without allocating memory
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 3)
    PersistentContainer ( const GridType &grid, const int codim)
    : BaseType( grid, codim, grid.leafIndexSet(), 1.1 )
#else
    PersistentContainer ( const GridType &grid, const int codim, const Allocator &allocator = Allocator() )
    : BaseType( grid, codim, grid.leafIndexSet(), 1.1, allocator )
#endif
    {}
  };

} // end namespace Dune

#endif // end DUNE_CPGRID_PERSISTENTCONTAINER_HH
