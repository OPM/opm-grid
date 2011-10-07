#ifndef DUNE_CPGRID_PERSISTENTCONTAINER_HH
#define DUNE_CPGRID_PERSISTENTCONTAINER_HH

#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/grid/CpGrid.hpp>

namespace Dune
{
  // PersistentContainer for CpGrid
  // -------------------------------
  
  template< class Data, class Allocator >
  class PersistentContainer< CpGrid, Data, Allocator >
  : public PersistentContainerVector< CpGrid, 
                                      typename CpGrid::Traits::LeafIndexSet,
                                      std::vector<Data,Allocator> >
  {
  public:
    typedef CpGrid  GridType;
  private:
    typedef PersistentContainerVector< GridType, typename GridType::Traits::LeafIndexSet, std::vector<Data,Allocator> > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor 
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const GridType &grid, const int codim, const Allocator &allocator = Allocator() )
    : BaseType( grid, codim, grid.leafIndexSet(), 1.1, allocator )
    {}
  };

} // end namespace Dune

#endif // end DUNE_CPGRID_PERSISTENTCONTAINER_HH
