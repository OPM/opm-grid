#ifndef DUNE_CPGRID_PERSISTENTCONTAINER_HH
#define DUNE_CPGRID_PERSISTENTCONTAINER_HH

#include <dune/common/version.hh>

#include <dune/grid/utility/persistentcontainer.hh>
#include <opm/grid/CpGrid.hpp>

#include <dune/grid/utility/persistentcontainervector.hh>

namespace Dune
{
  // PersistentContainer for CpGrid
  // -------------------------------
  template< class Data >
  class PersistentContainer< CpGrid, Data >
  : public PersistentContainerVector< CpGrid,
                                      typename CpGrid::Traits::LeafIndexSet,
                                      std::vector<Data> >
  {
  public:
    typedef CpGrid  GridType;
    typedef typename std::vector<Data>::allocator_type Allocator;

  private:
    typedef PersistentContainerVector< GridType, typename GridType::Traits::LeafIndexSet, std::vector<Data> > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const GridType &grid, const int codim, const Data& data = Data(), const Allocator &allocator = Allocator() )
    : BaseType( grid.leafIndexSet(), codim, data, allocator )
    {}
  };

} // end namespace Dune

#endif // end DUNE_CPGRID_PERSISTENTCONTAINER_HH
