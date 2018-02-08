// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_PERSISTENTCONTAINER_HH
#define DUNE_POLYHEDRALGRID_PERSISTENTCONTAINER_HH

#include <dune/common/version.hh>

#include <dune/grid/utility/persistentcontainer.hh>
#include <opm/grid/polyhedralgrid/grid.hh>

#include <dune/grid/utility/persistentcontainervector.hh>

namespace Dune
{
  // PersistentContainer for CpGrid
  // -------------------------------
  template< int dim, int dimworld, class Data >
  class PersistentContainer< PolyhedralGrid< dim, dimworld >, Data >
  : public PersistentContainerVector< PolyhedralGrid< dim, dimworld >,
                                      typename PolyhedralGrid< dim, dimworld >::Traits::LeafIndexSet,
                                      std::vector<Data> >
  {
  public:
    typedef PolyhedralGrid< dim, dimworld >  GridType;
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

#endif
