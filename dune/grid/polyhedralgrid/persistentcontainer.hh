// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_PERSISTENTCONTAINER_HH
#define DUNE_POLYHEDRALGRID_PERSISTENTCONTAINER_HH

#include <dune/common/version.hh>

#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/grid/polyhedralgrid/grid.hh>

#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 3)
#include <dune/grid/utility/persistentcontainervector.hh>
#endif

namespace Dune
{
  // PersistentContainer for CpGrid
  // -------------------------------
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 3)
  template< int dim, int dimworld, class Data >
  class PersistentContainer< PolyhedralGrid< dim, dimworld >, Data >
  : public PersistentContainerVector< PolyhedralGrid< dim, dimworld >,
                                      typename PolyhedralGrid< dim, dimworld >::Traits::LeafIndexSet,
                                      std::vector<Data> >
#else
  template< int dim, int dimworld, class Data, class Allocator >
  class PersistentContainer< PolyhedralGrid< dim, dimworld >, Data, Allocator >
  : public PersistentContainerVector< PolyhedralGrid< dim, dimworld >,
                                      typename PolyhedralGrid< dim, dimworld >::Traits::LeafIndexSet,
                                      std::vector<Data,Allocator> >
#endif
  {
  public:
    typedef PolyhedralGrid< dim, dimworld >  GridType;
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 3)
    typedef typename std::vector<Data>::allocator_type Allocator;
#endif

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
    PersistentContainer ( const GridType &grid, const int codim, const Data& data = Data(), const Allocator &allocator = Allocator() )
    : BaseType( grid.leafIndexSet(), codim, data, allocator )
#else
    PersistentContainer ( const GridType &grid, const int codim, const Allocator &allocator = Allocator() )
    : BaseType( grid, codim, grid.leafIndexSet(), 1.1, allocator )
#endif
    {}
  };

} // end namespace Dune

#endif
