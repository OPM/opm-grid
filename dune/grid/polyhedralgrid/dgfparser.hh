// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_DGFPARSER_HH
#define DUNE_POLYHEDRALGRID_DGFPARSER_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include <dune/grid/polyhedralgrid/grid.hh>

namespace Dune
{

#warning TODO: non-trivial DGFGridFactory

  // DGFGridFactory for PolyhedralGrid
  // -------------------------

  template< int dim, int dimworld >
  struct DGFGridFactory< PolyhedralGrid< dim, dimworld > >
  {
    typedef PolyhedralGrid< dim, dimworld > Grid;

    const static int dimension = Grid::dimension;
    typedef MPIHelper::MPICommunicator MPICommunicator;
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<dimension>::Entity Vertex;

    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicator comm = MPIHelper::getCommunicator() )
    : grid_( nullptr )
    {
    }

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicator comm = MPIHelper::getCommunicator() )
    : grid_( nullptr )
    {
    }

    Grid *grid () const
    {
      return grid_;
    }

    template< class Intersection >
    bool wasInserted ( const Intersection &intersection ) const
    {
        return false;
    }

    template< class Intersection >
    int boundaryId ( const Intersection &intersection ) const
    {
        return false;
    }

    bool haveBoundaryParameters () const
    {
        return false;
    }

    template< int codim >
    int numParameters () const
    {
        return 0;
    }

    template< class Intersection >
    const typename DGFBoundaryParameter::type &
    boundaryParameter ( const Intersection &intersection ) const
    {
        return DGFBoundaryParameter::defaultValue();
    }

    template< class Entity >
    std::vector< double > &parameter ( const Entity &entity )
    {
        static std::vector<double> dummy;
        return dummy;
    }

  private:
    Grid *grid_;
  };



  // DGFGridInfo for PolyhedralGrid
  // ----------------------

  template< int dim, int dimworld >
  struct DGFGridInfo< PolyhedralGrid< dim, dimworld > >
  {
    static int refineStepsForHalf ()
    {
        return 0;
    }

    static double refineWeight ()
    {
        return 0;
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_DGFPARSER_HH
