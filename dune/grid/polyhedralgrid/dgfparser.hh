#ifndef DUNE_POLYHEDRALGRID_DGFPARSER_HH
#define DUNE_POLYHEDRALGRID_DGFPARSER_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include <dune/grid/polyhedralgrid/grid.hh>
#include <dune/grid/polyhedralgrid/hostgridaccess.hh>

namespace Dune
{

  // DGFGridFactory for PolyhedralGrid
  // -------------------------

  template< class HostGrid >
  struct DGFGridFactory< PolyhedralGrid< HostGrid > >
  {
    typedef PolyhedralGrid< HostGrid > Grid;

    const static int dimension = Grid::dimension;
    typedef MPIHelper::MPICommunicator MPICommunicator;
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<dimension>::Entity Vertex;

    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicator comm = MPIHelper::getCommunicator() )
    : dgfHostFactory_( input, comm ),
      grid_( 0 )
    {
      HostGrid *hostGrid = dgfHostFactory_.grid();
      assert( hostGrid != 0 );
      grid_ = new Grid( *hostGrid );
    }

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicator comm = MPIHelper::getCommunicator() )
    : dgfHostFactory_( filename, comm ),
      grid_( 0 )
    {
      HostGrid *hostGrid = dgfHostFactory_.grid();
      assert( hostGrid != 0 );
      std::ifstream input( filename.c_str() );
      grid_ = new Grid( *hostGrid );
    }

    Grid *grid () const
    {
      return grid_;
    }

    template< class Intersection >
    bool wasInserted ( const Intersection &intersection ) const
    {
      return dgfHostFactory_.wasInserted( HostGridAccess< Grid >::hostIntersection( intersection ) );
    }

    template< class Intersection >
    int boundaryId ( const Intersection &intersection ) const
    {
      return dgfHostFactory_.boundaryId( HostGridAccess< Grid >::hostIntersection( intersection ) );
    }

    bool haveBoundaryParameters () const
    {
      return dgfHostFactory_.haveBoundaryParameters();
    }

    template< int codim >
    int numParameters () const
    {
      return dgfHostFactory_.template numParameters< codim >();
    }

    template< class Intersection >
    const typename DGFBoundaryParameter::type &
    boundaryParameter ( const Intersection &intersection ) const
    {
      return dgfHostFactory_.boundaryParameter( HostGridAccess< Grid >::hostIntersection( intersection ) );
    }

    template< class Entity >
    std::vector< double > &parameter ( const Entity &entity )
    {
      return dgfHostFactory_.parameter( HostGridAccess< Grid >::hostEntity( entity ) );
    }

  private:
    DGFGridFactory< HostGrid > dgfHostFactory_;
    Grid *grid_;
  };



  // DGFGridInfo for PolyhedralGrid
  // ----------------------

  template< class HostGrid >
  struct DGFGridInfo< PolyhedralGrid< HostGrid > >
  {
    static int refineStepsForHalf ()
    {
      return DGFGridInfo< HostGrid >::refineStepsForHalf();
    }

    static double refineWeight ()
    {
      return DGFGridInfo< HostGrid >::refineWeight();
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_DGFPARSER_HH
