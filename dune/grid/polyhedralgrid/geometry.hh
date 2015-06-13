#ifndef DUNE_POLYHEDRALGRID_GEOMETRY_HH
#define DUNE_POLYHEDRALGRID_GEOMETRY_HH

#include <dune/grid/common/geometry.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< int, int, class > class PolyhedralGridGeometry;
  template< int, int, class > class PolyhedralGridLocalGeometry;

  // PolyhedralGridBasicGeometry
  // -------------------

  template< int mydim, int cdim, class Grid >
  struct PolyhedralGridBasicGeometry
  {
    typedef typename remove_const< Grid >::type::HostGrid HostGrid;

    static const int dimension = HostGrid::dimension;
    static const int mydimension = mydim;
    static const int codimension = dimension - mydimension;
    typedef typename HostGrid::template Codim< codimension >::Geometry HostGeometry;

    static const int coorddimension = HostGeometry::coorddimension;
    static const int dimensionworld = HostGeometry::dimensionworld;

    typedef typename HostGeometry::ctype ctype;
    typedef typename HostGeometry::LocalCoordinate  LocalCoordinate;
    typedef typename HostGeometry::GlobalCoordinate GlobalCoordinate;

    typedef typename HostGeometry::JacobianTransposed JacobianTransposed;
    typedef typename HostGeometry::JacobianInverseTransposed JacobianInverseTransposed;

    PolyhedralGridBasicGeometry ( const HostGeometry &hostGeometry )
    : hostGeometry_( hostGeometry )
    {}

    GeometryType type () const { return hostGeometry().type(); }
    bool affine () const { return hostGeometry().affine(); }

    int corners () const { return hostGeometry().corners(); }
    GlobalCoordinate corner ( const int i ) const { return hostGeometry().corner( i ); }
    GlobalCoordinate center () const { return hostGeometry().center(); }

    GlobalCoordinate global ( const LocalCoordinate &local   ) const { return hostGeometry().global( local ); }
    LocalCoordinate  local  ( const GlobalCoordinate &global ) const { return hostGeometry().local( global ); }

    ctype integrationElement ( const LocalCoordinate &local ) const { return hostGeometry().integrationElement( local ); }
    ctype volume () const { return hostGeometry().volume(); }

    JacobianTransposed jacobianTransposed ( const LocalCoordinate &local ) const
    {
      return hostGeometry().jacobianTransposed( local );
    }

    JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local ) const
    {
      return hostGeometry().jacobianInverseTransposed( local );
    }

  protected:
    const HostGeometry &hostGeometry () const { return hostGeometry_; }

    // host geometry object
    HostGeometry hostGeometry_;
  };


  // PolyhedralGridGeometry
  // --------------

  template< int mydim, int cdim, class Grid >
  class PolyhedralGridGeometry
  : public PolyhedralGridBasicGeometry< mydim, cdim, Grid >
  {
    typedef PolyhedralGridBasicGeometry< mydim, cdim, Grid > Base;

  public:
    typedef typename Base::HostGeometry HostGeometry;

    PolyhedralGridGeometry ()
    {}

    PolyhedralGridGeometry ( const HostGeometry &hostGeometry )
    : Base( hostGeometry )
    {}
  };


  // PolyhedralGridLocalGeometry
  // -------------------

  template< int mydim, int cdim, class Grid >
  class PolyhedralGridLocalGeometry
  : public PolyhedralGridBasicGeometry< mydim, cdim, Grid >
  {
    typedef PolyhedralGridBasicGeometry< mydim, cdim, Grid >  Base ;

  public:
    typedef typename Base::HostGeometry HostGeometry;

    PolyhedralGridLocalGeometry ()
    {}

    PolyhedralGridLocalGeometry ( const HostGeometry &hostGeometry )
    : Base( hostGeometry )
    {}
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_GEOMETRY_HH
