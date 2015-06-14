// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_GEOMETRY_HH
#define DUNE_POLYHEDRALGRID_GEOMETRY_HH

#include <dune/common/fmatrix.hh>
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
    static const int dimension = Grid::dimension;
    static const int mydimension = mydim;
    static const int codimension = dimension - mydimension;

    static const int dimensionworld = Grid::dimensionworld;
    static const int coorddimension = dimensionworld;

    typedef typename Grid::ctype ctype;
    typedef Dune::FieldVector< ctype, coorddimension > GlobalCoordinate;
    typedef Dune::FieldVector< ctype, mydimension >    LocalCoordinate;

    //! type of jacobian inverse transposed
    typedef FieldMatrix<ctype,cdim,mydim> JacobianInverseTransposed;

    //! type of jacobian transposed
    typedef FieldMatrix< ctype, mydim, cdim > JacobianTransposed;

    typedef typename Grid::Traits::ExtraData  ExtraData;
    typedef typename Grid::Traits::template Codim<codimension>::EntitySeed EntitySeed;

    explicit PolyhedralGridBasicGeometry ( ExtraData data )
    : data_( data ),
      seed_( )
    {}

    PolyhedralGridBasicGeometry ( ExtraData data, const EntitySeed& seed )
    : data_( data ),
      seed_( seed )
    {}

    GeometryType type () const { return GeometryType( GeometryType::cube, mydimension ); }
    bool affine () const { return false; }

    int corners () const { return data()->corners( seed_ ); }
    GlobalCoordinate corner ( const int i ) const { return data()->corner( seed_, i ); }
    GlobalCoordinate center () const { return data()->centroids( seed_ ); }

    GlobalCoordinate global ( const LocalCoordinate &local   ) const
    {
      DUNE_THROW(NotImplemented,"global not implemented");
      return GlobalCoordinate( 0 );
    }

    LocalCoordinate  local  ( const GlobalCoordinate &global ) const
    {
      DUNE_THROW(NotImplemented,"local not implemented");
      return LocalCoordinate( 0 );
    }

    ctype integrationElement ( const LocalCoordinate &local ) const { return volume(); }
    ctype volume () const { return data()->volumes( seed_ ); }

    JacobianTransposed jacobianTransposed ( const LocalCoordinate &local ) const
    {
      DUNE_THROW(NotImplemented,"jacobianTransposed not implemented");
      return JacobianTransposed( 0 );
    }

    JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate &local ) const
    {
      DUNE_THROW(NotImplemented,"jacobianInverseTransposed not implemented");
      return JacobianInverseTransposed( 0 );
    }

    ExtraData data() const { return data_; }

  protected:
    ExtraData  data_;
    // host geometry object
    EntitySeed seed_;
  };


  // PolyhedralGridGeometry
  // --------------

  template< int mydim, int cdim, class Grid >
  class PolyhedralGridGeometry
  : public PolyhedralGridBasicGeometry< mydim, cdim, Grid >
  {
    typedef PolyhedralGridBasicGeometry< mydim, cdim, Grid > Base;

  public:
    typedef typename Base::ExtraData  ExtraData;
    typedef typename Base::EntitySeed EntitySeed;

    explicit PolyhedralGridGeometry ( ExtraData data )
    : Base( data )
    {}

    PolyhedralGridGeometry ( ExtraData data, const EntitySeed& seed )
    : Base( data, seed )
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
    typedef typename Base::ExtraData  ExtraData;
    typedef typename Base::EntitySeed EntitySeed;

    explicit PolyhedralGridLocalGeometry ( ExtraData data )
    : Base( data )
    {}

    PolyhedralGridLocalGeometry ( ExtraData data, const EntitySeed& seed )
    : Base( data, seed )
    {}
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_GEOMETRY_HH