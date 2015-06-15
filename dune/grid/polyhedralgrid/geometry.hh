// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_GEOMETRY_HH
#define DUNE_POLYHEDRALGRID_GEOMETRY_HH

#include <dune/common/version.hh>
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

    GeometryType type () const { return data()->geomTypes(codimension)[0]; }
    bool affine () const { return false; }

    int corners () const { return data()->corners( seed_ ); }
    GlobalCoordinate corner ( const int i ) const { return data()->corner( seed_, i ); }
    GlobalCoordinate center () const { return data()->centroids( seed_ ); }

    GlobalCoordinate global(const LocalCoordinate& local) const
    {
      assert( mydimension == 3 );
      assert( coorddimension == 3 );
      // uvw = { (1-u, 1-v, 1-w), (u, v, w) }
      LocalCoordinate uvw[2] = { LocalCoordinate(1.0), local };
      uvw[0] -= local;
      // Access pattern for uvw matching ordering of corners.
      const int pat[8][3] = { { 0, 0, 0 },
                              { 1, 0, 0 },
                              { 0, 1, 0 },
                              { 1, 1, 0 },
                              { 0, 0, 1 },
                              { 1, 0, 1 },
                              { 0, 1, 1 },
                              { 1, 1, 1 } };
      GlobalCoordinate xyz(0.0);
      for (int i = 0; i < 8; ++i) {
        GlobalCoordinate corner_contrib = corner(i);
        double factor = 1.0;
        for (int j = 0; j < 3; ++j) {
          factor *= uvw[pat[i][j]][j];
        }
        corner_contrib *= factor;
        xyz += corner_contrib;
      }
      return xyz;
    }

    /// Mapping from the cell to the reference domain.
    /// May be slow.
    LocalCoordinate local(const GlobalCoordinate& y) const
    {
      /*
      assert( mydimension == 3 );
      assert( coorddimension == 3 );
      // This code is modified from dune/grid/genericgeometry/mapping.hh
      // \todo: Implement direct computation.
      const ctype epsilon = 1e-12;
      const ReferenceElement< ctype , mydimension > & refElement =
        ReferenceElements< ctype, mydimension >::general(type());

      LocalCoordinate x = refElement.position(0,0);
      LocalCoordinate dx;
      do {
        using namespace GenericGeometry;
        // DF^n dx^n = F^n, x^{n+1} -= dx^n
        JacobianTransposed JT = jacobianTransposed(x);
        GlobalCoordinate z = global(x);
        z -= y;
        MatrixHelper<DuneCoordTraits<double> >::template xTRightInvA<3, 3>(JT, z, dx );
        x -= dx;
      } while (dx.two_norm2() > epsilon*epsilon);
      return x;
      */
      return LocalCoordinate( 0 );
    }

    ctype integrationElement ( const LocalCoordinate &local ) const { return volume(); }
    ctype volume () const { return data()->volumes( seed_ ); }

#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
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
#else
    const JacobianTransposed& jacobianTransposed ( const LocalCoordinate &local ) const
    {
      DUNE_THROW(NotImplemented,"jacobianTransposed not implemented");
      static const JacobianTransposed jac( 0 );
      return jac;
    }

    const JacobianInverseTransposed& jacobianInverseTransposed ( const LocalCoordinate &local ) const
    {
      DUNE_THROW(NotImplemented,"jacobianInverseTransposed not implemented");
      static const JacobianInverseTransposed jac( 0 );
      return jac;
    }
#endif

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

  template< int mydim, int cdim, class Grid >
  class PolyhedralGridLocalGeometry
  : public PolyhedralGridBasicGeometry< mydim, cdim, Grid >
  {
    typedef PolyhedralGridBasicGeometry< mydim, cdim, Grid > Base;

  public:
    typedef typename Base::ExtraData  ExtraData;

    explicit PolyhedralGridLocalGeometry ( ExtraData data )
    : Base( data )
    {}
  };


#if ! DUNE_VERSION_NEWER(DUNE_GRID,2,4)
  namespace FacadeOptions
  {

    //! \brief Traits class determining whether the Dune::Geometry facade
    //!        class stores the implementation object by reference or by value
    template< int mydim, int cdim, class GridImp >
    struct StoreGeometryReference< mydim, cdim, GridImp, PolyhedralGridGeometry >
    {
      //! Whether to store by reference.
      static const bool v = false;
    };

    template< int mydim, int cdim, class GridImp >
    struct StoreGeometryReference< mydim, cdim, GridImp, PolyhedralGridLocalGeometry >
    {
      //! Whether to store by reference.
      static const bool v = false;
    };

  }
#endif


} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_GEOMETRY_HH
