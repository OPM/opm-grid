// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_GEOMETRY_HH
#define DUNE_POLYHEDRALGRID_GEOMETRY_HH

#include <dune/common/version.hh>
#include <dune/common/fmatrix.hh>
#include <dune/grid/common/geometry.hh>

#if DUNE_VERSION_NEWER(DUNE_GEOMETRY, 2, 5 )
#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>
#else
#include <dune/geometry/genericgeometry/geometrytraits.hh>
#include <dune/geometry/genericgeometry/matrixhelper.hh>
#endif


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

    typedef typename Grid::Traits::ExtraData  ExtraData;
    typedef typename Grid::Traits::template Codim<codimension>::EntitySeed EntitySeed;

    template< class ct >
    struct PolyhedralMultiLinearGeometryTraits
      : public Dune::MultiLinearGeometryTraits< ct >
    {
      struct Storage
      {
        struct Iterator
          : public Dune::ForwardIteratorFacade< Iterator, GlobalCoordinate, GlobalCoordinate >
        {
          const Storage* data_;
          int count_;
          explicit Iterator( const Storage* ptr, int count ) : data_( ptr ), count_( count ) {}

          GlobalCoordinate dereference() const { return data_->corner( count_ ); }
          void increment() { ++count_; }

          bool equals( const Iterator& other ) const { return count_ == other.count_; }
        };

        ExtraData  data_;
        // host geometry object
        EntitySeed seed_;

        Storage( ExtraData data, EntitySeed seed )
          : data_( data ), seed_( seed )
        {}

        Storage( ExtraData data )
          : data_( data ), seed_()
        {}

        ExtraData data() const { return data_; }
        bool isValid () const { return seed_.isValid(); }

        GlobalCoordinate operator [] (const int i) const { return corner( i ); }

        Iterator begin() const { return Iterator(this, 0); }
        Iterator end ()  const { return Iterator(this, corners()); }

        int corners () const { return data()->corners( seed_ ); }
        GlobalCoordinate corner ( const int i ) const { return data()->corner( seed_, i ); }
        GlobalCoordinate center () const { return data()->centroids( seed_ ); }

        ctype volume() const { return data()->volumes( seed_ ); }
      };

      template <int mdim, int cordim>
      struct CornerStorage
      {
        typedef Storage Type;
      };
    };

    typedef Dune::MultiLinearGeometry< ctype, mydimension, coorddimension, PolyhedralMultiLinearGeometryTraits<ctype> >
      MultiLinearGeometryType;

    typedef typename PolyhedralMultiLinearGeometryTraits< ctype > ::template
      CornerStorage<mydimension, coorddimension >::Type CornerStorageType;

    //! type of jacobian inverse transposed
    typedef FieldMatrix< ctype, cdim, mydim > JacobianInverseTransposed;

    //! type of jacobian transposed
    typedef FieldMatrix< ctype, mydim, cdim > JacobianTransposed;


#if DUNE_VERSION_NEWER(DUNE_GRID,2,5)
    typedef Dune::Impl::FieldMatrixHelper< ctype >  MatrixHelperType;
#else
    typedef Dune::GenericGeometry::MatrixHelper< Dune::GenericGeometry::DuneCoordTraits< ctype > >  MatrixHelperType;
#endif

    explicit PolyhedralGridBasicGeometry ( ExtraData data )
    : storage_( data )
    {}

    PolyhedralGridBasicGeometry ( ExtraData data, const EntitySeed& seed )
    : storage_( data, seed )
    {
      GeometryType myType = type();
      if( ! myType.isNone() && storage_.isValid() )
      {
        if( myType.isLine() && storage_.corners() != 2 )
        {
          std::cout << myType << "  " << storage_.corners() << std::endl;
          std::abort();
        }

        //std::cout << myType << "  " << storage_.corners() << std::endl;
        geometryImpl_.reset( new MultiLinearGeometryType(myType, storage_) );
      }
    }

    GeometryType type () const { return data()->geomTypes(codimension)[0]; }
    bool affine () const { return (geometryImpl_) ? geometryImpl_->affine() : false; }

    int corners () const { return storage_.corners(); }
    GlobalCoordinate corner ( const int i ) const { return storage_.corner( i ); }
    GlobalCoordinate center () const { return storage_.center(); }

    GlobalCoordinate global(const LocalCoordinate& local) const
    {
      if( geometryImpl_ )
      {
        return geometryImpl_->global( local );
      }

      return center();
    }

    /// Mapping from the cell to the reference domain.
    /// May be slow.
    LocalCoordinate local(const GlobalCoordinate& global) const
    {
      if( geometryImpl_ )
      {
        return geometryImpl_->local( global );
      }

      // if no geometry type return a vector filled with 1
      return LocalCoordinate( 1 );
    }

    ctype integrationElement ( const LocalCoordinate &local ) const
    {
      if( geometryImpl_ )
      {
        return geometryImpl_->integrationElement( local );
      }
      return volume();
    }

    ctype volume () const
    {
      if( geometryImpl_ )
      {
        return geometryImpl_->volume();
      }
      return storage_.volume();
    }

#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
    JacobianTransposed jacobianTransposed ( const LocalCoordinate & local ) const
    {
      if( geometryImpl_ )
      {
        return geometryImpl_->jacobianTransposed( local );
      }

      DUNE_THROW(NotImplemented,"jacobianTransposed not implemented");
      return JacobianTransposed( 0 );
    }

    JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate & local ) const
    {
      if( geometryImpl_ )
      {
        return geometryImpl_->jacobianInverseTransposed( local );
      }

      DUNE_THROW(NotImplemented,"jacobianInverseTransposed not implemented");
      return JacobianInverseTransposed( 0 );
    }
#else
    const JacobianTransposed& jacobianTransposed ( const LocalCoordinate &local ) const
    {
      if( geometryImpl_ )
      {
        return geometryImpl_->jacobianTransposed( local );
      }

      DUNE_THROW(NotImplemented,"jacobianTransposed not implemented");
      static const JacobianTransposed jac( 0 );
      return jac;
    }

    const JacobianInverseTransposed& jacobianInverseTransposed ( const LocalCoordinate &local ) const
    {
      if( geometryImpl_ )
      {
        return geometryImpl_->jacobianInverseTransposed( local );
      }

      DUNE_THROW(NotImplemented,"jacobianInverseTransposed not implemented");
      static const JacobianInverseTransposed jac( 0 );
      return jac;
    }
#endif

    ExtraData data() const { return storage_.data(); }

  protected:
    CornerStorageType storage_;
    std::shared_ptr< MultiLinearGeometryType > geometryImpl_;
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
