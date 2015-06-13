// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_INTERSECTION_HH
#define DUNE_POLYHEDRALGRID_INTERSECTION_HH

//- dune-common includes
#include <dune/common/nullptr.hh>

//- local includes
#include <dune/grid/polyhedralgrid/declaration.hh>

namespace Dune
{

  // PolyhedralGridIntersection
  // ------------------

  template< class Grid >
  class PolyhedralGridIntersection
  {
  protected:
    typedef typename remove_const< Grid >::type::Traits Traits;

    typedef typename Traits :: ExtraDataType ExtraData ;
    typedef typename Traits::template Codim< 0 >::EntityPointerImpl EntityPointerImpl;
  public:
    typedef typename Traits::ctype ctype;

    static const int dimension = Traits::dimension;
    static const int dimensionworld = Traits::dimensionworld;

    typedef typename Traits::template Codim< 0 >::Entity Entity;
    typedef typename Traits::template Codim< 0 >::EntityPointer EntityPointer;
    typedef typename Traits::template Codim< 1 >::Geometry Geometry;
    typedef typename Traits::template Codim< 1 >::LocalGeometry LocalGeometry;

  public:
    explicit PolyhedralGridIntersection ( ExtraData data )
    : data_( data )
    {}

    const EntityImpl inside () const
    {
        return EntityImpl( data(), seed_ );
    }

    EntityImpl outside () const
    {
        return data()->neighbor(seed_, intersectionIdx_);
    }

    bool boundary () const { return !neighbor(); }

    bool conforming () const { return false; }

    bool neighbor () const { return !data()->neighbor(seed_, intersectionIdx_).isValid(); }

    int boundaryId () const { return 1; }

    size_t boundarySegmentIndex () const
    {
        // not implementable
        return -1;
    }

    LocalGeometry geometryInInside () const
    {
      return LocalGeometry( data() );
    }

    LocalGeometry geometryInOutside () const
    {
        return LocalGeometryImpl(data());
    }

    Geometry geometry () const
    {
        return Geometry( data(), data()->template subEntitySeed<1>(seed_, intersectionIdx_) );
    }

    GeometryType type () const { return GeometryType::none; }

    int indexInInside () const
    { return faceTagToSubentityNum(data().faceTag(seed_, intersectionIdx_)); }

    int indexInInside () const
    {
        return data().indexInInside(seed_, intersectionIdx_));
    }

    FieldVector< ctype, dimensionworld >
    integrationOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
    {
        return data().indexInOutside(seed_, intersectionIdx_));
    }

    FieldVector< ctype, dimensionworld >
    outerNormal ( const FieldVector< ctype, dimension-1 > &local ) const
    { return data().outerNormal(seed_, intersectionIdx_); }

    FieldVector< ctype, dimensionworld >
    unitOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
    {
      return data().unitOuterNormal(seed_, intersectionIdx_); }
    }

    FieldVector< ctype, dimensionworld > centerUnitOuterNormal () const
    {
      return data().unitOuterNormal(seed_, intersectionIdx_); }
    }

    ExtraData data() const { return data_; }

  protected:
    ExtraData  data_;
    EntitySeed seed_;
    int intersectionIdx_; // the element-local index
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_INTERSECTION_HH
