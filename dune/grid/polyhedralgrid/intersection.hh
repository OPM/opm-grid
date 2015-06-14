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
    typedef PolyhedralGridIntersection< Grid > This;
  protected:
    typedef typename Grid :: Traits Traits;

    typedef typename Traits :: ExtraData ExtraData ;

  public:
    typedef typename Traits::ctype ctype;
    typedef typename Traits::GlobalCoordinate GlobalCoordinate;

    static const int dimension = Traits::dimension;
    static const int dimensionworld = Traits::dimensionworld;

    typedef typename Traits::template Codim< 0 >::Entity Entity;
    typedef typename Traits::template Codim< 0 >::EntityPointer EntityPointer;
    typedef typename Traits::template Codim< 0 >::EntitySeed    EntitySeed;
    typedef typename Traits::template Codim< 1 >::Geometry Geometry;
    typedef typename Traits::template Codim< 1 >::LocalGeometry LocalGeometry;

  protected:
    typedef typename Traits::template Codim< 0 >::EntityImpl EntityImpl;
    typedef typename Traits::template Codim< 1 >::GeometryImpl GeometryImpl;
    //typedef typename Traits::ExtraData  ExtraData;

  public:
    explicit PolyhedralGridIntersection ( ExtraData data )
    : data_( data ),
      seed_(),
      intersectionIdx_( -1 )
    {}

    explicit PolyhedralGridIntersection ( ExtraData data, const EntitySeed& seed, const int intersectionIdx )
    : data_( data ),
      seed_( seed ),
      intersectionIdx_( intersectionIdx )
    {}

    const EntityImpl inside () const
    {
        return EntityImpl( data(), seed_ );
    }

    EntityImpl outside () const
    {
        return data()->neighbor(seed_, intersectionIdx_);
    }

    bool operator == ( const This& other ) const
    {
      return (seed_ == other.seed_) &&
             (intersectionIdx_ == other.intersectionIdx_);
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
      return Geometry(GeometryImpl(data(), data()->template subEntitySeed<1>(seed_, intersectionIdx_)));
    }

    GeometryType type () const { return GeometryType(GeometryType::cube, dimension); }

    int indexInInside () const
    {
        return data().indexInInside(seed_, intersectionIdx_);
    }

    GlobalCoordinate
    integrationOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
    {
      return outerNormal( local );
    }

    GlobalCoordinate
    outerNormal ( const FieldVector< ctype, dimension-1 > &local ) const
    { return outerNormal(); }

    GlobalCoordinate outerNormal () const
    { return data().outerNormal(seed_, intersectionIdx_); }

    GlobalCoordinate
    unitOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
    {
      return centerUnitOuterNormal();
    }

    GlobalCoordinate
    centerUnitOuterNormal () const
    {
      GlobalCoordinate normal( outerNormal() );
      normal /= normal.two_norm();
      return std::move( normal );
    }

    ExtraData data() const { return data_; }

    bool equals(const This& other) const
    {
      return seed_.equals(other.seed_) && intersectionIdx_ == other.intersectionIdx_;
    }

  protected:
    ExtraData  data_;
    EntitySeed seed_;
  public:
    int intersectionIdx_; // the element-local index
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_INTERSECTION_HH
