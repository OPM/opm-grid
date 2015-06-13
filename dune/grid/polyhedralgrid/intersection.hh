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

  template< class Grid, class HostIntersection >
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
    : hostIntersection_( nullptr ),
      data_( data )
    {}

    PolyhedralGridIntersection ( ExtraData data, const HostIntersection &hostIntersection )
    : hostIntersection_( &hostIntersection ),
      data_( data )
    {}

    operator bool () const { return bool( hostIntersection_ ); }

    const EntityPointer inside () const
    {
      return EntityPointer( EntityPointerImpl( data(), hostIntersection().inside() ) );
    }

    EntityPointer outside () const
    {
      return EntityPointer( EntityPointerImpl( data(), hostIntersection().outside() ) );
    }

    bool boundary () const { return hostIntersection().boundary(); }

    bool conforming () const { return hostIntersection().conforming(); }

    bool neighbor () const { return hostIntersection().neighbor(); }

    int boundaryId () const { return hostIntersection().boundaryId(); }

    size_t boundarySegmentIndex () const
    {
      return hostIntersection().boundarySegmentIndex();
    }

    LocalGeometry geometryInInside () const
    {
      return LocalGeometry( hostIntersection().geometryInInside() );
    }

    LocalGeometry geometryInOutside () const
    {
      return LocalGeometry( hostIntersection().geometryInOutside() );
    }

    Geometry geometry () const
    {
      return Geometry( hostIntersection().geometry() );
    }

    GeometryType type () const { return hostIntersection().type(); }

    int indexInInside () const { return hostIntersection().indexInInside(); }
    int indexInOutside () const { return hostIntersection().indexInOutside(); }

    FieldVector< ctype, dimensionworld >
    integrationOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
    {
      return hostIntersection().integrationOuterNormal( local );
    }

    FieldVector< ctype, dimensionworld >
    outerNormal ( const FieldVector< ctype, dimension-1 > &local ) const
    {
      return hostIntersection().outerNormal( local );
    }

    FieldVector< ctype, dimensionworld >
    unitOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
    {
      return hostIntersection().unitOuterNormal( local );
    }

    FieldVector< ctype, dimensionworld > centerUnitOuterNormal () const
    {
      return hostIntersection().centerUnitOuterNormal();
    }

    const HostIntersection &hostIntersection () const
    {
      assert( *this );
      return *hostIntersection_;
    }

    ExtraData data() const { return data_; }

  protected:
    const HostIntersection *hostIntersection_;
    ExtraData  data_;
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_INTERSECTION_HH
