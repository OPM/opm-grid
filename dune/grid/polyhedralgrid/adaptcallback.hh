#ifndef DUNE_POLYHEDRALGRID_ADAPTCALLBACK_HH
#define DUNE_POLYHEDRALGRID_ADAPTCALLBACK_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/common/adaptcallback.hh>

namespace Dune
{

  // PolyhedralGridAdaptDataHandle
  // ---------------------

  template< class Grid, class WrappedHandle >
  class PolyhedralGridAdaptDataHandle
    : public AdaptDataHandle< typename Grid::HostGrid, PolyhedralGridAdaptDataHandle< Grid, WrappedHandle > >
  {
  protected:
    typedef PolyhedralGridAdaptDataHandle< Grid, WrappedHandle > This;

    typedef typename remove_const< Grid >::type::Traits Traits;
    typedef typename Traits :: ExtraDataType  ExtraData;

    typedef typename Traits::HostGrid::template Codim< 0 >::Entity HostEntity;

    typedef typename Traits::template Codim< 0 >::EntityImpl  EntityImpl;

  public:
    typedef typename Traits::template Codim< 0 >::Entity      Entity;

  private:
    PolyhedralGridAdaptDataHandle ( const This & );
    This &operator= ( const This & );

  public:
    PolyhedralGridAdaptDataHandle( ExtraData data, WrappedHandle &handle )
    : wrappedHandle_( handle ),
      data_( data )
    {}

    void preAdapt ( unsigned int estimateAdditionalElements )
    {
      wrappedHandle_.preAdapt( estimateAdditionalElements );
    }

    void postAdapt () { wrappedHandle_.postAdapt(); }

    void preCoarsening ( const HostEntity &hostFather ) const
    {
      const Entity father( EntityImpl( data_, hostFather ) );
      wrappedHandle_.preCoarsening( father );
    }

    void postRefinement ( const HostEntity &hostFather ) const
    {
      const Entity father( EntityImpl( data_, hostFather ) );
      wrappedHandle_.postRefinement( father );
    }

  protected:
    WrappedHandle &wrappedHandle_;
    ExtraData data_;
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_ADAPTCALLBACK_HH
