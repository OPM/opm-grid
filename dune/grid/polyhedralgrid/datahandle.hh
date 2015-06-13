#ifndef DUNE_POLYHEDRALGRID_DATAHANDLE_HH
#define DUNE_POLYHEDRALGRID_DATAHANDLE_HH

//- dune-common includes
#include <dune/common/typetraits.hh>

//- dune-grid includes
#include <dune/grid/common/datahandleif.hh>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/common/ldbhandleif.hh>
#endif

namespace Dune
{

  // PolyhedralGridDataHandleBase
  // ---------------------------

#if HAVE_DUNE_ALUGRID
  template< class WrappedHandle, bool hasCompressAndReserve = Conversion< WrappedHandle, LoadBalanceHandleWithReserveAndCompress >::exists >
  class PolyhedralGridDataHandleBase;
#else // #if HAVE_DUNE_ALUGRID
  template< class WrappedHandle, bool hasCompressAndReserve = false >
  class PolyhedralGridDataHandleBase;
#endif // #else // #if HAVE_DUNE_ALUGRID

#if HAVE_DUNE_ALUGRID
  template< class WrappedHandle >
  struct PolyhedralGridDataHandleBase< WrappedHandle, true >
    : public LoadBalanceHandleWithReserveAndCompress
  {
    explicit PolyhedralGridDataHandleBase ( WrappedHandle &handle )
      : wrappedHandle_( handle )
    {}

    void reserveMemory ( std::size_t newElements )
    {
      return wrappedHandle_.reserveMemory( newElements );
    }

    void compress () { return wrappedHandle_.compress(); }

  protected:
    WrappedHandle &wrappedHandle_;
  };
#endif // #if HAVE_DUNE_ALUGRID

  template< class WrappedHandle >
  struct PolyhedralGridDataHandleBase< WrappedHandle, false >
  {
    explicit PolyhedralGridDataHandleBase ( WrappedHandle &handle )
      : wrappedHandle_( handle )
    {}

  protected:
    WrappedHandle &wrappedHandle_;
  };


  // PolyhedralGridDataHandle
  // ----------------

  template< class WrappedHandle, class Grid >
  class PolyhedralGridDataHandle
  : public PolyhedralGridDataHandleBase< WrappedHandle >,
    public CommDataHandleIF< PolyhedralGridDataHandle< WrappedHandle, Grid >, typename WrappedHandle::DataType >
  {
  protected:
    typedef PolyhedralGridDataHandle< WrappedHandle, Grid > This;
    typedef PolyhedralGridDataHandleBase< WrappedHandle >   Base;

    using Base :: wrappedHandle_ ;

    // type of traits
    typedef typename remove_const< Grid >::type::Traits Traits;

    typedef typename Traits :: ExtraDataType ExtraData ;

    template< int codim >
    struct Codim
    {
      // type of entity
      typedef typename Traits::template Codim< codim >::Entity Entity;
    };

  public:
    // type of data to be communicated
    typedef typename WrappedHandle::DataType DataType;

    typedef CommDataHandleIF< This, DataType > DataHandleIF;

  private:
    // prohibit copying
    PolyhedralGridDataHandle ( const This & );

  public:
    PolyhedralGridDataHandle ( ExtraData data, WrappedHandle &wrappedHandle )
    : Base( wrappedHandle ),
      data_( data )
    {}

    bool contains ( int dim, int codim ) const
    {
      return wrappedHandle_.contains( dim, codim );
    }

    bool fixedsize ( int dim, int codim ) const
    {
      return wrappedHandle_.fixedsize( dim, codim );
    }

    template< class HostEntity >
    size_t size ( const HostEntity &hostEntity ) const
    {
      typedef typename Codim< HostEntity::codimension >::Entity Entity;
      const Entity entity( typename Entity::Implementation( data(), hostEntity ) );
      return wrappedHandle_.size( entity );
    }

    template< class MessageBuffer, class HostEntity >
    void gather ( MessageBuffer &buffer, const HostEntity &hostEntity ) const
    {
      typedef typename Codim< HostEntity::codimension >::Entity Entity;
      const Entity entity( typename Entity::Implementation( data(), hostEntity ) );
      wrappedHandle_.gather( buffer, entity );
    }

    template< class MessageBuffer, class HostEntity >
    void scatter ( MessageBuffer &buffer, const HostEntity &hostEntity, size_t size )
    {
      typedef typename Codim< HostEntity::codimension >::Entity Entity;
      const Entity entity( typename Entity::Implementation( data(), hostEntity ) );
      wrappedHandle_.scatter( buffer, entity, size );
    }

    ExtraData data() const { return data_; }

  protected:
    ExtraData data_;
  };



  // PolyhedralGridWrappedDofManager
  // -----------------------

  template< class WrappedDofManager, class Grid >
  class PolyhedralGridWrappedDofManager
  {
    typedef PolyhedralGridWrappedDofManager< WrappedDofManager, Grid > This;

    typedef typename remove_const< Grid >::type::Traits Traits;

    // type of element (i.e., entity of codimension 0)
    typedef typename Traits::template Codim< 0 >::Entity Element;

    // type of host element (i.e., host entity of codimension 0)
    typedef typename Traits::HostGrid::template Codim< 0 >::Entity HostElement;

  private:
    // prohibit copy constructor
    PolyhedralGridWrappedDofManager ( const This & );

  public:
    explicit PolyhedralGridWrappedDofManager ( WrappedDofManager &wrappedDofManager )
    : wrappedDofManager_( wrappedDofManager )
    {}

    template< class MessageBuffer >
    void inlineData ( MessageBuffer &buffer, const HostElement &hostElement )
    {
      const Element element( (typename Element::Implementation)( hostElement ) );
      wrappedDofManager_.inlineData( buffer, element );
    }

    template< class MessageBuffer >
    void xtractData ( MessageBuffer &buffer, const HostElement &hostElement, std::size_t newElements )
    {
      const Element element( (typename Element::Implementation)( hostElement ) );
      wrappedDofManager_.xtractData( buffer, element, newElements );
    }

    void compress () { wrappedDofManager_.compress(); }

  private:
    WrappedDofManager &wrappedDofManager_;
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_DATAHANDLE_HH
