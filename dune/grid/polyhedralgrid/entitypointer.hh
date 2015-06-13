#ifndef DUNE_POLYHEDRALGRID_ENTITYPOINTER_HH
#define DUNE_POLYHEDRALGRID_ENTITYPOINTER_HH

//- dune-grid includes
#include <dune/grid/common/grid.hh>

//- dune-metagrid includes
#include <dune/grid/polyhedralgrid/declaration.hh>

namespace Dune
{
  // PolyhedralGridEntityPointer
  // -------------------

  template< class Grid, class HostIterator >
  class PolyhedralGridEntityPointer
  {
    typedef PolyhedralGridEntityPointer< Grid, HostIterator > This;

  protected:
    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    /** \brief grid dimension */
    static const int dimension = remove_const< Grid >::type::dimension;
    /** \brief world dimension */
    static const int codimension = HostIterator::codimension;

    /** \brief type of entity */
    typedef typename Traits::template Codim< codimension >::Entity Entity;

  protected:
    typedef typename Traits::HostGrid::template Codim< codimension >::EntityPointer HostEntityPointer;

    typedef typename Traits::ExtraDataType ExtraData;

    typedef typename Traits::template Codim< codimension > :: EntityImpl EntityImpl;

  private:
    friend class PolyhedralGridEntityPointer< Grid, HostEntityPointer >;

  public:
    PolyhedralGridEntityPointer ( ExtraData data, const HostIterator &hostIterator )
    : entity_( EntityImpl( data ) ),
      hostIterator_( hostIterator )
    {}

    explicit PolyhedralGridEntityPointer ( const EntityImpl &entity )
    : entity_( EntityImpl( entity.data() ) ),
      hostIterator_( entity.hostEntity() )
    {}

    PolyhedralGridEntityPointer ( const This &other )
    : entity_( EntityImpl( other.data() ) ),
      hostIterator_( other.hostIterator_ )
    {}

    template< class HI >
    explicit PolyhedralGridEntityPointer ( const PolyhedralGridEntityPointer< Grid, HI > &other )
    : entity_( EntityImpl( other.data() ) ),
      hostIterator_( other.hostIterator_ )
    {}

    const This &operator= ( const This &other )
    {
      entityImpl() = EntityImpl( other.data() );
      hostIterator_ = other.hostIterator_;
      return *this;
    }

    template< class HI >
    const This &operator= ( const PolyhedralGridEntityPointer< Grid, HI > &other )
    {
      entityImpl() = EntityImpl( other.data() );
      hostIterator_ = other.hostIterator_;
      return *this;
    }

    /** \brief check for equality */
    template< class HI >
    bool equals ( const PolyhedralGridEntityPointer< Grid, HI > &other ) const
    {
      return (hostIterator() == other.hostIterator());
    }

    /** \brief dereference entity */
    Entity &dereference () const
    {
      if( !entityImpl() )
        entityImpl() = EntityImpl( data(), *hostIterator() );
      return entity_;
    }

    /** \brief obtain level */
    int level () const { return hostIterator().level(); }

    /** \brief obtain host iterator */
    const HostIterator &hostIterator() const { return hostIterator_; }

  protected:
    void releaseEntity () { entityImpl() = EntityImpl( data() ); }

    EntityImpl &entityImpl () const
    {
      return Grid::getRealImplementation( entity_ );
    }

    ExtraData data () const { return entityImpl().data(); }

  protected:
    mutable Entity entity_;
    HostIterator hostIterator_;
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_ENTITYPOINTER_HH
