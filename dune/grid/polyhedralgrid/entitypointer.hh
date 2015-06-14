// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
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

  template< class Grid, int codim >
  class PolyhedralGridEntityPointer
  {
    typedef PolyhedralGridEntityPointer< Grid, codim > This;

  protected:
    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    /** \brief grid dimension */
    static const int dimension = remove_const< Grid >::type::dimension;
    /** \brief world dimension */
    static const int codimension = codim;

    /** \brief type of entity */
    typedef typename Traits::template Codim< codimension >::Entity Entity;

  protected:
    typedef typename Traits::ExtraData ExtraData;

    typedef typename Traits::template Codim< codimension > :: EntityImpl EntityImpl;

  public:
    PolyhedralGridEntityPointer ( ExtraData data )
    : entity_( EntityImpl( data ) )
    {}

    explicit PolyhedralGridEntityPointer ( const EntityImpl &entity )
    : entity_( EntityImpl( entity ) )
    {}

    PolyhedralGridEntityPointer ( const This &other )
    : entity_( other.entity_ )
    {}

    const This &operator= ( const This &other )
    {
      entity_ = other.entity_;
      return *this;
    }

    /** \brief check for equality */
    bool equals ( const This &other ) const
    {
      return (entity_ == other.entity_);
    }

    /** \brief dereference entity */
    Entity& dereference () const
    {
      return entity_;
    }

    /** \brief obtain level */
    int level () const { return entity_.level(); }

  protected:
    EntityImpl &entityImpl () const
    {
      return Grid::getRealImplementation( entity_ );
    }

    ExtraData data () const { return entityImpl().data(); }

  protected:
    mutable Entity entity_;
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_ENTITYPOINTER_HH
