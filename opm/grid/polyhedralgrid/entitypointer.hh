// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_ENTITYPOINTER_HH
#define DUNE_POLYHEDRALGRID_ENTITYPOINTER_HH

//- dune-grid includes
#include <dune/grid/common/grid.hh>

//- dune-metagrid includes
#include <opm/grid/polyhedralgrid/declaration.hh>

namespace Dune
{
  // PolyhedralGridEntityPointer
  // -------------------

  template< int codim, class Grid >
  class PolyhedralGridEntityPointer
  {
    typedef PolyhedralGridEntityPointer< codim, Grid > This;

  protected:
    typedef typename Grid::Traits Traits;

  public:
    /** \brief grid dimension */
    static const int dimension = Grid::dimension;
    /** \brief world dimension */
    static const int codimension = codim;

    /** \brief type of entity */
    typedef typename Traits::template Codim< codimension >::Entity Entity;

  protected:
    typedef typename Traits::ExtraData ExtraData;

    typedef typename Traits::template Codim< codimension > :: EntityImpl EntityImpl;

  public:
    explicit PolyhedralGridEntityPointer ( ExtraData data )
    : entity_( EntityImpl( data ) )
    {}

    explicit PolyhedralGridEntityPointer ( const EntityImpl &entity )
    : entity_( EntityImpl( entity ) )
    {}

    PolyhedralGridEntityPointer ( const This &other )
    : entity_( EntityImpl( other.entityImpl() ) )
    {}

    const This &operator= ( const This &other )
    {
      entityImpl() = other.entityImpl();
      return *this;
    }

    /** \brief check for equality */
    bool equals ( const This &other ) const
    {
      return entityImpl().equals( other.entityImpl() );
    }

    /** \brief dereference entity */
    Entity& dereference () const
    {
      return entity_;
    }

    operator const Entity& () const { return entity_; }
    operator       Entity& ()       { return entity_; }

    /** \brief obtain level */
    int level () const { return entity_.level(); }

  protected:
    EntityImpl &entityImpl () const
    {
      return entity_.impl();
    }

    ExtraData data () const { return entityImpl().data(); }

  protected:
    mutable Entity entity_;
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_ENTITYPOINTER_HH
