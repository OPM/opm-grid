#ifndef DUNE_POLYHEDRALGRID_ENTITY_HH
#define DUNE_POLYHEDRALGRID_ENTITY_HH

//- dune-common includes
#include <dune/common/nullptr.hh>
#include <dune/common/typetraits.hh>

//- dune-grid includes
#include <dune/grid/common/entity.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class, class > class PolyhedralGridIntersection;
  template< class, class > class PolyhedralGridIntersectionIterator;

  template< class, class > class PolyhedralGridIterator;


  // PolyhedralGridEntityBasic
  // -----------------

  /** \copydoc PolyhedralGridEntity
   *
   *  \nosubgrouping
   */
  template< int codim, int dim, class Grid >
  class PolyhedralGridEntityBasic
  {
  protected:
    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    /** \name Attributes
     *  \{ */

    //! codimensioon of the entity
    static const int codimension = codim;
    //! dimension of the grid
    static const int dimension = Traits::dimension;
    //! dimension of the entity
    static const int mydimension = dimension - codimension;
    //! dimension of the world
    static const int dimensionworld = Traits::dimensionworld;

    /** \} */

    /** \name Types Required by DUNE
     *  \{ */

    //! coordinate type of the grid
    typedef typename Traits::ctype ctype;

    //! type of corresponding entity seed
    typedef typename Grid::template Codim< codimension >::EntitySeed EntitySeed;
    //! type of corresponding geometry
    typedef typename Traits::template Codim< codimension >::Geometry Geometry;

    /** \} */

  protected:
    // type of the host grid
    typedef typename Traits::HostGrid  HostGrid;

    // type of extra data, e.g. a pointer to grid (here empty)
    typedef typename Traits::ExtraDataType ExtraData;

  public:
    /** \name Construction, Initialization and Destruction
     *  \{ */

    /** \brief construct a null entity */
    explicit PolyhedralGridEntityBasic ( ExtraData data )
    : seed_()
      data_( data )
    {}

    /** \brief construct an initialized entity
     *
     *  \param[in]  hostEntity  corresponding entity in the host grid
     *
     *  \note The reference to the host entity must remain valid  as long as
     *        this entity is in use.
     */
    PolyhedralGridEntityBasic ( ExtraData data, const EntitySeed& seed )
    : seed_( seed )
    , data_( data )
    {}

    /** \} */

    /** \brief return true if entity hold a vaild host entity */
    operator bool () const { return seed_.isValid(); }

    /** \name Methods Shared by Entities of All Codimensions
     *  \{ */

    /** \brief obtain the name of the corresponding reference element
     *
     *  This type can be used to access the DUNE reference element.
     */
    GeometryType type () const
    {
      return GeometryType();
    }

    /** \brief obtain the level of this entity */
    int level () const
    {
      return 0;
    }

    /** \brief obtain the partition type of this entity */
    PartitionType partitionType () const
    {
      return data->partitionType( *this );
    }

    /** obtain the geometry of this entity */
    Geometry geometry () const
    {
      return Geometry( seed_ );
    }

    /** \brief return EntitySeed of host grid entity */
    EntitySeed seed () const { return seed_; }

    /** \} */


    /** \name Methods Supporting the Grid Implementation
     *  \{ */
    ExtraData data() const { return data_; }

    /** \} */

  protected:
    EntitySeed  seed_;
    ExtraData   data_;
  };



  // PolyhedralGridEntity
  // ------------

  template< int codim, int dim, class Grid >
  class PolyhedralGridEntity : public PolyhedralGridEntityBasic< codim, dim, Grid >
  {
    typedef PolyhedralGridEntityBasic< codim, dim, Grid > Base ;
  protected:
    typedef typename remove_const< Grid >::type::Traits Traits;

  protected:
    // type of the host grid
    typedef typename Traits::HostGrid  HostGrid;

    // type of extra data, e.g. a pointer to grid (here empty)
    typedef typename Traits::ExtraDataType ExtraData;

  public:
    using Base :: codimension ;

    /** \name Host Types
     *  \{ */

    //! type of corresponding host entity
    typedef typename HostGrid::template Codim< codimension >::Entity HostEntity;
    /** \} */

    explicit PolyhedralGridEntity ( ExtraData data )
    : Base( data )
    {}

    PolyhedralGridEntity ( ExtraData data, const EntitySeed& seed )
    : Base( data, seed )
    {}
  };


  // PolyhedralGridEntity for codimension 0
  // ----------------------------------

  /** \copydoc PolyhedralGridEntity
   *
   *  \nosubgrouping
   */
  template< int dim, class Grid >
  class PolyhedralGridEntity< 0, dim, Grid > : public PolyhedralGridEntityBasic< 0, dim, Grid >
  {
    typedef PolyhedralGridEntityBasic< 0, dim, Grid > Base ;
  protected:
    typedef typename Base::Traits Traits;
    typedef typename Base::HostGrid HostGrid;

    // type of extra data, e.g. a pointer to grid (here empty)
    typedef typename Base::ExtraData ExtraData;

  public:
    using Base::codimension ;
    using Base::data ;
    using Base::hostEntity ;
    /** \name Host Types
     *  \{ */

    //! type of corresponding host entity
    typedef typename HostGrid::template Codim< codimension >::Entity HostEntity;
    /** \} */

  protected:
    typedef typename Traits :: LeafIntersectionIteratorImpl  LeafIntersectionIteratorImpl;
    typedef typename Traits :: LevelIntersectionIteratorImpl LevelIntersectionIteratorImpl;

    typedef typename Traits :: LeafIntersectionImpl          LeafIntersectionImpl;
    typedef typename Traits :: LevelIntersectionImpl         LevelIntersectionImpl;

  public:
    /** \name Types Required by DUNE
     *  \{ */

    //! type of corresponding local geometry
    typedef typename Traits::template Codim< codimension >::LocalGeometry LocalGeometry;
    //! type of corresponding entity pointer
    typedef typename Traits::template Codim< codimension >::EntityPointer EntityPointer;

    //! type of hierarchic iterator
    typedef typename Traits::HierarchicIterator        HierarchicIterator;
    //! type of leaf intersection iterator
    typedef typename Traits::LeafIntersectionIterator  LeafIntersectionIterator;
    //! type of level intersection iterator
    typedef typename Traits::LevelIntersectionIterator LevelIntersectionIterator;

    /** \} */

    /** \name Construction, Initialization and Destruction
     *  \{ */

    /** \brief construct a null entity */
    explicit PolyhedralGridEntity ( ExtraData data )
    : Base( data )
    {}

    /** \brief construct an initialized entity
     *
     *  \param[in]  hostEntity  corresponding entity in the host grid
     *
     *  \note The reference to the host entity must remain valid as long as
     *        this entity is in use.
     */
    PolyhedralGridEntity ( ExtraData data, const EntitySeed& seed )
    : Base( data, seed )
    {}

    /** \} */

    unsigned int subEntities( const unsigned int codim ) const
    {
      if( codim == 0 )
        return 1;
      else
        return data->subEntities( seed_ );
    }

    template< int codim >
    int count () const
    {
      return subEntities( codim );
    }

    template< int codim >
    typename Grid::template Codim< codim >::EntityPointer
    subEntity ( int i ) const
    {
      typedef typename Traits::template Codim< codim >::EntityPointerImpl EntityPointerImpl;
      return EntityPointerImpl( data(), data->subEntity( seed_, codim, i ) );
    }

    bool hasBoundaryIntersections () const
    {
      const int faces = subEntities( 1 );
      for( int i=0; i<faces; ++i )
      {
        if( ! this->template subEntity< 1 >( i ) )
          return true;
      }
      return false;
    }

    bool isLeaf () const
    {
      return true;
    }

    EntityPointer father () const
    {
      DUNE_THROW(InvalidStateException,"no father available");
      typedef typename Traits::template Codim< 0 >::EntityPointerImpl EntityPointerImpl;
      return EntityPointerImpl( data() );
    }

    bool hasFather () const
    {
      return false;
    }

    LocalGeometry geometryInFather () const
    {
      DUNE_THROW(InvalidStateException,"no father available");
      return LocalGeometry( hostEntity().geometryInFather() );
    }

    HierarchicIterator hbegin ( int maxLevel ) const
    {
      typedef typename Traits :: HierarchicIteratorImpl HierarchicIteratorImpl ;
      return HierarchicIteratorImpl( data() );
    }

    HierarchicIterator hend ( int maxLevel ) const
    {
      typedef typename Traits :: HierarchicIteratorImpl HierarchicIteratorImpl ;
      return HierarchicIteratorImpl( data() );
    }

    bool isRegular () const
    {
      return true;
    }

    bool isNew () const
    {
      return false;
    }

    bool mightVanish () const
    {
      return false;
    }

    /** \} */

  };


} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_ENTITY_HH
