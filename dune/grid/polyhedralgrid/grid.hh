#ifndef DUNE_POLYHEDRALGRID_GRID_HH
#define DUNE_POLYHEDRALGRID_GRID_HH

//- dune-common includes
#include <dune/common/nullptr.hh>

//- dune-grid includes
#include <dune/grid/common/grid.hh>

//- dune-metagrid includes
#include <dune/grid/common/hostgridinoutstreams.hh>
#include <dune/grid/polyhedralgrid/adaptcallback.hh>
#include <dune/grid/polyhedralgrid/backuprestore.hh>
#include <dune/grid/polyhedralgrid/capabilities.hh>
#include <dune/grid/polyhedralgrid/datahandle.hh>
#include <dune/grid/polyhedralgrid/declaration.hh>
#include <dune/grid/polyhedralgrid/entity.hh>
#include <dune/grid/polyhedralgrid/entitypointer.hh>
#include <dune/grid/polyhedralgrid/entityseed.hh>
#include <dune/grid/polyhedralgrid/geometry.hh>
#include <dune/grid/polyhedralgrid/gridview.hh>
#include <dune/grid/polyhedralgrid/idset.hh>

#include <opm/core/grid.h>
#include <opm/core/grid/GridManager.hpp>

namespace Dune
{

  // PolyhedralGridExportParams
  // ------------------

  template< class HG >
  struct PolyhedralGridExportParams
  : public HostGridHasInOutStreams< HG, Conversion< HG, HasObjectStream >::exists >
  {
    typedef HG HostGrid;
  };



  // PolyhedralGridFamily
  // ------------

  template< class HostGrid >
  struct PolyhedralGridFamily
  {
    struct Traits
    : public PolyhedralGridExportParams< HostGrid >
    {
      typedef PolyhedralGrid< HostGrid > Grid;

      typedef typename HostGrid::ctype ctype;

      struct EmptyData{};

      // type of data passed to entities, intersections, and iterators
      // for PolyhedralGrid this is just an empty place holder
      typedef EmptyData ExtraDataType;

      static const int dimension = HostGrid::dimension;
      static const int dimensionworld = HostGrid::dimensionworld;

      typedef PolyhedralGridIntersection< const Grid, typename HostGrid::LeafIntersection > LeafIntersectionImpl;
      typedef Dune::Intersection< const Grid, LeafIntersectionImpl > LeafIntersection;

      typedef PolyhedralGridIntersection< const Grid, typename HostGrid::LevelIntersection > LevelIntersectionImpl;
      typedef Dune::Intersection< const Grid, LevelIntersectionImpl > LevelIntersection;

      typedef PolyhedralGridIntersectionIterator< const Grid, typename HostGrid::LeafIntersectionIterator > LeafIntersectionIteratorImpl;
      typedef Dune::IntersectionIterator< const Grid, LeafIntersectionIteratorImpl, LeafIntersectionImpl > LeafIntersectionIterator;

      typedef PolyhedralGridIntersectionIterator< const Grid, typename HostGrid::LevelIntersectionIterator > LevelIntersectionIteratorImpl;
      typedef Dune::IntersectionIterator< const Grid, LevelIntersectionIteratorImpl, LevelIntersectionImpl > LevelIntersectionIterator;

      typedef PolyhedralGridIterator< const Grid, typename HostGrid::HierarchicIterator > HierarchicIteratorImpl;
      typedef Dune::EntityIterator< 0, const Grid, HierarchicIteratorImpl > HierarchicIterator;

      template< int codim >
      struct Codim
      {
        typedef Dune::Geometry< dimension-codim, dimensionworld, const Grid, PolyhedralGridGeometry > Geometry;
        typedef Dune::Geometry< dimension-codim, dimension, const Grid, PolyhedralGridLocalGeometry > LocalGeometry;

        typedef PolyhedralGridEntityPointer< const Grid, typename HostGrid::template Codim< codim >::EntityPointer > EntityPointerImpl;
        typedef Dune::EntityPointer< const Grid, EntityPointerImpl > EntityPointer;

        typedef PolyhedralGridEntity< codim, dimension, const Grid > EntityImpl;
        typedef Dune::Entity< codim, dimension, const Grid, PolyhedralGridEntity > Entity;

        typedef Dune::EntitySeed< const Grid, PolyhedralGridEntitySeed< codim, const Grid > > EntitySeed;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef typename HostGrid::template Codim< codim >::template Partition< pitype > HostPartition;

          typedef PolyhedralGridIterator< const Grid, typename HostPartition::LeafIterator > LeafIteratorImpl;
          typedef Dune::EntityIterator< codim, const Grid, LeafIteratorImpl > LeafIterator;

          typedef PolyhedralGridIterator< const Grid, typename HostPartition::LevelIterator > LevelIteratorImpl;
          typedef Dune::EntityIterator< codim, const Grid, LevelIteratorImpl > LevelIterator;
        };

        typedef typename Partition< All_Partition >::LeafIterator LeafIterator;
        typedef typename Partition< All_Partition >::LevelIterator LevelIterator;
      };

      typedef PolyhedralGridIndexSet< const Grid, typename HostGrid::Traits::LeafIndexSet > LeafIndexSet;
      typedef PolyhedralGridIndexSet< const Grid, typename HostGrid::Traits::LevelIndexSet > LevelIndexSet;

      typedef PolyhedralGridIdSet< const Grid, typename HostGrid::Traits::GlobalIdSet > GlobalIdSet;
      typedef PolyhedralGridIdSet< const Grid, typename HostGrid::Traits::LocalIdSet > LocalIdSet;

      typedef typename HostGrid::Traits::CollectiveCommunication CollectiveCommunication;

      template< PartitionIteratorType pitype >
      struct Partition
      {
        typedef Dune::GridView< PolyhedralGridViewTraits< typename HostGrid::template Partition< pitype >::LeafGridView, pitype > > LeafGridView;
        typedef Dune::GridView< PolyhedralGridViewTraits< typename HostGrid::template Partition< pitype >::LevelGridView, pitype > > LevelGridView;
      };
    };
  };



  // PolyhedralGrid
  // ------

  /** \class PolyhedralGrid
   *  \brief identical grid wrapper
   *  \ingroup PolyhedralGrid
   *
   *  \tparam  HostGrid   DUNE grid to be wrapped (called host grid)
   *
   *  \nosubgrouping
   */
  template < int dim, int dimworld >
  class PolyhedralGrid
  /** \cond */
  : public GridDefaultImplementation
      < HostGrid::dimension, HostGrid::dimensionworld, typename HostGrid::ctype, PolyhedralGridFamily< HostGrid > >,
    public PolyhedralGridExportParams< HostGrid >,
    public PolyhedralGridBackupRestoreFacilities< PolyhedralGrid< HostGrid > >
  /** \endcond */
  {
    typedef PolyhedralGrid< HostGrid > Grid;

    typedef GridDefaultImplementation
      < HostGrid::dimension, HostGrid::dimensionworld, typename HostGrid::ctype, PolyhedralGridFamily< HostGrid > >
      Base;

    template< int, int, class > friend class PolyhedralGridEntity;
    template< class, class > friend class PolyhedralGridEntityPointer;
    template< class, class > friend class PolyhedralGridIntersection;
    template< class, class > friend class PolyhedralGridIntersectionIterator;
    template< class, class > friend class PolyhedralGridIdSet;
    template< class, class > friend class PolyhedralGridIndexSet;
    template< class > friend class HostGridAccess;

  public:
    /** \cond */
    typedef PolyhedralGridFamily< HostGrid > GridFamily;
    /** \endcond */

    /** \name Traits
     *  \{ */

    //! type of the grid traits
    typedef typename GridFamily::Traits Traits;

    /** \brief traits structure containing types for a codimension
     *
     *  \tparam codim  codimension
     *
     *  \nosubgrouping
     */
    template< int codim >
    struct Codim;

    /** \} */

    /** \name Iterator Types
     *  \{ */

    //! iterator over the grid hierarchy
    typedef typename Traits::HierarchicIterator HierarchicIterator;
    //! iterator over intersections with other entities on the leaf level
    typedef typename Traits::LeafIntersectionIterator LeafIntersectionIterator;
    //! iterator over intersections with other entities on the same level
    typedef typename Traits::LevelIntersectionIterator LevelIntersectionIterator;

    /** \} */

    /** \name Grid View Types
     *  \{ */

    /** \brief Types for GridView */
    template< PartitionIteratorType pitype >
    struct Partition
    {
      typedef typename GridFamily::Traits::template Partition< pitype >::LevelGridView LevelGridView;
      typedef typename GridFamily::Traits::template Partition< pitype >::LeafGridView LeafGridView;
    };

    /** \brief View types for All_Partition */
    typedef typename Partition< All_Partition >::LevelGridView LevelGridView;
    typedef typename Partition< All_Partition >::LeafGridView LeafGridView;

    /** \} */

    /** \name Index and Id Set Types
     *  \{ */

    /** \brief type of leaf index set
     *
     *  The index set assigns consecutive indices to the entities of the
     *  leaf grid. The indices are of integral type and can be used to access
     *  arrays.
     *
     *  The leaf index set is a model of Dune::IndexSet.
     */
    typedef typename Traits::LeafIndexSet LeafIndexSet;

    /** \brief type of level index set
     *
     *  The index set assigns consecutive indices to the entities of a grid
     *  level. The indices are of integral type and can be used to access
     *  arrays.
     *
     *  The level index set is a model of Dune::IndexSet.
     */
    typedef typename Traits::LevelIndexSet LevelIndexSet;

    /** \brief type of global id set
     *
     *  The id set assigns a unique identifier to each entity within the
     *  grid. This identifier is unique over all processes sharing this grid.
     *
     *  \note Id's are neither consecutive nor necessarily of an integral
     *        type.
     *
     *  The global id set is a model of Dune::IdSet.
     */
    typedef typename Traits::GlobalIdSet GlobalIdSet;

    /** \brief type of local id set
     *
     *  The id set assigns a unique identifier to each entity within the
     *  grid. This identifier needs only to be unique over this process.
     *
     *  Though the local id set may be identical to the global id set, it is
     *  often implemented more efficiently.
     *
     *  \note Ids are neither consecutive nor necessarily of an integral
     *        type.
     *  \note Local ids need not be compatible with global ids. Also, no
     *        mapping from local ids to global ones needs to exist.
     *
     *  The global id set is a model of Dune::IdSet.
     */
    typedef typename Traits::LocalIdSet LocalIdSet;

    /** \} */

    /** \name Miscellaneous Types
     * \{ */

    //! type of vector coordinates (e.g., double)
    typedef typename Traits::ctype ctype;

    //! communicator with all other processes having some part of the grid
    typedef typename Traits::CollectiveCommunication CollectiveCommunication;

    /** \} */

    /** \name Construction and Destruction
     *  \{ */

    /** \brief constructor
     *
     *  The references to host grid and coordinate function are stored in the
     *  grid. Therefore, they must remain valid until the grid is destroyed.
     *
     *  \param[in]  hostGrid       reference to the grid to wrap
     */
    explicit PolyhedralGrid ( Opm::DeckConstPtr deck,
                              const  std::vector<double>& poreVolumes = std::vector<double> ())
    : grid_( createGrid( deck, poreVolumes )
      // levelIndexSets_( hostGrid.maxLevel()+1, nullptr )
    {}

    /** \brief destructor
     */
    ~PolyhedralGrid ()
    {
      for( unsigned int i = 0; i < levelIndexSets_.size(); ++i )
      {
        if( levelIndexSets_[ i ] )
          delete( levelIndexSets_[ i ] );
      }
    }

    /** \} */

    /** \name Size Methods
     *  \{ */

    /** \brief obtain maximal grid level
     *
     *  Grid levels are numbered 0, ..., L, where L is the value returned by
     *  this method.
     *
     *  \returns maximal grid level
     */
    int maxLevel () const
    {
      return hostGrid().maxLevel();
    }

    /** \brief obtain number of entites on a level
     *
     *  \param[in]  level  level to consider
     *  \param[in]  codim  codimension to consider
     *
     *  \returns number of entities of codimension \em codim on grid level
     *           \em level.
     */
    int size ( int level, int codim ) const
    {
      return hostGrid().size( level, codim );
    }

    /** \brief obtain number of leaf entities
     *
     *  \param[in]  codim  codimension to consider
     *
     *  \returns number of leaf entities of codimension \em codim
     */
    int size ( int codim ) const
    {
      return hostGrid().size( codim );
    }

    /** \brief obtain number of entites on a level
     *
     *  \param[in]  level  level to consider
     *  \param[in]  type   geometry type to consider
     *
     *  \returns number of entities with a geometry of type \em type on grid
     *           level \em level.
     */
    int size ( int level, GeometryType type ) const
    {
      return hostGrid().size( level, type );
    }

    /** \brief returns the number of boundary segments within the macro grid
     *
     *  \returns number of boundary segments within the macro grid
     */
    int size ( GeometryType type ) const
    {
      return hostGrid().size( type );
    }

    /** \brief obtain number of leaf entities
     *
     *  \param[in]  type   geometry type to consider
     *
     *  \returns number of leaf entities with a geometry of type \em type
     */
    size_t numBoundarySegments () const
    {
      return hostGrid().numBoundarySegments( );
    }
    /** \} */

    template< int codim >
    typename Codim< codim >::LevelIterator lbegin ( int level ) const
    {
      return lbegin< codim, All_Partition >( level );
    }

    template< int codim >
    typename Codim< codim >::LevelIterator lend ( int level ) const
    {
      return lend< codim, All_Partition >( level );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim >::template Partition< pitype >::LevelIterator
    lbegin ( int level ) const
    {
      typedef typename Traits::template Codim< codim >::template Partition< pitype >::LevelIteratorImpl Impl;
      return Impl( extraData(), hostGrid().template lbegin< codim, pitype >( level ) );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim >::template Partition< pitype >::LevelIterator
    lend ( int level ) const
    {
      typedef typename Traits::template Codim< codim >::template Partition< pitype >::LevelIteratorImpl Impl;
      return Impl( extraData(), hostGrid().template lend< codim, pitype >( level ) );
    }

    template< int codim >
    typename Codim< codim >::LeafIterator leafbegin () const
    {
      return leafbegin< codim, All_Partition >();
    }

    template< int codim >
    typename Codim< codim >::LeafIterator leafend () const
    {
      return leafend< codim, All_Partition >();
    }

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim >::template Partition< pitype >::LeafIterator
    leafbegin () const
    {
      typedef typename Traits::template Codim< codim >::template Partition< pitype >::LeafIteratorImpl Impl;
      return Impl( extraData(), hostGrid().template leafbegin< codim, pitype >() );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim >::template Partition< pitype >::LeafIterator
    leafend () const
    {
      typedef typename Traits::template Codim< codim >::template Partition< pitype >::LeafIteratorImpl Impl;
      return Impl( extraData(), hostGrid().template leafend< codim, pitype >() );
    }

    const GlobalIdSet &globalIdSet () const
    {
      if( !globalIdSet_ )
        globalIdSet_ = GlobalIdSet( hostGrid().globalIdSet() );
      assert( globalIdSet_ );
      return globalIdSet_;
    }

    const LocalIdSet &localIdSet () const
    {
      if( !localIdSet_ )
        localIdSet_ = LocalIdSet( hostGrid().localIdSet() );
      assert( localIdSet_ );
      return localIdSet_;
    }

    const LevelIndexSet &levelIndexSet ( int level ) const
    {
      assert( levelIndexSets_.size() == (size_t)(maxLevel()+1) );
      if( (level < 0) || (level > maxLevel()) )
      {
        DUNE_THROW( GridError, "LevelIndexSet for nonexisting level " << level
                               << " requested." );
      }

      LevelIndexSet *&levelIndexSet = levelIndexSets_[ level ];
      if( !levelIndexSet )
        levelIndexSet = new LevelIndexSet( hostGrid().levelIndexSet( level ) );
      assert( levelIndexSet );
      return *levelIndexSet;
    }

    const LeafIndexSet &leafIndexSet () const
    {
      if( !leafIndexSet_ )
        leafIndexSet_ = LeafIndexSet( hostGrid().leafIndexSet() );
      assert( leafIndexSet_ );
      return leafIndexSet_;
    }

    void globalRefine ( int refCount )
    {
      hostGrid().globalRefine( refCount );
      // update overall status
      update();
    }

    bool mark ( int refCount, const typename Codim< 0 >::Entity &entity )
    {
      return hostGrid().mark( refCount, getHostEntity< 0 >( entity ) );
    }

    int getMark ( const typename Codim< 0 >::Entity &entity ) const
    {
      return hostGrid().getMark( getHostEntity< 0 >( entity ) );
    }

    /** \brief  @copydoc Dune::Grid::preAdapt() */
    bool preAdapt ()
    {
      return false;
    }

    /** \brief  @copydoc Dune::Grid::adapt() */
    bool adapt ()
    {
      return false ;
    }

    /** \brief  @copydoc Dune::Grid::adapt()
        \param handle handler for restriction and prolongation operations
        which is a Model of the AdaptDataHandleInterface class.
    */
    template< class GridImp, class DataHandle >
    bool adapt ( AdaptDataHandleInterface< GridImp, DataHandle > &datahandle )
    {
      return false;
    }

    /** \brief  @copydoc Dune::Grid::postAdapt() */
    void postAdapt ()
    {
    }

    /** \name Parallel Data Distribution and Communication Methods
     *  \{ */

    /** \brief obtain size of overlap region for the leaf grid
     *
     *  \param[in]  codim  codimension for with the information is desired
     */
    int overlapSize ( int codim ) const
    {
      return 0;
    }

    /** \brief obtain size of ghost region for the leaf grid
     *
     *  \param[in]  codim  codimension for with the information is desired
     */
    int ghostSize( int codim ) const
    {
      return (codim == 0 ) ? 1 : 0;
    }

    /** \brief obtain size of overlap region for a grid level
     *
     *  \param[in]  level  grid level (0, ..., maxLevel())
     *  \param[in]  codim  codimension (0, ..., dimension)
     */
    int overlapSize ( int level, int codim ) const
    {
      return 0;
    }

    /** \brief obtain size of ghost region for a grid level
     *
     *  \param[in]  level  grid level (0, ..., maxLevel())
     *  \param[in]  codim  codimension (0, ..., dimension)
     */
    int ghostSize ( int level, int codim ) const
    {
      return ghostSize( codim );
    }

    /** \brief communicate information on a grid level
     *
     *  \param      dataHandle  communication data handle (user defined)
     *  \param[in]  interface   communication interface (one of
     *                          InteriorBorder_InteriorBorder_Interface,
     *                          InteriorBorder_All_Interface,
     *                          Overlap_OverlapFront_Interface,
     *                          Overlap_All_Interface,
     *                          All_All_Interface)
     *  \param[in]  direction   communication direction (one of
     *                          ForwardCommunication or BackwardCommunication)
     *  \param[in]  level       grid level to communicate
     */
    template< class DataHandle, class Data >
    void communicate ( CommDataHandleIF< DataHandle, Data > &dataHandle,
                       InterfaceType interface,
                       CommunicationDirection direction,
                       int level ) const
    {
       //levelGridView( level ).communicate( dataHandle, interface, direction );
    }

    /** \brief communicate information on leaf entities
     *
     *  \param      dataHandle  communication data handle (user defined)
     *  \param[in]  interface   communication interface (one of
     *                          InteriorBorder_InteriorBorder_Interface,
     *                          InteriorBorder_All_Interface,
     *                          Overlap_OverlapFront_Interface,
     *                          Overlap_All_Interface,
     *                          All_All_Interface)
     *  \param[in]  direction   communication direction (one of
     *                          ForwardCommunication, BackwardCommunication)
     */
    template< class DataHandle, class Data >
    void communicate ( CommDataHandleIF< DataHandle, Data > &dataHandle,
                       InterfaceType interface,
                       CommunicationDirection direction ) const
    {
      //leafGridView().communicate( dataHandle, interface, direction );
    }

    /** \brief obtain CollectiveCommunication object
     *
     *  The CollectiveCommunication object should be used to globally
     *  communicate information between all processes sharing this grid.
     *
     *  \note The CollectiveCommunication object returned is identical to the
     *        one returned by the host grid.
     */
    const CollectiveCommunication &comm () const
    {
      return comm_;
    }

    // data handle interface different between geo and interface

    /** \brief rebalance the load each process has to handle
     *
     *  A parallel grid is redistributed such that each process has about
     *  the same load (e.g., the same number of leaf entites).
     *
     *  \note DUNE does not specify, how the load is measured.
     *
     *  \returns \b true, if the grid has changed.
     */
    bool loadBalance ()
    {
      return false ;
    }

    /** \brief rebalance the load each process has to handle
     *
     *  A parallel grid is redistributed such that each process has about
     *  the same load (e.g., the same number of leaf entites).
     *
     *  The data handle is used to communicate the data associated with
     *  entities that move from one process to another.
     *
     *  \note DUNE does not specify, how the load is measured.
     *
     *  \param  datahandle  communication data handle (user defined)
     *
     *  \returns \b true, if the grid has changed.
     */

    template< class DataHandle, class Data >
    bool loadBalance ( CommDataHandleIF< DataHandle, Data > &datahandle )
    {
      return false;
    }

    /** \brief rebalance the load each process has to handle
     *
     *  A parallel grid is redistributed such that each process has about
     *  the same load (e.g., the same number of leaf entites).
     *
     *  The data handle is used to communicate the data associated with
     *  entities that move from one process to another.
     *
     *  \note DUNE does not specify, how the load is measured.
     *
     *  \param  dataHandle  data handle following the ALUGrid interface
     *
     *  \returns \b true, if the grid has changed.
     */
    template< class DofManager >
    bool loadBalance ( DofManager &dofManager )
    {
      return false;
    }

    /** \brief View for a grid level */
    template< PartitionIteratorType pitype >
    typename Partition< pitype >::LevelGridView levelGridView ( int level ) const
    {
      typedef typename Partition< pitype >::LevelGridView View;
      typedef typename View::GridViewImp ViewImp;
      return View( ViewImp( *this ) );
    }

    /** \brief View for the leaf grid */
    template< PartitionIteratorType pitype >
    typename Partition< pitype >::LeafGridView leafGridView () const
    {
      typedef typename Traits::template Partition< pitype >::LeafGridView View;
      typedef typename View::GridViewImp ViewImp;
      return View( ViewImp( *this ) );
    }

    /** \brief View for a grid level for All_Partition */
    LevelGridView levelGridView ( int level ) const
    {
      typedef typename LevelGridView::GridViewImp ViewImp;
      return LevelGridView( ViewImp( *this ) );
    }

    /** \brief View for the leaf grid for All_Partition */
    LeafGridView leafGridView () const
    {
      typedef typename LeafGridView::GridViewImp ViewImp;
      return LeafGridView( ViewImp( *this ) );
    }

    /** \brief obtain EntityPointer from EntitySeed. */
    template< class EntitySeed >
    typename Traits::template Codim< EntitySeed::codimension >::EntityPointer
    entityPointer ( const EntitySeed &seed ) const
    {
      typedef typename Traits::template Codim< EntitySeed::codimension >::EntityPointerImpl EntityPointerImpl;
      return EntityPointerImpl( extraData(), seed ) ;
    }

    /** \} */

    /** \name Miscellaneous Methods
     *  \{ */

    /** \brief update grid caches
     *
     *  This method has to be called whenever the underlying host grid changes.
     *
     *  \note If you adapt the host grid through this geometry grid's
     *        adaptation or load balancing methods, update is automatically
     *        called.
     */
    void update ()
    {
    }

    /** \} */

  protected:
    UnstructuredGridType* createGrid( Opm::DeckConstPtr deck, const std::vector< double >& poreVolumes ) const
    {
        auto eclipseGrid = std::make_shared<const Opm::EclipseGrid>(deck);
        struct Opm::grdecl g;
        std::vector<int> actnum;
        std::vector<double> coord;
        std::vector<double> zcorn;
        std::vector<double> mapaxes;

        g.dims[0] = eclipseGrid->getNX();
        g.dims[1] = eclipseGrid->getNY();
        g.dims[2] = eclipseGrid->getNZ();

        eclipseGrid->exportMAPAXES( mapaxes );
        eclipseGrid->exportCOORD( coord );
        eclipseGrid->exportZCORN( zcorn );
        eclipseGrid->exportACTNUM( actnum );

        g.coord = coord.data();
        g.zcorn = zcorn.data();
        g.actnum = actnum.data();
        g.mapaxes = mapaxes.data();

        if (!poreVolumes.empty() && (eclipseGrid->getMinpvMode() != MinpvMode::ModeEnum::Inactive)) {
            MinpvProcessor mp(g.dims[0], g.dims[1], g.dims[2]);
            const double minpv_value  = eclipseGrid->getMinpvValue();
            mp.process(poreVolumes, minpv_value, actnum, zcorn.data());
        }

        const double z_tolerance = eclipseGrid->isPinchActive() ?
            eclipseGrid->getPinchThresholdThickness() : 0.0;
        UnstructuredGridType* cgrid = create_grid_cornerpoint(&g, z_tolerance);
        if (!cgrid) {
            OPM_THROW(std::runtime_error, "Failed to construct grid.");
        }
      return cgrid;
    }


  public:
    using Base::getRealImplementation;

    template< int codim >
    static const typename HostGrid::template Codim< codim >::Entity &
    getHostEntity( const typename Codim< codim >::Entity &entity )
    {
      return getRealImplementation( entity ).hostEntity();
    }

    typedef typename Traits :: ExtraDataType ExtraData;
    ExtraData extraData () const  { return ExtraData(); }

    template <class EntitySeed>
    int faces( const EntitySeed&

  protected:
    std::unique_ptr< UnstructuredGridType > grid_;
    mutable std::vector< LevelIndexSet * > levelIndexSets_;
    mutable LeafIndexSet leafIndexSet_;
    mutable GlobalIdSet globalIdSet_;
    mutable LocalIdSet localIdSet_;
  };



  // PolyhedralGrid::Codim
  // -------------

  template< class HostGrid >
  template< int codim >
  struct PolyhedralGrid< HostGrid >::Codim
  : public Base::template Codim< codim >
  {
    /** \name Entity and Entity Pointer Types
     *  \{ */

    /** \brief type of entity
     *
     *  The entity is a model of Dune::Entity.
     */
    typedef typename Traits::template Codim< codim >::Entity Entity;

    /** \brief type of entity pointer
     *
     *  The entity pointer is a model of Dune::EntityPointer.
     */
    typedef typename Traits::template Codim< codim >::EntityPointer EntityPointer;

    /** \} */

    /** \name Geometry Types
     *  \{ */

    /** \brief type of world geometry
     *
     *  Models the geomtry mapping of the entity, i.e., the mapping from the
     *  reference element into world coordinates.
     *
     *  The geometry is a model of Dune::Geometry, implemented through the
     *  generic geometries provided by dune-grid.
     */
    typedef typename Traits::template Codim< codim >::Geometry Geometry;

    /** \brief type of local geometry
     *
     *  Models the geomtry mapping into the reference element of dimension
     *  \em dimension.
     *
     *  The local geometry is a model of Dune::Geometry, implemented through
     *  the generic geometries provided by dune-grid.
     */
    typedef typename Traits::template Codim< codim >::LocalGeometry LocalGeometry;

    /** \} */

    /** \name Iterator Types
     *  \{ */

    template< PartitionIteratorType pitype >
    struct Partition
    {
      typedef typename Traits::template Codim< codim >
        ::template Partition< pitype >::LeafIterator
        LeafIterator;
      typedef typename Traits::template Codim< codim >
        ::template Partition< pitype >::LevelIterator
        LevelIterator;
    };

    /** \brief type of level iterator
     *
     *  This iterator enumerates the entites of codimension \em codim of a
     *  grid level.
     *
     *  The level iterator is a model of Dune::LevelIterator.
     */
    typedef typename Partition< All_Partition >::LeafIterator LeafIterator;

    /** \brief type of leaf iterator
     *
     *  This iterator enumerates the entites of codimension \em codim of the
     *  leaf grid.
     *
     *  The leaf iterator is a model of Dune::LeafIterator.
     */
    typedef typename Partition< All_Partition >::LevelIterator LevelIterator;

    /** \} */
  };

} // namespace Dune

#include <dune/grid/polyhedralgrid/persistentcontainer.hh>
#include <dune/grid/polyhedralgrid/twistutility.hh>

#endif // #ifndef DUNE_POLYHEDRALGRID_GRID_HH
