// $Id: gridcheck.cc 10840 2009-06-26 07:14:23Z atgeirr $
#ifndef GRIDCHECK_CC
#define GRIDCHECK_CC

/**

  Implements a generic grid check

  
  \todo check return types
  
*/

#include <dune/grid/common/capabilities.hh>
#include <dune/common/helpertemplates.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>
#include <dune/grid/common/referenceelements.hh>

#include <dune/grid/genericgeometry/conversion.hh>

#include "checkindexset.cc"

#include <limits>

// machine epsilon is multiplied by this factor 
static double factorEpsilon = 1.e8;

class CheckError : public Dune::Exception {};

#if 0
// --- compile-time check of element-interface

template <class Geometry, bool doCheck>
struct JacobianInverse
{
  static void check(const Geometry &e) 
    {
      typedef typename Geometry::ctype ctype;
      Dune::FieldVector<ctype, Geometry::mydimension> v;
      e.jacobianInverseTransposed(v);
    }
  JacobianInverse() 
    {
      c = check;
    };
  void (*c)(const Geometry&);
};

template <class Geometry>
struct JacobianInverse<Geometry, false>
{
  static void check(const Geometry &e) 
    {
    }
  JacobianInverse() 
    {
      c = check;
    };
  void (*c)(const Geometry&);
};
#endif

template <class Geometry, int codim, int dim>
struct GeometryInterface 
{
  static void check ( const Geometry &geo )
  {
    IsTrue<dim-codim == Geometry::mydimension>::yes();
    IsTrue<dim == Geometry::dimension>::yes();
    
    typedef typename Geometry::ctype ctype;
    
    geo.type();
    geo.corners();
    geo.corner( 0 );
    
    Dune::FieldVector<ctype, Geometry::mydimension> v;
    geo.global(v);
    Dune::FieldVector<ctype, Geometry::coorddimension> g;
    geo.local(g);
    geo.integrationElement(v);
    geo.jacobianTransposed( v );
    geo.jacobianInverseTransposed( v );
#if 0
    JacobianInverse<Geometry,
      (int)Geometry::coorddimension == (int)Geometry::mydimension>();
#endif
  }

  GeometryInterface() 
    {
      c = check;
    };
  void (*c)(const Geometry&);
};

// reduced test on vertices
template <class Geometry, int dim>
struct GeometryInterface <Geometry, dim, dim>
{
  static void check(const Geometry &e) 
    {
      IsTrue<0 == Geometry::mydimension>::yes();
      IsTrue<dim == Geometry::dimension>::yes();
      
      // vertices have only a subset of functionality
      e.type();
      e.corners();
      e[0];
    }
  GeometryInterface() 
    {
      c = check;
    };
  void (*c)(const Geometry&);
};

// --- compile-time check of entity-interface

// tests that should work on entities of all codimensions
template <class Entity>
void DoEntityInterfaceCheck (Entity &e) 
{
  // exported types
  typedef typename Entity::ctype ctype;
  
  // methods on each entity
  e.level();
  e.partitionType();
  e.geometry();
  
  // check interface of attached element-interface
  GeometryInterface<typename Entity::Geometry, Entity::codimension, Entity::dimension>();
}

// recursive check of codim-0-entity methods count(), entity()
template <class Grid, int cd, bool hasEntity>
struct ZeroEntityMethodCheck 
{
  typedef typename Grid::template Codim<0>::Entity Entity;
  static void check(Entity &e)
    {
      // check types
      typedef typename Entity::IntersectionIterator IntersectionIterator;
      typedef typename Entity::HierarchicIterator HierarchicIterator;
      typedef typename Entity::EntityPointer EntityPointer;

      e.template count<cd>();
      e.template entity<cd>(0);

      // recursively check on
      ZeroEntityMethodCheck<Grid, cd - 1,
        Dune::Capabilities::hasEntity<Grid, cd - 1>::v >();
    }
  ZeroEntityMethodCheck () 
    {
      c = check;
    }
  void (*c)(Entity &e);
};

// just the recursion if the grid does not know about this codim-entity
template<class Grid, int cd>
struct ZeroEntityMethodCheck<Grid, cd, false>
{
  typedef typename Grid::template Codim<0>::Entity Entity;
  static void check(Entity &e)
    {
      // check types
      typedef typename Entity::IntersectionIterator IntersectionIterator;
      typedef typename Entity::HierarchicIterator HierarchicIterator;
      typedef typename Entity::EntityPointer EntityPointer;
      
      // recursively check on
      ZeroEntityMethodCheck<Grid, cd - 1,
         Dune::Capabilities::hasEntity<Grid, cd - 1>::v >();
    }
  ZeroEntityMethodCheck () 
    {
      c = check;
    }
  void (*c)(Entity &e);
};

// end recursive checking
template <class Grid>
struct ZeroEntityMethodCheck<Grid, 0, true>
{
  typedef typename Grid::template Codim<0>::Entity Entity;
  static void check(Entity &e)
    {
      // check types
      typedef typename Entity::IntersectionIterator IntersectionIterator;
      typedef typename Entity::HierarchicIterator HierarchicIterator;
      typedef typename Entity::EntityPointer EntityPointer;
      
      e.template count<0>();
      e.template entity<0>(0);

    }
  ZeroEntityMethodCheck () 
    {
      c = check;
    }
  void (*c)(Entity &e);
};

// end recursive checking - same as true
// ... codim 0 is always needed
template <class Grid>
struct ZeroEntityMethodCheck<Grid, 0, false>
{
  typedef typename Grid::template Codim<0>::Entity Entity;
  static void check(Entity &e)
    {
      // check types
      typedef typename Entity::IntersectionIterator IntersectionIterator;
      typedef typename Entity::HierarchicIterator HierarchicIterator;
      typedef typename Entity::EntityPointer EntityPointer;
      
      e.template count<0>();
      e.template entity<0>(0);
    }
  ZeroEntityMethodCheck () 
    {
      c = check;
    }
  void (*c)(Entity &e);
};

// IntersectionIterator interface check
template <class Grid>
struct IntersectionIteratorInterface
{
  typedef typename Grid::template Codim<0>::IntersectionIterator IntersectionIterator;
  enum { dim = Grid::dimension };
  typedef typename Grid::ctype ct;
  
  static void check (IntersectionIterator &i)
    {
      // increment / equality / ...
      IntersectionIterator j = i;
      j++;
      i == j;
      i != j;
      j = i;

      // state
      i.boundary();
      i.neighbor();

      // id of boundary segment
      i.boundaryId();

      // neighbouring elements
      i.inside();
      if(i.neighbor()) i.outside();

      // geometry
      i.intersectionSelfLocal();
      if(i.neighbor()) i.intersectionNeighborLocal();
      i.intersectionGlobal();

      i.numberInSelf();
      if(i.neighbor()) i.numberInNeighbor();

      Dune::FieldVector<ct, dim-1> v(0);
      i.outerNormal(v);
      i.integrationOuterNormal(v);
      i.unitOuterNormal(v);
    }
  IntersectionIteratorInterface ()
    {
      c = check;  
    }
  void (*c)(IntersectionIterator&);    
};

// check codim-entity and pass on to codim + 1 
template <class Grid, int codim, int dim, bool hasEntity>
struct EntityInterface
{
  typedef typename Grid::template Codim<codim>::Entity Entity;
  
  static void check (Entity &e)
    {
      // consistent?
      IsTrue<codim == Entity::codimension>::yes();
      IsTrue<dim == Entity::dimension>::yes();      

      // do the checking
      DoEntityInterfaceCheck(e);
      
      // recursively check sub-entities
      EntityInterface<Grid, codim + 1, dim,
        Dune::Capabilities::hasEntity<Grid, codim + 1>::v >();
    }
  EntityInterface ()
    {
      c = check;  
    }
  void (*c)(Entity&);    
};

// just the recursion if the grid does not know about this codim-entity
template <class Grid, int codim, int dim>
struct EntityInterface<Grid, codim, dim, false>
{
  typedef typename Grid::template Codim<codim>::Entity Entity;
  
  static void check (Entity &e)
    {
      // recursively check sub-entities
      EntityInterface<Grid, codim + 1, dim,
        Dune::Capabilities::hasEntity<Grid, codim + 1>::v >();
    }
  EntityInterface ()
    {
      c = check;
    }
  void (*c)(Entity&);    
};

// codim-0 entities have different interface
template <class Grid, int dim>
struct EntityInterface<Grid, 0, dim, true>
{
  typedef typename Grid::template Codim<0>::Entity Entity;
  
  static void check (Entity &e,bool checkLevelIter=true) 
    {
      // consistent?
      IsTrue<0 == Entity::codimension>::yes();
      IsTrue<dim == Entity::dimension>::yes();      

      // do the common checking
      DoEntityInterfaceCheck(e);

      // special codim-0-entity methods which are parametrized by a codimension
      ZeroEntityMethodCheck
        <Grid, dim, Dune::Capabilities::hasEntity<Grid, dim>::v >();

      // grid hierarchy
      e.father();
      e.geometryInFather();


      // intersection iterator
      if (checkLevelIter) {
        e.ilevelbegin();
        e.ilevelend();
        IntersectionIteratorInterface<Grid>(e.ilevelbegin());
      }
      e.ileafbegin();
      e.ileafend();

      if(e.isLeaf())
        IntersectionIteratorInterface<Grid>(e.ileafbegin());
      
      // hierarchic iterator
      e.hbegin(0);
      e.hend(0);

      // adaption
      e.state();

      // recursively check sub-entities
      EntityInterface<Grid, 1, dim,
        Dune::Capabilities::hasEntity<Grid, 1>::v >();      
    }
  EntityInterface ()
    {
      c = check;  
    }
  void (*c)(Entity&);    
};

// non existinng codim-0 entity
template <class Grid, int dim>
struct EntityInterface<Grid, 0, dim, false>
{
  typedef typename Grid::template Codim<0>::Entity Entity;

  static void check (Entity &e)
    {
      // recursively check sub-entities
      EntityInterface<Grid, 1, dim,
        Dune::Capabilities::hasEntity<Grid, 1>::v >();      
    }
  EntityInterface ()
    {
      c = check;
    }
  void (*c)(Entity&);    
};

// end the recursion over entity-codimensions
template <class Grid, int dim>
struct EntityInterface<Grid, dim, dim, true>
{
  typedef typename Grid::template Codim<dim>::Entity Entity;
  
  // end recursion
  static void check (Entity &e)
    {
      // consistent?
      IsTrue<dim == Entity::codimension>::yes();
      IsTrue<dim == Entity::dimension>::yes();
      
      // run common test
      DoEntityInterfaceCheck(e);      
    }
  
  EntityInterface()
    {
      c = check;  
    }
  void (*c)(Entity&);    
};

// end the recursion over entity-codimensions
// ... codim dim entity does not exist
template <class Grid, int dim>
struct EntityInterface<Grid, dim, dim, false>
{
  typedef typename Grid::template Codim<dim>::Entity Entity;
  
  // end recursion
  static void check (Entity &e)
    {
    }
  
  EntityInterface()
    {
      c = check;  
    }
  void (*c)(Entity&);    
};

template<class Grid>
struct LeafInterface
{
  static void check(Grid &g)
    {
      g.template leafbegin<0>();
      g.template leafend<0>();
    }
  LeafInterface()
    {
      c = check;  
    }
  void (*c)(Grid&);  
};

template <class Grid>
struct GridInterface
{
  static void check (const Grid &g)
    {
      // check for exported types
      typedef typename Grid::ctype ctype;
      typedef typename Grid::template Codim<0>::LevelIterator LevelIterator;
      typedef typename Grid::template Codim<0>::EntityPointer EntityPointer;
      typedef typename Grid::template Codim<0>::LeafIterator LeafIterator;
      
      // check for member functions
      g.maxLevel();
      // number of grid entities of a given codim on a given level
      g.size(0,0);
      // number of leaf entities per codim in this process
      g.size(0);
      // number of entities per level and geometry type in this process
      g.size(0, Dune::GeometryType(Dune::GeometryType::cube,Grid::dimension));
      // number of leaf entities per geometry type in this process
      g.size(Dune::GeometryType(Dune::GeometryType::cube,Grid::dimension));

      // Check overlap and ghost size on level 0
      g.overlapSize(0,0);
      g.ghostSize(0,0);

      // Check overlap and ghost size on the leaf level
      g.overlapSize(0);
      g.ghostSize(0);

      // check for iterator functions
      g.template lbegin<0>(0);
      g.template lend<0>(0);

      LeafInterface< Grid>();
      
      // Check for index sets
      typedef typename Grid::template Codim<0>::LevelIndexSet LevelIndexSet;
      typedef typename Grid::template Codim<0>::LeafIndexSet LeafIndexSet;
      typedef typename Grid::template Codim<0>::LocalIdSet LocalIdSet;
      typedef typename Grid::template Codim<0>::GlobalIdSet GlobalIdSet;
      
      g.levelIndexSet(0);
      if (g.template lbegin<0>(0) !=g.template lend<0>(0) ) {
	// Instantiate all methods of LevelIndexSet
	g.levelIndexSet(0).index(*g.template lbegin<0>(0));
	/** \todo Test for subindex is missing, because I don't know yet
	    how to test for the existence of certain codims */
      }
      g.levelIndexSet(0).
	size(Dune::GeometryType(Dune::GeometryType::simplex,Grid::dimension));
      for (int codim = 0; codim < Grid::dimension; codim++)
        g.levelIndexSet(0).geomTypes(codim);

      if (g.template lbegin<0>(0) !=g.template lend<0>(0) ) {
	// Instantiate all methods of LeafIndexSet
	g.leafIndexSet().index(*g.template leafbegin<0>());
      }
      /** \todo Test for subindex is missing, because I don't know yet
       how to test for the existence of certain codims */
      g.leafIndexSet().size(Dune::GeometryType(Dune::GeometryType::simplex,Grid::dimension));
      for (int codim = 0; codim < Grid::dimension; codim++)
        g.leafIndexSet().geomTypes(codim);

      if (g.template lbegin<0>(0) !=g.template lend<0>(0) ) {
	// Instantiate all methods of LocalIdSet
	/** \todo Test for subindex is missing, because I don't know yet
	    how to test for the existence of certain codims */
	g.localIdSet().id(*g.template lbegin<0>(0));
	// Instantiate all methods of GlobalIdSet
	/** \todo Test for subindex is missing, because I don't know yet
	    how to test for the existence of certain codims */
	g.globalIdSet().id(*g.template lbegin<0>(0));
      }
      // recursively check entity-interface
      // ... we only allow grids with codim 0 zero entites
      dune_static_assert((Dune::Capabilities::hasEntity<Grid, 0>::v),"Grid must have codim 0 entities");
      dune_static_assert((Dune::Capabilities::hasEntity<const Grid, 0>::v),"Grid must have codim 0 entities");
      
      //EntityInterface<Grid, 0, Grid::dimension,
      //  Dune::Capabilities::hasEntity<Grid, 0>::v >();

      // !!! check for parallel grid?
      /*
      g.template lbegin<0, Dune::Ghost_Partition>(0);
      g.template lend<0, Dune::Ghost_Partition>(0);
      */
    }
  GridInterface()
    {
      c = check;
    }
  // member just to avoid "unused variable"-warning in constructor
  void (*c)(const Grid&);
};

// check
// Entity::geometry()[c] == Entity::entity<dim>.geometry()[0]
// for codim=cd
template <int cd, class Grid, class Entity, bool doCheck>
struct subIndexCheck
{
  subIndexCheck ( const Grid &g, const Entity &e )
  {
    typedef typename Grid::template Codim< cd >::EntityPointer EntityPointer;
    const int imax = e.template count<cd>();
    for( int i = 0; i < imax; ++i )
    {
      // check construction of entity pointers   
      EntityPointer ep( *(e.template subEntity< cd >( i ) ) );
      assert( ep == e.template subEntity< cd >( i ) );

      // test compactify 
      ep.compactify();

      const typename Grid::LevelIndexSet &levelIndexSet = g.levelIndexSet( e.level() );

      if( !levelIndexSet.contains( *ep ) )
      {
        std::cerr << "Error: Level index set does not contain all subentities." << std::endl;
        assert( false );
      }
      if( levelIndexSet.index( *ep ) != levelIndexSet.subIndex( e, i, cd ) )
      {
        int id_e = levelIndexSet.index( e );
        int id_e_i = levelIndexSet.index( *ep );
        int subid_e_i = levelIndexSet.subIndex( e, i, cd );
        std::cerr << "Error: levelIndexSet.index( *(e.template subEntity< cd >( i ) ) ) "
                  << "!= levelIndexSet.subIndex( e, i, cd )  "
                  << "[with cd=" << cd << ", i=" << i << "]" << std::endl;
        std::cerr << "       ... index( e ) = " << id_e << std::endl;
        std::cerr << "       ... index( e.subEntity< cd >( i ) ) = " << id_e_i << std::endl;
        std::cerr << "       ... subIndex( e, i, cd ) = " << subid_e_i << std::endl;
        assert( false );
      }

#if 0
      typedef Dune::GenericGeometry::MapNumberingProvider< Entity::dimension > Numbering;
      const unsigned int tid = Dune::GenericGeometry::topologyId( e.type() );
      const int oldi = Numbering::template generic2dune< cd >( tid, i );

      if( levelIndexSet.subIndex( e, i, cd ) != levelIndexSet.template subIndex< cd >( e, oldi ) )
      {
        int id_e = levelIndexSet.index( e );
        int subid_e_i = levelIndexSet.template subIndex< cd >( e, oldi );
        int subid_e_i_cd = levelIndexSet.subIndex( e, i, cd );

        std::cerr << "Error: levelIndexSet.subIndex( e, i, cd ) "
                  << "!= levelIndexSet.subIndex< cd >( e, generic2dune( i ) )  "
                  << "[with cd=" << cd << ", i=" << i << "]" << std::endl;
        std::cerr << "       ... index(e)=" << id_e << std::endl;
        std::cerr << "       ... subIndex<cd>(e,i)=" << subid_e_i << std::endl;
        std::cerr << "       ... subIndex(e,dune2generic(i),cd)=" << subid_e_i_cd << std::endl;
        assert( false );
      }
#endif
    }

    subIndexCheck< cd-1, Grid, Entity, Dune::Capabilities::hasEntity< Grid, cd-1 >::v > sick( g, e );
  }
};


// end recursion of subIndexCheck
template <class Grid, class Entity, bool doCheck>
struct subIndexCheck<-1, Grid, Entity, doCheck>
{
  subIndexCheck (const Grid & g, const Entity & e)
    {
      return;
    }
};
// do nothing if doCheck==false
template <int cd, class Grid, class Entity>
struct subIndexCheck<cd, Grid, Entity, false>
{
  subIndexCheck (const Grid & g, const Entity & e)
    {
      subIndexCheck<cd-1,Grid,Entity,
        Dune::Capabilities::hasEntity<Grid,cd-1>::v> sick(g,e);
    }
};
template <class Grid, class Entity>
struct subIndexCheck<-1, Grid, Entity, false>
{
  subIndexCheck (const Grid & g, const Entity & e)
    {
      return;
    }
};

// name says all
template <class Grid>
void zeroEntityConsistency (Grid &g)
{
  typedef typename Grid::ctype ctype;
  const int dimGrid = Grid::dimension;
  const int dimWorld = Grid::dimensionworld;

  typedef typename Grid::template Codim<0>::LevelIterator LevelIterator;
  typedef typename Grid::template Codim<0>::Geometry Geometry;
  typedef typename Grid::template Codim<0>::Entity Entity;
  typedef typename Grid::template Codim<0>::EntityPointer EntityPointer;

  const typename Grid::LevelIndexSet &levelIndexSet = g.levelIndexSet( g.maxLevel() );
  const typename Grid::LeafIndexSet &leafIndexSet = g.leafIndexSet();
  const typename Grid::LocalIdSet &localIdSet = g.localIdSet();
  const typename Grid::GlobalIdSet &globalIdSet = g.globalIdSet();

  LevelIterator it = g.template lbegin<0>(g.maxLevel());
  const LevelIterator endit = g.template lend<0>(g.maxLevel());
  
  for (; it!=endit; ++it)
  {
    // check construction from entity
    EntityPointer ep( *it ) ;
    assert( ep == it );

    // check compactify of entity pointer 
    ep.compactify();

    // Entity::subEntity<0>(0) == Entity
    EntityPointer subEntity = it->template subEntity< 0 >( 0 );
    if( levelIndexSet.index( *subEntity ) != levelIndexSet.index( *it ) )
    {
      std::cerr << "Error: Level index for subEntity< 0 >( 0 ) differs." << std::endl;
      assert( false );
    }
    if( leafIndexSet.index( *subEntity ) != leafIndexSet.index( *it ) )
    {
      std::cerr << "Error: Leaf index for subEntity< 0 >( 0 ) differs." << std::endl;
      assert( false );
    }
    if( globalIdSet.id( *subEntity ) != globalIdSet.id( *it ) )
    {
      std::cerr << "Error: Global id for subEntity< 0 >( 0 ) differs." << std::endl;
      assert( false );
    }
    if( localIdSet.id( *subEntity ) != localIdSet.id( *it ) )
    {
      std::cerr << "Error: Local id for subEntity< 0 >( 0 ) differs." << std::endl;
      assert( false );
    }

    if( subEntity->level() != it->level() )
    {
      std::cerr << "Error: Level for subEntity< 0 >( 0 ) differs." << std::endl;
      assert( false );
    }

    // Entity::count<dim>() == Entity::geometry().corners();
    // Entity::geometry()[c] == Entity::entity<dim>.geometry()[0];
    const int numCorners = it->template count< dimGrid >();
    if( numCorners != it->geometry().corners() )
    {
      std::cerr << "Error: Entity::count< dimGrid >() != Entity::geometry().corners()." << std::endl;
      assert( false );
    }
    for( int c = 0; c < numCorners; ++c )
    {
      Dune::FieldVector< ctype, dimWorld > c1( it->geometry().corner( c ) );
      Dune::FieldVector< ctype, dimWorld > c2( it->template subEntity< dimGrid >( c )->geometry().corner( 0 ) );
      if( (c2-c1).two_norm() > 10*std::numeric_limits< ctype >::epsilon() )
      {
        std::cerr << "Error: | geometry.corner( c ) - subEntity< dimGrid >( c ) | = | "
                  << c1 << " - " << c2 << " | = " << (c2-c1).two_norm()
                  << " != 0  [with c = " << c << " ]" << std::endl;
        assert( false );
      }
    }

    // check that corner coordinates are distinct 
    const int corners= it->geometry().corners();
    for (int c=0; c<corners; ++c)
    {
      typedef Dune::FieldVector<typename Grid::ctype, Grid::dimensionworld> CoordinateType;
      for(int d=c+1; d<corners; ++d) 
      {
        const CoordinateType c1 = it->geometry().corner( c );
        const CoordinateType c2 = it->geometry().corner( d );
        if( (c2-c1).two_norm() < 10 * std::numeric_limits<typename Grid::ctype>::epsilon() )
        {
          DUNE_THROW(CheckError, "geometry["<<c<<"] == geometry["<<d<<"]");
        }
      }
    }
    subIndexCheck<Grid::dimension, Grid, Entity, true> sick(g,*it);
  }
}

/*
 * search the LevelIterator for each IntersectionIterator
 */
template <class Grid>
void assertNeighbor (Grid &g)
{
  typedef typename Grid::template Codim<0>::LevelIterator LevelIterator;
  typedef typename Grid::template Codim<0>::LevelIntersectionIterator IntersectionIterator;
  enum { dim = Grid::dimension };
  typedef typename Grid::ctype ct;

  typedef typename Grid::template Codim<0>::GlobalIdSet GlobalIdSet;
  const GlobalIdSet & globalid = g.globalIdSet();

  typedef typename Grid::template Codim< 0 >::Entity Entity;
  typedef typename Grid::template Codim< 0 >::EntityPointer EntityPointer;
  
  LevelIterator e = g.template lbegin<0>(0);
  const LevelIterator eend = g.template lend<0>(0);

  LevelIterator next = e; 
  if (next != eend)
  {
    ++next;

    // do nothing for AlbertaGrid
    if( g.name() == "AlbertaGrid" )
    {
      std::cerr << "WARNING: skip indices test using LevelIntersectionIterator "
                << "for AlbertaGrid!" << std :: endl;
      return;
    }

    // iterate over level
    for (;e != eend; ++e)
    {
      const Entity &entity = *e;

      const bool isGhost = (entity.partitionType() == Dune::GhostEntity);

      // call global id 
      globalid.id( entity );

      if( !entity.isLeaf() )
      {
        if( entity.ileafbegin() != entity.ileafend() )
        {
          DUNE_THROW(CheckError, "On non-leaf entities ileafbegin should be equal to ileafend!");
        }
      }

      const int numFaces = entity.template count< 1 >();
      bool hasFaces = (numFaces != 0);
      // flag vector for elements faces
      std::vector< bool > visited( numFaces, false );

      // loop over intersections
      IntersectionIterator endit = entity.ilevelend();
      IntersectionIterator it = entity.ilevelbegin();

      if( it == endit )
      {
        // ALU has no intersection iterators on ghosts, so make this case non-fatal
        if( !isGhost )
        {
          std::cerr << "Error: Entity without intersections encountered." << std::endl;
          assert( false );
        }
        else
          std::cerr << "Warning: Ghost Entity without intersections encountered." << std :: endl;
      }

      // for all intersections
      for(; it != endit; ++it)
      {
        // numbering
	if (hasFaces) {
	  const int numberInSelf = it->indexInInside();
	  if( (numberInSelf < 0) || (numberInSelf >= numFaces) )
	    {
	      std :: cout << "Error: Invalid numberInSelf: " << numberInSelf
			  << " (should be between 0 and " << (numFaces-1) << ")"
			  << std :: endl;
	      assert( false );
	    }
	  else
	    {
	      // mark face visited
		visited[ numberInSelf ] = true;
	    }
	}

        // id of boundary segment
        if( it->boundary() )
          it->boundaryId();

        // check id
        assert( globalid.id(*(it->inside())) == globalid.id( entity ) );


        // geometry
        it->geometryInInside();
        it->geometry();
        
        // normal vectors
        Dune::FieldVector<ct, dim-1> v(0);
        it->outerNormal(v);
        it->integrationOuterNormal(v);
        it->unitOuterNormal(v);

        if( isGhost && !it->neighbor() )
        {
          std::cerr << "Error: On ghosts, all intersections must have neighbors." << std::endl;
          assert( false );
        }
       
        if( it->neighbor() )
        {
          const EntityPointer outsidePtr = it->outside();
          const Entity &outside = *outsidePtr;

          if( isGhost && (outside.partitionType() == Dune::GhostEntity) )
          {
            std::cerr << "Error: Intersections between ghosts shall not be returned." << std::endl;
            assert( false );
          }
          
          // geometry
          it->geometryInOutside();

          // numbering
          const int numFaces = outside.template count< 1 >();
	  bool hasFaces = (numFaces != 0);
	  if (hasFaces) {
            const int numberInNeighbor = it->indexInOutside();
	    if( (numberInNeighbor < 0) || (numberInNeighbor >= numFaces) )
	      {
		std :: cout << "Error: Invalid numberInNeighbor: " << numberInNeighbor
			    << " (should be between 0 and " << (numFaces-1) << ")"
			    << std :: endl;
		assert( false );
	      }
	  }

          // search neighbouring cell
          if( outsidePtr->partitionType() == Dune::InteriorEntity )
          {
            assert( globalid.id( outside ) != globalid.id( entity ) );
            const Dune::PartitionIteratorType pitype = Dune::InteriorBorder_Partition;

            typedef typename Grid::template Codim< 0 >
              ::template Partition< pitype >::LevelIterator
              LevelIterator;

            const int level = entity.level();
            bool foundNeighbor = false;
            LevelIterator nit = g.template lbegin< 0, pitype >( level );
            const LevelIterator nend = g.template lend< 0, pitype > ( level );
            for( ; nit != nend; ++nit ) 
            {
              if( nit->partitionType() != Dune::InteriorEntity )
              {
                std::cerr << "Error: LevelIterator for InteriorBorder_Partition "
                          << "stops on non-interior entity." << std :: endl;
                assert( false );
              }

              if( nit != outsidePtr )
                assert( globalid.id( outside ) != globalid.id( *nit ) );
              else
                foundNeighbor = true;
            }
            if( !foundNeighbor )
            {
              std :: cerr << "Error: Interior neighbor returned by "
                          << "LevelIntersectionIterator not found on that level."
                          << std :: endl;
              assert( false );
            }
          }
        }
      }

      // check that all faces were visited
      // note: This check is wrong on ghosts, where only intersections with
      //       the domain are allowed.
      if( !isGhost )
      {
        for(size_t i=0; i<visited.size(); i++) assert(visited[i] == true);
      }
    }
  }
}

template <class GridType, bool c> 
struct CheckMark
{
  template <class IteratorType>
  static void check(GridType & grid, IteratorType & it)
  {
      // last marker is 0, so the grid is not changed after this check
      const int refCount[4] = {1,0,-1,0};
      for(int k=0; k<4; ++k) 
      {
        // mark entity 
        bool marked = grid.mark( refCount[k] , *it );
        // if element was marked, check that the marker was set correctly 
        if(marked)
        {
          // now getMark should return the mark we just set, otherwise error
          if( grid.getMark( *it ) != refCount[k] )
            DUNE_THROW(CheckError,"mark/getMark method not working correctly!");
        }
      }
  }
};

template <class GridType> 
struct CheckMark<GridType,false>
{
  template <class IteratorType>
  static void check(const GridType & grid, IteratorType & it )
  {
  }
};

/*
 * Iterate over the grid und do some runtime checks
 */

template <bool checkMark , class Grid> 
void iterate(Grid &g)
{
  typedef typename Grid::template Codim<0>::LevelIterator LevelIterator;
  typedef typename Grid::template Codim<0>::EntityPointer EntityPointer;
  typedef typename Grid::template Codim<0>::HierarchicIterator HierarchicIterator;
  typedef typename Grid::template Codim<0>::Geometry Geometry;
  int l = g.maxLevel();
  LevelIterator it = g.template lbegin<0>(l);
  const LevelIterator endit = g.template lend<0>(l);

  Dune::FieldVector<typename Grid::ctype, Grid::dimension> origin(1);
  Dune::FieldVector<typename Grid::ctype, Grid::dimension> result;

  for (;it != endit; ++it)
  {
    LevelIterator l1 = it;
    LevelIterator l2 = l1; ++l1;
    assert( (l2 == it) && (it == l2) );
    assert( (l1 != it) && (it != l1) );
    ++l2;
    assert( (l1 == l2) && (l2 == l1) );

    const Geometry &geo = it->geometry();
    
    result = geo.local( geo.global( origin ) );
    typename Grid::ctype error = (result-origin).two_norm();
    if(!it->type().isSingular() && error >= factorEpsilon * std::numeric_limits<typename Grid::ctype>::epsilon())
      {
        DUNE_THROW(CheckError, "|| geom.local(geom.global(" << origin
                   << ")) - origin || != 0 ( || " << result << " - origin || ) = " << error);
      };
    geo.integrationElement( origin );

    if(!it->type().isSingular() && (int)Geometry::coorddimension == (int)Geometry::mydimension)
      geo.jacobianInverseTransposed( origin );

    if( geo.type() != it->type() )
    {
      std::cerr << "inconsistent geometry type: entity.type() = " << it->type()
                << ", geometry.type() = " << geo.type() << "." << std::endl;
      assert( false );
    }
    if(geo.corners() > 0) {
	geo.corner( 0 );
    }
  }
  
  typedef typename Grid::template Codim<0>::LeafIterator LeafIterator;
  LeafIterator lit = g.template leafbegin<0>();
  const LeafIterator lend = g.template leafend<0>();

  // if empty grid, do nothing 
  if(lit == lend) return;

  for (;lit != lend; ++lit)
  {
    LeafIterator l1 = lit;
    LeafIterator l2 = l1; ++l1;
    assert(l2 == lit);
    assert(l1 != lit);
    ++l2;
    assert(l1 == l2);

    // leaf check 
    if( !lit->isLeaf() ) 
      DUNE_THROW(CheckError,"LeafIterator gives non-leaf entity!");

    // check adaptation mark for leaf entity mark 
    CheckMark<Grid,checkMark>::check(g,lit);

    result = lit->geometry().local(lit->geometry().global(origin));
    typename Grid::ctype error = (result-origin).two_norm();
    if(!lit->type().isSingular() && error >= factorEpsilon * std::numeric_limits<typename Grid::ctype>::epsilon())
      {
        DUNE_THROW(CheckError, "|| geom.local(geom.global(" << origin
                   << ")) - origin || != 0 ( || " << result << " - origin || ) = " << error);
      };
    lit->geometry().integrationElement(origin);
    if(!lit->type().isSingular() && (int)Geometry::coorddimension == (int)Geometry::mydimension)
      lit->geometry().jacobianInverseTransposed(origin);

    lit->geometry().type();
    lit->type();
    if (lit->geometry().corners() > 0) {
	lit->geometry().corner( 0 );
    }
  }

}

template <class Grid>
void iteratorEquals (Grid &g)
{
  typedef typename Grid::template Codim<0>::LevelIterator LevelIterator;
  typedef typename Grid::template Codim<0>::LeafIterator LeafIterator;
  typedef typename Grid::template Codim<0>::HierarchicIterator HierarchicIterator;
  typedef typename Grid::template Codim<0>::LeafIntersectionIterator IntersectionIterator;
  typedef typename Grid::template Codim<0>::EntityPointer EntityPointer;

  // assignment tests
  LevelIterator l1 = g.template lbegin<0>(0);
  LevelIterator l2 = g.template lbegin<0>(0);
  LeafIterator L1 = g.template leafbegin<0>();
  LeafIterator L2 = g.template leafbegin<0>();

  // if grid empty, leave 
  if (l1 == g.template lend<0>(0)) 
    return;

  // check '==' consistency
  EntityPointer a = g.template levelView<Dune::All_Partition>(0).template begin<0>();
  EntityPointer i = g.template levelView<Dune::Interior_Partition>(0).template begin<0>();
  
  assert(
    (g.levelIndexSet(0).index(*a) != g.levelIndexSet(0).index(*i)) // index equal
    == (a != i) // entitpointer equal
    );
  assert(
    (g.levelIndexSet(0).index(*a) == g.levelIndexSet(0).index(*i)) // index equal
    == (a == i) // entitpointer equal
    );

  HierarchicIterator h1 = l1->hbegin(99);
  HierarchicIterator h2 = l2->hbegin(99);
  IntersectionIterator i1 = l1->ileafbegin();
  IntersectionIterator i2 = l2->ileafbegin();
  EntityPointer e1 = l1;
  EntityPointer e2 = h2;

  // assign
  l1 = l2;
  L1 = L2;
  h1 = h2;
  i1 = i2;
  e1 = e2;

  // equals
  #define TestEquals(i) { \
      i == e2; \
      i == l2; \
      i == h2; \
      i == L2; \
      if (i2 != l2->ileafend()) i == i2->inside(); \
      if (i2 != l2->ileafend() && i2->neighbor()) i == i2->outside(); \
  }
  TestEquals(e1);
  TestEquals(l1);
  TestEquals(h1);
  TestEquals(L1);
  if (i1 != l1->ileafend()) TestEquals(i1->inside());
  if (i1 != l1->ileafend() && i1->neighbor()) TestEquals(i1->outside());
}

template <class Grid>
void gridcheck (Grid &g)
{
  /*
   * first do the compile-test: this will not produce any code but
   * fails if an interface-component is missing
   */
  GridInterface<Grid>();

  enum { dim      = Grid :: dimension };
  enum { dimworld = Grid :: dimensionworld };
  typedef typename Grid  :: ctype ctype;
  typedef typename Grid  :: GridFamily GridFamily;

  // type of GridInterface == GridDefaultImplementation 
  typedef Dune::GridDefaultImplementation<dim,dimworld,ctype,GridFamily> GridIF;
  const GridIF & gridIF = g;
  // check functionality when grid is interpreted as reference to interface
  GridInterface<GridIF>::check(gridIF);
  /*
   * now the runtime-tests
   */
  const Grid & cg = g;
  iteratorEquals(g);
  iteratorEquals(cg);
  iterate<true>(g);
  iterate<false>(cg);
  zeroEntityConsistency(g);
  zeroEntityConsistency(cg);
  assertNeighbor(g);
  assertNeighbor(cg);
  // note that for some grid this might fail
  // then un comment this test 
  Dune :: checkIndexSet( g, g.leafView(), Dune :: dvverb );
  for( int level = 0; level <= g.maxLevel(); ++level )
    Dune :: checkIndexSet( g, g.levelView( level ), Dune :: dvverb, true );
}

#endif // #ifndef GRIDCHECK_CC
