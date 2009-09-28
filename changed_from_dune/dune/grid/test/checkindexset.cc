#ifndef DUNE_CHECKINDEXSET_CC
#define DUNE_CHECKINDEXSET_CC

#include <iostream>
#include <algorithm>

#include <dune/common/fvector.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/common/genericreferenceelements.hh>

#include <map>
#include <set>
#include <vector>
#include <limits>

/** @file
  @author Robert Kloefkorn
  @brief Provides a check of the grids index set.
*/

namespace Dune
{

// compare 2 FieldVectors 
template <typename ctype, int dim> 
bool compareVec(const FieldVector<ctype,dim> & vx1 , const FieldVector<ctype,dim> & vx2 ) 
{
  const ctype eps = 1e5 * std::numeric_limits<ctype>::epsilon();
  bool comp = true;
  for(int i=0;i<dim; i++)
  {
    if(std::abs( vx1[i] - vx2[i] ) > eps ) comp = false;
  }
  return comp; 
}

// check som functionality of grid 
template <int codim, class GridType, class EntityType, 
  class IndexSetType, class OutputStreamImp,
  class MapType1 , class MapType2 , class MapType3 >
void checkSubEntity ( const GridType & grid, 
    const EntityType &en , const IndexSetType & lset,
    OutputStreamImp & sout , MapType1 & subEntities , MapType2 & vertices ,
    MapType3 & vertexCoordsMap )
{
  enum { dim = EntityType::dimension };
  const int dimworld = GridType::dimensionworld;
  typedef typename EntityType::ctype coordType;

  const GeometryType type = en.type();
  assert( type == en.geometry().type() );

#if !defined DUNE_ENABLE_OLD_NUMBERING || defined NEW_SUBENTITY_NUMBERING
    const GenericReferenceElement< coordType, dim > &refElem
      = GenericReferenceElements< coordType, dim >::general( type );
#else
    const ReferenceElement< coordType, dim > &refElem
      = ReferenceElements< coordType, dim >::general( type );
#endif

      // check all subEntities of codimension  codim
      if(en.template count<codim>() != refElem.size(0,0,codim))
      {
        std::cerr << "entity index = " << lset.index(en)
                  << ", type = " << type
                  << std::endl
                  << "codim = " << codim
                  << std::endl
                  << "count<codim>() = " << en.template count<codim>()
                  << std::endl
                  << "refElem.size(codim) = " << refElem.size(0,0,codim)
                  << std::endl;
        DUNE_THROW(GridError,
                   "wrong number of subEntities of codim " << codim);
      }

      for( int subEntity = 0; subEntity < refElem.size( 0, 0, codim ); ++subEntity )
      {
        typedef std::pair< int, GeometryType > SubEntityKeyType;
        typedef Dune::GenericGeometry::MapNumberingProvider< dim > MapNumbering;

#if defined DUNE_ENABLE_OLD_NUMBERING && !defined NEW_SUBENTITY_NUMBERING
        const unsigned int topologyId = Dune::GenericGeometry::topologyId( type );
        const int duneSubEntity = MapNumbering::generic2dune( topologyId, subEntity, codim );
#endif

        {
          int numSubEntities = refElem.size( subEntity, codim, dim );

          // every entity have at least one vertex
          assert( numSubEntities > 0 );
          
          // create vectors of number of vertices on sub entity  
          std::vector<int> local (numSubEntities,-1);
          std::vector<int> global(numSubEntities,-1);

          for( int j = 0; j < numSubEntities; ++j )
          {
#if !defined DUNE_ENABLE_OLD_NUMBERING || defined NEW_SUBENTITY_NUMBERING
            local[ j ] = refElem.subEntity ( subEntity, codim, j, dim );
#else
            const int k = refElem.subEntity( duneSubEntity, codim, j, dim );
            local[ j ] = MapNumbering::dune2generic( topologyId, k, dim );
#endif
          }

          sout << numSubEntities << " vertices on subEntity< codim = " << codim << " >" << std::endl;
          sout << "check subentity [" << local[ 0 ];
          for( int j = 1; j < numSubEntities; ++j )
            sout << ", " << local[ j ];
          sout << "]" << std::endl;
         
          for( int j = 0; j < numSubEntities; ++j )
            global[ j ] = lset.subIndex( en, local[ j ], dim );
          
          typedef typename GridType::template Codim< codim >::EntityPointer SubEntityPointer;
          const SubEntityPointer subEntityPtr = en.template subEntity< codim >( subEntity );
          if( lset.subIndex( en, subEntity, codim ) != lset.index( *subEntityPtr) )
          {
            std::cerr << "Index for subEntity does not match." << std::endl;
            assert( lset.subIndex( en, subEntity, codim ) == lset.index( *subEntityPtr) );
          }

          SubEntityKeyType globalSubEntity( lset.index( *subEntityPtr ), subEntityPtr->type() );
          assert( globalSubEntity.first >= 0 );
          sout << "local subentity " << subEntity << " consider subentity with global key (" << globalSubEntity.first << "," << globalSubEntity.second << ") on en = " << lset.index(en) << std::endl;

          if( subEntityPtr->type() != subEntityPtr->geometry().type() )
          {
            std::cerr << "Geometry types for subEntity don't match." << std::endl;
            assert( subEntityPtr->type() == subEntityPtr->geometry().type() );
          }

          // assert that all sub entities have the same level 
          // otherwise one of the theoretical conditions is violated 
          assert( subEntityPtr.level() == en.level() );

          sout << "Found global numbers of entity [ ";
          for( int j = 0; j < numSubEntities; ++j )
            sout << global[ j ] << " ";
          sout << "]" << std::endl;

          for( int j = 0; j < numSubEntities; ++j )
          {
#if !defined DUNE_ENABLE_OLD_NUMBERING || defined NEW_SUBENTITY_NUMBERING
            const int gj = j;
#else
            const int tid = Dune::GenericGeometry::topologyId( refElem.type( duneSubEntity, codim ) );
            const int gj = Dune::GenericGeometry::MapNumberingProvider< dim-codim >::template dune2generic< dim-codim >( tid, j );
#endif

            {
              // get entity pointer of sub entity codim=dim (Vertex)
              typedef typename GridType::template Codim< dim >::EntityPointer VertexPointer;
              VertexPointer vxp = en.template subEntity< dim >( local[ j ] );

              FieldVector< coordType, dimworld > vx = vxp->geometry().corner( 0 );
              if(vertexCoordsMap.find(global[j]) != vertexCoordsMap.end())
              {
                FieldVector<coordType,dimworld> vxcheck ( vertexCoordsMap[global[j]] );
                if( ! compareVec( vxcheck, vx ) ) 
                {
                  std::cerr << "ERROR map global vertex [" << global[j] << "] vx " << vxcheck << " is not " << vx << "\n";
                  assert( compareVec( vxcheck, vx ) );
                }
              }
            }
            
            FieldVector< coordType, dimworld > vx = subEntityPtr->geometry().corner( gj );
            if(vertexCoordsMap.find(global[j]) != vertexCoordsMap.end())
            {
              FieldVector<coordType,dimworld> vxcheck ( vertexCoordsMap[global[j]] );
              if( ! compareVec( vxcheck, vx ) ) 
              {
                std::cerr << "Error map global vertex [" << global[j] << "] vx " << vxcheck << " is not " << vx << "\n";
                //assert( compareVec( vxcheck, vx ) );
              }
            }
            sout << "vx[" << global[j] << "] = "  << vx << "\n";
          }
      	  sout << "sort vector of global vertex\n";
         
          // sort vector of global vertex number for storage in map 
          // the smallest entry is the first entry 
          std::sort( global.begin(), global.end() );

          // check whether vertex key is already stored in map
          if(vertices.find(global) == vertices.end())
          {
              vertices[global] = globalSubEntity;
          }
          else 
          {
              SubEntityKeyType otherSubEntity = vertices[global];
              assert( globalSubEntity == otherSubEntity );
          }

          // check whether subEntity is already stored in map 
          if(subEntities.find(globalSubEntity) == subEntities.end() )
          {
            subEntities[globalSubEntity] = global;
          }
          else 
          {
            std::vector<int> globalcheck = subEntities[globalSubEntity];
            if(! (global == globalcheck ))
            {
              std::cerr << "For subEntity key (" << globalSubEntity.first << "," << globalSubEntity.second << ") \n";
              std::cerr << "Got ";
              for(int j=0 ;j<numSubEntities; j++ )
              {
                std::cerr << global[j] << " "; 
              }
              std::cerr << "\n";
              std::cerr << "Found ";
              for(int j=0 ;j<numSubEntities; j++ )
              {
                std::cerr << globalcheck [j] << " "; 
              }
              std::cerr << "\n";
              DUNE_THROW(Dune::GridError, "global != globalcheck");
            }
          }
        }
      } // end check sub entities 
      sout << "end check sub entities\n";
}


// check some functionality of grid 
template< int codim, class Grid, class GridView, class OutputStream >
void checkIndexSetForCodim ( const Grid &grid, const GridView &view,
                             OutputStream &sout, bool levelIndex )
{
  enum { dim = Grid :: dimension };
  enum { dimworld = Grid :: dimensionworld };
  
  typedef typename Grid :: ctype coordType;

  //typedef typename GridView :: template Codim< 0 > :: Entity EntityCodim0Type;
  typedef typename GridView :: IndexSet IndexSetType;
  typedef typename GridView :: template Codim< codim > :: Iterator IteratorType;

  const IndexSetType &lset = view.indexSet();
  
  sout <<"\n\nStart consistency check of index set \n\n";

  // ////////////////////////////////////////////////////////////////
  //   Check whether geomTypes() returns correct result
  // ////////////////////////////////////////////////////////////////

  std :: set< GeometryType > geometryTypes;
  
  const IteratorType endit = view.template end< codim >();
  IteratorType it = view.template begin< codim >();

  if (it == endit) return;

  for (; it!=endit; ++it) {
      // while we're here: check whether the GeometryTypes returned by the entity
      // and the Geometry match
      assert(it->type()==it->type());
      geometryTypes.insert(it->type());
  }

  bool geomTypesError = false;
  // Check whether all entries in the official geometry types list are contained in our self-computed one
  for (size_t i=0; i<lset.geomTypes(codim).size(); i++) 
      if (geometryTypes.find(lset.geomTypes(codim)[i])==geometryTypes.end())
          geomTypesError = true;
          

  // And vice versa
  for (std::set<GeometryType>::iterator it = geometryTypes.begin(); it!=geometryTypes.end(); ++it) {
      bool found = false;
      for (size_t i=0; i<lset.geomTypes(codim).size(); i++) 
          if (*it == lset.geomTypes(codim)[i]) {
              found = true;
              break;
          }
      
      if (!found)
          geomTypesError = true;

  }

  if (geomTypesError) {

      std::cerr << "There is a mismatch in the list of geometry types of codim " << codim << "." << std::endl;
      std::cerr << "Geometry types present in the grid are:" << std::endl;
      for (std::set<GeometryType>::iterator it = geometryTypes.begin(); it!=geometryTypes.end(); ++it)
          std::cerr << "  " << *it << std::endl;
      
      std::cerr << std::endl << "but the method geomTypes() returned:" << std::endl;
      for (size_t j=0; j<lset.geomTypes(codim).size(); j++)
          std::cerr << "  " << lset.geomTypes(codim)[j] << std::endl;
      
      DUNE_THROW(GridError, "!");
  }

  //*****************************************************************
  // check size of index set 
  int gridsize = 0;
  {
    typedef typename GridView :: template Codim< codim > :: Iterator IteratorType;

    int count = 0;
    const IteratorType endit = view.template end< codim >();
    for( IteratorType it = view.template begin< codim >(); it != endit; ++it )
      ++count;

    int lsetsize = lset.size(codim);
    if( count != lsetsize)
    {
      derr << "WARNING: walk = "<< count << " entities | set = "
        << lsetsize << " for codim " << codim << std::endl;
    }
    gridsize = count;
    // lsetsize should be at least the size of iterated entities 
    assert( count <= gridsize );
  }

  { 
    typedef typename GridView :: template Codim< 0 > :: Iterator Iterator;
    typedef typename Grid :: Traits :: LocalIdSet LocalIdSetType;
    typedef typename LocalIdSetType :: IdType IdType;

    std::set< IdType > entityfound;

    const Iterator endit = view.template end< 0 >();
    Iterator it = view.template begin< 0 >();
    if( it == endit )
      return;

    const LocalIdSetType &localIdSet = grid.localIdSet();
    for( ; it != endit; ++it )
    {
      const typename Iterator::Entity &entity = *it;
      if( !lset.contains( entity ) )
      {
        std::cerr << "Error: IndexSet does not contain all entities." << std::endl;
        assert( false );
      }
      const int subcount = entity.template count< codim >();
      for( int i = 0; i < subcount; ++i )
      {
        const IdType id = localIdSet.id( *(entity.template subEntity< codim >( i ) ) );
        entityfound.insert( id );
      }
    }

    if( (size_t)gridsize != entityfound.size() )
    {
      std::cerr << "Warning: gridsize = " << gridsize << " entities"
               << ", set of entities = " << entityfound.size()
               << " [codim " << codim << "]" << std::endl;
    }

    // gridsize should be at least the size of found entities 
    //assert( gridsize <= (int) entityfound.size() );
  }

  //******************************************************************
 
  typedef std::pair < int , GeometryType > SubEntityKeyType; 
  typedef std::map < int , std::pair<int,int> > subEntitymapType;
  std::map < SubEntityKeyType , std::vector<int> > subEntities;
  std::map < std::vector<int> , SubEntityKeyType > vertices;
  
  std::map < int , FieldVector<coordType,dimworld> > vertexCoordsMap;
  // setup vertex map , store vertex coords for vertex number 
  {
    unsigned int count = 0;
    typedef typename GridView :: template Codim< dim > :: Iterator VxIterator;
    const VxIterator end = view.template end< dim >();
    for( VxIterator it = view.template begin< dim >(); it != end; ++it )
    {
      ++count;
      // get coordinates of vertex 
      FieldVector< coordType, dimworld > vx ( it->geometry().corner( 0 ) );

      // get index of vertex 
      sout << "Vertex " << vx << "\n";
      assert( lset.contains ( *it ) );
      int idx = lset.index( *it );
      
      sout << "Vertex " << idx << " = [" << vx << "]\n";
      
      // if vertex not in map insert it 
      if( vertexCoordsMap.find(idx) == vertexCoordsMap.end())
        vertexCoordsMap[idx] = vx;
    }
    sout << "Found " << vertexCoordsMap.size() << " vertices for that index set!\n\n";

    // check whether size of map equals all found vertices 
    assert( vertexCoordsMap.size() == count );
    
    // check whether size of vertices of set equals all found vertices 
    sout << "Checking size of vertices " 
	 << count 
	 << " equals all found vertices " 
	 << (unsigned int)lset.size(Dune::GeometryType(0))
	 << "\n";
    // assertion goes wrong for parallel grid since no iteration over ghost
    // subentities
    assert( count == (unsigned int)lset.size(Dune::GeometryType(0)) );
  }

  {
    typedef typename GridView :: template Codim< 0 > :: Iterator Iterator;
    // choose the right reference element 
    const Iterator refend = view.template end< 0 >();
    Iterator refit = view.template begin< 0 >();
    assert( refit != refend );
      
    GeometryType type = refit->type();
    
#if !defined DUNE_ENABLE_OLD_NUMBERING || defined NEW_SUBENTITY_NUMBERING
    const GenericReferenceElement< coordType, dim > &refElem
      = GenericReferenceElements< coordType, dim >::general( type );
#else
    const ReferenceElement< coordType, dim > &refElem
      = ReferenceElements< coordType, dim >::general( type );
#endif

    // print dune reference element 
    sout << "Dune reference element provides: \n";
    for(int i = 0; i < refElem.size( codim ); ++i )
    {
      sout << i << " subEntity [";
      int s = refElem.size(i,codim,dim);
      for(int j=0; j<s; j++)
      {
        sout << refElem.subEntity(i , codim , j , dim );
        if(j != s-1) sout << ",";
      }
      sout << "]\n";
    }
  }

  {
    typedef typename GridView :: template Codim< 0 > :: Iterator Iterator;
    const Iterator endit = view.template end< 0 >();
    for( Iterator it = view.template begin< 0 >(); it != endit; ++it )
    {
      // if (it->partitionType()==4) continue;
      sout << "****************************************\n";
      sout << "Element = " << lset.index(*it) << " on level " << it->level () << "\n";
      sout << "Vertices      = [";
      int svx = it->template count<dim>();

      // print all vertex numbers 
      for( int i = 0; i < svx; ++i )
      {
        const typename IndexSetType::IndexType idx = lset.subIndex( *it, i, dim );
        sout << idx << (i < svx-1 ? ", " : "]\n");
      }

      // print all vertex coordinates 
      sout << "Vertex Coords = ["; 
      for( int i = 0; i < svx; ++i )
      {
        // get entity pointer of sub entity codim=dim (Vertex)
        typedef typename Grid::template Codim< dim >::EntityPointer VertexPointer;
        VertexPointer vxp = it->template subEntity< dim >( i );
       
        // get coordinates of entity pointer 
        FieldVector< coordType, dimworld > vx( vxp->geometry().corner( 0 ) );

        // output vertex coordinates 
        sout << vx << (i < svx-1 ? ", " : "]\n");
        
        const typename IndexSetType::IndexType vxidx = lset.subIndex( *it, i, dim );
        
        // the subIndex and the index for subEntity must be the same
        if( vxidx != lset.index( *vxp ) )
        {
          std::cerr << "Error: index( *subEntity< dim >( i ) ) != subIndex( entity, i, dim )" << std::endl;
          assert( vxidx == lset.index( *vxp ) );
        }
        
#if 0
        typedef GenericGeometry::MapNumberingProvider< dim > Numbering;
        const unsigned int tid = GenericGeometry::topologyId( it->type() );
        const int di = Numbering::template generic2dune< dim >( tid, i );

        // static and dynamic method must yield the same result
        if( vxidx != lset.template subIndex< dim >( *it, di ) )
        {
          std::cerr << "Error: subIndex< dim >( entity, generic2dune( i ) ) != subIndex( entity, i, dim )" << std::endl;
          assert( vxidx == lset.template subIndex< dim >( *it, dim ) );
        }
#endif
          
        // check whether the coordinates are the same 
      	assert(vertexCoordsMap.find(vxidx)!=vertexCoordsMap.end());
        FieldVector<coordType,dimworld> vxcheck ( vertexCoordsMap[vxidx] );
        if( !compareVec( vxcheck, vx ) ) 
        {
          std::cerr << "Error: Inconsistent map of global vertex " << vxidx
                    << ": " << vxcheck << " != " << vx
                    << "  (type = " << it->partitionType() << ")" << std::endl;
          assert( compareVec( vxcheck, vx ) );
        }
      }

      ////////////////////////////////////////////////////////////
      // check sub entities 
      //////////////////////////////////////////////////////////// 
//       checkSubEntity< codim >( grid, *it, lset, sout,
//                                subEntities, vertices, vertexCoordsMap );

      // check neighbors 
      if( codim == 1 )
      {
        typedef typename GridView :: IntersectionIterator IntersectionIterator;

        const std :: string name = grid.name();
        if( !levelIndex || (name != "AlbertaGrid") )
        {
          const IntersectionIterator endnit = view.iend( *it );
          for( IntersectionIterator nit = view.ibegin( *it ); nit != endnit; ++nit )
          {
            if( !nit->neighbor() )
              continue;

            checkSubEntity< codim >( grid, *(nit->outside()), lset, sout,
                                     subEntities, vertices, vertexCoordsMap );
          }
        }
        else
        {
          static bool called = false;
          if( !called ) 
          {
            std::cerr << "WARNING: skip indices test using LevelIntersectionIterator for AlbertaGrid!\n";
            called = true;
          }
        }
      }
    }
  }
}


  template< class Grid, class GridView, class OutputStream, int codim, bool hasCodim >
  struct CheckIndexSet
  {
    static void checkIndexSet ( const Grid &grid, const GridView &view,
                                OutputStream &sout, bool levelIndex )
    {
      checkIndexSetForCodim< codim >( grid, view, sout, levelIndex );
      typedef Dune :: Capabilities :: hasEntity< Grid, codim-1 > hasNextCodim;
      CheckIndexSet< Grid, GridView, OutputStream, codim-1, hasNextCodim :: v >
        :: checkIndexSet( grid, view, sout, levelIndex );
    }
  };

  template< class Grid, class GridView, class OutputStream, int codim >
  struct CheckIndexSet< Grid, GridView, OutputStream, codim, false >
  {
    static void checkIndexSet ( const Grid &grid, const GridView &view,
                                OutputStream &sout, bool levelIndex )
    {
      derr << "WARNING: Entities for codim " << codim
           << " are not being tested!" << std::endl;
      typedef Dune :: Capabilities :: hasEntity< Grid, codim-1 > hasNextCodim;
      CheckIndexSet< Grid, GridView, OutputStream, codim-1, hasNextCodim :: v >
        :: checkIndexSet( grid, view, sout, levelIndex );
    }
  };

  template< class Grid, class GridView, class OutputStream >
  struct CheckIndexSet< Grid, GridView, OutputStream, 0, true >
  {
    static void checkIndexSet ( const Grid &grid, const GridView &view,
                                OutputStream &sout, bool levelIndex )
    {
      checkIndexSetForCodim< 0 >( grid, view, sout, levelIndex );
    }
  };

  template< class Grid, class GridView, class OutputStream >
  void checkIndexSet ( const Grid &grid, const GridView &view,
                       OutputStream &sout,  bool levelIndex = false )
  {
    CheckIndexSet< Grid, GridView, OutputStream, Grid :: dimension, true >
      :: checkIndexSet ( grid, view, sout, levelIndex );
  }

} // end namespace Dune 

#endif
