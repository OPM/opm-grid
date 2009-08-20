#ifndef DUNE_CHECK_INTERSECTIONITERATOR_CC
#define DUNE_CHECK_INTERSECTIONITERATOR_CC

#include <cmath>

#include <dune/grid/common/quadraturerules.hh>
#include <dune/grid/common/genericreferenceelements.hh>

/** \file
    \brief Tests for the IntersectionIterator
*/

struct CheckIntersectionIteratorErrorState
{
  unsigned int sumNormalsNonZero;

  CheckIntersectionIteratorErrorState ()
  : sumNormalsNonZero( 0 )
  {}
};



/** \brief Helper routine: Test a general geometry
 */
template <class GeometryImp>
void checkGeometry(const GeometryImp& geometry)
{
    using namespace Dune;

    // Get the dimensions and types
    const int dim      = GeometryImp::mydimension;
    const int dimworld = GeometryImp::coorddimension;

    typedef typename GeometryImp::ctype ctype;
#if NEW_SUBENTITY_NUMBERING
    // Get the corresponding reference element
    const GenericReferenceElement< double, dim > &genericRefElement
      = GenericReferenceElements< double, dim >::general( geometry.type() );

    // Check whether the number of corners is correct
    if (geometry.corners() != genericRefElement.size(dim))
        DUNE_THROW(GridError, "Geometry has wrong number of corners!");

    // check consistency between operator[] and global()
    for( int i = 0; i < genericRefElement.size( dim ); ++i )
    {
      FieldVector< double, dim > localPos = genericRefElement.position( i, dim );
      if( (geometry.corner( i ) - geometry.global( localPos )).infinity_norm() > 1e-6 )
        DUNE_THROW( GridError, "Methods 'corner' and 'global' are inconsistent." );
    }
#else
    // Get the corresponding reference element
    const ReferenceElement< double, dim > &refElement
      = ReferenceElements< double, dim >::general( geometry.type() );

    // Check whether the number of corners is correct
    if (geometry.corners() != refElement.size(dim))
        DUNE_THROW(GridError, "Geometry has wrong number of corners!");

    // check consistency between operator[] and global()
    for( int i = 0; i < refElement.size( dim ); ++i )
    {
      FieldVector< double, dim > localPos = refElement.position( i, dim );
      if( (geometry.corner( i ) - geometry.global( localPos )).infinity_norm() > 1e-6 )
        DUNE_THROW( GridError, "Methods 'corner' and 'global' are inconsistent." );
    }
#endif
    // Use a quadrature rule to create a few test points for the following checks
    const QuadratureRule<double, dim>& quad 
        = QuadratureRules<double, dim>::rule(geometry.type(), 2);

    for( size_t i = 0; i < quad.size(); ++i )
    {

        const FieldVector< double, dim > &testPoint = quad[i].position();

#if NEW_SUBENTITY_NUMBERING
        // Check whether point is within the intersection
        if( !genericRefElement.checkInside( testPoint ) )
        {
          std::cerr << "Test point (" << testPoint
                    << ") not within reference element." << std::endl;
        }
#endif    
        // Transform to global coordinates
        FieldVector< ctype, dimworld > global = geometry.global( testPoint );

        // The back to local coordinates
        FieldVector< ctype, dim > local = geometry.local( global );
    
        // check for correctness
        if ((testPoint-local).infinity_norm() > 1e-6)
            DUNE_THROW(GridError, "local() and global() are not inverse to each other!");
    
        // The integration element at the element center
        ctype intElement = geometry.integrationElement(testPoint);
        if (intElement <=0)
            DUNE_THROW(GridError, "nonpositive integration element found!");
    
#if 0
    // This method exists in the interface, but it is not expected to work
    // unless dim==dimworld
    const FieldMatrix<ctype, dim, dim> jacobi
        = intersectionGlobal.jacobianInverseTransposed(testPoint);
#endif
    }
}

// Check that normal and normal2 pointing in the same direction
template< class ctype, int dimworld, class String >
inline void checkParallel ( const Dune :: FieldVector< ctype, dimworld > &normal,
                            const Dune :: FieldVector< ctype, dimworld > &normal2,
                            const String & name )
{
  if( std :: abs( normal * normal2 - normal.two_norm() * normal2.two_norm() ) > 1e-8 )
  {
    std :: cerr << "Error: " << name
                << " does not point in the direction of outerNormal."
                << std :: endl;
    std :: cerr << "       " << name << " = " << normal2
                << ", normal = " << normal << std :: endl;
    assert( false );
  }
}


// Check whether the normal is orthogonal to the intersection, i.e.,
// whether (J^-1 * n) = 0. Here J is the jacobian of the intersection
// geometry (intersectionGlobal) and n is a normal.
template< class ctype, int dimworld, int facedim, class String >
inline void checkJIn ( const Dune :: FieldVector< ctype, dimworld > &normal,
                       const Dune :: FieldMatrix< ctype, dimworld, facedim > &jit,
                       const String & name )
{
  Dune :: FieldVector< ctype, facedim > x( ctype( 0 ) );
  jit.umtv( normal, x );
  if (x.infinity_norm() > 1e-8)
  {
    std :: cerr << "Error:  (J^-1 * n) != 0." << std :: endl;
    std :: cerr << "       " << name << " = " << normal
                << std :: endl;
    std :: cerr << "       J^-1^T = \n" << jit << std::endl;
    assert( false );
  }
}


/** \brief Test the IntersectionIterator
*/
template <class GridViewType, class ErrorState >
void checkIntersectionIterator(const GridViewType& view,
                               const typename GridViewType::template Codim<0>::Iterator& eIt,
                               ErrorState &errorState )
{
  using namespace Dune;

  typedef typename GridViewType::Grid GridType;
  typedef typename GridViewType::IntersectionIterator IntersectionIterator;

  const GridType& grid = view.grid();
  const bool checkOutside = (grid.name() != "AlbertaGrid");
  const typename GridViewType::IndexSet& indexSet = view.indexSet();

  typedef typename GridType :: ctype ctype;
  
  const int dim      = GridType :: dimension;
  const int dimworld = GridType :: dimensionworld;

  FieldVector< ctype,dimworld > sumNormal( ctype( 0 ) ); 

  // /////////////////////////////////////////////////////////
  //   Check the types defined by the iterator
  // /////////////////////////////////////////////////////////
  dune_static_assert((is_same<
      typename IntersectionIterator::ctype,
          typename GridType::ctype>::value),
      "IntersectionIterator has wrong ctype");

  dune_static_assert((is_same<
      typename IntersectionIterator::Intersection,
          typename GridViewType::Intersection>::value),
      "IntersectionIterator has wrong Intersection type");

  dune_static_assert((static_cast<int>(IntersectionIterator::dimension)
    == static_cast<int>(GridType::dimension)),"IntersectionIterator has wrong dimension");

  dune_static_assert((static_cast<int>(IntersectionIterator::dimensionworld)
    == static_cast<int>(GridType::dimensionworld)),"IntersectionIterator has wrong dimensionworld");

  IntersectionIterator iIt    = view.ibegin(*eIt);
  IntersectionIterator iEndIt = view.iend(*eIt);
  
  bool hasBoundaryIntersection = false;

  for (;iIt!=iEndIt; ++iIt) 
  { 
      // //////////////////////////////////////////////////////////////////////
      //   Compute the integral of the outer normal over the whole element.
      //   This has to be zero.
      // //////////////////////////////////////////////////////////////////////
      const int interDim = IntersectionIterator::LocalGeometry::mydimension;
      const QuadratureRule< double, interDim > &quad
        = QuadratureRules< double, interDim >::rule( iIt->type(), 3 );
      //const QuadratureRule<double, interDim>& quad 
      //    = QuadratureRules<double, interDim>::rule(iIt->intersectionSelfLocal().type(), interDim);

      typedef typename IntersectionIterator::Entity EntityType; 
      typedef typename EntityType::EntityPointer EntityPointer;

      assert(eIt == iIt->inside());

      // check that boundary id has positive value and that intersection is
      // conform
      if( iIt->boundary() )
      {
        // entity has boundary intersections 
        hasBoundaryIntersection = true;
        
        if( iIt->boundaryId() < 0 )
        {
          DUNE_THROW(GridError, "boundary id has negative value (" << iIt->boundaryId() << ") !");
        }
        if( ! iIt->conforming() )
        {
          DUNE_THROW(GridError, "Boundary intersection should be conforming!");
        }
      }

      // //////////////////////////////////////////////////////////////////////
      //   Check whether the 'has-intersection-with'-relation is symmetric
      // //////////////////////////////////////////////////////////////////////

      if (iIt->neighbor() && checkOutside ) 
      {
          EntityPointer outside = iIt->outside();
          bool insideFound = false;

          IntersectionIterator outsideIIt    = view.ibegin(*outside);
          IntersectionIterator outsideIEndIt = view.iend(*outside);

          for (; outsideIIt!=outsideIEndIt; ++outsideIIt) {

              if (outsideIIt->neighbor() && outsideIIt->outside() == iIt->inside()) {

                  if (outsideIIt->indexInInside() != iIt->indexInOutside())
                      DUNE_THROW(GridError, "outside()->outside() == inside(), but with incorrect numbering!");
                  else
                      insideFound = true;

              }

          }

          if (!insideFound)
              DUNE_THROW(GridError, "Could not find inside() through intersection iterator of outside()!");

      }
      else if (!checkOutside) 
      {
        static bool called = false;
        if(!called)
        {
          derr << "WARNING: skip reverse intersection iterator test for " << grid.name() << "!"<< std::endl;
          called = true;
        }
      }

      // Check if conforming() methods is compatible with static
      // information on GridView
      if ( GridViewType::conforming && !iIt->conforming()) 
      {
        DUNE_THROW(GridError, "GridView says conforming but intersection is not conforming!");
      }

      // /////////////////////////////////////////////////////////////
      //   Check the consistency of numberInSelf, numberInNeighbor
      //   and the indices of the subface between.
      // /////////////////////////////////////////////////////////////
      if ( iIt->conforming() && iIt->neighbor() ) 
      {
        EntityPointer outside = iIt->outside();
        const int indexInInside  = iIt->indexInInside();
        const int indexInOutside = iIt->indexInOutside();

        if( indexSet.subIndex( *eIt, indexInInside, 1 ) != indexSet.subIndex( *outside, indexInOutside, 1 ) )
        {
          std::cerr << "Error: Index of conforming intersection differs when "
                    << "obtained from inside and outside." << std::endl;
          std::cerr << "       inside index = " << indexSet.subIndex( *eIt, indexInInside, 1 )
                    << ", outside index = " << indexSet.subIndex( *outside, indexInOutside, 1 ) << std::endl;
          assert( false );
        }

        const typename GridType::LocalIdSet &localIdSet = grid.localIdSet();
        if( localIdSet.subId( *eIt, indexInInside, 1 ) != localIdSet.subId( *outside, indexInOutside, 1 ) )
        {
          std::cerr << "Error: Local id of conforming intersection differs when "
                    << "obtained from inside and outside." << std::endl;
          std::cerr << "       inside id = " << localIdSet.subId( *eIt, indexInInside, 1 )
                    << ", outside id = " << localIdSet.subId( *outside, indexInOutside, 1 ) << std::endl;
          assert( false );
        }

        const typename GridType::GlobalIdSet &globalIdSet = grid.globalIdSet();
        if( globalIdSet.subId( *eIt, indexInInside, 1 ) != globalIdSet.subId( *outside, indexInOutside, 1 ) )
        {
          std::cerr << "Error: Global id of conforming intersection differs when "
                    << "obtained from inside and outside." << std::endl;
          std::cerr << "       inside id = " << globalIdSet.subId( *eIt, indexInInside, 1 )
                    << ", outside id = " << globalIdSet.subId( *outside, indexInOutside, 1 ) << std::endl;
          assert( false );
        }
      }

      // //////////////////////////////////////////////////////////
      //   Check the geometry returned by intersectionGlobal()
      // //////////////////////////////////////////////////////////
      typedef typename IntersectionIterator::Geometry Geometry;
      const Geometry &intersectionGlobal = iIt->geometry();

      checkGeometry(intersectionGlobal);

      if( intersectionGlobal.type() != iIt->type() )
      {
        std::cerr << "Error: Reference geometry type returned by intersection "
                  << "differs from the one returned by the global geometry." << std::endl;
        assert( false );
      }

      // //////////////////////////////////////////////////////////
      //   Check the geometry returned by intersectionSelfLocal()
      // //////////////////////////////////////////////////////////

      const typename IntersectionIterator::LocalGeometry &intersectionSelfLocal = iIt->geometryInInside();
      checkGeometry(intersectionSelfLocal);

      //  Check the consistency of intersectionSelfLocal() and intersectionGlobal
      
      if (intersectionSelfLocal.corners() != intersectionGlobal.corners())
          DUNE_THROW(GridError, "Geometry of intersection is inconsistent from left hand side and global view!");
      
      // Use a quadrature rule as a set of test points
      for (size_t i=0; i<quad.size(); i++) 
      {
          const FieldVector< ctype, dim-1 > &pt = quad[ i ].position();

          // Check outer normal
          const FieldVector< ctype, dimworld > normal = iIt->outerNormal( pt );

          // Check normal vector is orthogonal to all vectors connecting
          // the vertices
          for (int c=1; c<intersectionGlobal.corners(); c++)
          {
            Dune::FieldVector< ctype, dimworld > x = intersectionGlobal.corner( c-1 );
            x -= intersectionGlobal.corner( c );
            if( x*normal >= 10*std::numeric_limits< ctype >::epsilon() )
            {
              std::cerr << "outerNormal not orthogonal to line between corner "
                        << (c-1) << " and corner " << c << "." << std::endl;
              std::cerr << "Note: This is ok for curved faces, though." << std::endl;
              assert( false );
            }
          }
          
          // Check integration outer normal
          const FieldVector< ctype, dimworld > intNormal =
            iIt->integrationOuterNormal( pt );
          sumNormal.axpy( quad[ i ].weight(), intNormal );

          const ctype det = intersectionGlobal.integrationElement( pt );
          if( std :: abs( det - intNormal.two_norm() ) > 1e-8 )
          {
            std :: cerr << "Error: integrationOuterNormal yields wrong length."
                        << std :: endl;
            std :: cerr << "       |integrationOuterNormal| = " << intNormal.two_norm()
                        << ", integrationElement = " << det << std :: endl;
            assert( false );
          }

          checkParallel ( intNormal, normal, "integrationOuterNormal");

          // Check unit outer normal
          const FieldVector< ctype, dimworld > unitNormal = iIt->unitOuterNormal( pt );

          if( std :: abs( ctype( 1 ) - unitNormal.two_norm() ) > 1e-8 )
          {
            std :: cerr << "Error: unitOuterNormal yields wrong length." << std :: endl;
            std :: cerr << "       |unitOuterNormal| = " << unitNormal.two_norm()
                        << std :: endl;
            assert( false );
          }

          checkParallel ( unitNormal, normal, "unitOuterNormal");

	  if (!iIt->type().isSingular()) {
	      // Check JacobianInverseTransposed
	      // (J^-1 * n) = 0. Here J is the jacobian of the intersection
	      // geometry (intersectionGlobal) and n is a normal.
	      const FieldMatrix< ctype, dimworld, dim-1 > &jit
		  = intersectionGlobal.jacobianInverseTransposed( pt );
	      checkJIn( normal, jit, "outerNormal" );
	      checkJIn( intNormal, jit, "integrationOuterNormal" );
	      checkJIn( unitNormal, jit, "unitOuterNormal" );
	  }

	  // check intersectionSelfLocal
	  FieldVector<double,dimworld> globalPos = intersectionGlobal.global(quad[i].position());
	  FieldVector<double,dimworld> localPos  = eIt->geometry().global(intersectionSelfLocal.global(quad[i].position()));
          if( !iIt->type().isSingular() &&  (globalPos - localPos).infinity_norm() > 1e-6 )
            DUNE_THROW( GridError, "entity.geometry().global( intersection.geometryInInside().global( x ) ) != intersection.geometry().global( x ) at x = " << quad[ i ].position() << "." );
      }

      // ////////////////////////////////////////////////////////////////
      //   Check the geometry returned by intersectionNeighborLocal()
      // ////////////////////////////////////////////////////////////////

      if (iIt->neighbor() ) 
      {
          const typename IntersectionIterator::LocalGeometry &intersectionNeighborLocal = iIt->geometryInOutside();
          
          checkGeometry(intersectionNeighborLocal);

          if (intersectionSelfLocal.corners() != intersectionNeighborLocal.corners())
              DUNE_THROW(GridError, "Geometry of intersection is inconsistent from left and right hand side!");

          // (Ab)use a quadrature rule as a set of test points
          const int interDim = IntersectionIterator::LocalGeometry::mydimension;
          const QuadratureRule<double, interDim>& quad 
              = QuadratureRules<double, interDim>::rule(intersectionNeighborLocal.type(), 2);

          for (size_t i=0; i<quad.size(); i++) 
          {

              FieldVector<double,dimworld> globalPos = intersectionGlobal.global(quad[i].position());
              FieldVector<double,dimworld> localPos  = iIt->outside()->geometry().global(intersectionNeighborLocal.global(quad[i].position()));

              if ( !iIt->type().isSingular() && (globalPos - localPos).infinity_norm() > 1e-6)
                  DUNE_THROW(GridError, "global( intersectionNeighborLocal(global() ) is not the same as intersectionGlobal.global() at " << quad[i].position() << "!");
              
          }

      }

  }

  // check implementation of hasBoundaryIntersections
  if( hasBoundaryIntersection != eIt->hasBoundaryIntersections() )
  {
    DUNE_THROW(GridError,"Entity::hasBoundaryIntersections implemented wrong! \n");
  }

  // Chech whether integral over the outer normal is zero
  // 
  // Note: This is wrong on curved surfaces (take, e.g., the upper half sphere).
  //       Therefore we only issue a warning.
  if( (sumNormal.two_norm() > 1e-8) && (eIt->partitionType() != Dune::GhostEntity) )
  {
    std :: cout << "Warning: Integral over outer normals is nonzero: "
                << sumNormal << std :: endl;
    ++errorState.sumNormalsNonZero;
  }
}

/** \brief Test both IntersectionIterators 
 */
template <class GridViewType>
void checkViewIntersectionIterator(const GridViewType& view) {

    using namespace Dune;

    typedef typename GridViewType ::template Codim<0>::Iterator 
	    ElementIterator;
    ElementIterator eIt    = view.template begin<0>();
    ElementIterator eEndIt = view.template end<0>();
    
    CheckIntersectionIteratorErrorState errorState;
    for (; eIt!=eEndIt; ++eIt) 
      checkIntersectionIterator( view, eIt, errorState );

    if( errorState.sumNormalsNonZero > 0 )
    {
      std :: cerr << "Warning: Integral over outer normals is not always zero."
                  << std :: endl;
      std :: cerr << "         This behaviour may be correct for entities with"
                  << " nonzero curvature." << std :: endl;;
    }
}

template <class GridType>
void checkIntersectionIterator(const GridType& grid, bool skipLevelIntersectionTest = false) {

    using namespace Dune;

    // Loop over all levels
    if(skipLevelIntersectionTest) 
    {
      std::cerr<<"WARNING: skip test of LevelIntersectionIterator! \n";
    }
    else   
    {
      for (int i=0; i<=grid.maxLevel(); i++) 
	checkViewIntersectionIterator(grid.levelView(i));
    }

    // test leaf intersection iterator 
    {
      checkViewIntersectionIterator(grid.leafView());
    }

}

#endif
