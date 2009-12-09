#ifndef DUNE_GRID_INTERSECTION_HH
#define DUNE_GRID_INTERSECTION_HH

#include <dune/grid/common/grid.hh>

#include <dune/grid/genericgeometry/conversion.hh>

namespace Dune
{

/** \brief %Intersection of a mesh entities of codimension 0 ("elements")
    with a "neighboring" element or with the domain
    boundary.  

   Template parameters are:

   - <tt>GridImp</tt> Type that is a model of Dune::Grid
   - <tt>IntersectionImp</tt> Class template that is a model of 
   Dune::Intersection
    
   <h2>Overview</h2>
   
   Intersections are codimension 1 objects. These
   intersections are accessed via an Intersection. This allows
   the implementation of non-matching grids, as one face can now
   consist of several intersections.
   In a conforming mesh such an intersection corresponds to an entity of
   codimension 1 but in the general non-conforming case there will be no entity
   in the mesh that directly corresponds to the intersection. Thus, the
   Intersection describes these intersections implicitly.

   <H2>Engine Concept</H2>

   The Intersection class template wraps an object of type IntersectionImp
   and forwards all member 
   function calls to corresponding members of this class. In that sense Intersection
   defines the interface and IntersectionImp supplies the implementation.

   <h2>Methods neighbor and boundary </h2>
   
   The following holds for both the level and the leaf intersection
   :
   The %intersection  is started on a codimension 0 entity of the grid.
   If this entity belongs to the interior or the overlap region 
   (see. ???) then the union of all intersections is identical to the
   boundary of the entity. On ghost elements the  only stops
   on the border of the domain, i.e., only on the intersections with
   entities in the interior region. Depending on the boolean values returned by 
   the methods %boundary() and %neighbor() 
   one can detect the position of an intersection 
   relative to the boundary. In general
   the method boundary() returns true if and only if the intersection is
   part of the physical boundary of the domain. The method neighbor() returns
   true only if the method outside() has a well defined return value. 

  The following cases are possible if the intersection  is 
  started on an entity in the interior or overlap region. More
  details are given below:
   
  <table>
  <tr>
  <td></td><td>intersection</td><td>neighbor()</td><td>boundary()</td><td>outside()</td></tr>
  <tr>
  <td>1</td><td>with inner, overlap <br>
                or ghost entity</td>
  <td>true</td><td>false</td>
  <td>the neighbor entity</td></tr>
  <tr>
  <td>2</td><td>on domain boundary</td>
  <td>false</td><td>true</td><td><em>undefined</em></td></tr>
  <tr>
  <td>3</td><td>on periodic boundary</td>
  <td>true</td><td>true</td><td>Ghost-/Overlap cell<br>(with transformed geometry)</td></tr>
  <tr>
  <td>4</td><td>on processor boundary</td>
  <td>false <em>if grid has no ghosts</em><br>true <em>otherwise</em></td><td>false </td>
  <td>ghost entity <em>(if it exists)</em></td></tr>
  </table>

  <dl>
  <dt>Inner Intersections: </dt>
  <dd>
     The type of the neighboring entity can be determined through
     methods defined on the outside entity.
  </dd>
  <dt>Handling physical boundaries: </dt>
  <dd>
     Different types of physical boundaries can be modeled using either
     the global coordinates of the intersection or by using the
     boundaryID method. On some grids (AluGrid, AlbertaGrid) this
     method returns an integer value which can be individually assigned
     to each boundary intersection of the macro grid and which is
     prolonged to higher levels during grid refinement.<br>
     A more general concept will be included in latter releases along the
     following guidelines:
     - We require differently constructed geometries outside the domain
     - The kind of construction depends on the discrete problem 
     - Therefor these constructions can't be part of the Grid interface
     - Utility classes are required to do this construction
     - The utility classes must be parameterized with the intersection (in our 
       case the Intersection)
     - The utility classes return a suitable transformation of the inner() 
       entitys geometry (with respect to the intersection), e.g.,
       reflection at the intersection
       point reflection
       reflection combined with translation...
     .
  </dd>     
  <dt>Handling periodic boundaries: </dt>
  <dd>
     - The Intersection stops at periodic boundaries
     - periodic grids are handled in correspondence to parallel grids
     - %At the periodic boundary one can adjust an overlap- or ghost-layer.
     - outer() returns a ghost or overlap cell (for ghost and overlap look into 
       the documentation of the parallel grid interface)
     - outer() cell has a periodically transformed geometry (so that one does 
       not see a jump or something like that)
     - outer() cell has its own index
     - outer() cell has the same id as the corresponding "original" cell
     .
  </dd>     
  <dt>Processor boundaries: </dt>
  <dd>
     At processor boundaries, i.e. when an element has an intersection with 
     another element 
     in the sequential grid but this element is only stored in other processors 
     the intersection  stops but neither 
     neighbor() nor boundary()
     are true.
  </dd>
  </dl> 


   <h2>Geometry of an intersection</h2>

   The method intersectionGlobal returns a geometry mapping the intersection
   as a codim one structure to global coordinates. The methods 
   intersectionSelfLocal and intersectionNeighborLocal return geometries
   mapping the intersection into the reference elements of the 
   originating entity and the neighboring entity, respectively. 
   The numberInSelf and numberInNeighbor methods return the codim one
   subentities which contain the intersection. 
     

   @ingroup GIIntersectionIterator
 */
template<class GridImp, template<class> class IntersectionImp>
class Intersection
{
  IntersectionImp<const GridImp> real;

  enum { dim=GridImp::dimension };
  enum { dimworld=GridImp::dimensionworld };

public:
  
  // type of real implementation 
  typedef IntersectionImp<const GridImp> ImplementationType;
  
    /** \brief Type of entity that this Intersection belongs to */
  typedef typename GridImp::template Codim<0>::Entity Entity;

    /** \brief Pointer to the type of entities that this Intersection belongs to */
  typedef typename GridImp::template Codim<0>::EntityPointer EntityPointer;

    /** \brief Codim 1 geometry returned by intersectionGlobal() */
  typedef typename GridImp::template Codim<1>::Geometry Geometry;

    /** \brief local coordinate type used as parameter for the normals */
  typedef typename Geometry::LocalCoordinate LocalCoordinate;

    /** \brief global coordinate type used as parameter for the normals */
  typedef typename Geometry::GlobalCoordinate GlobalCoordinate;

    /** \brief Codim 1 geometry returned by intersectionSelfLocal() 
        and intersectionNeighborLocal() */
  typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;

  //! @brief Export grid dimension
  enum { dimension=dim /*!< grid dimension */ };

  //! @brief Export dimension of world
  enum { dimensionworld=dimworld /*!< dimension of world */ };

  //! define type used for coordinates in grid module
    typedef typename GridImp::ctype ctype;

  //! return true if intersection is with interior or exterior boundary (see the cases above)
    bool boundary () const
    {
      return this->real.boundary();
    }

  /**
     \brief Identifier for boundary segment from macro grid.

     One can attach a boundary Id to a boundary segment on the macro
     grid. This Id will also be used for all fragments of these
     boundary segments.

     The numbering is defined as:
     - Id==0 for all intersections without boundary()==false
     - Id>=0 for all intersections without boundary()==true
     
     The way the Identifiers are attached to the grid may differ
     between the different grid implementations.

   */
  int boundaryId () const
  {
    return this->real.boundaryId();
  } 

  //! @brief return true if intersection is shared with another element.
  bool neighbor () const 
    {
      return this->real.neighbor();
    }

  /*! @brief return EntityPointer to the Entity on the inside of this
    intersection. That is the Entity where we started this .
   */
  EntityPointer inside() const
    {
      return this->real.inside();
    }

  /*! @brief return EntityPointer to the Entity on the outside of this
    intersection. That is the neighboring Entity.

    @warning Don't call this method if there is no neighboring Entity
    (neighbor() returns false). In this case the result is undefined.
   */
  EntityPointer outside() const
    {
      return this->real.outside();
    }
  
  /*! @brief return true if intersection is conform.
  
      This method returns true, if 
      @code
      inside()->entity<1>(numberInSelf()) == 
          outside()->entity<1>(numberInNeighbor()) ||
      boundary()
      @endcode
      holds.
  */ 
  bool conforming () const 
    {
      return this->real.conforming();
    }

  /** \brief geometrical information about this intersection in local
   *         coordinates of the inside() entity.
   *
   *  This method returns a Geometry object that provides a mapping from
   *  local coordinates of the intersection to local coordinates of the
   *  inside() entity.
   */
  const LocalGeometry &geometryInInside () const
  {
    return this->real.geometryInInside();
  }

  /**
     \deprecated use geometryInInside (renamed to improve consistency)
     \brief please read the details
   */
  const LocalGeometry& intersectionSelfLocal () const DUNE_DEPRECATED
  {
    return geometryInInside();
  }

  /** \brief geometrical information about this intersection in local
   *         coordinates of the outside() entity.
   *
   * This method returns a Geometry object that provides a mapping from
   * local coordinates of the intersection to local coordinates of the
   * outside() entity.
   */
  const LocalGeometry &geometryInOutside () const
  {
    return this->real.geometryInOutside();
  }

  /**
     \deprecated use geometryInOutside (renamed to improve consistency)
     \brief please read the details
   */
  const LocalGeometry& intersectionNeighborLocal () const DUNE_DEPRECATED
  {
    return geometryInOutside();
  }

  /** \brief geometrical information about the intersection in global coordinates.
   *
   *  This method returns a Geometry object that provides a mapping from
   *  local coordinates of the intersection to global (world) coordinates.
   */
  const Geometry &geometry () const
  {
    return this->real.geometry();
  }

  /**
     \deprecated use geometry (renamed to improve consistency)
     \brief please read the details
   */
  const Geometry& intersectionGlobal () const DUNE_DEPRECATED
  {
    return geometry();
  }

  /** \brief obtain the type of reference element for this intersection */
  GeometryType type () const
  {
    return this->real.type();
  }

  /** \brief Local index of codim 1 entity in the inside() entity where
   *         intersection is contained in
   *
   *  \note This method returns the face number with respect to the generic
   *        reference element.
   *
   *  \returns the index of the inside entity's face containing this
   *           intersection (with respect to the generic reference element)
   */
  int indexInInside () const
  {
    return this->real.indexInInside();
  }

#ifdef DUNE_ENABLE_OLD_NUMBERING
  /**
     \deprecated use indexInInside (renamed to improve consistency)
     \brief please read the details
     
     \warning \{
     this method uses the numbering of the old referenceelements (see GridReferenceElements)
     the output of indexInInside will differ as it uses the numbering of the generic referenceelements
     (see GridGenericReferenceElements).
     \}
     
   */
  deprecated_int numberInSelf () const DUNE_DEPRECATED
  {
    const int number = indexInInside();

    typedef GenericGeometry::MapNumberingProvider< dimension > Numbering;
    const unsigned int tid = GenericGeometry::topologyId( inside()->type() );
    return deprecated_int(Numbering::template generic2dune< 1 >( tid, number ));
  }
#endif

  /** \brief Local index of codim 1 entity in outside() entity where
   *         intersection is contained in
   *
   *  \note This method returns the face number with respect to the generic
   *        reference element.
   *
   *  \returns the index of the outside entity's face containing this
   *           intersection (with respect to the generic reference element)
   */
  int indexInOutside () const
  {
    return this->real.indexInOutside();
  }

#ifdef DUNE_ENABLE_OLD_NUMBERING
  /**
     \deprecated use indexInOutside (renamed to improve consistency)
     \brief please read the details
     
     \warning \{
     this method uses the numbering of the old referenceelements (see GridReferenceElements)
     the output of indexInOutside will differ as it uses the numbering of the generic referenceelements
     (see GridGenericReferenceElements).
     \}
     
   */
  deprecated_int numberInNeighbor () const
  {
    const int number = indexInOutside();

    typedef GenericGeometry::MapNumberingProvider< dimension > Numbering;
    const unsigned int tid = GenericGeometry::topologyId( outside()->type() );
    return deprecated_int(Numbering::template generic2dune< 1 >( tid, number ));
  }
#endif

  /*! @brief Return an outer normal (length not necessarily 1)

    The returned vector may depend on local position within the intersection.
  */
  GlobalCoordinate outerNormal (const LocalCoordinate& local) const
    {
      return this->real.outerNormal(local);
    }

  /*! @brief return outer normal scaled with the integration element
	@copydoc Dune::Intersection::outerNormal
    The normal is scaled with the integration element of the intersection. This
	method is redundant but it may be more efficent to use this function
	rather than computing the integration element via intersectionGlobal().
  */
  GlobalCoordinate integrationOuterNormal (const LocalCoordinate& local) const
    {
      return this->real.integrationOuterNormal(local);
    }

  /*! @brief Return unit outer normal (length == 1)

  The returned vector may depend on the local position within the intersection.
  It is scaled to have unit length.
  */
  GlobalCoordinate unitOuterNormal (const LocalCoordinate& local) const
    {
      return this->real.unitOuterNormal(local);
    }

  /*! @brief Return unit outer normal (length == 1)

  The returned vector is the normal at the center() of the
  intersection's geometry.
  It is scaled to have unit length.
  */
  GlobalCoordinate centerUnitOuterNormal () const
    {
      return this->real.centerUnitOuterNormal();
    }

  //===========================================================
  /** @name Implementor interface
   */
  //@{
  //===========================================================

  /** Copy Constructor from IntersectionImp */
  Intersection(const IntersectionImp<const GridImp> & i) :
    real(i) {};

  /** Copy constructor */
  Intersection(const Intersection& i) :
    real(i.real) {}
  //@}

  typedef typename remove_const<GridImp>::type mutableGridImp;
protected:
  //! give the GridDefaultImplementation class access to the realImp 
  friend class GridDefaultImplementation<
            GridImp::dimension, GridImp::dimensionworld,
            typename GridImp::ctype,
            typename GridImp::GridFamily> ;

  //! give the pseudo IntersectionIterator class access to the realImp
  //! \todo cleanup this hack
  friend class IntersectionIterator<GridImp, IntersectionImp, IntersectionImp>;
 
  //! return reference to the real implementation 
  ImplementationType & getRealImp() { return real; }
  //! return reference to the real implementation 
  const ImplementationType & getRealImp() const { return real; }

};

//**********************************************************************
/**
   @brief Default Implementations of integrationOuterNormal and unitOuterNormal for IntersectionImp

   @ingroup GridDevel
*/
template<class GridImp, template<class> class IntersectionImp>
class IntersectionDefaultNormalVectors
{
  enum { dim=GridImp::dimension };
  enum { dimworld=GridImp::dimensionworld };
  typedef typename GridImp::ctype ct;
public:
  //! return unit outer normal, this should be dependent on
  //! local coordinates for higher order boundary
  //! the normal is scaled with the integration element of the intersection.
  FieldVector<ct, dimworld> integrationOuterNormal (const FieldVector<ct, dim-1>& local) const
    {
        FieldVector<ct, dimworld> n = asImp().unitOuterNormal(local);
        n *= asImp().intersectionGlobal().integrationElement(local);
        return n;
    }
  //! return unit outer normal
  FieldVector<ct, dimworld> unitOuterNormal (const FieldVector<ct, dim-1>& local) const
    {
      FieldVector<ct, dimworld> n = asImp().outerNormal(local);
      n /= n.two_norm();
      return n;
    }
  //! return unit outer normal at center of intersection geometry
  FieldVector<ct, dimworld> centerUnitOuterNormal () const
    {
        // For now, we do this...
        GeometryType type = asImp().geometry().type();
        const GenericReferenceElement<ct, dim-1> & refElement =
            GenericReferenceElements<ct, dim-1>::general(type);
        return asImp().unitOuterNormal(refElement.position(0,0));
        // But later, if we change the meaning of center(), we may have to change to this...
        // return asImp().unitOuterNormal(asImp().local(asImp().center()));
    }

private:
  //!  Barton-Nackman trick
  IntersectionImp<GridImp>& asImp ()
    {return static_cast<IntersectionImp<GridImp>&>(*this);}
  const IntersectionImp<GridImp>& asImp () const
    {return static_cast<const IntersectionImp<GridImp>&>(*this);}
};

}

#endif // DUNE_GRID_INTERSECTION_HH
