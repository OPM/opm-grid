#ifndef DUNE_GRID_GEOMETRY_HH
#define DUNE_GRID_GEOMETRY_HH

/** \file
    \brief Wrapper and interface classes for element geometries
*/

#include <cassert>

#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/genericgeometry/conversion.hh>

namespace Dune
{

// External Forward Declarations
// -----------------------------

template< int dim, int dimworld, class ct, class GridFamily >
class GridDefaultImplementation;



//*****************************************************************************
//
// Geometry
// forwards the interface to the implementation
//
//*****************************************************************************

  /** 
   @brief Wrapper class for geometries


   Template parameters are:

   - <tt>mydim</tt> Dimension of the domain
   - <tt>cdim</tt> Dimension of the range
   - <tt>GridImp</tt> Type that is a model of Dune::Grid
   - <tt>GeometryImp</tt> Class template that is a model of Dune::Geometry
 
   <H3>Maps</H3>

   A Geometry defines a map \f[ g : D \to W\f] where
   \f$D\subseteq\mathbf{R}^\textrm{mydim}\f$ and
   \f$W\subseteq\mathbf{R}^\textrm{cdim}\f$. 
   The domain \f$D\f$ is one of a set of predefined convex polytopes, the
   so-called reference elements (\see Dune::GenericReferenceElement). The dimensionality
   of \f$D\f$ is <tt>mydim</tt>. 
   In general \f$\textrm{mydim}\leq\textrm{cdim}\f$, i.e.
   the convex polytope may be mapped to a manifold. Moreover, we require that 
   \f$ g\in \left( C^1(D) \right)^\textrm{cdim}\f$ and one-to-one.


   <H3>Engine Concept</H3>

   The Geometry class template wraps an object of type GeometryImp and forwards all member 
   function calls to corresponding members of this class. In that sense Geometry
   defines the interface and GeometryImp supplies the implementation.



    \ingroup GIGeometry
  */
template<int mydim, int cdim, class GridImp, template<int,int,class> class GeometryImp>  
class Geometry {

protected:
  GeometryImp<mydim,cdim,GridImp> realGeometry;
  
  // type of underlying implementation, for internal use only 
  typedef GeometryImp<mydim,cdim,GridImp> ImplementationType;
public:
  //! @brief export grid dimension
  enum { dimension=GridImp::dimension /*!< grid dimension */ };
  //! @brief export geometry dimension
  enum { mydimension=mydim /*!< geometry dimension */ };
  //! @brief export coordinate dimension
  enum { coorddimension=cdim /*!< dimension of embedding coordsystem */ };

  //! @brief export dimension of world
  enum { dimensionworld=GridImp::dimensionworld /*!< dimension of world */ };
  //! define type used for coordinates in grid module
  typedef typename GridImp::ctype ctype;

  //! type of local coordinates
  typedef FieldVector<ctype, mydim> LocalCoordinate;

  //! type of the global coordinates
  typedef FieldVector< ctype, cdim > GlobalCoordinate;
  
  //! type of jacobian (also of jacobian inverse transposed)
  typedef FieldMatrix<ctype,cdim,mydim> Jacobian;

  //! type of jacobian transposed
  typedef FieldMatrix< ctype, mydim, cdim > JacobianTransposed;

  /** \brief Return the name of the reference element. The type can
    be used to access the Dune::GenericReferenceElement.
   */
  GeometryType type () const { return realGeometry.type(); }

  /** \brief Return true if the geometry mapping is affine and false otherwise */
  bool affine() const { return realGeometry.affine(); }

  /** \brief Return the number of corners of the reference element. Since
    this is a convex polytope the number of corners is a well-defined concept.
    The method is redundant because this information is also available
    via the reference element. It is here for efficiency and ease of use.
   */
  int corners () const { return realGeometry.corners(); }

#ifdef DUNE_ENABLE_OLD_NUMBERING
  /** \brief Access to corners of the geometry. 
    \param[in] i The number of the corner
    \return const reference to a vector containing the position \f$g(c_i)\f$ where
    \f$c_i\f$ is the position of the i'th corner of the reference element.
   */
  const GlobalCoordinate operator[] ( deprecated_int i ) const DUNE_DEPRECATED
    {
        typedef GenericGeometry::MapNumberingProvider< dimension > Numbering;
        const unsigned int tid = GenericGeometry::topologyId( type() );
        const int j = Numbering::template dune2generic< dimension >( tid, i.value() );
        return realGeometry.corner(j);
    }
#endif

  /** \brief Obtain a corner of the geometry
   *
   *  This method is for convenient access to the corners of the geometry. The
   *  same result could be achieved by by calling
   *  \code
   *  global( genericReferenceElement.position( i, mydimension )
   *  \endcode
   *
   *  \note This method replaces the deprecated operator[]. Yet, there are two
   *        differences:
   *        - operator[] returns a reference instead of an object
   *        - operator[] uses the old sub entity numbering
   *
   *  \param[in]  i  number of the corner (with respect to the generic reference
   *                 element)
   *  \returns position of the i'th corner
   */
  GlobalCoordinate corner ( int i ) const
  {
    return realGeometry.corner( i );
  }

  /** \brief Evaluate the map \f$ g\f$.
    \param[in] local Position in the reference element \f$D\f$
    \return Position in \f$W\f$
  */
  GlobalCoordinate global (const LocalCoordinate& local) const
    {
      return realGeometry.global(local);
    }

  /** \brief Evaluate the inverse map \f$ g^{-1}\f$
    \param[in] global Position in \f$W\f$
    \return Position in \f$D\f$ that maps to global
  */
  LocalCoordinate local (const GlobalCoordinate& global) const
    {
      return realGeometry.local(global);
    }

    /** \brief Checks whether a point is in the reference element \f$D\f$
     *
     *  \param[in]  local  point \f$x\f$ to verify
     *  \returns true, if \f$x\in D\f$ of.
     *
     *  \deprecated Use the corresponding method in GenericReferenceElement
     */
    bool checkInside ( const LocalCoordinate& local ) const DUNE_DEPRECATED
    {
      const GenericReferenceElement< ctype, mydim > &refElement
        = GenericReferenceElements< ctype, mydim >::general( type() );
      return refElement.checkInside( local );
    }

  /** \brief Return the factor appearing in the integral transformation formula

    Let \f$ g : D \to W\f$ denote the transformation described by the Geometry.
    Then the jacobian of the transformation is defined as the 
    \f$\textrm{cdim}\times\textrm{mydim}\f$ matrix
    \f[ J_g(x) = \left( \begin{array}{ccc} \frac{\partial g_0}{\partial x_0} &
    \cdots & \frac{\partial g_0}{\partial x_{n-1}} \\
    \vdots & \ddots & \vdots \\ \frac{\partial g_{m-1}}{\partial x_0} &
    \cdots & \frac{\partial g_{m-1}}{\partial x_{n-1}}
    \end{array} \right).\f]
    Here we abbreviated \f$m=\textrm{cdim}\f$ and \f$n=\textrm{mydim}\f$ for ease of
    readability.

    The integration element \f$\mu(x)\f$ for any \f$x\in D\f$ is then defined as
    \f[ \mu(x) = \sqrt{|\det J_g^T(x)J_g(x)|}.\f]

    \param[in] local Position \f$x\in D\f$
    \return    integration element \f$\mu(x)\f$

    \note Each implementation computes the integration element with optimal
    efficieny. For example in an equidistant structured mesh it may be as
    simple as \f$h^\textrm{mydim}\f$.
  */
  ctype integrationElement (const LocalCoordinate& local) const
  {
    return realGeometry.integrationElement(local);
  }
 
  /** \brief return volume of geometry */
  ctype volume () const
  {
    return realGeometry.volume();
  }

  /** \brief return center of geometry 
   *  Note that this method is still subject to a change of name and semantics.
   *  At the moment, the center is not required to be the centroid of the
   *  geometry, or even the centroid of its corners. This makes the current
   *  default implementation acceptable, which maps the centroid of the
   *  reference element to the geometry.
   *  We may change the name (and semantic) of the method to centroid() if we
   *  find reasonably efficient ways to implement it properly.
   */
  GlobalCoordinate center () const
  {
    return realGeometry.center();
  }

  /** \brief Return the transposed of the Jacobian
   *
   *  The Jacobian is defined in the documentation of
   *  \ref Dune::Geometry::integrationElement "integrationElement".
   *
   *  \param[in]  local  position \f$x\in D\f$
   *
   *  \return \f$J_g^T(x)\f$
   */
  const JacobianTransposed &
  jacobianTransposed ( const LocalCoordinate& local ) const
  {
    return realGeometry.jacobianTransposed( local );
  }

  /** \brief Return inverse of transposed of Jacobian
   *
   *  The Jacobian is defined in the documentation of
   *  \ref Dune::Geometry::integrationElement "integrationElement".
   *
   *  \param[in]  local  position \f$x\in D\f$
   *  \return \f$J_g^{-T}(x)\f$
   *
   *  The use of this function is to compute the gradient of some function
   *  \f$f : W \to \textbf{R}\f$ at some position \f$y=g(x)\f$, where
   *  \f$x\in D\f$ and \f$g\f$ the transformation of the Geometry.
   *  When we set \f$\hat{f}(x) = f(g(x))\f$ and apply the chain rule we obtain
   *  \f[\nabla f(g(x)) = J_g^{-T}(x) \nabla \hat{f}(x).\f]
   *
   *  \note In the non-symmetric case \f$\textrm{cdim} \neq \textrm{mydim}\f$, the
   *        pseudoinverse of \f$J_g^T(x)\f$ is returned.
   *        This means that it is inverse for all tangential vectors in
   *        \f$g(x)\f$ while mapping all normal vectors to zero.
   */
  const Jacobian& jacobianInverseTransposed (const LocalCoordinate& local) const
  {
    return realGeometry.jacobianInverseTransposed(local);
  }

public:    
  //! copy constructor from GeometryImp
  explicit Geometry(const GeometryImp<mydim,cdim,GridImp> & e) : realGeometry(e) {};
  
protected:
  // give the GridDefaultImplementation class access to the realImp 
  friend class GridDefaultImplementation< 
            GridImp::dimension, GridImp::dimensionworld,
            typename GridImp::ctype,
            typename GridImp::GridFamily> ;

  //! return reference to the real implementation 
  GeometryImp<mydim,cdim,GridImp> & getRealImp() { return realGeometry; }
  //! return reference to the real implementation 
  const GeometryImp<mydim,cdim,GridImp> & getRealImp() const { return realGeometry; }
 
protected:  
  /** hide copy constructor */
  Geometry(const Geometry& rhs) : realGeometry(rhs.realGeometry) {};
  /** hide assignment operator */
  Geometry & operator = (const Geometry& rhs) {
    realGeometry = rhs.realGeometry;
    return *this;
  };
};


//************************************************************************
// GEOMETRY Default Implementations
//*************************************************************************
//
// --GeometryDefault 
// 
//! Default implementation for class Geometry
template<int mydim, int cdim, class GridImp, template<int,int,class> class GeometryImp>  
class GeometryDefaultImplementation
{
public:
  // save typing
  typedef typename GridImp::ctype ctype;
 
  //! return volume of the geometry 
  ctype volume () const 
  {
    GeometryType type = asImp().type();

    // get corresponding reference element
    const GenericReferenceElement< ctype , mydim > & refElement =
      GenericReferenceElements< ctype, mydim >::general(type);

    FieldVector<ctype,mydim> localBaryCenter (0.0);
    // calculate local bary center 
    const int corners = refElement.size(0,0,mydim);
    for(int i=0; i<corners; ++i) localBaryCenter += refElement.position(i,mydim);
    localBaryCenter *= (ctype) (1.0/corners);

    // volume is volume of reference element times integrationElement
    return refElement.volume() * asImp().integrationElement(localBaryCenter);
  }
  
  //! return center of the geometry 
  FieldVector<ctype, cdim> center () const 
  {
    GeometryType type = asImp().type();

    // get corresponding reference element
    const GenericReferenceElement< ctype , mydim > & refElement =
      GenericReferenceElements< ctype, mydim >::general(type);

    // center is (for now) the centroid of the reference element mapped to
    // this geometry.
    return asImp().global(refElement.position(0,0));
  }
  
private:
  //!  Barton-Nackman trick 
  GeometryImp<mydim,cdim,GridImp>& asImp () {return static_cast<GeometryImp<mydim,cdim,GridImp>&>(*this);}
  const GeometryImp<mydim,cdim,GridImp>& asImp () const {return static_cast<const GeometryImp<mydim,cdim,GridImp>&>(*this);}
}; // end GeometryDefault 

template<int cdim, class GridImp, template<int,int,class> class GeometryImp>  
class GeometryDefaultImplementation<0,cdim,GridImp,GeometryImp>
{
  // my dimension 
  enum { mydim = 0 };
public:
  // save typing
  typedef typename GridImp::ctype ctype;

  //! return the only coordinate 
  FieldVector<ctype, cdim> global (const FieldVector<ctype, mydim>& local) const
  {
    return asImp()[0];
  }

  //! return empty vector 
  FieldVector<ctype, mydim> local (const FieldVector<ctype, cdim>& ) const
  {
    return FieldVector<ctype, mydim>();
  }
 
  //! checkInside here returns true  
  bool checkInside (const FieldVector<ctype, mydim>& ) const
  {
    return true; 
  }
 
  //! return volume of the geometry 
  ctype volume () const { return 1.0; }

  //! return center of the geometry 
  FieldVector<ctype, cdim> center () const 
  {
      return asImp()[0];
  }
  
private:
  //!  Barton-Nackman trick 
  GeometryImp<mydim,cdim,GridImp>& asImp () {return static_cast<GeometryImp<mydim,cdim,GridImp>&>(*this);}
  const GeometryImp<mydim,cdim,GridImp>& asImp () const {return static_cast<const GeometryImp<mydim,cdim,GridImp>&>(*this);}
}; // end GeometryDefault 

}
#endif // DUNE_GRID_GEOMETRY_HH
