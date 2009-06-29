#ifndef DUNE_QUADRATURERULES_HH
#define DUNE_QUADRATURERULES_HH

#include<iostream>
#include<vector>
#include<map>

#include<dune/common/fvector.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/stdstreams.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/referenceelements.hh>

/**
   \file
   Interface for quadrature points and rules
*/

namespace Dune {

  /** \brief Exception thrown if an desired QuadratureRule is not available,
     because the order is to high
     \ingroup Quadrature
  */
  class QuadratureOrderOutOfRange : public NotImplemented {};
  
  /** \brief Single evaluation point in a quadrature rule
      \ingroup Quadrature
  */
  template<typename ct, int dim>
  class QuadraturePoint {
  public:
	// compile time parameters
	enum { d=dim };
	typedef ct CoordType;

	//! set up quadrature of given order in d dimensions
	QuadraturePoint (const FieldVector<ct, dim>& x, double w) : local(x)
      {
        wght = w;
      }

	//! return local coordinates of integration point i 
	const FieldVector<ct, dim>& position () const
      {
        return local;
      }

	//! return weight associated with integration point i
	double weight () const
      {
        return wght;
      }
    virtual ~QuadraturePoint(){}
    
  protected:
	FieldVector<ct, dim> local;
	double wght;
  };

  /** \brief Defines an \p enum for currently available quadrature rules.
      \ingroup Quadrature
  */
  namespace QuadratureType {
    enum Enum {
      Gauss      = 0,
      
      Jacobian_1_0 = 1,
      Jacobian_2_0 = 2,
      
      Simpson    = 3,
      Trap       = 4,
      Grid       = 5,
      
      Clough     = 21,
      
      Invalid_Rule = 127
    };
  }

  /** \brief Abstract base class for quadrature rules
      \ingroup Quadrature
  */
  template<typename ct, int dim>
  class QuadratureRule : public std::vector<QuadraturePoint<ct,dim> >
  {
  public:

      /** \brief Default constructor */
    QuadratureRule() : delivered_order(-1) {}

      /** \brief Constructor for a given geometry type.  Leaves the quadrature order invalid  */
    QuadratureRule(GeometryType t) : geometry_type(t), delivered_order(-1) {}

      /** \brief Constructor for a given geometry type and a given quadrature order */
    QuadratureRule(GeometryType t, int order) : geometry_type(t), delivered_order(order) {}

      /** \brief The space dimension */
	enum { d=dim };

      /** \brief The type used for coordinates */
	typedef ct CoordType;

	//! return order
	virtual int order () const { return delivered_order; }

	//! return type of element
	virtual GeometryType type () const { return geometry_type; }
    virtual ~QuadratureRule(){}

    //! this container is always a const container,
    //! therefore iterator is the same as const_iterator
    typedef typename std::vector<QuadraturePoint<ct,dim> >::const_iterator iterator;

  protected:
    GeometryType geometry_type;
	int delivered_order;
    
    void tensor_product (const QuadratureRule<ct,1> & q1d)
      {
        // fill in all the gauss points
        int m = q1d.size();
        int n = power(m,dim);
        for (int i=0; i<n; i++)
		{
		  // compute multi index for Gauss point
		  int x[dim];
		  int z = i;
		  for (int k=0; k<dim; ++k)
          {
            x[k] = z%m;
            z = z/m;
          }

		  // compute coordinates and weight
		  double weight = 1.0;
		  FieldVector<ct, dim> local;
		  for (int k=0; k<dim; k++) 
          {
            local[k] = q1d[x[k]].position()[0];
            weight *= q1d[x[k]].weight();
          }

		  // put in container
		  push_back(QuadraturePoint<ct,dim>(local,weight));
		}
      }

	int power (int y, int d)
      {
        int m=1;
        for (int i=0; i<d; i++) m *= y;
        return m;
      }
  };

  // Forward declaration of the factory class,
  // needed internally by the QuadratureRules container class.
  template<typename ctype, int dim> struct QuadratureRuleFactory;

  /** \brief A container for all quadrature rules of dimension <tt>dim</tt>
      \ingroup Quadrature
  */
  template<typename ctype, int dim>
  class QuadratureRules {

    /** \todo Please doc me! */
    typedef std::pair<GeometryType,int> QuadratureRuleKey;

    /** \todo Please doc me! */
    typedef Dune::QuadratureRule<ctype, dim> QuadratureRule;

    //! real rule creator
    const QuadratureRule& _rule(const GeometryType& t, int p, QuadratureType::Enum qt=QuadratureType::Gauss)
      {
        static std::map<QuadratureRuleKey, QuadratureRule> _quadratureMap;
        QuadratureRuleKey key(t,p);
        if (_quadratureMap.find(key) == _quadratureMap.end()) {
          /*
            The rule must be acquired before we can store it.
            If we write this in one command, an invalid rule
            would get stored in case of an exception.
          */
          QuadratureRule rule =
            QuadratureRuleFactory<ctype,dim>::rule(t,p,qt);
          _quadratureMap[key] = rule;
        }
        return _quadratureMap[key];
      }
    //! singleton provider
    static QuadratureRules& instance()
      {
        static QuadratureRules instance;
        return instance;
      }
    //! private constructor
    QuadratureRules () {};
  public:
	//! select the appropriate QuadratureRule for GeometryType t and order p
    static const QuadratureRule& rule(const GeometryType& t, int p, QuadratureType::Enum qt=QuadratureType::Gauss)
      {
        return instance()._rule(t,p,qt);
      }
	//! @copydoc rule
    static const QuadratureRule& rule(const GeometryType::BasicType t, int p, QuadratureType::Enum qt=QuadratureType::Gauss)
      {
        GeometryType gt(t,dim);
        return instance()._rule(gt,p,qt);
      }
  };

  /***********************************************************/

  /** \brief Gauss quadrature for the n-dimensional cube
      \ingroup Quadrature
  */
  template<typename ct, int dim>
  class CubeQuadratureRule :
    public QuadratureRule<ct,dim>
  {
  public:
      /** \brief Space dimension */
	enum { d=dim };

      /** \brief The highest possible quadrature order */
	enum { highest_order=CubeQuadratureRule<ct,1>::highest_order };

      /** \brief The type used for coordinates */
	typedef ct CoordType;

      /** \todo Please doc me! */
	typedef CubeQuadratureRule value_type;

    ~CubeQuadratureRule(){}
  private:
    friend class QuadratureRuleFactory<ct,dim>;
	//! set up quadrature of given order in d dimensions
	CubeQuadratureRule (int p) : QuadratureRule<ct,dim>(GeometryType(GeometryType::cube, d))
      {
        QuadratureRule<ct,1> q1D = QuadratureRules<ct,1>::rule(GeometryType::cube, p);
        tensor_product( q1D );
        this->delivered_order = q1D.order();
      }

  };

  //! @copydoc CubeQuadratureRule
  //! Specialization for 0D.
  template<typename ct>
  class CubeQuadratureRule<ct,0> :
    public QuadratureRule<ct,0>
  {
  public:
	// compile time parameters
	enum { d=0 };
	enum { dim=0 };
	enum { highest_order=61 };
	typedef ct CoordType;
	typedef CubeQuadratureRule value_type;

    ~CubeQuadratureRule(){}
  private:
    friend class QuadratureRuleFactory<ct,dim>;
    CubeQuadratureRule (int p):
      QuadratureRule<ct,0>(GeometryType(GeometryType::cube, 0))
    {
      FieldVector<ct, dim> point(0.0);

      if (p > highest_order)
        DUNE_THROW(QuadratureOrderOutOfRange, "Quadrature rule " << p << " not supported!");

      this->delivered_order = highest_order;
      this->push_back(QuadraturePoint<ct,dim>(point, 1.0));
    }
  };
 
  //! @copydoc CubeQuadratureRule
  //! Specialization for 1D.
  template<typename ct>
  class CubeQuadratureRule<ct,1> :
    public QuadratureRule<ct,1>
  {
  public:
	// compile time parameters
	enum { d=1 };
	enum { dim=1 };
	enum { highest_order=61 };
	typedef ct CoordType;
	typedef CubeQuadratureRule value_type;

    ~CubeQuadratureRule(){}
  private:
    void init(int p,
	      std::vector< FieldVector<ct, dim> > & _points,
	      std::vector< double > & _weight,
	      int & delivered_order);
    friend class QuadratureRuleFactory<ct,dim>;
    CubeQuadratureRule (int p)
      : QuadratureRule<ct,1>(GeometryType(GeometryType::cube, 1))
    {
      //! set up quadrature of given order in d dimensions
      std::vector< FieldVector<ct, dim> > _points;
      std::vector< double > _weight;
      
      init(p, _points, _weight, this->delivered_order);
      
      assert(_points.size() == _weight.size());
      for (size_t i = 0; i < _points.size(); i++)
        this->push_back(QuadraturePoint<ct,dim>(_points[i], _weight[i]));
    }
  };
 
    /** \brief Jacobi-Gauss quadrature for alpha=1, beta=0
        \ingroup Quadrature
    */
  template<typename ct, int dim>
  class Jacobi1QuadratureRule;
  
    /** \brief Jacobi-Gauss quadrature for alpha=1, beta=0
        \ingroup Quadrature
    */
  template<typename ct>
  class Jacobi1QuadratureRule<ct,1> :
    public QuadratureRule<ct,1>
  {
  public:
      /** \brief The space dimension */
	enum { d=1 };
      /** \brief The space dimension
          \todo Do we have it twice? */
	enum { dim=1 };

      /** \brief The highest quadrature order available */
	enum { highest_order=61 };

      /** \brief The type used for coordinates */
	typedef ct CoordType;

      /** \todo Please doc me! */
	typedef Jacobi1QuadratureRule value_type;

    ~Jacobi1QuadratureRule(){}
  private:
    void init(int p,
	      std::vector< FieldVector<ct, dim> > & _points,
	      std::vector< double > & _weight,
	      int & delivered_order);
    friend class QuadratureRuleFactory<ct,dim>;
    Jacobi1QuadratureRule (int p)
      : QuadratureRule<ct,1>(GeometryType(GeometryType::cube, 1))
    {
      //! set up quadrature of given order in d dimensions
      std::vector< FieldVector<ct, dim> > _points;
      std::vector< double > _weight;
      
      int delivered_order;
      
      init(p, _points, _weight, delivered_order);
      this->delivered_order = delivered_order;
      assert(_points.size() == _weight.size());
      for (size_t i = 0; i < _points.size(); i++)
        this->push_back(QuadraturePoint<ct,dim>(_points[i], _weight[i]));
    }
  };
 
    /** \brief Jacobi-Gauss quadrature for alpha=2, beta=0
        \ingroup Quadrature
    */
  template<typename ct, int dim>
  class Jacobi2QuadratureRule;
  
    /** \brief Jacobi-Gauss quadrature for alpha=2, beta=0
        \ingroup Quadrature
    */
  template<typename ct>
  class Jacobi2QuadratureRule<ct,1> :
    public QuadratureRule<ct,1>
  {
  public:
      /** \brief The space dimension */
	enum { d=1 };

      /** \brief The space dimension
          \todo Do we have it twice? */
	enum { dim=1 };

      /** \brief The highest quadrature order available */
	enum { highest_order=61 };

      /** \brief The type used for coordinates */
	typedef ct CoordType;

      /** \todo Please doc me! */
	typedef Jacobi2QuadratureRule value_type;

    ~Jacobi2QuadratureRule(){}
  private:
    void init(int p,
	      std::vector< FieldVector<ct, dim> > & _points,
	      std::vector< double > & _weight,
	      int & delivered_order);
    friend class QuadratureRuleFactory<ct,dim>;
    Jacobi2QuadratureRule (int p)
      : QuadratureRule<ct,1>(GeometryType(GeometryType::cube, 1))
    {
      //! set up quadrature of given order in d dimensions
      std::vector< FieldVector<ct, dim> > _points;
      std::vector< double > _weight;
      
      int delivered_order;
      
      init(p, _points, _weight, delivered_order);
      
      this->delivered_order = delivered_order;
      assert(_points.size() == _weight.size());
      for (size_t i = 0; i < _points.size(); i++)
        this->push_back(QuadraturePoint<ct,dim>(_points[i], _weight[i]));
    }
  };
 
  /************************************************
   * Quadraturerule for Simplices/Triangle
   *************************************************/

    /** \brief Quadrature rules for simplices
        \ingroup Quadrature
    */
  template<typename ct, int dim>
  class SimplexQuadratureRule;

  /** \brief Quadrature rules for triangles
      \ingroup Quadrature
  */
  template<typename ct>
  class SimplexQuadratureRule<ct,2> : public QuadratureRule<ct,2>
  {
  public:

    /** \brief The space dimension */
    enum{d=2};

    /** \brief The highest quadrature order available */
    enum { highest_order=CubeQuadratureRule<ct,1>::highest_order -1 };

    /** \brief The type used for coordinates */
    typedef ct CoordType;

    /** export my type */
    typedef SimplexQuadratureRule value_type;
    ~SimplexQuadratureRule(){}
  private:
    friend class QuadratureRuleFactory<ct,d>;
    SimplexQuadratureRule (int p);
  };

  /** \brief Quadrature rules for tetrahedra
      \ingroup Quadrature
  */
  template<typename ct>
  class SimplexQuadratureRule<ct,3> : public QuadratureRule<ct,3>
  {
  public:

      /** \brief The space dimension */
    enum{d=3};

    /** \brief The highest quadrature order available */
    enum { highest_order=CubeQuadratureRule<ct,1>::highest_order -2 };

      /** \brief The type used for coordinates */
    typedef ct CoordType;

      /** \todo Please doc me! */
    typedef SimplexQuadratureRule<ct,3> value_type;
    ~SimplexQuadratureRule(){}
  private:
    friend class QuadratureRuleFactory<ct,d>;
    SimplexQuadratureRule (int p);
  };

/***********************************
 * quadrature for Prism
 **********************************/

    /** \todo Please doc me! */
  template<int dim>
  class PrismQuadraturePoints;

    /** \todo Please doc me! */
  template<>
  class PrismQuadraturePoints<3>
  {
  public:
	enum { MAXP=6};
	enum { highest_order=2 };

	//! initialize quadrature points on the interval for all orders
	PrismQuadraturePoints ()
      {
        int m = 0;
        O[m] = 0;
          
        // polynom degree 0  ???
        m = 6;
        G[m][0][0] = 0.0;
        G[m][0][1] = 0.0;
        G[m][0][2] = 0.0;

        G[m][1][0] = 1.0;
        G[m][1][1] = 0.0;
        G[m][1][2] = 0.0;

        G[m][2][0] = 0.0;
        G[m][2][1] = 1.0;
        G[m][2][2] = 0.0;

        G[m][3][0] = 0.0;
        G[m][3][1] = 0.0;
        G[m][3][2] = 1.0;
          
        G[m][4][0] = 1.0;
        G[m][4][1] = 0.0;
        G[m][4][2] = 1.0;

        G[m][5][0] = 0.0;
        G[m][5][1] = 0.1;
        G[m][5][2] = 1.0;

        W[m][0] = 0.16666666666666666 / 2.0;
        W[m][1] = 0.16666666666666666 / 2.0;
        W[m][2] = 0.16666666666666666 / 2.0;
        W[m][3] = 0.16666666666666666 / 2.0;
        W[m][4] = 0.16666666666666666 / 2.0;
        W[m][5] = 0.16666666666666666 / 2.0;
	  
        O[m] = 0;// verify ????????
	  

        // polynom degree 2  ???
        m = 6;
        G[m][0][0] =0.66666666666666666 ;
        G[m][0][1] =0.16666666666666666 ;
        G[m][0][2] =0.211324865405187 ;

        G[m][1][0] = 0.16666666666666666;
        G[m][1][1] =0.66666666666666666 ;
        G[m][1][2] = 0.211324865405187;

        G[m][2][0] = 0.16666666666666666;
        G[m][2][1] = 0.16666666666666666;
        G[m][2][2] = 0.211324865405187;

        G[m][3][0] = 0.66666666666666666;
        G[m][3][1] = 0.16666666666666666;
        G[m][3][2] = 0.788675134594813;
          
        G[m][4][0] = 0.16666666666666666;
        G[m][4][1] = 0.66666666666666666;
        G[m][4][2] = 0.788675134594813;

        G[m][5][0] = 0.16666666666666666;
        G[m][5][1] = 0.16666666666666666;
        G[m][5][2] = 0.788675134594813;

        W[m][0] = 0.16666666666666666 / 2.0;
        W[m][1] = 0.16666666666666666 / 2.0;
        W[m][2] = 0.16666666666666666 / 2.0;
        W[m][3] = 0.16666666666666666 / 2.0;
        W[m][4] = 0.16666666666666666 / 2.0;
        W[m][5] = 0.16666666666666666 / 2.0;
	  
        O[m] = 2;// verify ????????
	 
      }

      /** \todo Please doc me! */
      FieldVector<double, 3> point(int m, int i)
      {
        return G[m][i];
      }

      /** \todo Please doc me! */
	double weight (int m, int i)
      {
        return W[m][i];
      }

      /** \todo Please doc me! */
	int order (int m)
      {
        return O[m];
      }

  private:
    FieldVector<double, 3> G[MAXP+1][MAXP];//positions
    
	double W[MAXP+1][MAXP]; // weights associated with points       
	int    O[MAXP+1];       // order of the rule
  };


  /** \brief Singleton holding the Prism Quadrature points 
   \ingroup Quadrature
  */
  template<int dim>  
  struct PrismQuadraturePointsSingleton {
	static PrismQuadraturePoints<3> prqp;
  };

    /** \brief Singleton holding the Prism Quadrature points 
   \ingroup Quadrature
  */
  template<>  
  struct PrismQuadraturePointsSingleton<3> {
	static PrismQuadraturePoints<3> prqp;
  };

    /** \brief Quadrature rules for prisms
        \ingroup Quadrature
    */
  template<typename ct, int dim>
  class PrismQuadratureRule;

    /** \brief Quadrature rules for prisms
        \ingroup Quadrature
    */
  template<typename ct>
  class PrismQuadratureRule<ct,3> : public QuadratureRule<ct,3>
  {
  public:

      /** \brief The space dimension */
    enum{ d=3 };

      /** \brief The highest quadrature order available */
    enum{
      /* min(Line::order, Triangle::order) */
      highest_order =
        (int)CubeQuadratureRule<ct,1>::highest_order
        < (int)SimplexQuadratureRule<ct,2>::highest_order
        ? (int)CubeQuadratureRule<ct,1>::highest_order
        : (int)SimplexQuadratureRule<ct,2>::highest_order
        };

      /** \brief The type used for coordinates */
    typedef ct CoordType;

      /** \todo Please doc me! */
    typedef PrismQuadratureRule<ct,3> value_type;

    ~PrismQuadratureRule(){}
  private:
    friend class QuadratureRuleFactory<ct,d>;
    PrismQuadratureRule(int p) : QuadratureRule<ct,3>(GeometryType(GeometryType::prism, d))
      {
        if (p>highest_order)
          DUNE_THROW(QuadratureOrderOutOfRange,
                     "QuadratureRule for order " << p << " and GeometryType "
                     << this->type() << " not available");
          
        if (p<=2) {
          int m=6;
          this->delivered_order = PrismQuadraturePointsSingleton<3>::prqp.order(m);
          for(int i=0;i<m;++i)
          { 
            FieldVector<ct,3> local;
            for (int k=0; k<d; k++)
              local[k] = PrismQuadraturePointsSingleton<3>::prqp.point(m,i)[k];
            double weight =
              PrismQuadraturePointsSingleton<3>::prqp.weight(m,i);
            // put in container
            push_back(QuadraturePoint<ct,d>(local,weight));
          }
        }
        else {
          const QuadratureRule<ct,2> & triangle = QuadratureRules<ct,2>::rule(GeometryType::simplex, p);
          const QuadratureRule<ct,1> & line = QuadratureRules<ct,1>::rule(GeometryType::cube, p);
          
          this->delivered_order = std::min(triangle.order(),line.order());

          for (typename QuadratureRule<ct,1>::const_iterator
                 lit = line.begin(); lit != line.end(); ++lit)
          {
            for (typename QuadratureRule<ct,2>::const_iterator
                   tit = triangle.begin(); tit != triangle.end(); ++tit)
            {
              FieldVector<ct, d> local;
              local[0] = tit->position()[0];
              local[1] = tit->position()[1];
              local[2] = lit->position()[0];
              
              double weight = tit->weight() * lit->weight();

              // put in container
              push_back(QuadraturePoint<ct,d>(local,weight));
            }
          }
        }
      }
  };

    /** \todo Please doc me! */
  class PyramidQuadraturePoints
  {
  public:
	enum { MAXP=8};
	enum { highest_order=2 };

	//! initialize quadrature points on the interval for all orders
	PyramidQuadraturePoints()
      {
        int m = 0;
        O[m] = 0;
         

        // polynom degree 2  ???
        m = 8;
        G[m][0][0] =0.58541020;
        G[m][0][1] =0.72819660;
        G[m][0][2] =0.13819660;

        G[m][1][0] =0.13819660;
        G[m][1][1] =0.72819660;
        G[m][1][2] =0.13819660;

        G[m][2][0] =0.13819660;
        G[m][2][1] =0.27630920;
        G[m][2][2] =0.58541020;

        G[m][3][0] =0.13819660;
        G[m][3][1] =0.27630920;
        G[m][3][2] =0.13819660;
          
        G[m][4][0] =0.72819660;
        G[m][4][1] =0.13819660;
        G[m][4][2] =0.13819660;

        G[m][5][0] =0.72819660;
        G[m][5][1] =0.58541020;
        G[m][5][2] =0.13819660;

        G[m][6][0] =0.27630920;
        G[m][6][1] =0.13819660;
        G[m][6][2] =0.58541020;

        G[m][7][0] =0.27630920;
        G[m][7][1] =0.13819660;
        G[m][7][2] =0.13819660;

        W[m][0] = 0.125 / 3.0;
        W[m][1] = 0.125 / 3.0;
        W[m][2] = 0.125 / 3.0;
        W[m][3] = 0.125 / 3.0;
        W[m][4] = 0.125 / 3.0;
        W[m][5] = 0.125 / 3.0;
        W[m][6] = 0.125 / 3.0;
        W[m][7] = 0.125 / 3.0;
	  
        O[m] = 2;// verify ????????
	 
      }

    /** \todo Please doc me! */
    FieldVector<double, 3> point(int m, int i)
      {
        return G[m][i];
      }

    /** \todo Please doc me! */
	double weight (int m, int i)
      {
        return W[m][i];
      }

    /** \todo Please doc me! */
	int order (int m)
      {
        return O[m];
      }

  private:
    FieldVector<double, 3> G[MAXP+1][MAXP];
	double W[MAXP+1][MAXP]; // weights associated with points       
	int    O[MAXP+1];       // order of the rule
  };

    /** \brief Singleton holding the pyramid quadrature points 
        \ingroup Quadrature
    */
  template<int dim>
  struct PyramidQuadraturePointsSingleton {};

    /** \brief Singleton holding the pyramid quadrature points
        \ingroup Quadrature
    */
  template<>
  struct PyramidQuadraturePointsSingleton<3> {
	static PyramidQuadraturePoints pyqp;
  };

    /** \brief Quadrature rules for Pyramids
        \ingroup Quadrature
    */
  template<typename ct, int dim>
  class PyramidQuadratureRule; 

    /** \brief Quadrature rules for Pyramids
        \ingroup Quadrature
    */
  template<typename ct>
  class PyramidQuadratureRule<ct,3> : public QuadratureRule<ct,3>
  {
  public:

      /** \brief The space dimension */
    enum{d=3};

      /** \brief The highest quadrature order available */
    enum{highest_order=CubeQuadratureRule<ct,1>::highest_order};

      /** \brief The type used for coordinates */
    typedef ct CoordType;
      
      /** \todo Please doc me! */
    typedef PyramidQuadratureRule<ct,3> value_type;

    ~PyramidQuadratureRule(){}
  private:
    friend class QuadratureRuleFactory<ct,d>;
    PyramidQuadratureRule(int p) : QuadratureRule<ct,3>(GeometryType(GeometryType::pyramid, d))
      {
        int m;
        
        if (p>highest_order)
          DUNE_THROW(QuadratureOrderOutOfRange,
                     "QuadratureRule for order " << p << " and GeometryType "
                     << this->type() << " not available");
          
        if(false) {
//        if(p<=2) {
          m=8;
          this->delivered_order = PyramidQuadraturePointsSingleton<3>::pyqp.order(m);
          FieldVector<ct, d> local;
          double weight;
          for(int i=0;i<m;++i)
          { 
            for(int k=0;k<d;++k)
              local[k]=PyramidQuadraturePointsSingleton<3>::pyqp.point(m,i)[k];
            weight=PyramidQuadraturePointsSingleton<3>::pyqp.weight(m,i);
            // put in container
            push_back(QuadraturePoint<ct,d>(local,weight));	
          }
        }
        else
        {
          // Define the quadrature rules...
          QuadratureRule<ct,3> simplex =
            QuadratureRules<ct,3>::rule(GeometryType::simplex,p);

          for (typename QuadratureRule<ct,3>::const_iterator
                 it=simplex.begin(); it != simplex.end(); ++it)
          {
            FieldVector<ct,3> local = it->position();
            ct weight = it->weight();
            // Simplex 1
            // x:=x+y
            local[0] = local[0]+local[1];
            push_back(QuadraturePoint<ct,d>(local,weight));
            // Simplex 2
            // y:=x+y
            local[0] = it->position()[0];
            local[1] = local[0]+local[1];
            push_back(QuadraturePoint<ct,d>(local,weight));
          }

          this->delivered_order = simplex.order();
        }
      }
  };


  //! @copydoc Quadrature rule for the singular referece element.
    template<typename ct, int dim>
  class SingularQuadratureRule :
    public QuadratureRule<ct,dim>
  {
  public:
	// compile time parameters
	enum { d=dim };
	enum { highest_order=61 };
	typedef ct CoordType;
	typedef SingularQuadratureRule value_type;

    ~SingularQuadratureRule(){}
  private:
    friend class QuadratureRuleFactory<ct,dim>;
      SingularQuadratureRule (int /*p*/)
      : QuadratureRule<ct,dim>(GeometryType(GeometryType::singular, dim))
    {
	FieldVector<ct, dim> origin(0.0);
	this->push_back(QuadraturePoint<ct,dim>(origin, 1.0));
    }
  };


    /** \brief Factory class for creation of quadrature rules,
        depending on GeometryType, order and QuadratureType.

        The whole class is private and can only be accessed
        by the singleton container class QuadratureRules.
  */
  template<typename ctype, int dim>
  class QuadratureRuleFactory {
  private:
    friend class QuadratureRules<ctype, dim>;
    static QuadratureRule<ctype, dim> rule(const GeometryType& t, int p, QuadratureType::Enum qt)
      {
        if (t.isSingular())
        {
          return SingularQuadratureRule<ctype,dim>(p);
        }
        if (t.isCube())
        {
          return CubeQuadratureRule<ctype,dim>(p);
        }
        if (t.isSimplex())
        {
          return SimplexQuadratureRule<ctype,dim>(p);
        }
        DUNE_THROW(Exception, "Unknown GeometryType");
      }
  };

  template<typename ctype>
  class QuadratureRuleFactory<ctype, 0> {
  private:
    enum { dim = 0 };
    friend class QuadratureRules<ctype, dim>;
    static QuadratureRule<ctype, dim> rule(const GeometryType& t, int p, QuadratureType::Enum qt)
      {
        if (t.isSingular())
        {
          return SingularQuadratureRule<ctype,dim>(p);
        }
        if (t.isVertex())
        {
	  return CubeQuadratureRule<ctype,dim>(p);
        }
        DUNE_THROW(Exception, "Unknown GeometryType");
      }
  };

  template<typename ctype>
  class QuadratureRuleFactory<ctype, 1> {
  private:
    enum { dim = 1 };
    friend class QuadratureRules<ctype, dim>;
    static QuadratureRule<ctype, dim> rule(const GeometryType& t, int p, QuadratureType::Enum qt)
      {
        if (t.isSingular())
        {
          return SingularQuadratureRule<ctype,dim>(p);
        }
        if (t.isLine())
        {
          switch (qt) {
          case QuadratureType::Gauss:
            return CubeQuadratureRule<ctype,dim>(p);
          case QuadratureType::Jacobian_1_0:
            return Jacobi1QuadratureRule<ctype,dim>(p);
          case QuadratureType::Jacobian_2_0:
            return Jacobi2QuadratureRule<ctype,dim>(p);
          default:
            DUNE_THROW(Exception, "Unknown QuadratureType");
          }
        }
        DUNE_THROW(Exception, "Unknown GeometryType");
      }
  };

  template<typename ctype>
  class QuadratureRuleFactory<ctype, 3> {
  private:
    enum { dim = 3 };
    friend class QuadratureRules<ctype, dim>;
    static QuadratureRule<ctype, dim> rule(const GeometryType& t, int p, QuadratureType::Enum qt)
      {
        if (t.isSingular())
        {
          return SingularQuadratureRule<ctype,dim>(p);
        }
        if (t.isCube())
        {
          return CubeQuadratureRule<ctype,dim>(p);
        }
        if (t.isSimplex())
        {
          return SimplexQuadratureRule<ctype,dim>(p);
        }
        if (t.isPrism())
	    {
          return PrismQuadratureRule<ctype,dim>(p);
	    }
        if (t.isPyramid())
	    {
          return PyramidQuadratureRule<ctype,dim>(p);
	    }
        DUNE_THROW(Exception, "Unknown GeometryType");
      }
  };

} // end namespace

#endif
