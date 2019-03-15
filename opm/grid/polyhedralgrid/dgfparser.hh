// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_DGFPARSER_HH
#define DUNE_POLYHEDRALGRID_DGFPARSER_HH

#include <algorithm>
#include <numeric>

#include <dune/common/typetraits.hh>
#include <dune/common/version.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#if DUNE_VERSION_NEWER(DUNE_GRID,2,5)
#include <dune/grid/io/file/dgfparser/blocks/polyhedron.hh>
#endif

#include <opm/grid/polyhedralgrid/gridfactory.hh>

#if HAVE_OPM_PARSER
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#endif

namespace Dune
{

  namespace dgf
  {

#if ! DUNE_VERSION_NEWER(DUNE_GRID,2,5)
    namespace PolyhedralGrid
    {

      // PolygonBlock
      // ------------

      struct PolygonBlock
        : public BasicBlock
      {
        PolygonBlock ( std::istream &in, int numVtx, int vtxOfs )
          : BasicBlock( in, "Polygon" ), vtxBegin_( vtxOfs ), vtxEnd_( vtxOfs + numVtx )
        {}

        int get ( std::vector< std::vector< int > > &polygons )
        {
          reset();
          std::vector< int > polygon;
          while( getnextline() )
          {
            polygon.clear();
            for( int vtxIdx; getnextentry( vtxIdx ); )
            {
              if( (vtxBegin_ > vtxIdx) || (vtxIdx >= vtxEnd_) )
                DUNE_THROW( DGFException, "Error in " << *this << ": Invalid vertex index (" << vtxIdx << " not int [" << vtxBegin_ << ", " << vtxEnd_ << "[)" );
              polygon.push_back( vtxIdx - vtxBegin_ );
            }

            polygons.push_back( polygon );
          }
          return polygons.size();
        }

      private:
        int vtxBegin_, vtxEnd_;
      };



      // PolyhedronBlock
      // ---------------

      struct PolyhedronBlock
        : public BasicBlock
      {
        explicit PolyhedronBlock ( std::istream &in, int numPolys )
          : BasicBlock( in, "Polyhedron" ), numPolys_( numPolys )
        {}

        int get ( std::vector< std::vector< int > > &polyhedra )
        {
          reset();
          std::vector< int > polyhedron;
          int minPolyId = 1;
          while( getnextline() )
          {
            polyhedron.clear();
            for( int polyIdx; getnextentry( polyIdx ); )
            {
              if( (polyIdx < 0) || (polyIdx > numPolys_) )
                DUNE_THROW( DGFException, "Error in " << *this << ": Invalid polygon index (" << polyIdx << " not int [0, " << numPolys_ << "])" );

              minPolyId = std::min( minPolyId, polyIdx );
              polyhedron.push_back( polyIdx );
            }

            polyhedra.push_back( polyhedron );
          }

          // substract minimal number to have 0 starting numbering
          if( minPolyId > 0 )
          {
            const size_t polySize = polyhedra.size();
            for( size_t i=0; i<polySize; ++i )
            {
              const size_t pSize = polyhedra[ i ].size();
              for( size_t j=0; j<pSize; ++j )
              {
                polyhedra[ i ][ j ] -= minPolyId;
              }
            }
          }
          return polyhedra.size();
        }

      private:
        const int numPolys_;
      };

    } // namespace PolyhedralGrid

    using PolyhedralGrid :: PolygonBlock;
    using PolyhedralGrid :: PolyhedronBlock;
#endif

  } // namespace dgf



  // DGFGridFactory for PolyhedralGrid
  // ---------------------------------

  template< int dim, int dimworld >
  struct DGFGridFactory< PolyhedralGrid< dim, dimworld > >
  {
    typedef PolyhedralGrid< dim, dimworld > Grid;

    const static int dimension = Grid::dimension;
    typedef MPIHelper::MPICommunicator MPICommunicator;
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<dimension>::Entity Vertex;

    explicit DGFGridFactory ( std::istream &input, MPICommunicator comm = MPIHelper::getCommunicator() )
      : grid_()
    {
      input.clear();
      input.seekg( 0 );
      if( !input )
        DUNE_THROW( DGFException, "Error resetting input stream" );
      generate( input );
    }

    explicit DGFGridFactory ( const std::string &filename, MPICommunicator comm = MPIHelper::getCommunicator() )
      : grid_()
    {
      std::ifstream input( filename );
      if( !input )
        DUNE_THROW( DGFException, "Macrofile '" << filename << "' not found" );

#if HAVE_OPM_PARSER
      if( !DuneGridFormatParser::isDuneGridFormat( input ) )
      {
        Opm::Parser parser;
        Opm::ParseContext parseContext;
        const auto deck = parser.parseFile(filename, parseContext);
        std::vector<double> porv;

        grid_.reset( new Grid( deck, porv ) );
        return ;
      }
      else
#endif
      {
        generate( input );
      }
    }

    Grid *grid () const { return grid_.operator->(); }

    template< class Intersection >
    bool wasInserted ( const Intersection &intersection ) const
    {
      return false;
    }

    template< class Intersection >
    int boundaryId ( const Intersection &intersection ) const
    {
      return false;
    }

    bool haveBoundaryParameters () const { return false; }

    template< int codim >
    int numParameters () const
    {
      //return (codim == dimension ? numVtxParams_ : 0);;
      return 0;
    }

    template< class Intersection >
    const typename DGFBoundaryParameter::type &
    boundaryParameter ( const Intersection &intersection ) const
    {
      return DGFBoundaryParameter::defaultValue();;
    }

    template< class Entity >
    std::vector< double > &parameter ( const Entity &entity )
    {
      static std::vector< double > dummy;
      return dummy;
    }

  private:
    int readVertices ( std::istream &input, std::vector< std::vector< double > > &vertices )
    {
      int dimWorld = 3;
      dgf::VertexBlock vtxBlock( input, dimWorld );
      if( !vtxBlock.isactive() )
        DUNE_THROW( DGFException, "Vertex block not found" );

      vtxBlock.get( vertices, vtxParams_, numVtxParams_ );
      return vtxBlock.offset();
    }

    std::vector< std::vector< int > > readPolygons ( std::istream &input, int numVtx, int vtxOfs )
    {
      dgf::PolygonBlock polygonBlock( input, numVtx, vtxOfs );
      if( !polygonBlock.isactive() )
        DUNE_THROW( DGFException, "Polygon block not found" );

      std::vector< std::vector< int > > polygons;
      polygonBlock.get( polygons );
      return polygons;
    }

    std::vector< std::vector< int > > readPolyhedra ( std::istream &input, int numPolygons )
    {
      dgf::PolyhedronBlock polyhedronBlock( input, numPolygons );
      std::vector< std::vector< int > > polyhedra;
      if( polyhedronBlock.isactive() )
      {
        polyhedronBlock.get( polyhedra );
      }
      return polyhedra;
    }

    template< class Iterator >
    void copy ( Iterator begin, Iterator end, double *dest )
    {
      for( ; begin != end; ++begin )
        dest = std::copy( begin->begin(), begin->end(), dest );
    }

    template< class Iterator >
    void copy ( Iterator begin, Iterator end, int *dest, int *offset )
    {
      int size = 0;
      for( ; begin != end; ++begin )
      {
        *(offset++) = size;
        size += begin->size();
        dest = std::copy( begin->begin(), begin->end(), dest );
      }
      *offset = size;
    }

    void generate ( std::istream &input )
    {
      if( !DuneGridFormatParser::isDuneGridFormat( input ) )
      {
        DUNE_THROW( DGFException, "Not in DGF format" );
      }

      typedef std::vector< std::vector< double > > CoordinateVectorType;
      CoordinateVectorType nodes;
      const int vtxOfs = readVertices( input, nodes );

      typedef std::vector< std::vector< int > > IndexVectorType;
      IndexVectorType faces = readPolygons ( input, nodes.size(), vtxOfs );
      IndexVectorType cells = readPolyhedra( input, faces.size() );

      if( cells.empty() )
      {
        DUNE_THROW( DGFException, "Polyhedron block not found" );
      }

      typedef GridFactory< Grid > GridFactoryType;
      typedef typename GridFactoryType :: Coordinate Coordinate ;

      GridFactoryType gridFactory;

      const int nNodes = nodes.size();
      Coordinate node( 0 );
      for( int i=0; i<nNodes; ++i )
      {
        for( int d=0; d<Coordinate::dimension; ++d )
          node[ d ] = nodes[ i ][ d ];

        gridFactory.insertVertex( node );
      }
      //nodes.swap( CoordinateVectorType() );

      // insert faces with type none/dim-1
      GeometryType type;
      type.makeNone( Grid::dimension - 1 );
      std::vector< unsigned int > numbers;

      const int nFaces = faces.size();
      for(int i = 0; i < nFaces; ++ i )
      {
        // copy values into appropriate data type
        std::vector<int>& face = faces[ i ];
        numbers.resize( face.size() );
        std::copy( face.begin(), face.end(), numbers.begin() );
        gridFactory.insertElement( type, numbers );
      }

      // free memory
      //faces.swap( IndexVectorType() );

      // insert cells with type none/dim
      type.makeNone( Grid::dimension );

      const int nCells = cells.size();
      for(int i = 0; i < nCells; ++ i )
      {
        // copy values into appropriate data type
        std::vector<int>& cell = cells[ i ];
        numbers.resize( cell.size() );
        std::copy( cell.begin(), cell.end(), numbers.begin() );
        gridFactory.insertElement( type, numbers );
      }
      //cells.swap( IndexVectorType() );

      grid_ = gridFactory.createGrid();
    }

    std::unique_ptr< Grid > grid_;
    int numVtxParams_;
    std::vector< std::vector< double > > vtxParams_;
  };



  // DGFGridInfo for PolyhedralGrid
  // ------------------------------

  template< int dim, int dimworld >
  struct DGFGridInfo< PolyhedralGrid< dim, dimworld > >
  {
    static int refineStepsForHalf ()
    {
      return 0;
    }

    static double refineWeight ()
    {
      return 0;
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_DGFPARSER_HH
