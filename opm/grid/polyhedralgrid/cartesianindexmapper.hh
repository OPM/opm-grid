#ifndef OPM_POLYHEDRALCARTESIANINDEXMAPPER_HEADER
#define OPM_POLYHEDRALCARTESIANINDEXMAPPER_HEADER

#include <opm/grid/common/CartesianIndexMapper.hpp>
#include <opm/grid/polyhedralgrid.hh>

namespace Dune
{
    template< int dim, int dimworld, typename coord_t >
    class CartesianIndexMapper< PolyhedralGrid< dim, dimworld, coord_t > >
    {
        typedef PolyhedralGrid< dim, dimworld, coord_t >  Grid;

        const Grid& grid_;
        const int cartesianSize_;

        int computeCartesianSize() const
        {
            int size = cartesianDimensions()[ 0 ];
            for( int d=1; d<dim; ++d )
                size *= cartesianDimensions()[ d ];
            return size ;
        }
    public:
        static const int dimension = Grid :: dimension ;

        explicit CartesianIndexMapper( const Grid& grid )
          : grid_( grid ),
            cartesianSize_( computeCartesianSize() )
        {}

        const std::array<int, dimension>& cartesianDimensions() const
        {
          return grid_.logicalCartesianSize();
        }

        int cartesianSize() const
        {
            return cartesianSize_;
        }

        int compressedSize() const
        {
            return grid_.size( 0 );
        }

        // Only for unifying calls with CartesianIndexMapper<CpGrid> where levels are relevant.
        int compressedLevelZeroSize() const
        {
            return grid_.size( 0 );
        }

        int cartesianIndex( const int compressedElementIndex ) const
        {
            assert( compressedElementIndex >= 0 && compressedElementIndex < compressedSize() );
            return grid_.globalCell()[ compressedElementIndex ];
        }

        void cartesianCoordinate(const int compressedElementIndex, std::array<int,dimension>& coords) const
        {
          int gc = cartesianIndex( compressedElementIndex );
          if( dimension >=2 )
          {
              for( int d=0; d<dimension-2; ++d )
              {
                coords[d] = gc % cartesianDimensions()[d];  gc /= cartesianDimensions()[d];
              }

              coords[dimension-2] = gc % cartesianDimensions()[dimension-2];
              coords[dimension-1] = gc / cartesianDimensions()[dimension-1];
          }
          else
              coords[ 0 ] = gc ;
        }

        // Only for unifying calls with CartesianIndexMapper<CpGrid> where levels are relevant.
        void cartesianCoordinateLevel(const int compressedElementIndexOnLevel, std::array<int,dimension>& coordsOnLevel, int level) const
        {
            if (level) {
                throw std::invalid_argument("Invalid level.\n");
            }
            cartesianCoordinate(compressedElementIndexOnLevel, coordsOnLevel);
        }
    };

} // end namespace Opm
#endif
