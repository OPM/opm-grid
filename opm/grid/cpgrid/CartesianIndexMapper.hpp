#ifndef OPM_CPGRIDCARTESIANINDEXMAPPER_HEADER
#define OPM_CPGRIDCARTESIANINDEXMAPPER_HEADER

#include <array>
#include <cassert>
#include <stdexcept>

#include <opm/grid/common/CartesianIndexMapper.hpp>
#include <opm/grid/CpGrid.hpp>

namespace Dune
{
    template<>
    class CartesianIndexMapper< CpGrid >
    {
    public:
        static const int dimension = 3 ;
    protected:
        typedef CpGrid Grid;
        const Grid& grid_;
        const int cartesianSize_;

        int computeCartesianSize() const
        {
            int size = cartesianDimensions()[ 0 ];
            for( int d=1; d<dimension; ++d )
                size *= cartesianDimensions()[ d ];
            return size;
        }

    public:
        explicit CartesianIndexMapper( const Grid& grid )
            : grid_( grid ),
              cartesianSize_( computeCartesianSize() )
        {
        }

        const std::array<int, dimension>& cartesianDimensions() const
        {
            // For now, return the level-zero logical Cartesian size.
            // Note: grid_.logicalCartesianSize() can vary depending on how refinement was applied
            // (e.g., via addLgrsUpdateLeafView(...), adapt(), or globalRefine()).
            // This includes cases where all elements of the level-zero grid have been refined.
            return grid_.currentData().front()->logicalCartesianSize();
        }

        int cartesianSize() const
        {
            return cartesianSize_;
        }

        int compressedSize() const
        {
            return grid_.globalCell().size();
        }

        int cartesianIndex( const int compressedElementIndex ) const
        {
            assert(  compressedElementIndex >= 0 && compressedElementIndex < compressedSize() );
            return grid_.globalCell()[ compressedElementIndex ];
        }

        void cartesianCoordinate(const int compressedElementIndex, std::array<int,dimension>& coords) const
        {
            grid_.getIJK( compressedElementIndex, coords );
        }
    };

} // end namespace Opm
#endif
