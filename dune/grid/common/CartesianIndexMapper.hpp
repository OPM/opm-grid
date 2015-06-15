#ifndef OPM_CARTESIANINDEXMAPPER_HEADER
#define OPM_CARTESIANINDEXMAPPER_HEADER

#include <array>
#include <cassert>

#include <dune/common/exceptions.hh>

namespace Dune
{
    template< class Grid >
    class CartesianIndexMapper
    {
    public:
        static const int dimension = Grid :: dimension ;
        CartesianIndexMapper( const Grid& )
        {
            DUNE_THROW(InvalidStateException,"CartesianIndexMapper not specialized for given grid");
        }

        const std::array<int, dimension>& cartesianDimensions() const
        {
            static std::array<int, dimension> a;
            return a;
        }

        int cartesianSize() const
        {
            return 0;
        }

        int compressedSize() const
        {
            return 0;
        }

        int cartesianIndex( const int compressedElementIndex ) const
        {
            return 0;
        }

        void cartesianCoordinate(const int compressedElementIndex, std::array<int,dimension>& coords) const
        {
        }
    };

} // end namespace Opm
#endif
