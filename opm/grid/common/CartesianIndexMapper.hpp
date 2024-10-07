#ifndef OPM_CARTESIANINDEXMAPPER_HEADER
#define OPM_CARTESIANINDEXMAPPER_HEADER

#include <array>
#include <cassert>

#include <dune/common/exceptions.hh>

namespace Dune
{
    /** \brief Interface class to access the logical Cartesian grid as used in industry
               standard simulator decks.
               */
    template< class Grid >
    class CartesianIndexMapper
    {
    public:
        /** \brief dimension of the grid */
        static const int dimension = Grid :: dimension ;

        /** \brief constructor taking grid */
        explicit CartesianIndexMapper( const Grid& )
        {
            DUNE_THROW(InvalidStateException,"CartesianIndexMapper not specialized for given grid");
        }

        /** \brief return Cartesian dimensions, i.e. number of cells in each direction  */
        const std::array<int, dimension>& cartesianDimensions() const
        {
            static std::array<int, dimension> a;
            return a;
        }

        /** \brief return total number of cells in the logical Cartesian grid */
        int cartesianSize() const
        {
            return 0;
        }

        /** \brief return number of cells in the active grid */
        int compressedSize() const
        {
            return 0;
        }

        /** \brief return index of the cells in the logical Cartesian grid */
        int cartesianIndex( const int /* compressedElementIndex */) const
        {
            return 0;
        }

        /** \brief return Cartesian coordinate, i.e. IJK, for a given cell */
        void cartesianCoordinate(const int /* compressedElementIndex */, std::array<int,dimension>& /* coords */) const
        {
        }

        /// Only relevant for CpGrid specialization.
        /** \brief return number of cells in the active level zero grid. Only relevant for CpGrid specialization. */
        int compressedLevelZeroSize() const
        {
            return 0;
        }

        /** \brief return Cartesian coordinate, i.e. IJK, for a given cell. Only relevant for CpGrid specialization.*/
        void cartesianCoordinateLevel(const int /* compressedElementIndexOnLevel */,
                                      std::array<int,dimension>& /* coordsOnLevel */, int /*level*/) const
        {
        }

        /** \brief return index of the cells in the logical Cartesian grid, for refined level grids. Only relevant for CpGrid specialization. */
        int cartesianIndexLevel( const int /* compressedElementIndex */ , const int level) const
        {
            return 0;
        }
        /// Only relevant for CpGrid specialization. END
    };

} // end namespace Opm
#endif
