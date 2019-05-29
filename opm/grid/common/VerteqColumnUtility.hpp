#ifndef OPM_VERTEQCOLUMNUTILITY_HEADER
#define OPM_VERTEQCOLUMNUTILITY_HEADER

#include <array>
#include <cassert>

#include <dune/common/exceptions.hh>
#include <opm/grid/verteq/topsurf.hpp>


namespace Dune
{

    class ColumnCell
    {
        typedef Opm::TopSurf  TopSurfaceGridType;
        const TopSurfaceGridType& topSurf_;
        const int index_;

    public:
        ColumnCell( const TopSurfaceGridType& topSurf, const int index )
          : topSurf_( topSurf ),
            index_( index )
        {}

        int index () const { return index_; }

        double dz () const { return topSurf_.dz[ index_ ]; }
        double h () const  { return topSurf_.h [ index_ ]; }

        int maxVertRes() const { return topSurf_.max_vert_res; }

        // this can be used to access an entity in the original 3D grid
        int fineCellIndex() const { return topSurf_.fine_col[ index_ ]; }
    };


    /** \brief Interface class to access the logical Cartesian grid as used in industry
               standard simulator decks.
               */
    template< class Grid >
    class VerteqColumnUtility
    {
    public:
        /** \brief dimension of the grid */
        static const int dimension = Grid :: dimension ;

        /** \brief constructor taking grid */
        explicit VerteqColumnUtility( const Grid& )
        {
            DUNE_THROW(InvalidStateException,"CartesianIndexMapper not specialized for given grid");
        }
    };

} // end namespace Opm
#endif
