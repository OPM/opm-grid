#ifndef OPM_POLYHEDRALGRID_VERTEQCOLUMNUTILITY_HEADER
#define OPM_POLYHEDRALGRID_VERTEQCOLUMNUTILITY_HEADER

#include <array>
#include <cassert>

#include <dune/common/exceptions.hh>
#include <dune/common/iteratorfacades.hh>

#include <opm/grid/common/VerteqColumnUtility.hpp>

namespace Dune
{
  /** \brief Interface class to access the columns in a verteq grid
    */
  template< int dim, int dimworld, typename coord_t >
  class VerteqColumnUtility< PolyhedralGrid< dim, dimworld, coord_t > >
  {
    typedef PolyhedralGrid< dim, dimworld, coord_t >  Grid;
    typedef typename Grid :: TopSurfaceGridType TopSurfaceGridType;

    const TopSurfaceGridType* topSurfaceGrid_;

    // iterator through columns,
    // ForwardIteratorFacade provides the standard iterator methods
    class Iterator : public ForwardIteratorFacade< Iterator, ColumnCell, ColumnCell >
    {
      // see VerteqColumnUtility.hpp
      typedef ColumnCell ColumnCellType;

      const TopSurfaceGridType* topSurfaceGrid_;
      int idx_;

    public:
      Iterator( const TopSurfaceGridType* topSurfaceGrid, const int idx )
        : topSurfaceGrid_( topSurfaceGrid ),
          idx_( idx )
      {}

      ColumnCellType dereference() const
      {
        // -1 means non-valid iterator
        assert( idx_ >= 0 );
        // if iterator is valid topSurf must be available.
        assert( topSurfaceGrid_ );
        return ColumnCellType( *topSurfaceGrid_, idx_ );
      }

      void increment ()
      {
        // -1 means non-valid iterator
        assert( idx_ >= 0 );
        ++idx_;
      }

      bool equals ( const Iterator& other ) const
      {
        return idx_ == other.idx_;
      }
    };

  public:
    typedef typename Grid :: template Codim< 0 > :: Entity   Entity;

    typedef Iterator  IteratorType;

    /** \brief dimension of the grid */
    static const int dimension = Grid :: dimension ;

    /** \brief constructor taking grid */
    explicit VerteqColumnUtility( const Grid& grid )
      : topSurfaceGrid_( grid.topSurfaceGrid() )
    {
    }

    IteratorType begin( const Entity& entity ) const
    {
      if( topSurfaceGrid_ )
      {
        const int colIdx = topSurfaceGrid_->col_cells[ entity.impl().index() ];
        return IteratorType( topSurfaceGrid_,
                             topSurfaceGrid_->col_cellpos[ colIdx ] );
      }
      else
        return end( entity );
    }

    IteratorType end( const Entity& entity ) const
    {
      if( topSurfaceGrid_ )
      {
        const int colIdx = topSurfaceGrid_->col_cells[ entity.impl().index() ];
        return IteratorType( topSurfaceGrid_,
                             topSurfaceGrid_->col_cellpos[ colIdx+1 ] );
      }
      else
        return IteratorType( nullptr, -1 );
    }
  };

} // end namespace Dune
#endif
