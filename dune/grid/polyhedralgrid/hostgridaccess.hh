#ifndef DUNE_POLYHEDRALGRID_HOSTGRIDACCESS_HH
#define DUNE_POLYHEDRALGRID_HOSTGRIDACCESS_HH

#include <dune/grid/common/intersection.hh>

#include <dune/grid/polyhedralgrid/declaration.hh>
#include <dune/grid/polyhedralgrid/intersection.hh>

namespace Dune
{

  // HostGridAccess
  // --------------

  template< class Grid >
  struct HostGridAccess;



  /** \class HostGridAccess
   *  \brief provides access to host grid objects
   *
   *  \tparam  Grid  meta grid, whose host grid shall be accessed
   *
   *  \nosubgrouping
   */
  template< class HG >
  struct HostGridAccess< PolyhedralGrid< HG > >
  {
    /** \name Exported Types
     * \{ */

    typedef PolyhedralGrid< HG > Grid;

    //! type of HostGrid
    typedef typename Grid::HostGrid HostGrid;

     /** \} */

    /** \brief A Traits struct that collects return types of class member methods.
     *
     *  \tparam codim codimension
     */
    template< int codim >
    struct Codim
    {
      //! type of the PolyhedralGrid entity
      typedef typename Grid::template Codim< codim >::Entity Entity;
      //! type of the PolyhedralGrid entity pointer
      typedef typename Grid::template Codim< codim >::EntityPointer EntityPointer;

      //! type of the host entity
      typedef typename HostGrid::template Codim< codim >::Entity HostEntity;
      //! type of the host entity pointer
      typedef typename HostGrid::template Codim< codim >::EntityPointer HostEntityPointer;
    };

    //! type of the PolyhedralGrid leaf intersection
    typedef typename Grid::Traits::LeafIntersection LeafIntersection;
    //! type of the PolyhedralGrid level intersection
    typedef typename Grid::Traits::LevelIntersection LevelIntersection;

    //! type of the host leaf intersection
    typedef typename HostGrid::Traits::LeafIntersection HostLeafIntersection;
    //! type of the host level intersection
    typedef typename HostGrid::Traits::LevelIntersection HostLevelIntersection;

    /** \brief Get underlying HostGrid.
     *  \param[in] grid  PolyhedralGrid
     *  \returns HostGrid
     */
    static const HostGrid &hostGrid ( const Grid &grid )
    {
      return grid.hostGrid();
    }

    template< class Entity >
    static const typename Codim< Entity::codimension >::HostEntity &
    hostEntity ( const Entity &entity )
    {
      return hostEntity< Entity::codimension >( entity );
    }

    template< int codim >
    static const typename Codim< codim >::HostEntity &
    hostEntity ( const typename Codim< codim >::Entity &entity )
    {
      return Grid::getRealImplementation( entity ).hostEntity();
    }

    template< class HostIntersection >
    static const HostIntersection &
    hostIntersection ( const Intersection< const Grid, PolyhedralGridIntersection< const Grid, HostIntersection > > &intersection )
    {
      return Grid::getRealImplementation( intersection ).hostIntersection();
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_HOSTGRIDACCESS_HH
