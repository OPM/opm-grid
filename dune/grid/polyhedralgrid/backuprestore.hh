#ifndef DUNE_POLYHEDRALGRID_BACKUPRESTORE_HH
#define DUNE_POLYHEDRALGRID_BACKUPRESTORE_HH

#include <dune/grid/common/backuprestore.hh>
#include <dune/grid/utility/grapedataioformattypes.hh>

#include <dune/grid/polyhedralgrid/capabilities.hh>

namespace Dune
{

  // BackupRestoreFacilities
  // -----------------------

  template< class Grid, bool hasBackupRestoreFacilities = Capabilities::hasBackupRestoreFacilities< Grid > ::v >
  class PolyhedralGridBackupRestoreFacilities
  {};

  template< class Grid >
  class PolyhedralGridBackupRestoreFacilities< Grid, true >
  {
    typedef PolyhedralGridBackupRestoreFacilities< Grid, true > This;

  protected:
    PolyhedralGridBackupRestoreFacilities ()
    {}

  private:
    PolyhedralGridBackupRestoreFacilities ( const This & );
    This &operator= ( const This & );

  public:
    template< GrapeIOFileFormatType type >
    bool writeGrid ( const std::string &filename, double time ) const
    {
      return asImp().hostGrid().template writeGrid< type >( filename, time );
    }

    template< GrapeIOFileFormatType type >
    bool readGrid ( const std::string &filename, double &time )
    {
      const bool success
        = asImp().hostGrid().template readGrid< type >( filename, time );
      asImp().update();
      return success;
    }

  protected:
    const Grid &asImp () const
    {
      return static_cast< const Grid & >( *this );
    }

    Grid &asImp ()
    {
      return static_cast< Grid & >( *this );
    }
  };


  // BackupRestoreFacility for PolyhedralGrid
  // --------------------------------

  template< int dim, int dimworld >
  struct BackupRestoreFacility< PolyhedralGrid< dim, dimworld > >
  {
    typedef PolyhedralGrid< dim, dimworld > Grid;
    typedef BackupRestoreFacility< HostGrid > HostBackupRestoreFacility;

    static void backup ( const Grid &grid, const std::string &path, const std::string &fileprefix )
    {
      HostBackupRestoreFacility::backup( grid.hostGrid(), path, fileprefix );
    }

    static void backup ( const Grid &grid, const std::ostream &stream )
    {
      HostBackupRestoreFacility::backup( grid.hostGrid(), stream );
    }

    static Grid *restore ( const std::string &path, const std::string &fileprefix )
    {
      return HostBackupRestoreFacility::restore( path, fileprefix );
    }

    static Grid *restore ( const std::istream &stream )
    {
      return HostBackupRestoreFacility::restore( stream );
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_BACKUPRESTORE_HH
