/**
 * @brief Holds the implementation of the CpGrid Data.
 * @author Markus Blatt <markus@dr-blatt.de>
 *         Atgeirr F Rasmussen <atgeirr@sintef.no>
 *         Bård Skaflestad     <bard.skaflestad@sintef.no>
 */

#ifndef OPM_CPGRIDDATA_HEADER
#define OPM_CPGRIDDATA_HEADER

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif
#include <array>
#include "OrientedEntityTable.hpp"
#include "DefaultGeometryPolicy.hpp"
#include <opm/core/grid/cpgpreprocess/preprocess.h>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <dune/common/collectivecommunication.hh>

namespace Dune
{
class CpGrid;

namespace cpgrid
{

class IndexSet;
class IdSet;
class GlobalIdSet;
template<int,int> class Geometry;
template<int> class Entity;
template<int> class EntityRep;

/** 
 * @brief Struct that hods all the data needed to represent a 
 * Cpgrid.
 */
class CpGridData
{
private:
    CpGridData(const CpGridData& g)
    : ccobj_(g.ccobj_)
    {}
public:
    /// Constructor
    /// \param grid  The grid that we are the data of.
    explicit CpGridData(CpGrid& grid);
    /// Default constructor.
    CpGridData();
    /// Destructor
    ~CpGridData();
    /// number of leaf entities per codim in this process
    int size(int codim) const;
    
    /// number of leaf entities per geometry type in this process
    int size (GeometryType type) const
    {
        if (type.isCube()) {
            return size(3 - type.dim());
        } else {
            return 0;
        }
    }
    /// Read the Sintef legacy grid format ('topogeom').
    /// \param grid_prefix the grid name, such that topology is
    /// found in <grid_prefix>-topo.dat etc.
    void readSintefLegacyFormat(const std::string& grid_prefix);

    /// Write the Sintef legacy grid format ('topogeom').
    /// \param grid_prefix the grid name, such that topology will be
    /// found in <grid_prefix>-topo.dat etc.
    void writeSintefLegacyFormat(const std::string& grid_prefix) const;

    /// Read the Eclipse grid format ('grdecl').
    /// \param filename the name of the file to read.
    /// \param z_tolerance points along a pillar that are closer together in z
    ///        coordinate than this parameter, will be replaced by a single point.
    /// \param periodic_extension if true, the grid will be (possibly) refined, so that
    ///        intersections/faces along i and j boundaries will match those on the other
    ///        side. That is, i- faces will match i+ faces etc.
    void readEclipseFormat(const std::string& filename, double z_tolerance, bool periodic_extension, bool turn_normals = false);

    /// Read the Eclipse grid format ('grdecl').
    /// \param input_data the data contained in a parser object.
    /// \param z_tolerance points along a pillar that are closer together in z
    ///        coordinate than this parameter, will be replaced by a single point.
    /// \param periodic_extension if true, the grid will be (possibly) refined, so that
    ///        intersections/faces along i and j boundaries will match those on the other
    ///        side. That is, i- faces will match i+ faces etc.
    /// \param turn_normals if true, all normals will be turned. This is intended for handling inputs with wrong orientations.
    /// \param clip_z if true, the grid will be clipped so that the top and bottom will be planar.
    void processEclipseFormat(const Opm::EclipseGridParser& input_parser, double z_tolerance, bool periodic_extension, bool turn_normals = false, bool clip_z = false);

    /// Read the Eclipse grid format ('grdecl').
    /// \param input_data the data in grdecl format, declared in preprocess.h.
    /// \param z_tolerance points along a pillar that are closer together in z
    ///        coordinate than this parameter, will be replaced by a single point.
    /// \param remove_ij_boundary if true, will remove (i, j) boundaries. Used internally.
    void processEclipseFormat(const grdecl& input_data, double z_tolerance, bool remove_ij_boundary, bool turn_normals = false);


    /// @brief
    ///    Extract Cartesian index triplet (i,j,k) of an active cell.
    ///
    /// @param [in] c
    ///    Active cell index.
    ///
    /// @param [out] ijk  Cartesian index triplet
    void getIJK(int c, std::array<int,3>& ijk) const
    {
        int gc = global_cell_[c];
        ijk[0] = gc % logical_cartesian_size_[0];  gc /= logical_cartesian_size_[0];
        ijk[1] = gc % logical_cartesian_size_[1];
        ijk[2] = gc / logical_cartesian_size_[1];
    }

    // Make unique boundary ids for all intersections.
    void computeUniqueBoundaryIds();

    /// Is the grid currently using unique boundary ids?
    /// \return true if each boundary intersection has a unique id
    ///         false if we use the (default) 1-6 ids for i- i+ j- j+ k- k+ boundaries.
    bool uniqueBoundaryIds() const
    {
        return use_unique_boundary_ids_;
    }
    
    /// Set whether we want to have unique boundary ids.
    /// \param uids if true, each boundary intersection will have a unique boundary id.
    void setUniqueBoundaryIds(bool uids)
    {
        use_unique_boundary_ids_ = uids;
        if (use_unique_boundary_ids_ && unique_boundary_ids_.empty()) {
            computeUniqueBoundaryIds();
        }
    }
    /// Get the index set. This is the lead as well as th level index set.
    /// \return The index set.
    const IndexSet& indexSet() const
    {
        return *index_set_;
    }
    /// The logical cartesian size of the grid.
    /// This function is not part of the Dune grid interface,
    /// and should be used with caution.
    const std::array<int, 3>& logicalCartesianSize() const
    {
        return logical_cartesian_size_;
    }

    /// \brief Redistribute a global grid.
    ///
    /// The whole grid must be available on all processors.
    /// \return The distributed grid Data
    void distributeGlobalGrid(const CpGrid& grid,
                              const CpGridData& view_data,
                              const std::vector<int>& cell_part);
    
private:
    // Representing the topology
    /** @brief Container for lookup of the faces attached to each cell. */
    cpgrid::OrientedEntityTable<0, 1> cell_to_face_;
    /** @brief Container for the lookup attaching cells for each face. */
    cpgrid::OrientedEntityTable<1, 0> face_to_cell_;
    /** @brief Container for the lookup of the points for each face. */
    Opm::SparseTable<int>             face_to_point_;
    /** @brief Vector that contains an arrays of the points of each cell*/
    std::vector< array<int,8> >       cell_to_point_;
    /** @brief The size of the underlying logical cartesian grid. 
     *
     * In a Eclipse a cornerpoint grid has the same number of cells
     * in each pillar. Note that of these some may have no volume
     * and this be inactive.
     */
    std::array<int, 3>                logical_cartesian_size_;
    /** @brief vector with the gobal cell index for each cell.
     *
     * Note the size of this container is determined by the 
     * underlying cartesian grid.
     */
    std::vector<int>                  global_cell_;
    /** @brief The tag of the faces. */
    cpgrid::EntityVariable<enum face_tag, 1> face_tag_; // {LEFT, BACK, TOP}
    /** @brief The geometries representing the grid. */
    cpgrid::DefaultGeometryPolicy geometry_;
    /** @brief The type of a point in the grid. */
    typedef FieldVector<double, 3> PointType;
    /** @brief The face normals of the grid. */ 
    cpgrid::SignedEntityVariable<PointType, 1> face_normals_;
    /** @brief All corners of the grid. */
    std::vector<PointType> allcorners_; // Yes, this is already stored in the point geometries. \TODO Improve by removing it.
    /** @brief The boundary ids. */
    cpgrid::EntityVariable<int, 1> unique_boundary_ids_;
    /** @brief The index set of the grid (level). */
    cpgrid::IndexSet* index_set_;
    /** @brief The local id set. */
    const cpgrid::IdSet* local_id_set_;
    /** @brief The global id set. */
    GlobalIdSet* global_id_set_;
    /// \brief The type of the collective communication.
    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::CollectiveCommunication<MPICommunicator> CollectiveCommunication;
    /** @brief Object for collective communication operations. */
    CollectiveCommunication ccobj_;
    
    // Boundary information (optional).
    bool use_unique_boundary_ids_;

    // Return the geometry vector corresponding to the given codim.
    template <int codim>
    const EntityVariable<Geometry<3 - codim, 3>, codim>& geomVector() const
    {
        return geometry_.geomVector<codim>();
    }

    friend class Dune::CpGrid;
    template<int> friend class Entity;
    template<int> friend class EntityRep;
    template<int> friend class EntityPointer;
    friend class Intersection;
    
};
} // end namspace cpgrid
} // end namespace Dune

#endif
