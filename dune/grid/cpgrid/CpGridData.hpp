//===========================================================================
//
// File: CpGrid.hpp
//
// Created: Sep 17 21:11:41 2013
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
//            Markus Blatt        <markus@dr-blatt.de>
//
// Comment: Major parts of this file originated in dune/grid/CpGrid.hpp
//          and got transfered here during refactoring for the parallelization.
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010, 2013 Statoil ASA.
  Copyright 2013 Dr. Blatt - HPC-Simulation-Software & Services

  This file is part of The Open Porous Media project  (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
 * @brief Holds the implementation of the CpGrid as a pimple.
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
#include <dune/common/tuples.hh>
#include "OrientedEntityTable.hpp"
#include "DefaultGeometryPolicy.hpp"
#include <opm/core/grid/cpgpreprocess/preprocess.h>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <dune/common/collectivecommunication.hh>
#include <dune/common/parallel/indexset.hh>
#include <dune/common/parallel/interface.hh>
#include <dune/common/parallel/plocalindex.hh>
#include <dune/common/parallel/variablesizecommunicator.hh>
#include <dune/grid/common/gridenums.hh>
#include "Entity2IndexDataHandle.hpp"

namespace Dune
{
class CpGrid;

namespace cpgrid
{

class IndexSet;
class IdSet;
class GlobalIdSet;
class PartitionTypeIndicator;
template<int,int> class Geometry;
template<int> class Entity;
template<int> class EntityRep;

namespace mover
{
template<class T, int i> class Mover;
}

/** 
 * @brief Struct that hods all the data needed to represent a 
 * Cpgrid.
 */
class CpGridData
{
    template<class T, int i> friend class mover::Mover;
    
private:
    CpGridData(const CpGridData& g);
    
public:
    /// Constructor
    /// \param grid  The grid that we are the data of.
    explicit CpGridData(CpGrid& grid);
#if HAVE_MPI
    /// Constructor for parallel grid data.
    /// \param comm The MPI communicator
    /// Default constructor.
    CpGridData(MPI_Comm comm);
#endif
    /// Constructor
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
    
    /// \brief communicate objects for all codims on a given level
    /// \param data The data handle describing the data. Has to adhere to the 
    /// Dune::DataHandleIF interface.
    /// \param iftype The interface to use for the communication.
    /// \param dir The direction of the communication along the interface (forward or backward).
    template<class DataHandle>
    void communicate(DataHandle& data, InterfaceType iftype, CommunicationDirection dir);

    template<bool forward, class DataHandle>
    void moveData(DataHandle& data, CpGridData* global_data, 
                  CpGridData* distributed_data);
    
private:

    template<bool forward, int codim, class DataHandle>
    void moveCodimData(DataHandle& data, CpGridData* global_data, 
                          CpGridData* distributed_data);

    typedef VariableSizeCommunicator<>::InterfaceMap InterfaceMap;

    /// \brief Communicates data of a given codimension
    /// \tparam codim The codimension
    /// \tparam DataHandle The type of the data handle describing, gathering,
    ///  and gathering the data.
    /// \param DataHandle The data handle describing, gathering,
    ///  and gathering the data.
    /// \param dir The direction of the communication.
    /// \param interface The information about the communication interface
    template<int codim, class DataHandle>
    void communicateCodim(DataHandle& data, CommunicationDirection dir,
                          Interface& interface);

    /// \brief Communicates data of a given codimension
    /// \tparam codim The codimension
    /// \tparam DataHandle The type of the data handle describing, gathering,
    ///  and gathering the data.
    /// \param DataHandle The data handle describing, gathering,
    ///  and gathering the data.
    /// \param dir The direction of the communication.
    /// \param interface The information about the communication interface
    template<int codim, class DataHandle>
    void communicateCodim(DataHandle& data, CommunicationDirection dir,
                          InterfaceMap& interface);

    // Representing the topology
    /** @brief Container for lookup of the faces attached to each cell. */
    cpgrid::OrientedEntityTable<0, 1> cell_to_face_;
    /** 
     * @brief Container for the lookup of attaching cells for each face. 
     *
     * All faces have two neighbours except for those at the domain boundary.
     * @warn  Note that along the front partition there are invalid neighbours
     * marked with index std::numeric_limits<int>::max()
     */
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
    /** @brief The indicator of the partition type of the entities */
    PartitionTypeIndicator* partition_type_indicator_;
    
    /// \brief The type of the collective communication.
    typedef MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::CollectiveCommunication<MPICommunicator> CollectiveCommunication;
    /// \brief Object for collective communication operations.
    CollectiveCommunication ccobj_;
    
    // Boundary information (optional).
    bool use_unique_boundary_ids_;

    /// \brief The type of the set of the attributes
    enum AttributeSet{owner, overlap};

    /// \brief The type of the parallel index set
    typedef Dune::ParallelIndexSet<int,ParallelLocalIndex<AttributeSet> > ParallelIndexSet;

    /// \brief The parallel index set of the cells.
    ParallelIndexSet cell_indexset_;
    
    /// \brief The type of the remote indices information
    typedef Dune::RemoteIndices<ParallelIndexSet> RemoteIndices;
    /*
    /// \brief The remote index information for the cells.
    RemoteIndices cell_remote_indices;
    */
    /// \brief Communication interface for the cells.
    tuple<Interface,Interface,Interface,Interface,Interface> cell_interfaces_;
    
    /// \brief Interface from interior and border to interior and border for the faces.
    tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap> 
    face_interfaces_;

    /// \brief Interface from interior and border to interior and border for the faces.
    tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap> 
    point_interfaces_;
    
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
    friend class PartitionTypeIndicator;
};

namespace
{
/// \brief Get a value from a tuple according to the interface type.
/// \tparam T The type of the values in the tuple.
/// \param iftype The interface type.
/// \param interfaces A tuple with the values order by interface type.
template<class T>
T& getInterface(InterfaceType iftype, 
                Dune::tuple<T,T,T,T,T>& interfaces)
{
    switch(iftype)
    {
    case 0:
        return get<0>(interfaces);
    case 1:
        return get<1>(interfaces);
    case 2:
        return get<2>(interfaces);
    case 3:
        return get<3>(interfaces);
    case 4:
        return get<4>(interfaces);
    }
    OPM_THROW(std::runtime_error, "Invalid Interface type was used during communication");
}

/*
/// \brief Get the interface map associated with an interface type.
/// \param iftype The interface type.
/// \param interfaces A tuple with the interfaces order by interface type.

VariableSizeCommunicator<>::InterfaceMap&
getInterface(InterfaceType iftype, 
                           tuple<Interface,Interface,Interface,Interface,Interface>& interfaces)
 {
    return get(iftype, interfaces).interfaces();
 }

/// \brief Get the interface map associated with an interface type.
/// \param iftype The interface type.
/// \param interfaces A tuple with the interfaces order by interface type.
VariableSizeCommunicator<>::InterfaceMap& 
getInterface(InterfaceType iftype, 
                           tuple<VariableSizeCommunicator<>::InterfaceMap,
                                 VariableSizeCommunicator<>::InterfaceMap,
                                 VariableSizeCommunicator<>::InterfaceMap,
                                 VariableSizeCommunicator<>::InterfaceMap,
                                 VariableSizeCommunicator<>::InterfaceMap>& interfaces)
{
    return get(iftype, interfaces);
}
*/
} // end unnamed namespace

template<int codim, class DataHandle>
void CpGridData::communicateCodim(DataHandle& data, CommunicationDirection dir,
                                  Interface& interface)
{
    Entity2IndexDataHandle<DataHandle, codim> data_wrapper(*this, data);
    VariableSizeCommunicator<> comm(interface);
    if(dir==ForwardCommunication)
        comm.forward(data_wrapper);
    else
        comm.backward(data_wrapper);
}

template<int codim, class DataHandle>
void CpGridData::communicateCodim(DataHandle& data, CommunicationDirection dir,
                                  InterfaceMap& interface)
{
    Entity2IndexDataHandle<DataHandle, codim> data_wrapper(*this, data);
    VariableSizeCommunicator<> comm(ccobj_, interface);
    if(dir==ForwardCommunication)
        comm.forward(data_wrapper);
    else
        comm.backward(data_wrapper);
}

template<class DataHandle>
void CpGridData::communicate(DataHandle& data, InterfaceType iftype,
                             CommunicationDirection dir)
{
    if(data.contains(3,0))
        communicateCodim<0>(data, dir, getInterface(iftype, cell_interfaces_));
    if(data.contains(3,1))
        communicateCodim<1>(data, dir, getInterface(iftype, point_interfaces_));
    if(data.contains(3,3))
        communicateCodim<3>(data, dir, getInterface(iftype, point_interfaces_));
}

namespace
{

template<class T>
class MoveBuffer
{
public:
    void read(T& data)
    {
        data=buffer_[index_++];
    }
    void write(const T& data)
    {
        buffer_[index_++]=data;
    }
    void reset()
    {
        index_=0;
    }
    void resize(std::size_t size)
    {
        buffer_.resize(size);
        index_=0;
    }
private:
    std::vector<T> buffer_;
    typename std::vector<T>::size_type index_;
};

}

namespace mover
{
template<class DataHandle,int codim>
struct Mover
{
};

template<class DataHandle>
struct BaseMover
{
    BaseMover(DataHandle& data)
    : data_(data)
    {}
    template<class E>
    void moveData(E& from, E& to)
    {
        std::size_t size=data_.size(from);
        buffer.resize(size);
        data_.gather(buffer, from);
        buffer.reset();
        data_.scatter(buffer, to, size);
    }
    DataHandle& data_;
    MoveBuffer<typename DataHandle::DataType> buffer;
};


template<class DataHandle>
struct Mover<DataHandle,0> : BaseMover<DataHandle>
{
    Mover<DataHandle,0>(DataHandle& data, CpGridData* gatherView,
                        CpGridData* scatterView)
    : BaseMover<DataHandle>(data), gatherView_(gatherView), scatterView_(scatterView)
    {}

    void operator()(std::size_t from_cell_index,std::size_t to_cell_index)
    {
        Entity<0> from_entity=Entity<0>(gatherView_, from_cell_index, true);
        Entity<0> to_entity=Entity<0>(scatterView_, to_cell_index, true);
        this->moveData(from_entity, to_entity);
    }
    CpGridData* gatherView_;
    CpGridData* scatterView_;
};

template<class DataHandle>
struct Mover<DataHandle,1> : BaseMover<DataHandle>
{
    Mover<DataHandle,1>(DataHandle& data, CpGridData* gatherView,
                        CpGridData* scatterView)
    : BaseMover<DataHandle>(data), gatherView_(gatherView), scatterView_(scatterView)
    {}

    void operator()(std::size_t from_cell_index,std::size_t to_cell_index)
    {
        typedef typename OrientedEntityTable<0,1>::row_type::iterator Iter;
        EntityRep<0> from_cell=EntityRep<0>(from_cell_index, true);
        EntityRep<0> to_cell=Entity<0>(to_cell_index, true);
        OrientedEntityTable<0,1>& table = gatherView_->cell_to_face_;
        Iter from_face=table.operator[](from_cell).begin();
        Iter end_from_face=gatherView_->cell_to_face_[from_cell].end();
        Iter to_face=scatterView_->cell_to_face_[to_cell].begin();

        for(;from_face!=end_from_face; ++from_face, ++to_face)
        {
            //assert(to_face->valid());
            assert(to_face!=scatterView_->cell_to_face_[to_cell].end());
            this->moveData(*from_face, *to_face);
        }
    }
    CpGridData *gatherView_;
    CpGridData *scatterView_;
};

template<class DataHandle>
struct Mover<DataHandle,3> : BaseMover<DataHandle>
{
    Mover<DataHandle,3>(DataHandle& data, CpGridData* gatherView,
                        CpGridData* scatterView)
    : BaseMover<DataHandle>(data), gatherView_(gatherView), scatterView_(scatterView)
    {}
    void operator()(std::size_t from_cell_index,std::size_t to_cell_index)
    {
        const std::array<int,8>& from_cell_points=
            gatherView_->cell_to_point_[from_cell_index];
        const std::array<int,8>& to_cell_points=
            scatterView_->cell_to_point_[to_cell_index];
        for(std::size_t i=0; i<8; ++i)
        {
            this->moveData(Entity<1>(*gatherView_, from_cell_points[i], true),
                           Entity<1>(*scatterView_, to_cell_points[i], true));
        }
    }
    CpGridData* gatherView_;
    CpGridData* scatterView_;
};

} // end mover namespace
template<bool forward, class DataHandle>
void CpGridData::moveData(DataHandle& data, CpGridData* global_data, 
                          CpGridData* distributed_data)
{
    if(data.contains(3,0))
       moveCodimData<forward,0>(data, global_data, distributed_data);
    if(data.contains(3,1))
       moveCodimData<forward,1>(data, global_data, distributed_data);
    if(data.contains(3,3))
       moveCodimData<forward,3>(data, global_data, distributed_data);
}

template<bool forward, int codim, class DataHandle>
void CpGridData::moveCodimData(DataHandle& data, CpGridData* global_data, 
                          CpGridData* distributed_data)
{
    CpGridData *gather_view, *scatter_view;
    
    if(forward)
    {
        gather_view=global_data;
        scatter_view=distributed_data;
    }
    else
    {
        gather_view=distributed_data;
        scatter_view=global_data;
    }
    
    mover::Mover<DataHandle,codim> mover(data, gather_view, scatter_view);
    
    typedef typename ParallelIndexSet::const_iterator Iter;
    for(Iter index=distributed_data->cell_indexset_.begin(),
            end = distributed_data->cell_indexset_.end();
        index!=end; ++index)
    {
        std::size_t from=forward?index->global():index->local();
        std::size_t to=forward?index->local():index->global();
        mover(from,to);
    }
}
} // end namspace cpgrid
} // end namespace Dune

#endif
