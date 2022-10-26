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
  Copyright 2009, 2010, 2013, 2022 Equinor ASA.
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

// Warning suppression for Dune includes.
#include <opm/grid/utility/platform_dependent/disable_warnings.h>

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>
#ifdef HAVE_DUNE_ISTL
#include <dune/istl/owneroverlapcopy.hh>
#endif

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 7)
#include <dune/common/parallel/communication.hh>
#else
#include <dune/common/parallel/collectivecommunication.hh>
#endif
#include <dune/common/parallel/indexset.hh>
#include <dune/common/parallel/interface.hh>
#include <dune/common/parallel/plocalindex.hh>
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 7)
#include <dune/common/parallel/variablesizecommunicator.hh>
#else
#include <opm/grid/utility/VariableSizeCommunicator.hpp>
#endif
#include <dune/grid/common/gridenums.hh>
#include <dune/geometry/type.hh>

#include <opm/grid/utility/platform_dependent/reenable_warnings.h>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/Grid/NNC.hpp>
#endif

#include <array>
#include <tuple>
#include <algorithm>
#include <set>

#include "OrientedEntityTable.hpp"
#include "DefaultGeometryPolicy.hpp"
#include <opm/grid/cpgpreprocess/preprocess.h>

#include "Entity2IndexDataHandle.hpp"
#include "DataHandleWrappers.hpp"
#include "GlobalIdMapping.hpp"
#include "Geometry.hpp"

namespace Opm
{
class EclipseState;
}
namespace Dune
{
class CpGrid;

namespace cpgrid
{

class IndexSet;
class IdSet;
class LevelGlobalIdSet;
class PartitionTypeIndicator;
template<int,int> class Geometry;
template<int> class Entity;
template<int> class EntityRep;
}
}

void refine_and_check(const Dune::cpgrid::Geometry<3, 3>&,
                      const std::array<int, 3>&,
                    bool);
void refinePatch_and_check(const std::array<int,3>&,
                           const std::array<int,3>&,
                           const std::array<int,3>&);

void refinePatch_and_check(Dune::CpGrid&,
                           const std::array<int,3>&,
                           const std::array<int,3>&,
                           const std::array<int,3>&);

namespace Dune
{
namespace cpgrid
{
namespace mover
{
template<class T, int i> struct Mover;
}

/**
 * @brief Struct that hods all the data needed to represent a
 * Cpgrid.
 */
class CpGridData
{
    template<class T, int i> friend struct mover::Mover;

    friend class GlobalIdSet;
    
    friend
    void ::refine_and_check(const Dune::cpgrid::Geometry<3, 3>&,
                            const std::array<int, 3>&,
                            bool);
    friend
    void::refinePatch_and_check(const std::array<int,3>&,
                        const std::array<int,3>&,
                        const std::array<int,3>&);

    friend
    void ::refinePatch_and_check(Dune::CpGrid&,
                                 const std::array<int,3>&,
                                 const std::array<int,3>&,
                                 const std::array<int,3>&);

private:
    CpGridData(const CpGridData& g);

public:
    enum{
#ifndef MAX_DATA_COMMUNICATED_PER_ENTITY
        /// \brief The maximum data items allowed per cell (DUNE < 2.5.2)
        ///
        /// Due to a bug in DUNE < 2.5.2 we need to limit this when
        /// communicating. 1 is big enough for OPM as we always use
        /// one block for all unknowns, but some DUNE's grid checks
        /// actually need 2. So 2 it is.
        MAX_DATA_PER_CELL = 2
#else
        /// \brief The maximum data items allowed per cell (DUNE < 2.5.2)
        ///
        /// Due to a bug in DUNE < 2.5.2 we need to limit this when
        /// communicating. Uses the define MAX_DATA_COMMUNICATED_PER_ENTITY.
        MAX_DATA_PER_CELL = MAX_DATA_COMMUNICATED_PER_ENTITY
#endif
    };     

    /// Constructor for parallel grid data.
    /// \param comm The MPI communicator
    /// Default constructor.
    explicit CpGridData(MPIHelper::MPICommunicator comm);

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
    /// \param periodic_extension if true, the grid will be (possibly) refined, so that
    ///        intersections/faces along i and j boundaries will match those on the other
    ///        side. That is, i- faces will match i+ faces etc.
    void readEclipseFormat(const std::string& filename, bool periodic_extension, bool turn_normals = false);

#if HAVE_ECL_INPUT
    /// Read the Eclipse grid format ('grdecl').
    /// \param deck the parsed deck from opm-parser (which is a low-level object)
    /// \param periodic_extension if true, the grid will be (possibly) refined, so that
    ///        intersections/faces along i and j boundaries will match those on the other
    ///        side. That is, i- faces will match i+ faces etc.
    /// \param turn_normals if true, all normals will be turned. This is intended for handling inputs with wrong orientations.
    /// \param clip_z if true, the grid will be clipped so that the top and bottom will be planar.
    /// \param poreVolume pore volumes for use in MINPV processing, if asked for in deck
    void processEclipseFormat(const Opm::Deck& deck, bool periodic_extension, bool turn_normals = false, bool clip_z = false,
                              const std::vector<double>& poreVolume = std::vector<double>());

    /// Read the Eclipse grid format ('grdecl').
    /// \param ecl_grid the high-level object from opm-parser which represents the simulation's grid
    ///        In a parallel run this may be a nullptr on all ranks but rank zero.
    /// \param ecl_state the object from opm-parser provide information regarding to pore volume, NNC,
    ///        aquifer information when ecl_state is available. NNC and aquifer connection
    ///        information will also be updated during the function call when available and necessary.
    /// \param periodic_extension if true, the grid will be (possibly) refined, so that
    ///        intersections/faces along i and j boundaries will match those on the other
    ///        side. That is, i- faces will match i+ faces etc.
    /// \param turn_normals if true, all normals will be turned. This is intended for handling inputs with wrong orientations.
    /// \param clip_z if true, the grid will be clipped so that the top and bottom will be planar.
    std::vector<std::size_t> processEclipseFormat(const Opm::EclipseGrid* ecl_grid, Opm::EclipseState* ecl_state,
                                                  bool periodic_extension, bool turn_normals = false, bool clip_z = false, bool pinchActive = true);
#endif

    /// Read the Eclipse grid format ('grdecl').
    /// \param input_data the data in grdecl format, declared in preprocess.h.
    ///
    /// \param ecl_state the object from opm-parser provide information regarding to pore volume, NNC,
    ///        aquifer information when ecl_state is available. NNC and aquifer connection
    ///        information will also be updated during the function call when available and necessary.
    /// \param nnc is the non-neighboring connections
    /// \param remove_ij_boundary if true, will remove (i, j) boundaries. Used internally.
    /// \param pinchActive If true, we will add faces between vertical cells that have only inactive cells or cells
    ///            with zero volume between them. If false these cells will not be connected.
    void processEclipseFormat(const grdecl& input_data,
#if HAVE_ECL_INPUT
                              Opm::EclipseState* ecl_state,
#endif
                              std::array<std::set<std::pair<int, int>>, 2>& nnc,
                              bool remove_ij_boundary, bool turn_normals, bool pinchActive);

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

    // Given a start {i,j,k} and an end {i,j,k}, compute the dimension of the patch, i.e.
    // aomunt of cells in each direction. 
    const std::array<int,3> getPatchDim(const std::array<int,3>& start_ijk, const std::array<int,3>& end_ijk) const
    {
        return {end_ijk[0]-start_ijk[0], end_ijk[1]-start_ijk[1], end_ijk[2]-start_ijk[2]};
    } 
    
    const std::array<std::vector<int>,3> getPatchGeomIndices(const std::array<int,3>& start_ijk, const std::array<int,3>& end_ijk) const
    {
        // Get the patch dimension (total cells in each direction). Used to 'reserve vectors'.
        const std::array<int,3> patch_dim = getPatchDim(start_ijk, end_ijk);
        // Get grid dim
        const std::array<int,3> grid_dim = this -> logicalCartesianSize();
        // CORNERS
        std::vector<int> patch_corners;
        patch_corners.reserve((patch_dim[0]+1)*(patch_dim[1]+1)*(patch_dim[2]+1));
        for (int j = start_ijk[1]; j < end_ijk[1]+1; ++j) {
            for (int i = start_ijk[0]; i < end_ijk[0]+1; ++i) {
                for (int k = start_ijk[2]; k < end_ijk[2]+1; ++k) {
                    patch_corners.push_back((j*(grid_dim[0]+1)*(grid_dim[2]+1)) + (i*(grid_dim[2]+1))+k);
                } // end i-for-loop
            } // end j-for-loop
        } // end k-for-loop
        // FACES
        std::vector<int> patch_faces;
        // Integers to reserve patch_faces.
        int i_patch_faces = (patch_dim[0]+1)*patch_dim[1]*patch_dim[2];
        int j_patch_faces = patch_dim[0]*(patch_dim[1]+1)*patch_dim[2];
        int k_patch_faces = patch_dim[0]*patch_dim[1]*(patch_dim[2]+1);
        patch_faces.reserve(i_patch_faces +j_patch_faces +k_patch_faces);
        // Integers to compute face indices.
        int i_grid_faces = (grid_dim[0]+1)*grid_dim[1]*grid_dim[2];
        int j_grid_faces = grid_dim[0]*(grid_dim[1]+1)*grid_dim[2];
        // I_FACES
        for (int j = start_ijk[1]; j < end_ijk[1]; ++j) {
            for (int i = start_ijk[0]; i < end_ijk[0]+1; ++i) {
                for (int k = start_ijk[2]; k < end_ijk[2]; ++k) {
                    int face_idx = (j*grid_dim[1]*grid_dim[2]) + (i*grid_dim[2])+ k;
                    patch_faces.push_back(face_idx);          
                } // end k-for-loop
            } // end i-for-loop     
        } // end j-for-loop
        // J_FACES
        for (int j = start_ijk[1]; j < end_ijk[1]+1; ++j) {
            for (int i = start_ijk[0]; i < end_ijk[0]; ++i) {
                for (int k = start_ijk[2]; k < end_ijk[2]; ++k) {
                    int face_idx = i_grid_faces + (j*grid_dim[0]*grid_dim[2]) + (i*grid_dim[2])+ k; 
                    patch_faces.push_back(face_idx);          
                } // end k-for-loop
            } // end i-for-loop     
        } // end j-for-loop
         // K_FACES
        for (int j = start_ijk[1]; j < end_ijk[1]; ++j) {
            for (int i = start_ijk[0]; i < end_ijk[0]; ++i) {
                for (int k = start_ijk[2]; k < end_ijk[2]+1; ++k) {
                    int face_idx = j_grid_faces + i_grid_faces + (j*grid_dim[0]*(grid_dim[2]+1)) + (i*(grid_dim[2]+1))+ k;
                    patch_faces.push_back(face_idx);          
                } // end k-for-loop
            } // end i-for-loop     
        } // end j-for-loop

        // CELLS
        std::vector<int> patch_cells;
        patch_cells.reserve(patch_dim[0]*patch_dim[1]*patch_dim[2]);
        for (int k = start_ijk[2]; k < end_ijk[2]; ++k) {
            for (int j = start_ijk[1]; j < end_ijk[1]; ++j) {
                for (int i = start_ijk[0]; i < end_ijk[0]; ++i) {
                    patch_cells.push_back((k*grid_dim[0]*grid_dim[1]) + (j*grid_dim[1]) +i);
                } // end i-for-loop
            } // end j-for-loop
        } // end k-for-loop
        return {patch_corners, patch_faces, patch_cells};
    }

    // Construct a 'huge cell' out of a patch of connected (consecutive in each direction) cells.
    // CELL-FICATION OF A PATCH
    // This function takes a patch and build a cell out of it.
    // The function takes a connected patch formed by the product of consecutive cells in each direction, and
    // returns a Geometry<3,3> object, 'a patch cell', or 'a cellFIED patch'.
    // Idea: Select 8 corners of the boundary of the patch, compute center and volume.
    // @param patch_cells_indices           Indices of the cells from the grid that we want to refine, or, equivalently,
    //                                      indices of the cells that belong to the patch.
    // @param cellfiedPatch_to_point        Indices of the 8 corners of each parent cell.
    // @param cellfied_patch_geometry
    Geometry<3,3> cellfyPatch(const std::vector<int> patch_cells_indices)
    {
        if (patch_cells_indices.empty()){
            OPM_THROW(std::logic_error, "Empty patch. Cannot convert patch into cell.");
        }
        DefaultGeometryPolicy cellfied_patch_geometry;
        // Get the minimum and maximum of "patch_cells_indices"
        // to find the min_i, max_i, min_j, max_j, min_k, max_k,
        // to 'cell-fy' the patch (treating the patch as a 'huge cell')
        const std::array<int,2> min_max_indices = {
            *std::min_element(patch_cells_indices.begin(), patch_cells_indices.end()),
            *std::max_element(patch_cells_indices.begin(), patch_cells_indices.end())};
        // Get grid dim
        const std::array<int,3> grid_dim = this -> logicalCartesianSize();
        // Get min/max-ijk indices out of "min_max_indices"
        std::vector<std::array<int,3>> min_max_ijk_indices;
        min_max_ijk_indices.reserve(2);
        for (auto& idx : min_max_indices) {
            int i = idx/grid_dim[0]; // i
            int j = ((idx - i)/grid_dim[0])/grid_dim[1]; // j
            int k = (((idx - i)/grid_dim[0]) -j)/grid_dim[1]; // k
            min_max_ijk_indices.push_back({i,j,k});
        }
        // Get indices of (at most) 8 selected cells located on the boundary of the patch.
        /* std::array<int,8> selected_boundary_cell_indices = {
            // Index of the boundary cell from where corner '0' will be extracted.
            min_max_indices[0],
            // Index of the boundary cell from where corner '1' will be extracted: '{max_i, min_j, min_k}'
            (min_max_ijk_indices[0][2]*grid_dim[0]*grid_dim[1]) + (min_max_ijk_indices[0][1]*grid_dim[0])
            + min_max_ijk_indices[1][0],
            // Index of the bounday cell from where corner '2' will be extracted: '{min_i, max_j, min_k}'
            (min_max_ijk_indices[0][2]*grid_dim[0]*grid_dim[1]) + (min_max_ijk_indices[1][1]*grid_dim[0])
            + min_max_ijk_indices[0][0],
            // Index of the boundary cell from where corner '3' will be extracted: '{max_i, max_j, min_k}'
            (min_max_ijk_indices[0][2]*grid_dim[0]*grid_dim[1]) + (min_max_ijk_indices[1][1]*grid_dim[0])
            + min_max_ijk_indices[1][0],
            // Index of the bounday cell from where corner '4' will be extracted: '{min_i, min_j, max_k}'
            (min_max_ijk_indices[1][2]*grid_dim[0]*grid_dim[1]) + (min_max_ijk_indices[0][1]*grid_dim[0])
            + min_max_ijk_indices[0][0],
            // Index of the boundary cell from where corner '5' will be extracted: '{max_i, min_j, max_k}'
            (min_max_ijk_indices[1][2]*grid_dim[0]*grid_dim[1]) + (min_max_ijk_indices[0][1]*grid_dim[0])
            + min_max_ijk_indices[1][0],
            // Index of the bounday cell from where corner '6' will be extracted: '{min_i, max_j, max_k}'
            (min_max_ijk_indices[1][2]*grid_dim[0]*grid_dim[1]) + (min_max_ijk_indices[1][1]*grid_dim[0])
            + min_max_ijk_indices[0][0],
            // Index of the boundary cell from where corner '7' will be extracted.
            min_max_indices[1]}; */
        EntityVariableBase<cpgrid::Geometry<0,3>>& cellfiedPatch_corners =
                    cellfied_patch_geometry.geomVector(std::integral_constant<int,3>());
        cellfiedPatch_corners.resize(8);
        // Get the 8 corner indices of the 'cellFIED patch'
         std::array<int,8> cellfiedPatch_to_point = { //  corner indices: (J*(grid_dim[0]+1)*(grid_dim[2]+1)) + (I*(grid_dim[2]+1)) +K
            // Index of corner '0' {min_i, min_j, min_k}
            (min_max_ijk_indices[0][1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + (min_max_ijk_indices[0][0]*(grid_dim[2]+1))
            +min_max_ijk_indices[0][2],
            // Index of corner '1' '{max_i +1, min_j, min_k}'
            (min_max_ijk_indices[0][1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + ((min_max_ijk_indices[1][0] +1)*(grid_dim[2]+1))
            + min_max_ijk_indices[0][2],
            // Index of corner '2' '{min_i, max_j +1, min_k}'
            ((min_max_ijk_indices[1][1] +1)*(grid_dim[0]+1)*(grid_dim[2]+1)) + (min_max_ijk_indices[0][0]*(grid_dim[2]+1))
            + min_max_ijk_indices[0][2],
            // Index of corner '3' '{max_i +1, max_j +1, min_k}'
            ((min_max_ijk_indices[1][1] +1)*(grid_dim[0]+1)*(grid_dim[2]+1)) + ((min_max_ijk_indices[1][0] +1)*(grid_dim[2]+1))
             + min_max_ijk_indices[0][2],
             // Index of corner '4' '{min_i, min_j, max_k +1}'
             (min_max_ijk_indices[0][1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + (min_max_ijk_indices[0][0]*(grid_dim[2]+1))
             + min_max_ijk_indices[1][2],
             // Index of corner '5' '{max_i +1, min_j, max_k +1}'
             (min_max_ijk_indices[0][1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + ((min_max_ijk_indices[1][0] +1)*(grid_dim[2]+1))
             + min_max_ijk_indices[1][2] +1,
             // Index of corner '6' '{min_i, max_j +1, max_k +1}'
             ((min_max_ijk_indices[1][1] +1)*(grid_dim[0]+1)*(grid_dim[2]+1)) + (min_max_ijk_indices[0][0]*(grid_dim[2]+1))
             + min_max_ijk_indices[1][2] +1,
             // Index of corner '7' {max_i +1, max_j +1, max_k +1}
             ((min_max_ijk_indices[1][1]+1)*(grid_dim[0]+1)*(grid_dim[2]+1)) + ((min_max_ijk_indices[1][0]+1)*(grid_dim[2]+1))
             + min_max_ijk_indices[1][2] +1};
            // Center of the cell'fied' patch
            Geometry<0,3>::GlobalCoordinate cellfiedPatch_center = {0., 0.,0.};
            for (int corn; corn < 8; ++corn) {
                cellfiedPatch_center +=
                    (this -> geometry_.geomVector(std::integral_constant<int,3>()).get(cellfiedPatch_to_point[corn]).center())/8.;
                cellfiedPatch_corners[corn] =
                    this -> geometry_.geomVector(std::integral_constant<int,3>()).get(cellfiedPatch_to_point[corn]).center();
            }
            // Volume of the cell'fied' patch
            double cellfiedPatch_volume = 0.;
            for (auto idx : patch_cells_indices) {
                cellfiedPatch_volume += (this -> geometry_.geomVector(std::integral_constant<int,0>())
                                         [EntityRep<0>(idx, true)]).volume();
            }
            // Create a pointer to the first element of "cellfiedPatch_to_point"
            // (required as the fourth argement to construct a Geometry<3,3> type object).
            int* cellfiedPatch_indices_storage_ptr = &cellfiedPatch_to_point[0];
            // Construct (and return) the Geometry of the CEELfied PATCH.
            return Geometry<3,3>(cellfiedPatch_center, cellfiedPatch_volume,
                                 cellfied_patch_geometry.geomVector(std::integral_constant<int,3>()),
                                 cellfiedPatch_indices_storage_ptr);
        }
   
    
    
    /*  // AREA (via sum of 4 triangles) and CENTROID of a face given its 4 corners.
    // ----------- IN PROGRESS --------------
    std::tuple<double,Geometry<0,3>::GlobalCoordinate> getFaceAreaCentroid(const std::array<int,4> corners)
    {
        // AREA
        double face_area = 0.;
        // Face CENTROID.
        Geometry<0,3>::GlobalCoordinate face_centroid = {0.,0.,0.};
        for (auto& corner : corners)
        {
            face_centroid += (this -> geometry_.geomVector(std::integral_constant<int,1>()).get(corner).center())/4.;
        }
        
        return {face_area, face_centroid};
        
        }*/
    
    // VOLUME (via sum of 24 tetrahedra)  and CENTER of a hexaedron
    // given its corner and face indices.
    std::tuple<double,Geometry<0,3>::GlobalCoordinate>
    getHexaVolumeCenter(const std::array<int,8> corners, const std::array<int,6> faces)
    {
        // VOLUME HEXAHEDRON
        double hexa_volume = 0.0;
        // Hexa center.
        Geometry<0,3>::GlobalCoordinate hexa_center = {0.,0.,0.};
        for (auto& corner : corners)
        {
            hexa_center += (this -> geometry_.geomVector(std::integral_constant<int,3>()).get(corner).center())/8.;
        }
        // CENTROIDS of the faces of the hexahedron.
        // (one of the 6 corners of all 4 tetrahedra based on that face).
        std::vector<Geometry<0,3>::GlobalCoordinate> face_centroids;
        face_centroids.resize(6);
        for (auto& face : faces) {
            face_centroids.push_back(this -> geometry_.geomVector(std::integral_constant<int,1>())
                                     [Dune::cpgrid::EntityRep<1>(face, true)].center());
        }
        // Container with 6 entries, one per face. Each entry has the
        // 4 indices of the 4 corners of each face.
        std::vector<std::array<int,4>> hexa_face_to_point;
        hexa_face_to_point.reserve(6);
        for (int face = 0; face < 6;  ++face) {
            hexa_face_to_point.push_back(//this -> face_to_point_[faces[face]]);
            { this -> face_to_point_[faces[face]][0],
             this -> face_to_point_[faces[face]][1],
            this -> face_to_point_[faces[face]][2],
            this -> face_to_point_[faces[face]][3]});
        }
        // Container with indices of the edges of the 4 tetrahedra per face
        // [according to description above]
        std::vector<std::vector<std::array<int,2>>> tetra_edges;
        tetra_edges.reserve(6);
        for (auto& face4corners : hexa_face_to_point)
        {
            std::vector<std::array<int,2>> face4edges = {
                { face4corners[0], face4corners[1]}, // fake '{0,1}'/'{4,5}'
                { face4corners[0], face4corners[2]}, // fake '{0,2}'/'{4,6}'
                { face4corners[1], face4corners[3]}, // fake '{1,3}'/'{5,7}'
                { face4corners[2], face4corners[3]} }; // fake '{2,3}'/'{6,7}'
            tetra_edges.push_back(face4edges);
        }
        // Sum of the 24 volumes to get the volume of the hexahedron,
        // stored in "refined_cell_volume".
        // Calculate the volume of each hexahedron, by adding
        // the 4 tetrahedra at each face (4x6 = 24 tetrahedra).
        for (int face = 0; face < 6; ++face) {
            for (int edge = 0; edge < 4; ++edge) {
                // Construction of each tetrahedron based on "face" with one
                // of its edges equal to "edge".
                const Geometry<0, 3>::GlobalCoordinate tetra_corners[4] = {
                    this ->geometry_.geomVector(std::integral_constant<int,3>()).get(tetra_edges[face][edge][0]).center(),  
                    this ->geometry_.geomVector(std::integral_constant<int,3>()).get(tetra_edges[face][edge][1]).center(), 
                    face_centroids[face],
                    hexa_center};  
                hexa_volume += std::fabs(simplex_volume(tetra_corners));
            } // end edge-for-loop
        } // end face-for-loop
        return {hexa_volume, hexa_center};    
    }
    
    // Refine a single cell and return a shared pointer of CpGridData type.
    // REFINE A SINGLE CELL
    // @param cells_per_dim                 Number of sub-cells in each direction.
    // @param parent_idx                    Parent index.
    // @return refined_grid_ptr             Shared pointer pointing at refined_grid.
    //         parent_to_redined_corners    For each corner of the parent cell, we store the index of the
    //                                      refined corner that coincides with the old one.
    //                                      We assume they are ordered 0,1,..7
    //                                                              6---7
    //                                                      2---3   |   | TOP FACE
    //                                                      |   |   4---5
    //                                                      0---1 BOTTOM FACE
    //         parent_to_children_faces/cellFor each parent face/cell, we store its child face indices.
    //                                      {parent face/cell index in coarse level, {indices of its children in refined level}}
    //         child_to_parent_faces/cells  {child index, parent index}
    //         isParent_faces               True when the face got refined.
    //         isParent_cells               True when the cell got refined (here only "parent cell").
    std::tuple< const std::shared_ptr<CpGridData>,
                const std::vector<std::array<int,2>>,                // parent_to_refined_corners(~boundary_old_to_new_corners)
                const std::vector<std::tuple<int,std::vector<int>>>, // parent_to_children_faces (~boundary_old_to_new_faces)
                const std::tuple<int, std::vector<int>>,             // parent_to_children_cells
                const std::vector<std::array<int,2>>,                // child_to_parent_faces
                const std::vector<std::array<int,2>>,                // child_to_parent_cells
                const std::map<int,bool>,                            // isParent_faces
                const std::map<int,bool>>                            // isParent_cells
    refineSingleCell(const std::array<int,3>& cells_per_dim, const int& parent_idx)
    {
        std::shared_ptr<CpGridData> refined_grid_ptr = std::make_shared<CpGridData>(ccobj_);
        auto& refined_grid = *refined_grid_ptr;
        DefaultGeometryPolicy& refined_geometries = refined_grid.geometry_;
        std::vector<std::array<int,8>>& refined_cell_to_point = refined_grid.cell_to_point_;
        cpgrid::OrientedEntityTable<0,1>& refined_cell_to_face = refined_grid.cell_to_face_;
        Opm::SparseTable<int>& refined_face_to_point = refined_grid.face_to_point_;
        cpgrid::OrientedEntityTable<1,0>& refined_face_to_cell = refined_grid.face_to_cell_;
        cpgrid::EntityVariable<enum face_tag,1>& refined_face_tags = refined_grid.face_tag_;
        cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>,1>& refined_face_normals = refined_grid.face_normals_;
        // Get parent cell
        cpgrid::Geometry<3,3> parent_cell = geometry_.geomVector(std::integral_constant<int,0>())[EntityRep<0>(parent_idx, true)];
        // Refine parent cell
        parent_cell.refine(cells_per_dim,
                           refined_geometries,
                           refined_cell_to_point,
                           refined_cell_to_face,
                           refined_face_to_point,
                           refined_face_to_cell,
                           refined_face_tags,
                           refined_face_normals);

        // Get parent cell corners.
        std::array<int,8> parent_to_point = this->cell_to_point_[parent_idx];
        std::vector<std::array<int,2>> parent_to_refined_corners{
            // replacing parent-cell corner '0'
            {parent_to_point[0], 0},
            // replacing parent-cell corner '1'
            {parent_to_point[1], cells_per_dim[0]*(cells_per_dim[2]+1)},
            // replacing parent-cell corner '2'
            {parent_to_point[2], cells_per_dim[1]*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)},
            // replacing parent-cell corner '3'
            {parent_to_point[3], (cells_per_dim[1]*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (cells_per_dim[0]*(cells_per_dim[2]+1))},
            // replacing parent-cell corner '4'
            {parent_to_point[4], cells_per_dim[2]},
            // replacing parent-cell corner '5'
            {parent_to_point[5], (cells_per_dim[0]*(cells_per_dim[2]+1)) + cells_per_dim[2]},
            // replacing parent-cell corner '6'
            {parent_to_point[6], cells_per_dim[1]*(cells_per_dim[0]+1)*(cells_per_dim[2]+1) + cells_per_dim[2]},
            // replacing parent-cell corner '7'
            {parent_to_point[7], (cells_per_dim[1]*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (cells_per_dim[0]*(cells_per_dim[2]+1))
             + cells_per_dim[2]}};
        // Get relation old face -> new born faces (children faces)
        std::vector<std::tuple<int,std::vector<int>>>  parent_to_children_faces;
        // Get parent_cell_to_face = { {face, orientation}, {another face, its orientation}, ...}
        auto parent_cell_to_face = (this-> cell_to_face_[EntityRep<0>(parent_idx, true)]);
        //
        std::vector<std::array<int,2>> child_to_parent_faces;
        child_to_parent_faces.reserve(refined_face_to_cell.size());
        // Children faces are ordered following the criteria introduced in refine(), Geometry.hpp.
        // First K_FACES, then I_FACES, then J_FACES.
        int k_faces = cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1);
        int i_faces = (cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2];
        for (auto& face : parent_cell_to_face) {
            // Check face tag to identify the type of face (bottom, top, left, right, front, or back).
            auto parent_face_tag = (this-> face_tag_[Dune::cpgrid::EntityRep<1>(face.index(), true)]);
            // To store the new born faces for each face.
            std::vector<int> children_faces;
            // K_FACES
            if (parent_face_tag == face_tag::K_FACE) {
                children_faces.reserve(cells_per_dim[0]*cells_per_dim[1]);
                for (int j = 0; j < cells_per_dim[1]; ++j) {
                    for (int i = 0; i < cells_per_dim[0]; ++i) {
                        if (!face.orientation()) { // false -> BOTTOM FACE -> k=0
                            children_faces.push_back((j*cells_per_dim[0]) + i);
                            child_to_parent_faces.push_back({(j*cells_per_dim[0]) + i, face.index()});
                        }
                        else { // true -> TOP FACE -> k=cells_per_dim[2]
                            children_faces.push_back((cells_per_dim[2]*cells_per_dim[0]*cells_per_dim[1])
                                                     +(j*cells_per_dim[0]) + i);
                            child_to_parent_faces.push_back({(cells_per_dim[2]*cells_per_dim[0]*cells_per_dim[1])
                                    +(j*cells_per_dim[0]) + i, face.index()});
                        }
                    }
                }
            }
            // I_FACES
            if (parent_face_tag == face_tag::I_FACE) {
                children_faces.reserve(cells_per_dim[1]*cells_per_dim[2]);
                for (int k = 0; k < cells_per_dim[2]; ++k) {
                    for (int j = 0; j < cells_per_dim[1]; ++j) {
                        if (!face.orientation()) { // false -> LEFT FACE -> i=0
                            children_faces.push_back(k_faces + (k*cells_per_dim[1]) + j);
                            child_to_parent_faces.push_back({k_faces + (k*cells_per_dim[1]) + j, face.index()});
                        }
                        else { // true -> RIGHT FACE -> i=cells_per_dim[0]
                            children_faces.push_back(k_faces + (cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2])
                                                     + (k*cells_per_dim[1]) + j);
                            child_to_parent_faces.push_back({k_faces + (cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2])
                                    + (k*cells_per_dim[1]) + j, face.index()});
                        }
                    }
                }
            }
            // J_FACES
            if (parent_face_tag == face_tag::J_FACE) {
                children_faces.reserve(cells_per_dim[0]*cells_per_dim[2]);
                for (int i = 0; i < cells_per_dim[0]; ++i) {
                    for (int k = 0; k < cells_per_dim[2]; ++k) {
                        if (!face.orientation()) { // false -> FRONT FACE -> j=0
                            children_faces.push_back(k_faces + i_faces + (i*cells_per_dim[2]) + k);
                            child_to_parent_faces.push_back({k_faces + i_faces + (i*cells_per_dim[2]) + k, face.index()});
                        }
                        else { // true -> BACK FACE -> j=cells_per_dim[1]
                            children_faces.push_back(k_faces + i_faces  + (cells_per_dim[1]*cells_per_dim[0]*cells_per_dim[2])
                                                     + (i*cells_per_dim[2]) + k);
                            child_to_parent_faces.push_back({k_faces + i_faces  + (cells_per_dim[1]*cells_per_dim[0]*cells_per_dim[2])
                                    + (i*cells_per_dim[2]) + k, face.index()});
                        }
                    }
                }
            }
            std::tuple<int,std::vector<int>> aux_tuple = std::make_tuple(face.index(), children_faces);
            parent_to_children_faces.push_back(aux_tuple);
        }
        // PARENT TO CHILDREN CELLS
        std::vector<int> children_cells;
        // Child to parent cell
        std::vector<std::array<int, 2>> child_to_parent_cell;
        children_cells.reserve(cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]);
        child_to_parent_cell.reserve(cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]);
        for (int cell = 0; cell < cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]; ++cell) {
            children_cells.push_back(cell);
            child_to_parent_cell.push_back({cell, parent_idx});
        }


        std::tuple<int, std::vector<int>> parent_to_children_cells =
            std::make_tuple(parent_idx, children_cells);
        // IS PARENT - FACES
        // Map {index face, true/false}
        // true-> the face got refined (has children), false -> the face hasn't been refined (does not have children).
        std::map<int,bool> isParent_faces;
        for (int face = 0; face < this-> face_to_cell_.size(); ++face) {
            isParent_faces[face] = false;
        }
        // Rewrite the entries of the map for those faces that got refined.
        for (auto& face : parent_cell_to_face) {
            isParent_faces[face.index()] = true;
        }
        // IS PARENT - CELLS
        // Map {index cell, true/false}
        // true-> cell got refined (has children), false -> cell hasn't been refined (does not have children).
        std::map<int,bool> isParent_cells;
        for (int cell = 0; cell < this-> size(0); ++cell) {
            isParent_cells[cell] = false;
        }
        // Rewrite the entries of the map for those faces that got refined.
        isParent_cells[parent_idx] = true;

        return {refined_grid_ptr, parent_to_refined_corners,
            parent_to_children_faces, parent_to_children_cells, child_to_parent_faces, child_to_parent_cell,
            isParent_faces, isParent_cells};
    }
    // Refine a (connected block of cells) patch
    // REFINE A PATCH of CONNECTED (CONSECUTIVE in each direction) cells with 'uniform' regular intervals.
    // (meaning that the amount of children per cell is the same for all parent cells (cells of the patch)).
    // @param cells_per_dim                       Number of sub-cells in each direction (for each cell).
    // @param start_ijk                           Minimum values of i,j,k to construct the patch (start).
    // @param end_ijk                             Maximum values of i,j,k to construct the patch (end).
    // @return refined_grid_ptr                   Shared pointer of CpGridData type, pointing at the refined_grid
    //         boundary_old_to_new_corners/faces  Corners/faces on the boundary of the patch get replaced by new born refined one(s).
    //         parent_to_children_faces/cells     To store the indices of 'all the face/cell children' of each parent.
    //         child_to_parent_faces/cells        For each child-face, we store its parent face index.
    //         isParent_faces/cells               True for each face/cell that got refined.

    std::tuple<std::shared_ptr<CpGridData>,
               const std::vector<std::array<int,2>>,                // boundary_old_to_new_corners
               const std::vector<std::tuple<int,std::vector<int>>>, // boundary_old_to_new_faces
               const std::vector<std::tuple<int,std::vector<int>>>, // parent_to_children_faces
               const std::vector<std::tuple<int,std::vector<int>>>, // parent_to_children_cell
               const std::vector<std::array<int,2>>,                // child_to_parent_faces
               const std::vector<std::array<int,2>>,                // child_to_parent_cells
               const std::map<int,bool>,                            // isParent_faces
               const std::map<int,bool>>                            // isParent_cells
    refineBlockPatch(const std::array<int,3>& cells_per_dim,
                     const std::array<int,3>& start_ijk, const std::array<int,3>& end_ijk)
    {
        std::shared_ptr<CpGridData> refined_grid_ptr = std::make_shared<CpGridData>(ccobj_);
        auto& refined_grid = *refined_grid_ptr;
        DefaultGeometryPolicy& refined_geometries = refined_grid.geometry_;
        std::vector<std::array<int,8>>& refined_cell_to_point = refined_grid.cell_to_point_;
        cpgrid::OrientedEntityTable<0,1>& refined_cell_to_face = refined_grid.cell_to_face_;
        Opm::SparseTable<int>& refined_face_to_point = refined_grid.face_to_point_;
        cpgrid::OrientedEntityTable<1,0>& refined_face_to_cell = refined_grid.face_to_cell_;
        cpgrid::EntityVariable<enum face_tag,1>& refined_face_tags = refined_grid.face_tag_;
        cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>,1>& refined_face_normals = refined_grid.face_normals_;

        // Patch information (built from the grid).
        const auto& patch_dim = getPatchDim(start_ijk, end_ijk);
        // Coarse grid dimension
        const std::array<int,3>& grid_dim = this -> logicalCartesianSize();
        // If the patch contains only one cell:
        if ((patch_dim[0] == 1) && (patch_dim[1] == 1) && (patch_dim[2] == 1)){
            const int& parent_cell = (start_ijk[2]*grid_dim[0]*grid_dim[1]) + (start_ijk[1]*grid_dim[0]) +start_ijk[0];
            auto [refined_grid_ptr0, parent_to_refined_corners,
                  parent_to_children_faces, parent_to_children_cells,child_to_parent_faces, child_to_parent_cell,
                  isParent_faces, isParent_cells] = this->refineSingleCell(cells_per_dim, parent_cell);
            // When the patch is only one cell,
            // - boundary_old_to_new_corners == parent_to_refined_corners.
            // - boundary_old_to_new_faces == parent_to_children_faces.
            // Fix the type of parent_to_children_cells to return it correctly.
            std::vector<std::tuple<int, std::vector<int>>> parent_to_children_cells_vec = {parent_to_children_cells};
            return {refined_grid_ptr0, parent_to_refined_corners,
                parent_to_children_faces, parent_to_children_faces, parent_to_children_cells_vec,
                child_to_parent_faces, child_to_parent_cell,
                isParent_faces, isParent_cells};
        }
        // When the patch consists in more than one cell:
        else {
            const auto& [patch_corners, patch_faces, patch_cells] = getPatchGeomIndices(start_ijk, end_ijk); 
            // Construct the Geometry of the CEELfied PATCH.
            cpgrid::Geometry<3,3> cellfied_patch = this -> cellfyPatch(patch_cells);
            
            // Some integers to reduce notation later.
            int xfactor = cells_per_dim[0]*patch_dim[0];
            int yfactor = cells_per_dim[1]*patch_dim[1];
            int zfactor = cells_per_dim[2]*patch_dim[2];
            // Refine the cell "cellfied_patch"
            cellfied_patch.refine({xfactor, yfactor, zfactor},
                                  refined_geometries,
                                  refined_cell_to_point,
                                  refined_cell_to_face,
                                  refined_face_to_point,
                                  refined_face_to_cell,
                                  refined_face_tags,
                                  refined_face_normals);
            // For each cell of the patch, we select the 8 refined corners that coincide with the
            // original cell corners, and the refined child-faces that ('all' together) 'coincide' with
            // each of the 6 original faces of the cell.
            // To store the 8 new refined corners coinciding with the 8 original corners of each cell of the patch,
            // we choose a map, key = cell index in the coarse grid, value = {{old corner '0', new corner '0'}, ...}
            std::vector<std::array<int,2>> boundary_old_to_new_corners;
            // RESERVE?
            for (int j = start_ijk[1]; j < end_ijk[1]+1; ++j) {
                for (int i = start_ijk[0]; i < end_ijk[0]+1; ++i) {
                    for (int k = start_ijk[2]; k < end_ijk[2]+1; ++k) {
                        // WE ARE ASSUMING THE GRID FROM LEVEL 0 HAS THE DUNE NUMBERING FOR CORNERS.
                        if ( (j == start_ijk[1]) || (j == end_ijk[1]) ){
                            int old_corner_idx = (j*(grid_dim[2]+1)*(grid_dim[0]+1)) + (i*(grid_dim[2]+1)) + k;
                            int new_corner_idx = (cells_per_dim[1]*(j-start_ijk[1])*(xfactor +1)*(zfactor +1))
                                + (cells_per_dim[0]*(i-start_ijk[0])*(zfactor +1))
                                + (cells_per_dim[2]*(k-start_ijk[2]));
                            boundary_old_to_new_corners.push_back({old_corner_idx, new_corner_idx});
                        }
                        if ( (i == start_ijk[0]) || (i == end_ijk[0]) ) {
                            int old_corner_idx = (j*(grid_dim[2]+1)*(grid_dim[0]+1)) + (i*(grid_dim[2]+1)) + k;
                            int new_corner_idx = (cells_per_dim[1]*(j-start_ijk[1])*(xfactor +1)*(zfactor +1))
                                + (cells_per_dim[0]*(i-start_ijk[0])*(zfactor +1))
                                + (cells_per_dim[2]*(k-start_ijk[2]));
                            boundary_old_to_new_corners.push_back({old_corner_idx, new_corner_idx});
                        }
                        if ( (k == start_ijk[2]) || (k == end_ijk[2]) ) {
                            int old_corner_idx = (j*(grid_dim[2]+1)*(grid_dim[0]+1)) + (i*(grid_dim[2]+1)) + k;
                            int new_corner_idx = (cells_per_dim[1]*(j-start_ijk[1])*(xfactor +1)*(zfactor +1))
                                + (cells_per_dim[0]*(i-start_ijk[0])*(zfactor +1))
                                + (cells_per_dim[2]*(k-start_ijk[2]));
                            boundary_old_to_new_corners.push_back({old_corner_idx, new_corner_idx});
                        }
                    } // end k-for-loop
                } // end i-for-loop
            } // end j-for-loop

            // FACE STUFF
            // Boundary patch faces.
            std::vector<int> boundary_patch_faces;
            // Auxiliary integers to simplify notation.
            int bound_patch_faces = 2*patch_dim[1]*patch_dim[2] // left/right boundary
            +  patch_dim[0]*2*patch_dim[2] // front/back boundary
            + patch_dim[0]*patch_dim[1]*2; // bottom/top boundary
            boundary_patch_faces.reserve(bound_patch_faces);
            // Boundary old faces associated with its (new born) children faces.
            std::vector<std::tuple<int, std::vector<int>>> boundary_old_to_new_faces; // {face index, its children indices}
            boundary_old_to_new_faces.reserve(bound_patch_faces);
            // Auxiliary integers to simplify notation.
            int i_grid_faces =  (grid_dim[0]+1)*grid_dim[1]*grid_dim[2];
            int j_grid_faces =  grid_dim[0]*(grid_dim[1]+1)*grid_dim[2];
            // Parent to children faces (all faces of the patch, not only the ones on the boundary).
            std::vector<std::tuple<int,std::vector<int>>> parent_to_children_faces;
            parent_to_children_faces.reserve(patch_faces.size());
            // For each child, we associate it with its parent.
            std::vector<std::array<int,2>> child_to_parent_faces; // {child index (in 'level 1'), parent index (in 'level 0')}
            child_to_parent_faces.reserve(refined_face_to_cell.size());
            // I_FACES
            for (int j = start_ijk[1]; j < end_ijk[1]; ++j) {
                for (int i = start_ijk[0]; i < end_ijk[0]+1; ++i) {
                    for (int k = start_ijk[2]; k < end_ijk[2]; ++k) {
                        int face_idx = (j*grid_dim[1]*grid_dim[2]) + (i*grid_dim[2])+ k;
                        // Vector to store new born faces, per face.
                        // CHILDREN FACES ARE ORDERED AS IN refine(), Geometry.hpp
                        std::vector<int> children_list;  // I_FACE ikj (xzy-direction)
                        // l,m,n play the role of 'x,y,z-direction'
                        // lnm = fake ikj (how I_FACES are 'ordered' in refine())
                        for (int l = (i-start_ijk[0])*cells_per_dim[0]; l < (i-start_ijk[0]+1)*cells_per_dim[0]; ++l) {
                            for (int n = (k-start_ijk[2])*cells_per_dim[2];n < (k-start_ijk[2]+1)*cells_per_dim[2]; ++n) {
                                for (int m = (j-start_ijk[1])*cells_per_dim[1]; m < (j-start_ijk[1]+1)*cells_per_dim[1]; ++m) {
                                    children_list.push_back((xfactor*yfactor*(zfactor+1)) +(l*yfactor*zfactor) + (n*yfactor) + m);
                                    child_to_parent_faces.push_back({(xfactor*yfactor*(zfactor+1)) +(l*yfactor*zfactor)
                                            + (n*yfactor) + m, face_idx});
                                } // end m-for-loop
                            } // end n-for-loop
                        } // end l-for-loop
                        // Add parent information of each face to "parent_to_children_faces".
                        parent_to_children_faces.push_back(std::make_tuple(face_idx, children_list));
                        if ((i == start_ijk[0]) || (i == end_ijk[0])) {
                            boundary_patch_faces.push_back(face_idx);
                            // Associate each old face on the boundary of the patch with the new born ones.
                            boundary_old_to_new_faces.push_back(std::make_tuple(face_idx, children_list));
                        }
                    } // end k-for-loop
                } // end i-for-loop
            } // end j-for-loop
            // J_FACES
            for (int j = start_ijk[1]; j < end_ijk[1]+1; ++j) {
                for (int i = start_ijk[0]; i < end_ijk[0]; ++i) {
                    for (int k = start_ijk[2]; k < end_ijk[2]; ++k) {
                        int face_idx = i_grid_faces + (j*grid_dim[0]*grid_dim[2]) + (i*grid_dim[2])+ k;
                        // Vector to store new born faces, per face.
                        // CHILDREN FACES ARE ORDERED AS IN refine(), Geometry.hpp
                        std::vector<int> children_list;  // J_FACE jik (yxz-direction)
                        // l,m,n play the role of 'x,y,z-direction'
                        // mln = fake jik (how J_FACES are 'ordered' in refine())
                        for (int m = (j-start_ijk[1])*cells_per_dim[1]; m < (j-start_ijk[1]+1)*cells_per_dim[1]; ++m) {
                            for (int l = (i-start_ijk[0])*cells_per_dim[0]; l < (i-start_ijk[0]+1)*cells_per_dim[0]; ++l) {
                                for (int n = (k-start_ijk[2])*cells_per_dim[2]; n < (k-start_ijk[2]+1)*cells_per_dim[2]; ++n) {
                                    children_list.push_back((xfactor*yfactor*(zfactor+1)) + ((xfactor+1)*yfactor*zfactor)
                                                            + (m*xfactor*zfactor) + (l*zfactor)+n);
                                    child_to_parent_faces.push_back({(xfactor*yfactor*(zfactor+1)) + ((xfactor+1)*yfactor*zfactor)
                                            + (m*xfactor*zfactor) + (l*zfactor)+n, face_idx});
                                } // end n-for-loop
                            } // end l-for-loop
                        } // end m-for-loop
                        // Add parent information of each face to "parent_to_children_faces".
                        parent_to_children_faces.push_back(std::make_tuple(face_idx, children_list));
                        if ((j == start_ijk[1]) || (j == end_ijk[1])) {
                            boundary_patch_faces.push_back(face_idx);
                            // Associate each old face on the boundary of the patch with the new born ones.
                            boundary_old_to_new_faces.push_back(std::make_tuple(face_idx, children_list));
                        }
                    } // end k-for-loop
                } // end i-for-loop
            } // end j-for-loop
            // K_FACES
            for (int j = start_ijk[1]; j < end_ijk[1]; ++j) {
                for (int i = start_ijk[0]; i < end_ijk[0]; ++i) {
                    for (int k = start_ijk[2]; k < end_ijk[2]+1; ++k) {
                        int face_idx = i_grid_faces + j_grid_faces + (j*grid_dim[0]*grid_dim[2]) + (i*grid_dim[2])+ k;
                        // Vector to store new born faces, per face.
                        // CHILDREN FACES ARE ORDERED AS IN refine(), Geometry.hpp
                        std::vector<int> children_list;  // K_FACE kji (zyx-direction)
                        // l,m,n play the role of 'x,y,z-direction'
                        // nml = fake kji (how K_FACES are 'ordered' in refine())
                        for (int n = (k-start_ijk[2])*cells_per_dim[2]; n < (k-start_ijk[2]+1)*cells_per_dim[2]; ++n) {
                            for (int m = (j-start_ijk[1])*cells_per_dim[1]; m < (j-start_ijk[1]+1)*cells_per_dim[1]; ++m) {
                                for (int l = (i-start_ijk[0])*cells_per_dim[0]; l < (i-start_ijk[0]+1)*cells_per_dim[0]; ++l) {
                                    children_list.push_back((n*xfactor*yfactor) + (m*xfactor)+ l);
                                    child_to_parent_faces.push_back({(n*xfactor*yfactor) + (m*xfactor)+ l, face_idx});
                                } // end m-for-loop
                            } // end n-for-loop
                        } // end l-for-loop
                          // Add parent information of each face to "parent_to_children_faces".
                        parent_to_children_faces.push_back(std::make_tuple(face_idx, children_list));
                        if ((k == start_ijk[2]) || (k == end_ijk[2])) {
                            boundary_patch_faces.push_back(face_idx);
                            // Associate each old face on the boundary of the patch with the new born ones.
                            boundary_old_to_new_faces.push_back(std::make_tuple(face_idx, children_list));
                        }
                    } // end k-for-loop
                } // end i-for-loop
            } // end j-for-loop
            //--- END PARENT TO CHILDREN FACES / child_to_parent_faces ---

            // PARENT TO CHILDREN - CELLS  / child_to_parent_cells
            //
            // To store children indices for each parent. Each entry looks like
            // {parent index in the patch (coarse grid), index of one of its children in the refined grid}
            std::vector<std::tuple<int,std::vector<int>>> parent_to_children_cells;
            parent_to_children_cells.reserve(patch_dim[0]*patch_dim[1]*patch_dim[2]);
            // To store parent index for each child. The children are numbering
            // following the rule of moving first in the x-axes (from left to right),
            // then y-axes (from front to back), finally z-axes (from bottom to top).
            std::vector<std::array<int,2>> child_to_parent_cells; // {child index (in 'level 1'), parent cell index (in 'level 0')}
            child_to_parent_cells.reserve(xfactor*yfactor*zfactor);
            for (int k = 0; k < grid_dim[2]; ++k) {
                for (int j = 0; j < grid_dim[1]; ++j) {
                    for (int i = 0; i < grid_dim[0]; ++i) {
                        int cell_idx = (k*grid_dim[0]*grid_dim[1]) + (j*grid_dim[0]) +i;
                        std::vector<int> children_list;
                        if ( (i > start_ijk[0]-1) && (i < end_ijk[0]) && (j > start_ijk[1]-1) && (j < end_ijk[1])
                             && (k > start_ijk[2]-1) && (k < end_ijk[2])) {
                            for (int n = (k-start_ijk[2])*cells_per_dim[2]; n < (k-start_ijk[2]+1)*cells_per_dim[2]; ++n) {
                                for (int m = (j-start_ijk[1])*cells_per_dim[1]; m < (j-start_ijk[1]+1)*cells_per_dim[1]; ++m) {
                                    for (int l = (i-start_ijk[0])*cells_per_dim[0]; l < (i-start_ijk[0]+1)*cells_per_dim[0]; ++l) {
                                        children_list.push_back((n*xfactor*yfactor) + (m*xfactor) + l);
                                        child_to_parent_cells.push_back({(n*xfactor*yfactor) + (m*xfactor) + l, cell_idx});
                                    }// end l-for-loop
                                } // end m-for-loop
                            } // end n-for-loop
                            parent_to_children_cells.push_back(std::make_tuple(cell_idx, children_list));
                        }// end if 'patch cells'
                    } // end i-for-loop
                } // end j-for-loop
            } // end k-for-loop

            // IS PARENT - FACES
            // Map {index face, true/false}
            // true-> the face got refined (has children), false -> the face hasn't been refined (does not have children).
            std::map<int,bool> isParent_faces;
            for (int face = 0; face < this-> face_to_cell_.size(); ++face) {
                isParent_faces[face] = false;
            }
            // Rewrite the entries of the map for those faces that got refined.
            for (auto& face : patch_faces) {
                isParent_faces[face] = true;
            }
            // IS PARENT - CELLS
            // Map {index cell, true/false}
            // true-> cell got refined (has children), false -> cell hasn't been refined (does not have children).
            std::map<int,bool> isParent_cells;
            for (int cell = 0; cell < this-> size(0); ++cell) {
                isParent_cells[cell] = false;
            }
            // Rewrite the entries of the map for those faces that got refined.
            for (auto& cell : patch_cells) {
                isParent_cells[cell] = true;
            }
            return {refined_grid_ptr, boundary_old_to_new_corners, boundary_old_to_new_faces,
                parent_to_children_faces, parent_to_children_cells, child_to_parent_faces, child_to_parent_cells,
                isParent_faces, isParent_cells};
        }
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

    /// Return the internalized zcorn copy from the grid processing, if
    /// no cells were adjusted during the minpvprocessing this can be
    /// and empty vector.
    const std::vector<double>& zcornData() const {
        return zcorn;
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
    void distributeGlobalGrid(CpGrid& grid,
                              const CpGridData& view_data,
                              const std::vector<int>& cell_part);

    /// \brief communicate objects for all codims on a given level
    /// \param data The data handle describing the data. Has to adhere to the
    /// Dune::DataHandleIF interface.
    /// \param iftype The interface to use for the communication.
    /// \param dir The direction of the communication along the interface (forward or backward).
    template<class DataHandle>
    void communicate(DataHandle& data, InterfaceType iftype, CommunicationDirection dir);

#if HAVE_MPI
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 7)
    /// \brief The type of the  Communicator.
    using Communicator = VariableSizeCommunicator<>;
#else
    /// \brief The type of the Communicator.
    using Communicator = Opm::VariableSizeCommunicator<>;
#endif

    /// \brief The type of the map describing communication interfaces.
    using InterfaceMap = Communicator::InterfaceMap;

    /// \brief type of OwnerOverlap communication for cells
    using CommunicationType = Dune::OwnerOverlapCopyCommunication<int,int>;

    /// \brief The type of the parallel index set
    using  ParallelIndexSet = CommunicationType::ParallelIndexSet;

    /// \brief The type of the remote indices information
    using RemoteIndices = Dune::RemoteIndices<ParallelIndexSet>;

    /// \brief Get the owner-overlap-copy communication for cells
    ///
    /// Suitable e.g. for parallel linear algebra used by CCFV
    CommunicationType& cellCommunication()
    {
        return cell_comm_;
    }

    /// \brief Get the owner-overlap-copy communication for cells
    ///
    /// Suitable e.g. for parallel linear algebra used by CCFV
    const CommunicationType& cellCommunication() const
    {
        return cell_comm_;
    }

    ParallelIndexSet& cellIndexSet()
    {
        return cellCommunication().indexSet();
    }

    const ParallelIndexSet& cellIndexSet() const
    {
        return cellCommunication().indexSet();
    }

    RemoteIndices& cellRemoteIndices()
    {
        return cellCommunication().remoteIndices();
    }

    const RemoteIndices& cellRemoteIndices() const
    {
            return cellCommunication().remoteIndices();
    }
#endif

#ifdef HAVE_DUNE_ISTL
    /// \brief The type of the set of the attributes
    typedef Dune::OwnerOverlapCopyAttributeSet::AttributeSet AttributeSet;
#else
    /// \brief The type of the set of the attributes
    enum AttributeSet{owner, overlap, copy};
#endif

    /// \brief Get sorted active cell indices of numerical aquifer
    const std::vector<int>& sortedNumAquiferCells() const
    {
        return aquifer_cells_;
    }

private:

    /// \brief Adds entries to the parallel index set of the cells during grid construction
    void populateGlobalCellIndexSet();

#if HAVE_MPI

    /// \brief Gather data on a global grid representation.
    /// \param data A data handle for getting or setting the data
    /// \param global_view The view of the global grid (to gather the data on)
    /// \param distributed_view The view of the distributed grid.
    /// \tparam DataHandle The type of the data handle used.
    template<class DataHandle>
    void gatherData(DataHandle& data, CpGridData* global_view,
                    CpGridData* distributed_view);


    /// \brief Gather data specific to given codimension on a global grid representation.
    /// \param data A data handle for getting or setting the data
    /// \param global_view The view of the global grid (to gather the data on)
    /// \param distributed_view The view of the distributed grid.
    /// \tparam DataHandle The type of the data handle used.
    /// \tparam codim The codimension
    template<int codim, class DataHandle>
    void gatherCodimData(DataHandle& data, CpGridData* global_data,
                         CpGridData* distributed_data);

    /// \brief Scatter data from a global grid representation
    /// to a distributed representation of the same grid.
    /// \param data A data handle for getting or setting the data
    /// \param global_view The view of the global grid (to gather the data on)
    /// \param distributed_view The view of the distributed grid.
    /// \tparam DataHandle The type of the data handle used.
    template<class DataHandle>
    void scatterData(DataHandle& data, CpGridData* global_data,
                     CpGridData* distributed_data, const InterfaceMap& cell_inf,
                     const InterfaceMap& point_inf);

    /// \brief Scatter data specific to given codimension from a global grid representation
    /// to a distributed representation of the same grid.
    /// \param data A data handle for getting or setting the data
    /// \param global_view The view of the global grid (to gather the data on)
    /// \param distributed_view The view of the distributed grid.
    /// \tparam DataHandle The type of the data handle used.
    /// \tparam codim The codimension.
    template<int codim, class DataHandle>
    void scatterCodimData(DataHandle& data, CpGridData* global_data,
                          CpGridData* distributed_data);

    /// \brief Communicates data of a given codimension
    /// \tparam codim The codimension
    /// \tparam DataHandle The type of the data handle describing, gathering,
    ///  and gathering the data.
    /// \param DataHandle The data handle describing, gathering,
    ///  and gathering the data.
    /// \param dir The direction of the communication.
    /// \param interface The information about the communication interface
    template<int codim, class DataHandle>
    void communicateCodim(Entity2IndexDataHandle<DataHandle, codim>& data, CommunicationDirection dir,
                          const Interface& interface);

    /// \brief Communicates data of a given codimension
    /// \tparam codim The codimension
    /// \tparam DataHandle The type of the data handle describing, gathering,
    ///  and gathering the data.
    /// \param DataHandle The data handle describing, gathering,
    ///  and gathering the data.
    /// \param dir The direction of the communication.
    /// \param interface The information about the communication interface
    template<int codim, class DataHandle>
    void communicateCodim(Entity2IndexDataHandle<DataHandle, codim>& data, CommunicationDirection dir,
                          const InterfaceMap& interface);

#endif

    void computeGeometry(CpGrid& grid,
                         const DefaultGeometryPolicy&  globalGeometry,
                         const std::vector<int>& globalAquiferCells,
                         const OrientedEntityTable<0, 1>& globalCell2Faces,
                         DefaultGeometryPolicy& geometry,
                         std::vector<int>& aquiferCells,
                         const OrientedEntityTable<0, 1>& cell2Faces,
                         const std::vector< std::array<int,8> >& cell2Points);

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
    std::vector< std::array<int,8> >       cell_to_point_;
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
     * the number of cells present on the process and the content
     * by the mapping to the underlying global cartesian mesh..
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
    /** @brief The boundary ids. */
    cpgrid::EntityVariable<int, 1> unique_boundary_ids_;
    /** @brief The index set of the grid (level). */
    cpgrid::IndexSet* index_set_;
    /** @brief The internal local id set (not exported). */
    const cpgrid::IdSet* local_id_set_;
    /** @brief The global id set (used also as local id set). */
    LevelGlobalIdSet* global_id_set_;
    /** @brief The indicator of the partition type of the entities */
    PartitionTypeIndicator* partition_type_indicator_;   

    /// \brief The type of the collective communication.
    typedef MPIHelper::MPICommunicator MPICommunicator;
    #if DUNE_VERSION_NEWER(DUNE_GRID, 2, 7)
    using Communication = Dune::Communication<MPICommunicator>;
#else
    using CollectiveCommunication = Dune::CollectiveCommunication<MPICommunicator>;
    using Communication = Dune::CollectiveCommunication<MPICommunicator>;
#endif
    /// \brief Object for collective communication operations.
    Communication ccobj_;

    // Boundary information (optional).
    bool use_unique_boundary_ids_;

    /// This vector contains zcorn values from the initialization
    /// process where a CpGrid instance has been created from
    /// cornerpoint input zcorn and coord. During the initialization
    /// the zcorn values will typically be modified, and we retain a
    /// copy here to be able to create an EclipseGrid for output.
    std::vector<double> zcorn;

    /// \brief Sorted vector of aquifer cell indices.
    std::vector<int> aquifer_cells_;

#if HAVE_MPI

    /// \brief OwnerOverlap communication for cells
    CommunicationType cell_comm_;

    /// \brief Communication interface for the cells.
    std::tuple<Interface,Interface,Interface,Interface,Interface> cell_interfaces_;
    /*
    // code deactivated, because users cannot access face indices and therefore
    // communication on faces makes no sense!
    /// \brief Interface from interior and border to interior and border for the faces.
    std::tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>
    face_interfaces_;
    */
    /// \brief Interface from interior and border to interior and border for the faces.
    std::tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>
    point_interfaces_;

#endif

    // Return the geometry vector corresponding to the given codim.
    template <int codim>
    const EntityVariable<Geometry<3 - codim, 3>, codim>& geomVector() const
    {
        return geometry_.geomVector<codim>();
    }

    friend class Dune::CpGrid;
    template<int> friend class Entity;
    template<int> friend class EntityRep;
    friend class Intersection;
    friend class PartitionTypeIndicator;
};



#if HAVE_MPI

namespace
{
/// \brief Get a value from a tuple according to the interface type.
/// \tparam T The type of the values in the tuple.
/// \param iftype The interface type.
/// \param interfaces A tuple with the values order by interface type.
template<class T>
T& getInterface(InterfaceType iftype,
                std::tuple<T,T,T,T,T>& interfaces)
{
    switch(iftype)
    {
    case 0:
        return std::get<0>(interfaces);
    case 1:
        return std::get<1>(interfaces);
    case 2:
        return std::get<2>(interfaces);
    case 3:
        return std::get<3>(interfaces);
    case 4:
        return std::get<4>(interfaces);
    }
    OPM_THROW(std::runtime_error, "Invalid Interface type was used during communication");
}

} // end unnamed namespace

template<int codim, class DataHandle>
void CpGridData::communicateCodim(Entity2IndexDataHandle<DataHandle, codim>& data, CommunicationDirection dir,
                                  const Interface& interface)
{
    this->template communicateCodim<codim>(data, dir, interface.interfaces());
}

template<int codim, class DataHandle>
void CpGridData::communicateCodim(Entity2IndexDataHandle<DataHandle, codim>& data_wrapper, CommunicationDirection dir,
                                  const InterfaceMap& interface)
{
    Communicator comm(ccobj_, interface);

    if(dir==ForwardCommunication)
        comm.forward(data_wrapper);
    else
        comm.backward(data_wrapper);
}
#endif

template<class DataHandle>
void CpGridData::communicate(DataHandle& data, InterfaceType iftype,
                             CommunicationDirection dir)
{
#if HAVE_MPI
    if(data.contains(3,0))
    {
        Entity2IndexDataHandle<DataHandle, 0> data_wrapper(*this, data);
        communicateCodim<0>(data_wrapper, dir, getInterface(iftype, cell_interfaces_));
    }
    if(data.contains(3,3))
    {
        Entity2IndexDataHandle<DataHandle, 3> data_wrapper(*this, data);
        communicateCodim<3>(data_wrapper, dir, getInterface(iftype, point_interfaces_));
    }
#else
    // Suppress warnings for unused arguments.
    (void) data;
    (void) iftype;
    (void) dir;
#endif
}
}}

#if HAVE_MPI
#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/cpgrid/Indexsets.hpp>

namespace Dune {
namespace cpgrid {

namespace mover
{
template<class T>
class MoveBuffer
{
    friend class Dune::cpgrid::CpGridData;
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
    void moveData(const E& from, const E& to)
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
struct Mover<DataHandle,0> : public BaseMover<DataHandle>
{
    Mover<DataHandle,0>(DataHandle& data, CpGridData* gatherView,
                        CpGridData* scatterView)
    : BaseMover<DataHandle>(data), gatherView_(gatherView), scatterView_(scatterView)
    {}

    void operator()(std::size_t from_cell_index,std::size_t to_cell_index)
    {
        Entity<0> from_entity=Entity<0>(*gatherView_, from_cell_index, true);
        Entity<0> to_entity=Entity<0>(*scatterView_, to_cell_index, true);
        this->moveData(from_entity, to_entity);
    }
    CpGridData* gatherView_;
    CpGridData* scatterView_;
};

template<class DataHandle>
struct Mover<DataHandle,1> : public BaseMover<DataHandle>
{
    Mover<DataHandle,1>(DataHandle& data, CpGridData* gatherView,
                        CpGridData* scatterView)
    : BaseMover<DataHandle>(data), gatherView_(gatherView), scatterView_(scatterView)
    {}

    void operator()(std::size_t from_cell_index,std::size_t to_cell_index)
    {
        typedef typename OrientedEntityTable<0,1>::row_type row_type;
        EntityRep<0> from_cell=EntityRep<0>(from_cell_index, true);
        EntityRep<0> to_cell=EntityRep<0>(to_cell_index, true);
        OrientedEntityTable<0,1>& table = gatherView_->cell_to_face_;
        row_type from_faces=table.operator[](from_cell);
        row_type to_faces=scatterView_->cell_to_face_[to_cell];

        for(int i=0; i<from_faces.size(); ++i)
            this->moveData(from_faces[i], to_faces[i]);
    }
    CpGridData *gatherView_;
    CpGridData *scatterView_;
};

template<class DataHandle>
struct Mover<DataHandle,3> : public BaseMover<DataHandle>
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
            this->moveData(Entity<3>(*gatherView_, from_cell_points[i], true),
                           Entity<3>(*scatterView_, to_cell_points[i], true));
        }
    }
    CpGridData* gatherView_;
    CpGridData* scatterView_;
};

} // end mover namespace

template<class DataHandle>
void CpGridData::scatterData(DataHandle& data, CpGridData* global_data,
                             CpGridData* distributed_data, const InterfaceMap& cell_inf,
                             const InterfaceMap& point_inf)
{
#if HAVE_MPI
    if(data.contains(3,0))
    {
        Entity2IndexDataHandle<DataHandle, 0> data_wrapper(*global_data, *distributed_data, data);
        communicateCodim<0>(data_wrapper, ForwardCommunication, cell_inf);
    }
    if(data.contains(3,3))
    {
        Entity2IndexDataHandle<DataHandle, 3> data_wrapper(*global_data, *distributed_data, data);
        communicateCodim<3>(data_wrapper, ForwardCommunication, point_inf);
    }
#endif
}

template<int codim, class DataHandle>
void CpGridData::scatterCodimData(DataHandle& data, CpGridData* global_data,
                          CpGridData* distributed_data)
{
    CpGridData *gather_view, *scatter_view;
    gather_view=global_data;
    scatter_view=distributed_data;

    mover::Mover<DataHandle,codim> mover(data, gather_view, scatter_view);


    for(auto index=distributed_data->cellIndexSet().begin(),
            end = distributed_data->cellIndexSet().end();
        index!=end; ++index)
    {
        std::size_t from=index->global();
        std::size_t to=index->local();
        mover(from,to);
    }
}

namespace
{

template<int codim, class T, class F>
void visitInterior(CpGridData& distributed_data, T begin, T endit, F& func)
{
    for(T it=begin; it!=endit; ++it)
    {
        Entity<codim> entity(distributed_data, it-begin, true);
        PartitionType pt = entity.partitionType();
        if(pt==Dune::InteriorEntity)
        {
            func(*it, entity);
        }
    }
}

template<class DataHandle>
struct GlobalIndexSizeGatherer
{
    GlobalIndexSizeGatherer(DataHandle& data_,
                            std::vector<int>& ownedGlobalIndices_,
                            std::vector<int>& ownedSizes_)
        : data(data_), ownedGlobalIndices(ownedGlobalIndices_), ownedSizes(ownedSizes_)
    {}

    template<class T, class E>
    void operator()(T& i, E& entity)
    {
            ownedGlobalIndices.push_back(i);
            ownedSizes.push_back(data.size(entity));
    }
    DataHandle& data;
    std::vector<int>& ownedGlobalIndices;
    std::vector<int>& ownedSizes;
};

template<class DataHandle>
struct DataGatherer
{
    DataGatherer(mover::MoveBuffer<typename DataHandle::DataType>& buffer_,
                 DataHandle& data_)
        : buffer(buffer_), data(data_)
    {}

    template<class T, class E>
    void operator()(T& /* it */, E& entity)
    {
        data.gather(buffer, entity);
    }
    mover::MoveBuffer<typename DataHandle::DataType>& buffer;
    DataHandle& data;
};

}

template<class DataHandle>
void CpGridData::gatherData(DataHandle& data, CpGridData* global_data,
                            CpGridData* distributed_data)
{
#if HAVE_MPI
    if(data.contains(3,0))
       gatherCodimData<0>(data, global_data, distributed_data);
    if(data.contains(3,3))
       gatherCodimData<3>(data, global_data, distributed_data);
#endif
}

template<int codim, class DataHandle>
void CpGridData::gatherCodimData(DataHandle& data, CpGridData* global_data,
                                 CpGridData* distributed_data)
{
#if HAVE_MPI
    // Get the mapping to global index from  the global id set
    const std::vector<int>& mapping =
        distributed_data->global_id_set_->getMapping<codim>();

    // Get the global indices and data size for the entities whose data is
    // to be sent, i.e. the ones that we own.
    std::vector<int>         owned_global_indices;
    std::vector<int> owned_sizes;
    owned_global_indices.reserve(mapping.size());
    owned_sizes.reserve(mapping.size());

    GlobalIndexSizeGatherer<DataHandle> gisg(data, owned_global_indices, owned_sizes);
    visitInterior<codim>(*distributed_data, mapping.begin(), mapping.end(), gisg);

    // communicate the number of indices that each processor sends
    int no_indices=owned_sizes.size();
    // We will take the address of the first elemet for MPI_Allgather below.
    // Make sure the containers have such an element.
    if ( owned_global_indices.empty() )
        owned_global_indices.resize(1);
    if ( owned_sizes.empty() )
        owned_sizes.resize(1);
    std::vector<int> no_indices_to_recv(distributed_data->ccobj_.size());
    distributed_data->ccobj_.allgather(&no_indices, 1, &(no_indices_to_recv[0]));
    // compute size of the vector capable for receiving all indices
    // and allgather the global indices and the sizes.
    // calculate displacements
    std::vector<int> displ(distributed_data->ccobj_.size()+1, 0);
    std::transform(displ.begin(), displ.end()-1, no_indices_to_recv.begin(), displ.begin()+1,
                   std::plus<int>());
    int global_size=displ[displ.size()-1];//+no_indices_to_recv[displ.size()-1];
    std::vector<int>         global_indices(global_size);
    std::vector<int> global_sizes(global_size);
    MPI_Allgatherv(&(owned_global_indices[0]), no_indices, MPITraits<int>::getType(),
                   &(global_indices[0]), &(no_indices_to_recv[0]), &(displ[0]),
                   MPITraits<int>::getType(),
                   distributed_data->ccobj_);
    MPI_Allgatherv(&(owned_sizes[0]), no_indices, MPITraits<int>::getType(),
                   &(global_sizes[0]), &(no_indices_to_recv[0]), &(displ[0]),
                   MPITraits<int>::getType(),
                   distributed_data->ccobj_);
    std::vector<int>().swap(owned_global_indices); // free data for reuse.
    // Compute the number of data items to send
    std::vector<int> no_data_send(distributed_data->ccobj_.size());
    for(typename std::vector<int>::iterator begin=no_data_send.begin(),
            i=begin, end=no_data_send.end(); i!=end; ++i)
        *i = std::accumulate(global_sizes.begin()+displ[i-begin],
                            global_sizes.begin()+displ[i-begin+1], std::size_t());
    // free at least some memory that can be reused.
    std::vector<int>().swap(owned_sizes);
    // compute the displacements for receiving with allgatherv
    displ[0]=0;
    std::transform(displ.begin(), displ.end()-1, no_data_send.begin(), displ.begin()+1,
                   std::plus<std::size_t>());
    // Compute the number of data items we will receive
    int no_data_recv = displ[displ.size()-1];//+global_sizes[displ.size()-1];

    // Collect the data to send, gather it
    mover::MoveBuffer<typename DataHandle::DataType> local_data_buffer, global_data_buffer;
    if ( no_data_send[distributed_data->ccobj_.rank()] )
    {
        local_data_buffer.resize(no_data_send[distributed_data->ccobj_.rank()]);
    }
    else
    {
        local_data_buffer.resize(1);
    }
    global_data_buffer.resize(no_data_recv);

    DataGatherer<DataHandle> gatherer(local_data_buffer, data);
    visitInterior<codim>(*distributed_data, mapping.begin(), mapping.end(), gatherer);
    MPI_Allgatherv(&(local_data_buffer.buffer_[0]), no_data_send[distributed_data->ccobj_.rank()],
                   MPITraits<typename DataHandle::DataType>::getType(),
                   &(global_data_buffer.buffer_[0]), &(no_data_send[0]), &(displ[0]),
                   MPITraits<typename DataHandle::DataType>::getType(),
                   distributed_data->ccobj_);
    Entity2IndexDataHandle<DataHandle, codim> edata(*global_data, data);
    int offset=0;
    for(int i=0; i< codim; ++i)
        offset+=global_data->size(i);

    typename std::vector<int>::const_iterator s=global_sizes.begin();
    for(typename std::vector<int>::const_iterator i=global_indices.begin(),
            end=global_indices.end();
        i!=end; ++s, ++i)
    {
        edata.scatter(global_data_buffer, *i-offset, *s);
    }
#endif
}

} // end namespace cpgrid
} // end namespace Dune

#endif

#endif
