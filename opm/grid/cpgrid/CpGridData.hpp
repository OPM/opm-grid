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

void check_global_refine(const Dune::CpGrid&,
                         const Dune::CpGrid&);

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
    
    friend
    void ::check_global_refine(const Dune::CpGrid&,
                               const Dune::CpGrid&);

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

private:
    /// @brief Compute amount of cells in each direction of a patch of cells. (Cartesian grid required).
    ///
    /// @param [in]  startIJK  Cartesian triplet index where the patch starts.
    /// @param [in]  endIJK    Cartesian triplet index where the patch ends.
    ///                        Last cell part of the lgr will be {endijk[0]-1, ... endIJK[2]-1}.



    ///
    /// @return patch_dim Patch dimension {#cells in x-direction, #cells in y-direction, #cells in z-direction}.
    const std::array<int,3> getPatchDim(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK) const
    {
        return {endIJK[0]-startIJK[0], endIJK[1]-startIJK[1], endIJK[2]-startIJK[2]};
    }

    /// @brief Compute corner, face, and cell indices of a patch of cells. (Cartesian grid required).
    ///
    /// @param [in]  startIJK  Cartesian triplet index where the patch starts.
    /// @param [in]  endIJK    Cartesian triplet index where the patch ends.
    ///                        Last cell part of the lgr will be {endijk[0]-1, ... endIJK[2]-1}.
    ///
    /// @return {patch_corners, patch_faces, patch_cells} Indices of corners, faces, and cells of the patch of cells.
    const std::array<std::vector<int>,3> getPatchGeomIndices(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK) const
    {
        // Get the patch dimension (total cells in each direction). Used to 'reserve vectors'.
        const std::array<int,3>& patch_dim = getPatchDim(startIJK, endIJK);
        // Get grid dimension (total cells in each direction).
        const std::array<int,3>& grid_dim = this -> logicalCartesianSize();
        /// PATCH CORNERS
        std::vector<int> patch_corners;
        patch_corners.reserve((patch_dim[0]+1)*(patch_dim[1]+1)*(patch_dim[2]+1));
        for (int j = startIJK[1]; j < endIJK[1]+1; ++j) {
            for (int i = startIJK[0]; i < endIJK[0]+1; ++i) {
                for (int k = startIJK[2]; k < endIJK[2]+1; ++k) {
                    patch_corners.push_back((j*(grid_dim[0]+1)*(grid_dim[2]+1)) + (i*(grid_dim[2]+1))+k);
                } // end i-for-loop
            } // end j-for-loop
        } // end k-for-loop
        /// PATCH FACES
        std::vector<int> patch_faces;
        patch_faces.reserve(((patch_dim[0]+1)*patch_dim[1]*patch_dim[2])     // i_patch_faces
                            + (patch_dim[0]*(patch_dim[1]+1)*patch_dim[2])   // j_patch_faces
                            + (patch_dim[0]*patch_dim[1]*(patch_dim[2]+1))); // k_patch_faces
        // I_FACES
        for (int j = startIJK[1]; j < endIJK[1]; ++j) {
            for (int i = startIJK[0]; i < endIJK[0]+1; ++i) {
                for (int k = startIJK[2]; k < endIJK[2]; ++k) {
                    patch_faces.push_back((j*(grid_dim[0]+1)*grid_dim[2]) +(i*grid_dim[2]) + k);
                } // end k-for-loop
            } // end i-for-loop
        } // end j-for-loop
        // J_FACES
        for (int j = startIJK[1]; j < endIJK[1]+1; ++j) {
            for (int i = startIJK[0]; i < endIJK[0]; ++i) {
                for (int k = startIJK[2]; k < endIJK[2]; ++k) {
                    patch_faces.push_back(((grid_dim[0]+1)*grid_dim[1]*grid_dim[2]) // i_grid_faces
                                          + (j*grid_dim[0]*grid_dim[2]) + (i*grid_dim[2]) + k);
                } // end k-for-loop
            } // end i-for-loop
        } // end j-for-loop
        // K_FACES
        for (int j = startIJK[1]; j < endIJK[1]; ++j) {
            for (int i = startIJK[0]; i < endIJK[0]; ++i) {
                for (int k = startIJK[2]; k < endIJK[2]+1; ++k) {
                    patch_faces.push_back((grid_dim[0]*(grid_dim[1]+1)*grid_dim[2]) //j_grid_faces
                                          + ((grid_dim[0]+1)*grid_dim[1]*grid_dim[2])          // i_grid_faces
                                          + (j*grid_dim[0]*(grid_dim[2]+1)) + (i*(grid_dim[2]+1))+ k);
                } // end k-for-loop
            } // end i-for-loop
        } // end j-for-loop
        /// PATCH CELLS
        std::vector<int> patch_cells;
        patch_cells.reserve(patch_dim[0]*patch_dim[1]*patch_dim[2]);
        for (int k = startIJK[2]; k < endIJK[2]; ++k) {
            for (int j = startIJK[1]; j < endIJK[1]; ++j) {
                for (int i = startIJK[0]; i < endIJK[0]; ++i) {
                    patch_cells.push_back((k*grid_dim[0]*grid_dim[1]) + (j*grid_dim[0]) +i);
                } // end i-for-loop
            } // end j-for-loop
        } // end k-for-loop
        return {patch_corners, patch_faces, patch_cells};
    }

    /// @brief Construct a 'fake cell (Geometry<3,3> object)' out of a patch of cells.(Cartesian grid required).
    ///
    /// cellifyPatch() builds a Geometry<3,3> object, 'a celliFIED patch', from a connected patch formed
    /// by the product of consecutive cells in each direction; selecting 8 corners of the patch boundary,
    /// computing center and volume.
    ///
    /// @param [in] startIJK                   Cartesian triplet index where the patch starts.
    /// @param [in] endIJK                     Cartesian triplet index where the patch ends.
    /// @param [in] patch_cells                Cell indices from the block-shaped patch.
    /// @param [out] cellifiedPatch_geometry   Required as an argument when creating a Geomtry<3,3> object.
    /// @param [out] cellifiedPatch_to_point   To store the 8 corners of the created cellifiedPatch.
    /// @param [out] allcorners_cellifiedPatch Required to build a Geometry<3,3> object.
    ///
    /// @return 'cellifiedPatchCell'         Geometry<3,3> object.
    const Geometry<3,3> cellifyPatch(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK,
                                     const std::vector<int>& patch_cells, DefaultGeometryPolicy& cellifiedPatch_geometry,
                                     std::array<int,8>& cellifiedPatch_to_point,
                                     std::array<int,8>& allcorners_cellifiedPatch) const
    {
        if (patch_cells.empty()){
            OPM_THROW(std::logic_error, "Empty patch. Cannot convert patch into cell.");
        }
        else{
            // Get grid dimension.
            const std::array<int,3>& grid_dim = this -> logicalCartesianSize();
            // Select 8 corners of the patch boundary to be the 8 corners of the 'cellified patch'.
            cellifiedPatch_to_point = { // Corner-index: (J*(grid_dim[0]+1)*(grid_dim[2]+1)) + (I*(grid_dim[2]+1)) +K
                // Index of corner '0' {startI, startJ, startK}
                (startIJK[1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + (startIJK[0]*(grid_dim[2]+1)) + startIJK[2],
                // Index of corner '1' '{endI, startJ, startK}'
                (startIJK[1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + (endIJK[0]*(grid_dim[2]+1)) + startIJK[2],
                // Index of corner '2' '{startI, endJ, startK}'
                (endIJK[1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + (startIJK[0]*(grid_dim[2]+1)) + startIJK[2],
                // Index of corner '3' '{endI, endJ, startK}'
                (endIJK[1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + (endIJK[0]*(grid_dim[2]+1)) + startIJK[2],
                // Index of corner '4' '{startI, startJ, endK}'
                (startIJK[1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + (startIJK[0]*(grid_dim[2]+1))+ endIJK[2],
                // Index of corner '5' '{endI, startJ, endK}'
                (startIJK[1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + (endIJK[0]*(grid_dim[2]+1)) + endIJK[2],
                // Index of corner '6' '{startI, endJ, endK}'
                (endIJK[1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + (startIJK[0]*(grid_dim[2]+1)) + endIJK[2],
                // Index of corner '7' {endI, endJ, endK}
                (endIJK[1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + (endIJK[0]*(grid_dim[2]+1)) + endIJK[2]};
            EntityVariableBase<cpgrid::Geometry<0,3>>& cellifiedPatch_corners =
                cellifiedPatch_geometry.geomVector(std::integral_constant<int,3>());
            cellifiedPatch_corners.resize(8);
            // Compute the center of the 'cellified patch' and its corners.
            Geometry<0,3>::GlobalCoordinate cellifiedPatch_center = {0., 0.,0.};
            for (int corn = 0; corn < 8; ++corn) {
                cellifiedPatch_center +=
                    (this -> geometry_.geomVector(std::integral_constant<int,3>()).get(cellifiedPatch_to_point[corn]).center())/8.;
                cellifiedPatch_corners[corn] =
                    this -> geometry_.geomVector(std::integral_constant<int,3>()).get(cellifiedPatch_to_point[corn]);
            }
            // Compute the volume of the 'cellified patch'.
            double cellifiedPatch_volume = 0.;
            for (const auto& idx : patch_cells) {
                cellifiedPatch_volume += (this -> geometry_.geomVector(std::integral_constant<int,0>())
                                          [EntityRep<0>(idx, true)]).volume();
            }
            // Indices of 'all the corners', in this case, 0-7 (required to construct a Geometry<3,3> object).
            allcorners_cellifiedPatch = {0,1,2,3,4,5,6,7};
            // Create a pointer to the first element of "cellfiedPatch_to_point" (required to construct a Geometry<3,3> object).
            const int* cellifiedPatch_indices_storage_ptr = &allcorners_cellifiedPatch[0];
            // Construct (and return) the Geometry<3,3> of the 'cellified patch'.
            return Geometry<3,3>(cellifiedPatch_center, cellifiedPatch_volume,
                                 cellifiedPatch_geometry.geomVector(std::integral_constant<int,3>()), cellifiedPatch_indices_storage_ptr);
        }
    }

public:
    /// @brief Refine a single cell and return a shared pointer of CpGridData type.
    ///
    /// refineSingleCell() takes a cell and refines it in a chosen amount of cells (per direction); creating the
    /// geometries, topological relations, etc. Stored in a CpGridData object. Additionally, containers for
    /// parent-to-new-born entities are buil, as well as, new-born-to-parent. Maps(<int,bool>) to detect parent
    /// faces or cells are also provided. (Cell with 6 faces required).
    ///
    /// @param [in] cells_per_dim                 Number of (refined) cells in each direction that each parent cell should be refined to.
    /// @param [in] parent_idx                    Parent cell index, cell to be refined.
    ///
    /// @return refined_grid_ptr             Shared pointer pointing at refined_grid.
    /// @return parent_to_refined_corners    For each corner of the parent cell, we store the index of the
    ///                                           refined corner that coincides with the old one.
    ///                                           We assume they are ordered 0,1,..7
    ///                                                              6---7
    ///                                                      2---3   |   | TOP FACE
    ///                                                      |   |   4---5
    ///                                                      0---1 BOTTOM FACE
    /// @return parent_to_children_faces/cell For each parent face/cell, we store its child-face/cell indices.
    ///                                            {parent face/cell index in coarse level, {indices of its children in refined level}}
    /// @return child_to_parent_faces/cells   {child index, parent index}
    /// @return isParent_faces/cells          Map with all the face/cell indices. True when the face/cell got refined.
    const std::tuple< const std::shared_ptr<CpGridData>,
                      const std::vector<std::array<int,2>>,                // parent_to_refined_corners(~boundary_old_to_new_corners)
                      const std::vector<std::tuple<int,std::vector<int>>>, // parent_to_children_faces (~boundary_old_to_new_faces)
                      const std::tuple<int, std::vector<int>>,             // parent_to_children_cells
                      const std::vector<std::array<int,2>>,                // child_to_parent_faces
                      const std::vector<std::array<int,2>>,                // child_to_parent_cells
                      const std::vector<int>,                              // isParent_faces
                      const std::vector<int>>                               // isParent_cells
    refineSingleCell(const std::array<int,3>& cells_per_dim, const int& parent_idx) const
    {
        // To store the LGR/refined-grid.
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
        const cpgrid::Geometry<3,3>& parent_cell = geometry_.geomVector(std::integral_constant<int,0>())[EntityRep<0>(parent_idx, true)];
        // Get parent cell corners.
        const std::array<int,8>& parent_to_point = this->cell_to_point_[parent_idx];
        if (parent_to_point.size() != 8){
            OPM_THROW(std::logic_error, "Cell is not a hexahedron. Cannot be refined (yet).");
        }
        // Refine parent cell
        parent_cell.refine(cells_per_dim, refined_geometries, refined_cell_to_point, refined_cell_to_face,
                           refined_face_to_point, refined_face_to_cell, refined_face_tags, refined_face_normals);
        const std::vector<std::array<int,2>>& parent_to_refined_corners{
            // corIdx (J*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (I*(cells_per_dim[2]+1)) +K
            // replacing parent-cell corner '0' {0,0,0}
            {parent_to_point[0], 0},
                // replacing parent-cell corner '1' {cells_per_dim[0], 0, 0}
            {parent_to_point[1], cells_per_dim[0]*(cells_per_dim[2]+1)},
                // replacing parent-cell corner '2' {0, cells_perd_dim[1], 0}
            {parent_to_point[2], cells_per_dim[1]*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)},
                // replacing parent-cell corner '3' {cells_per_dim[0], cells_per_dim[1], 0}
            {parent_to_point[3], (cells_per_dim[1]*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (cells_per_dim[0]*(cells_per_dim[2]+1))},
                // replacing parent-cell corner '4' {0, 0, cells_per_dim[2]}
            {parent_to_point[4], cells_per_dim[2]},
                // replacing parent-cell corner '5' {cells_per_dim[0], 0, cells_per_dim[2]}
            {parent_to_point[5], (cells_per_dim[0]*(cells_per_dim[2]+1)) + cells_per_dim[2]},
                // replacing parent-cell corner '6' {0, cells_per_dim[1], cells_per_dim[2]}
            {parent_to_point[6], (cells_per_dim[1]*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + cells_per_dim[2]},
                // replacing parent-cell corner '7' {cells_per_dim[0], cells_per_dim[1], cells_per_dim[2]}
            {parent_to_point[7], (cells_per_dim[1]*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (cells_per_dim[0]*(cells_per_dim[2]+1))
                    + cells_per_dim[2]}};
        // Get parent_cell_to_face = { {face, orientation}, {another face, its orientation}, ...}.
        const auto& parent_cell_to_face = (this-> cell_to_face_[EntityRep<0>(parent_idx, true)]);
        // To store relation old-face to new-born-faces (children faces).
        std::vector<std::tuple<int,std::vector<int>>>  parent_to_children_faces;
        parent_to_children_faces.reserve(6);
        // To store child-to-parent-face relation. Child-faces ordered with the criteria introduced in refine()(Geometry.hpp)K,I,Jfaces.
        std::vector<std::array<int,2>> child_to_parent_faces;
        child_to_parent_faces.reserve(refined_face_to_cell.size());
        // Auxiliary integers to simplify new-born-face-index notation.
        const int& k_faces = cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1);
        const int& i_faces = (cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2];
        // Populate parent_to_children_faces and child_to_parent_faces.
        for (const auto& face : parent_cell_to_face) {
            // Check face tag to identify the type of face (bottom, top, left, right, front, or back).
            auto& parent_face_tag = (this-> face_tag_[Dune::cpgrid::EntityRep<1>(face.index(), true)]);
            // To store the new born faces for each face.
            std::vector<int> children_faces; // Cannot reserve/resize "now", it depends of the type of face.
            // K_FACES
            if (parent_face_tag == face_tag::K_FACE) {
                children_faces.reserve(cells_per_dim[0]*cells_per_dim[1]);
                for (int j = 0; j < cells_per_dim[1]; ++j) {
                    for (int i = 0; i < cells_per_dim[0]; ++i) {
                        int child_face;
                        if (!face.orientation()) // false -> BOTTOM FACE -> k=0
                            child_face = (j*cells_per_dim[0]) + i;
                        else // true -> TOP FACE -> k=cells_per_dim[2]
                            child_face = (cells_per_dim[2]*cells_per_dim[0]*cells_per_dim[1]) +(j*cells_per_dim[0]) + i;
                        children_faces.push_back(child_face);
                        child_to_parent_faces.push_back({child_face, face.index()});
                    } // i-for-lopp
                } //j-for-loop
            } // if-K_FACE
            // I_FACES
            if (parent_face_tag == face_tag::I_FACE) {
                children_faces.reserve(cells_per_dim[1]*cells_per_dim[2]);
                for (int k = 0; k < cells_per_dim[2]; ++k) {
                    for (int j = 0; j < cells_per_dim[1]; ++j) {
                        int child_face;
                        if (!face.orientation()) // false -> LEFT FACE -> i=0
                            child_face = k_faces + (k*cells_per_dim[1]) + j;
                        else // true -> RIGHT FACE -> i=cells_per_dim[0]
                            child_face = k_faces + (cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j;
                        children_faces.push_back(child_face);
                        child_to_parent_faces.push_back({child_face, face.index()});
                    } // j-for-loop
                } // k-for-loop
            } // if-I_FACE
            // J_FACES
            if (parent_face_tag == face_tag::J_FACE) {
                children_faces.reserve(cells_per_dim[0]*cells_per_dim[2]);
                for (int i = 0; i < cells_per_dim[0]; ++i) {
                    for (int k = 0; k < cells_per_dim[2]; ++k) {
                        int child_face;
                        if (!face.orientation()) // false -> FRONT FACE -> j=0
                            child_face = k_faces + i_faces + (i*cells_per_dim[2]) + k;
                        else  // true -> BACK FACE -> j=cells_per_dim[1]
                            child_face = k_faces + i_faces  + (cells_per_dim[1]*cells_per_dim[0]*cells_per_dim[2])
                                + (i*cells_per_dim[2]) + k;
                        children_faces.push_back(child_face);
                        child_to_parent_faces.push_back({child_face, face.index()});
                    } // k-for-loop
                } // i-for-loop
            } // if-J_FACE
            parent_to_children_faces.push_back(std::make_tuple(face.index(), children_faces));
        }
        std::tuple<int, std::vector<int>> parent_to_children_cells;
        auto& [ parent_index, children_cells ] = parent_to_children_cells;
        children_cells.reserve(cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]);
        // To store the child to parent cell relation.
        std::vector<std::array<int,2>> child_to_parent_cell;
        child_to_parent_cell.reserve(cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]);
        // Populate children_cells and child_to_parent_cell.
        for (int cell = 0; cell < cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]; ++cell) {
            children_cells.push_back(cell);
            child_to_parent_cell.push_back({cell, parent_idx});
        }
        // To identify a parent face {index face, true/false}. True when it got refined.
        std::vector<int> isParent_faces(face_to_cell_.size(), false);
        // Rewrite the entries of the map for those faces that got refined.
        for (const auto& face : parent_cell_to_face) {
            isParent_faces[face.index()] = true;
        }
        // To identify the parent cell (only one in this case). {index cell, true/false}. True when it got refined.
        std::vector<int> isParent_cells(this-> size(0), false);
        // Rewrite the entry of the map for the single parent cell that got refined.
        isParent_cells[parent_idx] = true;
        return {refined_grid_ptr, parent_to_refined_corners, parent_to_children_faces, parent_to_children_cells,
            child_to_parent_faces, child_to_parent_cell, isParent_faces, isParent_cells};
    }

    /// @brief Refine a (connected block-shaped) patch of cells. Based on the patch, a Geometry<3,3> object is created and refined.
    ///
    /// @param [in] cells_per_dim            Number of (refined) cells in each direction that each parent cell should be refined to.
    /// @param [in] startIJK                 Cartesian triplet index where the patch starts.
    /// @param [in] endIJK                   Cartesian triplet index where the patch ends.
    ///                                      Last cell part of the lgr will be {endijk[0]-1, ... endIJK[2]-1}.
    ///
    /// @return refined_grid_ptr                   Shared pointer of CpGridData type, pointing at the refined_grid
    /// @return boundary_old_to_new_corners/faces  Corner/face indices on the patch-boundary associated with new-born-entity indices.
    /// @return parent_to_children_faces/cell      For each parent face/cell, we store its child-face/cell indices.
    ///                                            {parent face/cell index in coarse level, {indices of its children in refined level}}
    /// @return child_to_parent_faces/cells        {child index, parent index}
    /// @return isParent_faces/cells               Map with all the face/cell indices. True when the face/cell got refined.
    const std::tuple<std::shared_ptr<CpGridData>,
                     const std::vector<std::array<int,2>>,                // boundary_old_to_new_corners
                     const std::vector<std::tuple<int,std::vector<int>>>, // boundary_old_to_new_faces
                     const std::vector<std::tuple<int,std::vector<int>>>, // parent_to_children_faces
                     const std::vector<std::tuple<int,std::vector<int>>>, // parent_to_children_cell
                     const std::vector<std::array<int,2>>,                // child_to_parent_faces
                     const std::vector<std::array<int,2>>,                // child_to_parent_cells
                     const std::vector<int>,                            // isParent_faces
                     const std::vector<int>>                            // isParent_cells
    refinePatch(const std::array<int,3>& cells_per_dim, const std::array<int,3>& startIJK, const std::array<int,3>& endIJK) const
    {
        // Coarse grid dimension (amount of cells in each direction).
        const std::array<int,3>& grid_dim = this -> logicalCartesianSize();
        // Check that the grid is a Cartesian one.
        long unsigned int gXYZ = grid_dim[0]*grid_dim[1]*grid_dim[2];
        if (global_cell_.size() != gXYZ){
            OPM_THROW(std::logic_error, "Grid is not Cartesian. Patch cannot be refined.");
        }
        // To store LGR/refined-grid.
        std::shared_ptr<CpGridData> refined_grid_ptr = std::make_shared<CpGridData>(ccobj_);
        auto& refined_grid = *refined_grid_ptr;
        DefaultGeometryPolicy& refined_geometries = refined_grid.geometry_;
        std::vector<std::array<int,8>>& refined_cell_to_point = refined_grid.cell_to_point_;
        cpgrid::OrientedEntityTable<0,1>& refined_cell_to_face = refined_grid.cell_to_face_;
        Opm::SparseTable<int>& refined_face_to_point = refined_grid.face_to_point_;
        cpgrid::OrientedEntityTable<1,0>& refined_face_to_cell = refined_grid.face_to_cell_;
        cpgrid::EntityVariable<enum face_tag,1>& refined_face_tags = refined_grid.face_tag_;
        cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>,1>& refined_face_normals = refined_grid.face_normals_;
        // Patch dimension (amount of cells in each direction).
        const auto& patch_dim = getPatchDim(startIJK, endIJK);
        // If the patch contains only one cell, use refineSingleCell() to avoid unnecessary computations.
        if ((patch_dim[0] == 1) && (patch_dim[1] == 1) && (patch_dim[2] == 1)){
            auto [refined_grid_ptr0, parent_to_refined_corners,
                  parent_to_children_faces, parent_to_children_cells,child_to_parent_faces, child_to_parent_cell,
                  isParent_faces, isParent_cells] =
                this->refineSingleCell(cells_per_dim,(startIJK[2]*grid_dim[0]*grid_dim[1]) + (startIJK[1]*grid_dim[0]) + startIJK[0]);
            // When the patch contains only one cell:
            // - boundary_old_to_new_corners == parent_to_refined_corners.
            // - boundary_old_to_new_faces == parent_to_children_faces.
            // Fix the type of parent_to_children_cells to return it correctly.
            const std::vector<std::tuple<int, std::vector<int>>> parent_to_children_cells_vec = {parent_to_children_cells};
            return {refined_grid_ptr0, parent_to_refined_corners, parent_to_children_faces, parent_to_children_faces,
                parent_to_children_cells_vec, child_to_parent_faces, child_to_parent_cell, isParent_faces, isParent_cells};
        }
        // When the patch consists in more than one cell:
        else {
            const auto& [patch_corners, patch_faces, patch_cells] = getPatchGeomIndices(startIJK, endIJK);
            // Construct the Geometry of the cellified patch.
            DefaultGeometryPolicy cellified_patch_geometry;
            std::array<int,8> cellifiedPatch_to_point;
            std::array<int,8> allcorners_cellifiedPatch;
            cpgrid::Geometry<3,3> cellified_patch = this -> cellifyPatch(startIJK, endIJK, patch_cells, cellified_patch_geometry,
                                                                         cellifiedPatch_to_point, allcorners_cellifiedPatch);
            // Some integers to reduce notation later.
            const int& xfactor = cells_per_dim[0]*patch_dim[0];
            const int& yfactor = cells_per_dim[1]*patch_dim[1];
            const int& zfactor = cells_per_dim[2]*patch_dim[2];
            // Refine the "cellified_patch".
            cellified_patch.refine({xfactor, yfactor, zfactor}, refined_geometries, refined_cell_to_point, refined_cell_to_face,
                                   refined_face_to_point, refined_face_to_cell, refined_face_tags, refined_face_normals);
            // To store the relation between old-corner-indices and the equivalent new-born ones (laying on the patch boundary).
            std::vector<std::array<int,2>> boundary_old_to_new_corners;
            boundary_old_to_new_corners.reserve((2*(cells_per_dim[0]+1)*(cells_per_dim[2]+1))
                                                + (2*(cells_per_dim[1]-1)*(cells_per_dim[2]+1))
                                                + (2*(cells_per_dim[0]-1)*(cells_per_dim[1]-1)));
            for (int j = startIJK[1]; j < endIJK[1]+1; ++j) {
                for (int i = startIJK[0]; i < endIJK[0]+1; ++i) {
                    for (int k = startIJK[2]; k < endIJK[2]+1; ++k) {
                        if ( (j == startIJK[1]) || (j == endIJK[1]) ){ // Corners in the front/back of the patch.
                            boundary_old_to_new_corners.push_back({
                                    // Old corner index
                                    (j*(grid_dim[2]+1)*(grid_dim[0]+1)) + (i*(grid_dim[2]+1)) + k,
                                    // New-born corner index (equivalent corner).
                                    (cells_per_dim[1]*(j-startIJK[1])*(xfactor +1)*(zfactor +1))
                                    + (cells_per_dim[0]*(i-startIJK[0])*(zfactor +1))
                                    + (cells_per_dim[2]*(k-startIJK[2])) });
                        }
                        if ( (i == startIJK[0]) || (i == endIJK[0]) ) { // Corners in the left/right of the patch.
                            boundary_old_to_new_corners.push_back({
                                    // Old corner index.
                                    (j*(grid_dim[2]+1)*(grid_dim[0]+1)) + (i*(grid_dim[2]+1)) + k,
                                    // New-born corner index (equivalent corner).
                                    (cells_per_dim[1]*(j-startIJK[1])*(xfactor +1)*(zfactor +1))
                                    + (cells_per_dim[0]*(i-startIJK[0])*(zfactor +1))
                                    + (cells_per_dim[2]*(k-startIJK[2]))});
                        }
                        if ( (k == startIJK[2]) || (k == endIJK[2]) ) { // Corners in the bottom/top of the patch.
                            boundary_old_to_new_corners.push_back({
                                    // Old corner index.
                                    (j*(grid_dim[2]+1)*(grid_dim[0]+1)) + (i*(grid_dim[2]+1)) + k,
                                    // New-born corner index (equivalent corner)
                                    (cells_per_dim[1]*(j-startIJK[1])*(xfactor +1)*(zfactor +1))
                                    + (cells_per_dim[0]*(i-startIJK[0])*(zfactor +1))
                                    + (cells_per_dim[2]*(k-startIJK[2]))});
                        }
                    } // end k-for-loop
                } // end i-for-loop
            } // end j-for-loop
            // To store face-indices of faces on the boundary of the patch.
            std::vector<int> boundary_patch_faces;
            // Auxiliary integers to simplify notation.
            const int& bound_patch_faces = (2*patch_dim[1]*patch_dim[2]) + (patch_dim[0]*2*patch_dim[2]) + (patch_dim[0]*patch_dim[1]*2);
            boundary_patch_faces.reserve(bound_patch_faces);
            // To store relation between old-face-index and its new-born-face indices.
            std::vector<std::tuple<int, std::vector<int>>> boundary_old_to_new_faces; // {face index, its children-indices}
            boundary_old_to_new_faces.reserve(bound_patch_faces);
            // Auxiliary integers to simplify notation.
            const int& i_grid_faces =  (grid_dim[0]+1)*grid_dim[1]*grid_dim[2];
            const int& j_grid_faces =  grid_dim[0]*(grid_dim[1]+1)*grid_dim[2];
            // To store relation bewteen parent face and its children (all faces of the patch, not only the ones on the boundary).
            std::vector<std::tuple<int,std::vector<int>>> parent_to_children_faces;
            parent_to_children_faces.reserve(patch_faces.size());
            // To store relation child-face-index and its parent-face-index.
            std::vector<std::array<int,2>> child_to_parent_faces; // {child index (in 'level 1'), parent index (in 'level 0')}
            child_to_parent_faces.reserve(refined_face_to_cell.size());
            // Populate child_to_parent_faces, parent_to_children_faces, boundary_old_to_new_faces, boundary_faces.
            // I_FACES
            for (int j = startIJK[1]; j < endIJK[1]; ++j) {
                for (int i = startIJK[0]; i < endIJK[0]+1; ++i) {
                    for (int k = startIJK[2]; k < endIJK[2]; ++k) {
                        int face_idx = (j*grid_dim[1]*grid_dim[2]) + (i*grid_dim[2])+ k;
                        // To store new born faces, per face. CHILDREN-FACES ARE ORDERED AS IN refine(), Geometry.hpp
                        std::vector<int> children_list;  // I_FACE ikj (xzy-direction)
                        // l,m,n play the role of 'x,y,z-direction', lnm = fake ikj (how I_FACES are 'ordered' in refine())
                        for (int l = (i-startIJK[0])*cells_per_dim[0]; l < (i-startIJK[0]+1)*cells_per_dim[0]; ++l) {
                            for (int n = (k-startIJK[2])*cells_per_dim[2];n < (k-startIJK[2]+1)*cells_per_dim[2]; ++n) {
                                for (int m = (j-startIJK[1])*cells_per_dim[1]; m < (j-startIJK[1]+1)*cells_per_dim[1]; ++m) {
                                    children_list.push_back((xfactor*yfactor*(zfactor+1)) +(l*yfactor*zfactor) + (n*yfactor) + m);
                                    child_to_parent_faces.push_back({(xfactor*yfactor*(zfactor+1)) +(l*yfactor*zfactor)
                                            + (n*yfactor) + m, face_idx});
                                } // end m-for-loop
                            } // end n-for-loop
                        } // end l-for-loop
                        // Add parent information of each face to "parent_to_children_faces".
                        parent_to_children_faces.push_back(std::make_tuple(face_idx, children_list));
                        if ((i == startIJK[0]) || (i == endIJK[0])) { // Detecting if the face is on the patch boundary.
                            boundary_patch_faces.push_back(face_idx);
                            // Associate each old face on the boundary of the patch with the new born ones.
                            boundary_old_to_new_faces.push_back(std::make_tuple(face_idx, children_list));
                        }
                    } // end k-for-loop
                } // end i-for-loop
            } // end j-for-loop
            // J_FACES
            for (int j = startIJK[1]; j < endIJK[1]+1; ++j) {
                for (int i = startIJK[0]; i < endIJK[0]; ++i) {
                    for (int k = startIJK[2]; k < endIJK[2]; ++k) {
                        int face_idx = i_grid_faces + (j*grid_dim[0]*grid_dim[2]) + (i*grid_dim[2])+ k;
                        // To store new born faces, per face. CHILDREN FACES ARE ORDERED AS IN refine(), Geometry.hpp
                        std::vector<int> children_list;  // J_FACE jik (yxz-direction)
                        // l,m,n play the role of 'x,y,z-direction', mln = fake jik (how J_FACES are 'ordered' in refine())
                        for (int m = (j-startIJK[1])*cells_per_dim[1]; m < (j-startIJK[1]+1)*cells_per_dim[1]; ++m) {
                            for (int l = (i-startIJK[0])*cells_per_dim[0]; l < (i-startIJK[0]+1)*cells_per_dim[0]; ++l) {
                                for (int n = (k-startIJK[2])*cells_per_dim[2]; n < (k-startIJK[2]+1)*cells_per_dim[2]; ++n) {
                                    children_list.push_back((xfactor*yfactor*(zfactor+1)) + ((xfactor+1)*yfactor*zfactor)
                                                            + (m*xfactor*zfactor) + (l*zfactor)+n);
                                    child_to_parent_faces.push_back({(xfactor*yfactor*(zfactor+1)) + ((xfactor+1)*yfactor*zfactor)
                                            + (m*xfactor*zfactor) + (l*zfactor)+n, face_idx});
                                } // end n-for-loop
                            } // end l-for-loop
                        } // end m-for-loop
                        // Add parent information of each face to "parent_to_children_faces".
                        parent_to_children_faces.push_back(std::make_tuple(face_idx, children_list));
                        if ((j == startIJK[1]) || (j == endIJK[1])) { // Detecting if face is on the patch boundary.
                            boundary_patch_faces.push_back(face_idx);
                            // Associate each old face on the boundary of the patch with the new born ones.
                            boundary_old_to_new_faces.push_back(std::make_tuple(face_idx, children_list));
                        }
                    } // end k-for-loop
                } // end i-for-loop
            } // end j-for-loop
            // K_FACES
            for (int j = startIJK[1]; j < endIJK[1]; ++j) {
                for (int i = startIJK[0]; i < endIJK[0]; ++i) {
                    for (int k = startIJK[2]; k < endIJK[2]+1; ++k) {
                        int face_idx = i_grid_faces + j_grid_faces + (j*grid_dim[0]*grid_dim[2]) + (i*grid_dim[2])+ k;
                        // To store new born faces, per face. CHILDREN FACES ARE ORDERED AS IN refine(), Geometry.hpp
                        std::vector<int> children_list;  // K_FACE kji (zyx-direction)
                        // l,m,n play the role of 'x,y,z-direction', nml = fake kji (how K_FACES are 'ordered' in refine())
                        for (int n = (k-startIJK[2])*cells_per_dim[2]; n < (k-startIJK[2]+1)*cells_per_dim[2]; ++n) {
                            for (int m = (j-startIJK[1])*cells_per_dim[1]; m < (j-startIJK[1]+1)*cells_per_dim[1]; ++m) {
                                for (int l = (i-startIJK[0])*cells_per_dim[0]; l < (i-startIJK[0]+1)*cells_per_dim[0]; ++l) {
                                    children_list.push_back((n*xfactor*yfactor) + (m*xfactor)+ l);
                                    child_to_parent_faces.push_back({(n*xfactor*yfactor) + (m*xfactor)+ l, face_idx});
                                } // end m-for-loop
                            } // end n-for-loop
                        } // end l-for-loop
                          // Add parent information of each face to "parent_to_children_faces".
                        parent_to_children_faces.push_back(std::make_tuple(face_idx, children_list));
                        if ((k == startIJK[2]) || (k == endIJK[2])) { // Detecting if the face is on the patch boundary.
                            boundary_patch_faces.push_back(face_idx);
                            // Associate each old face on the boundary of the patch with the new born ones.
                            boundary_old_to_new_faces.push_back(std::make_tuple(face_idx, children_list));
                        }
                    } // end k-for-loop
                } // end i-for-loop
            } // end j-for-loop
            // To store the relation between parent cell and its new-born-cells.
            // {parent index (coarse grid), {child 0 index, child 1 index, ... (refined grid)}}
            std::vector<std::tuple<int,std::vector<int>>> parent_to_children_cells;
            parent_to_children_cells.reserve(patch_dim[0]*patch_dim[1]*patch_dim[2]);
            // To store the relation between a new-born-cell and its parent cell.
            std::vector<std::array<int,2>> child_to_parent_cells; // {child index (refined grid), parent cell index (coarse grid)}
            child_to_parent_cells.reserve(xfactor*yfactor*zfactor);
            for (int k = 0; k < grid_dim[2]; ++k) {
                for (int j = 0; j < grid_dim[1]; ++j) {
                    for (int i = 0; i < grid_dim[0]; ++i) {
                        int cell_idx = (k*grid_dim[0]*grid_dim[1]) + (j*grid_dim[0]) +i;
                        std::vector<int> children_list;
                        if ( (i > startIJK[0]-1) && (i < endIJK[0]) && (j > startIJK[1]-1) && (j < endIJK[1])
                             && (k > startIJK[2]-1) && (k < endIJK[2])) {
                            for (int n = (k-startIJK[2])*cells_per_dim[2]; n < (k-startIJK[2]+1)*cells_per_dim[2]; ++n) {
                                for (int m = (j-startIJK[1])*cells_per_dim[1]; m < (j-startIJK[1]+1)*cells_per_dim[1]; ++m) {
                                    for (int l = (i-startIJK[0])*cells_per_dim[0]; l < (i-startIJK[0]+1)*cells_per_dim[0]; ++l) {
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
            // To identify parent faces. {index face, true/false}. True when a face got refined.
            std::vector<int> isParent_faces(this->face_to_cell_.size(), false);
            // Rewrite the entries of the map for those faces that got refined.
            for (const auto& face : patch_faces) {
                isParent_faces[face] = true;
            }
            // To identify parent cells. {index cell, true/false}. True when the cell got refined.
            std::vector<int> isParent_cells(this->size(0), false);
            // Rewrite the entries of the map for those cells that got refined.
            for (const auto& cell : patch_cells) {
                isParent_cells[cell] = true;
            }
            return {refined_grid_ptr, boundary_old_to_new_corners, boundary_old_to_new_faces, parent_to_children_faces,
                parent_to_children_cells, child_to_parent_faces, child_to_parent_cells, isParent_faces, isParent_cells};
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
