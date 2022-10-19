//===========================================================================
//
// File: Geometry.hpp
//
// Created: Fri May 29 23:29:24 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            B�rd Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010, 2011 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010, 2011, 2022 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

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

#ifndef OPM_GEOMETRY_HEADER
#define OPM_GEOMETRY_HEADER

#include <cmath>

// Warning suppression for Dune includes.
#include <opm/grid/utility/platform_dependent/disable_warnings.h>

#include <dune/common/version.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/geometry.hh>

#include <dune/geometry/type.hh>

#include <opm/grid/cpgrid/EntityRep.hpp>
#include <opm/grid/cpgrid/DefaultGeometryPolicy.hpp>
#include <opm/grid/cpgrid/OrientedEntityTable.hpp>
#include <opm/grid/common/Volumes.hpp>
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>
#include <opm/grid/utility/SparseTable.hpp>

#include <opm/grid/utility/ErrorMacros.hpp>

namespace Dune
{
    namespace cpgrid
    {

        /// This class encapsulates geometry for vertices,
        /// intersections, and cells. The main template is empty,
        /// the actual dim == 3 (cell), dim == 2 (intersection),
        /// and dim == 0 (vertex) cases have specializations.
        /// For vertices and cells we use the cube type, and provide
        /// constant (vertex) or trilinear (cell) mappings.
        /// For intersections, we use the singular geometry type
        /// (None), and provide no mappings.
        template <int mydim, int cdim>
        class Geometry
        {
        };




        /// Specialization for 0 dimensional geometries, i.e. vertices.
        template <int cdim> // GridImp arg never used
        class Geometry<0, cdim>
        {
            static_assert(cdim == 3, "");
        public:
            /// Dimension of underlying grid.
            enum { dimension = 3 };
            /// Dimension of domain space of \see global().
            enum { mydimension = 0};
            /// Dimension of range space of \see global().
            enum { coorddimension = cdim };
            /// World dimension of underlying grid.
            enum { dimensionworld = 3 };

            /// Coordinate element type.
            typedef double ctype;

            /// Domain type of \see global().
            typedef FieldVector<ctype, mydimension> LocalCoordinate;
            /// Range type of \see global().
            typedef FieldVector<ctype, coorddimension> GlobalCoordinate;

            /// Type of Jacobian matrix.
            typedef FieldMatrix< ctype, coorddimension, mydimension >         Jacobian;
            /// Type of inverse of Jacobian matrix.
            typedef FieldMatrix< ctype, coorddimension, mydimension >         JacobianInverse;
            /// Type of transposed Jacobian matrix.
            typedef FieldMatrix< ctype, mydimension, coorddimension >         JacobianTransposed;
            /// Type of the inverse of the transposed Jacobian matrix
            typedef FieldMatrix< ctype, coorddimension, mydimension >         JacobianInverseTransposed;


            /// @brief Construct from vertex position
            /// @param pos the position of the vertex
            Geometry(const GlobalCoordinate& pos)
                : pos_(pos)
            {
            }

            /// @brief Default constructor, giving a non-valid geometry.
            Geometry()
                : pos_(0.0)
            {
            }

            /// Returns the position of the vertex.
            const GlobalCoordinate& global(const LocalCoordinate&) const
            {
                return pos_;
            }

            /// Meaningless for the vertex geometry.
            LocalCoordinate local(const GlobalCoordinate&) const
            {
                // return 0 to make the geometry check happy.
                return LocalCoordinate(0.0);
            }

            /// Returns 1 for the vertex geometry.
            double integrationElement(const LocalCoordinate&) const
            {
                return volume();
            }

            /// Using the cube type for vertices.
            GeometryType type() const
            {
                return Dune::GeometryTypes::cube(mydimension);
            }

            /// A vertex is defined by a single corner.
            int corners() const
            {
                return 1;
            }

            /// Returns the single corner: the vertex itself.
            GlobalCoordinate corner(int cor) const
            {
                static_cast<void>(cor);
                assert(cor == 0);
                return pos_;
            }

            /// Volume of vertex is arbitrarily set to 1.
            ctype volume() const
            {
                return 1.0;
            }

            /// Returns the centroid of the geometry.
            const GlobalCoordinate& center() const
            {
                return pos_;
            }

            /// This method is meaningless for singular geometries.
            JacobianTransposed
            jacobianTransposed(const LocalCoordinate& /* local */) const
            {

                // Meaningless to call jacobianTransposed() on singular geometries. But we need to make DUNE happy.
                return {};
            }

            /// This method is meaningless for singular geometries.
            JacobianInverseTransposed
            jacobianInverseTransposed(const LocalCoordinate& /*local*/) const
            {
                // Meaningless to call jacobianInverseTransposed() on singular geometries. But we need to make DUNE happy.
                return {};
            }

            /// This method is meaningless for singular geometries.
            Jacobian
            jacobian(const LocalCoordinate& /*local*/) const
            {
                return {};
            }
            
            /// This method is meaningless for singular geometries.
            JacobianInverse
            jacobianInverse(const LocalCoordinate& /*local*/) const
            {
                return {};
            }
            
            /// The mapping implemented by this geometry is constant, therefore affine.
            bool affine() const
            {
                return true;
            }

        private:
            GlobalCoordinate pos_;
        };  // class Geometry<0,cdim>




        /// Specialization for 2 dimensional geometries, that is
        /// intersections (since codim 1 entities are not in CpGrid).
        template <int cdim> // GridImp arg never used
        class Geometry<2, cdim>
        {
            static_assert(cdim == 3, "");
        public:
            /// Dimension of underlying grid.
            enum { dimension = 3 };
            /// Dimension of domain space of \see global().
            enum { mydimension = 2 };
            /// Dimension of range space of \see global().
            enum { coorddimension = cdim };
            /// World dimension of underlying grid.
            enum { dimensionworld = 3 };

            /// Coordinate element type.
            typedef double ctype;

            /// Domain type of \see global().
            typedef FieldVector<ctype, mydimension> LocalCoordinate;
            /// Range type of \see global().
            typedef FieldVector<ctype, coorddimension> GlobalCoordinate;

            /// Type of Jacobian matrix.
            typedef FieldMatrix< ctype, coorddimension, mydimension >         Jacobian;
            /// Type of inverse of Jacobian matrix.
            typedef FieldMatrix< ctype, coorddimension, mydimension >         JacobianInverse;
            /// Type of transposed Jacobian matrix.
            typedef FieldMatrix< ctype, mydimension, coorddimension >         JacobianTransposed;
            /// Type of the inverse of the transposed Jacobian matrix
            typedef FieldMatrix< ctype, coorddimension, mydimension >         JacobianInverseTransposed;

            /// @brief Construct from centroid and volume (1- and 0-moments).
            /// @param pos the centroid of the entity
            /// @param vol the volume(area) of the entity
            Geometry(const GlobalCoordinate& pos,
                     ctype vol)
                : pos_(pos), vol_(vol)
            {
            }

            /// Default constructor, giving a non-valid geometry.
            Geometry()
                : pos_(0.0), vol_(0.0)
            {
            }

            /// This method is meaningless for singular geometries.
            const GlobalCoordinate& global(const LocalCoordinate&) const
            {
                OPM_THROW(std::runtime_error, "Geometry::global() meaningless on singular geometry.");
            }

            /// This method is meaningless for singular geometries.
            LocalCoordinate local(const GlobalCoordinate&) const
            {
                OPM_THROW(std::runtime_error, "Geometry::local() meaningless on singular geometry.");
            }

            /// For the singular geometry, we return a constant
            /// integration element equal to the volume.
            double integrationElement(const LocalCoordinate&) const
            {
                return vol_;
            }

            /// We use the singular type (None) for intersections.
            GeometryType type() const
            {
                return Dune::GeometryTypes::none(mydimension);
            }

            /// The number of corners of this convex polytope.
            /// Since this geometry is singular, we have no corners as such.
            int corners() const
            {
                return 0;
            }

            /// This method is meaningless for singular geometries.
            GlobalCoordinate corner(int /* cor */) const
            {
                // Meaningless call to cpgrid::Geometry::corner(int):
                //"singular geometry has no corners.
                // But the DUNE tests assume at least one corner.
                return GlobalCoordinate( 0.0 );
            }

            /// Volume (area, actually) of intersection.
            ctype volume() const
            {
                return vol_;
            }

            /// Returns the centroid of the geometry.
            const GlobalCoordinate& center() const
            {
                return pos_;
            }

            /// This method is meaningless for singular geometries.
            const FieldMatrix<ctype, mydimension, coorddimension>&
            jacobianTransposed(const LocalCoordinate& /* local */) const
            {
                OPM_THROW(std::runtime_error, "Meaningless to call jacobianTransposed() on singular geometries.");
            }

            /// This method is meaningless for singular geometries.
            const FieldMatrix<ctype, coorddimension, mydimension>&
            jacobianInverseTransposed(const LocalCoordinate& /*local*/) const
            {
                OPM_THROW(std::runtime_error, "Meaningless to call jacobianInverseTransposed() on singular geometries.");
            }

            /// @brief The jacobian.
            Jacobian
            jacobian(const LocalCoordinate& /*local*/) const
            {
                return jacobianTransposed({}).transposed();
            }
            
            /// @brief The inverse of the jacobian
            JacobianInverse
            jacobianInverse(const LocalCoordinate& /*local*/) const
            {
                return jacobianInverseTransposed({}).transposed();
            }

            /// Since integrationElement() is constant, returns true.
            bool affine() const
            {
                return true;
            }

        private:
            GlobalCoordinate pos_;
            ctype vol_;
        };





        /// Specialization for 3-dimensional geometries, i.e. cells.
        template <int cdim>
        class Geometry<3, cdim>
        {
            static_assert(cdim == 3, "");
        public:
            /// Dimension of underlying grid.
            enum { dimension = 3 };
            /// Dimension of domain space of \see global().
            enum { mydimension = 3 };
            /// Dimension of range space of \see global().
            enum { coorddimension = cdim };
            /// World dimension of underlying grid.
            enum { dimensionworld = 3 };

            /// Coordinate element type.
            typedef double ctype;

            /// Domain type of \see global().
            typedef FieldVector<ctype, mydimension> LocalCoordinate;
            /// Range type of \see global().
            typedef FieldVector<ctype, coorddimension> GlobalCoordinate;

            /// Type of Jacobian matrix.
            typedef FieldMatrix< ctype, coorddimension, mydimension >         Jacobian;
            /// Type of inverse of Jacobian matrix.
            typedef FieldMatrix< ctype, coorddimension, mydimension >         JacobianInverse;
            /// Type of transposed Jacobian matrix.
            typedef FieldMatrix< ctype, mydimension, coorddimension >         JacobianTransposed;
            /// Type of the inverse of the transposed Jacobian matrix
            typedef FieldMatrix< ctype, coorddimension, mydimension >         JacobianInverseTransposed;

            typedef Dune::Impl::FieldMatrixHelper< double >  MatrixHelperType;

            /// @brief Construct from center, volume (1- and 0-moments) and
            ///        corners.
            /// @param pos the centroid of the entity
            /// @param vol the volume(area) of the entity
            /// @param allcorners array of all corner positions in the grid
            /// @param corner_indices array of 8 indices into allcorners array. The
            ///                       indices must be given in lexicographical order
            ///                       by (kji), i.e. i running fastest.
            Geometry(const GlobalCoordinate& pos,
                     ctype vol,
                     const EntityVariable<cpgrid::Geometry<0, 3>, 3>& allcorners,
                     const int* corner_indices)
                : pos_(pos), vol_(vol), allcorners_(allcorners.data()), cor_idx_(corner_indices)
            {
                assert(allcorners_ && corner_indices);
            }

            /// @brief Construct from centroid and volume (1- and
            ///        0-moments).  Note that since corners are not
            ///        given, the geometry provides no mappings, and
            ///        some calls (corner(), global() etc.) will fail.
            ///        This possibly dangerous constructor is
            ///        available for the benefit of
            ///        CpGrid::readSintefLegacyFormat().
            /// @param pos the centroid of the entity
            /// @param vol the volume(area) of the entity
            Geometry(const GlobalCoordinate& pos,
                     ctype vol)
                : pos_(pos), vol_(vol)
            {
            }

            /// Default constructor, giving a non-valid geometry.
            Geometry()
                : pos_(0.0), vol_(0.0), allcorners_(0), cor_idx_(0)
            {
            }

            /// Provide a trilinear mapping.
            /// Note that this does not give a proper space-filling
            /// embedding of the cell complex in the general (faulted)
            /// case. We should therefore revisit this at some point.
            /// Map g from (local) reference domain to (global) cell
            GlobalCoordinate global(const LocalCoordinate& local_coord) const
            {
                static_assert(mydimension == 3, "");
                static_assert(coorddimension == 3, "");
                // uvw = { (1-u, 1-v, 1-w), (u, v, w) }
                LocalCoordinate uvw[2] = { LocalCoordinate(1.0), local_coord };
                uvw[0] -= local_coord;
                // Access pattern for uvw matching ordering of corners.
                const int pat[8][3] = { { 0, 0, 0 },
                                        { 1, 0, 0 },
                                        { 0, 1, 0 },
                                        { 1, 1, 0 },
                                        { 0, 0, 1 },
                                        { 1, 0, 1 },
                                        { 0, 1, 1 },
                                        { 1, 1, 1 } };
                GlobalCoordinate xyz(0.0);
                for (int i = 0; i < 8; ++i) {
                    GlobalCoordinate corner_contrib = corner(i);
                    double factor = 1.0;
                    for (int j = 0; j < 3; ++j) {
                        factor *= uvw[pat[i][j]][j];
                    }
                    corner_contrib *= factor;
                    xyz += corner_contrib;
                }
                return xyz;
            }

            /// Mapping from the cell to the reference domain.
            /// May be slow.
            LocalCoordinate local(const GlobalCoordinate& y) const
            {
                static_assert(mydimension == 3, "");
                static_assert(coorddimension == 3, "");
                // This code is modified from dune/grid/genericgeometry/mapping.hh
                // \todo: Implement direct computation.
                const ctype epsilon = 1e-12;
                auto refElement = Dune::ReferenceElements<ctype, 3>::cube();
                LocalCoordinate x = refElement.position(0,0);
                LocalCoordinate dx;
                do {
                    // DF^n dx^n = F^n, x^{n+1} -= dx^n
                    JacobianTransposed JT = jacobianTransposed(x);
                    GlobalCoordinate z = global(x);
                    z -= y;
                    MatrixHelperType::template xTRightInvA<3, 3>(JT, z, dx );
                    x -= dx;
                } while (dx.two_norm2() > epsilon*epsilon);
                return x;
            }

            /// Equal to \sqrt{\det{J^T J}} where J is the Jacobian.
            /// J_{ij} = (dg_i/du_j)
            /// where g is the mapping from the reference domain,
            /// and {u_j} are the reference coordinates.
            double integrationElement(const LocalCoordinate& local_coord) const
            {
                JacobianTransposed Jt = jacobianTransposed(local_coord);
                return MatrixHelperType::template sqrtDetAAT<3, 3>(Jt);
            }

            /// Using the cube type for all entities now (cells and vertices),
            /// but we use the singular type for intersections.
            GeometryType type() const
            {
                return Dune::GeometryTypes::cube(mydimension);
            }

            /// The number of corners of this convex polytope.
            /// Returning 8, since we treat all cells as hexahedral.
            int corners() const
            {
                return 8;
            }

            /// @brief Get the cor-th of 8 corners of the hexahedral base cell.
            GlobalCoordinate corner(int cor) const
            {
                assert(allcorners_ && cor_idx_);
                return allcorners_[cor_idx_[cor]].center();
            }

            /// Cell volume.
            ctype volume() const
            {
                return vol_;
            }

            void set_volume(ctype volume) {
                vol_ = volume;
            }

            /// Returns the centroid of the geometry.
            const GlobalCoordinate& center() const
            {
                return pos_;
            }

            /// @brief Jacobian transposed.
            /// J^T_{ij} = (dg_j/du_i)
            /// where g is the mapping from the reference domain,
            /// and {u_i} are the reference coordinates.
            /// g = g(u) = (g_1(u), g_2(u), g_3(u)), u=(u_1,u_2,u_3)
            /// g = map from (local) reference domain to global cell
            const JacobianTransposed
            jacobianTransposed(const LocalCoordinate& local_coord) const
            {
                static_assert(mydimension == 3, "");
                static_assert(coorddimension == 3, "");

                // uvw = { (1-u, 1-v, 1-w), (u, v, w) }
                LocalCoordinate uvw[2] = { LocalCoordinate(1.0), local_coord };
                uvw[0] -= local_coord;
                // Access pattern for uvw matching ordering of corners.
                const int pat[8][3] = { { 0, 0, 0 },
                                        { 1, 0, 0 },
                                        { 0, 1, 0 },
                                        { 1, 1, 0 },
                                        { 0, 0, 1 },
                                        { 1, 0, 1 },
                                        { 0, 1, 1 },
                                        { 1, 1, 1 } };
                JacobianTransposed  Jt(0.0);
                for (int i = 0; i < 8; ++i) {
                    for (int deriv = 0; deriv < 3; ++deriv) {
                        // This part contributing to dg/du_{deriv}
                        double factor = 1.0;
                        for (int j = 0; j < 3; ++j) {
                            factor *= (j != deriv) ? uvw[pat[i][j]][j]
                                : (pat[i][j] == 0 ? -1.0 : 1.0);
                        }
                        GlobalCoordinate corner_contrib = corner(i);
                        corner_contrib *= factor;
                        Jt[deriv] += corner_contrib; // using FieldMatrix row access.
                    }
                }
                return Jt;
            }

            /// @brief Inverse of Jacobian transposed. \see jacobianTransposed().
            const JacobianInverseTransposed
            jacobianInverseTransposed(const LocalCoordinate& local_coord) const
            {
                JacobianInverseTransposed Jti = jacobianTransposed(local_coord);
                Jti.invert();
                return Jti;
            }

            /// @brief The jacobian.
            Jacobian
            jacobian(const LocalCoordinate& local_coord) const
            {
                return jacobianTransposed(local_coord).transposed();
            }
            
            /// @brief The inverse of the jacobian
            JacobianInverse
            jacobianInverse(const LocalCoordinate& local_coord) const
            {
                return jacobianInverseTransposed(local_coord).transposed();
            }
            
            /// The mapping implemented by this geometry is not generally affine.
            bool affine() const
            {
                return false;
            }

            /**
             * @brief Refine a single cell with regular intervals.
             *
             * For each cell to be created, storage must be passed for its corners and the indices. That storage
             * must be externally managed, since the newly created geometry structures only store pointers and do
             * not free them on destruction.
             *
             * @param cells_per_dim                                The number of sub-cells in each direction,
             * @param all_geom                                     Geometry Policy for the refined geometries. Those will be added there.
             * @param cell_to_face                                 Mapping from cell to oriented faces.
             *                                                     later be used to construct inverse map
             *                                                     with makeInverseRelation.
             * @param face_to_point                                Map from face to its points.
             * @param face_to_cell                                 Map from face to its neighboring cells.
             *
             * @param global_refined_cell8corners_indices_storage  A vector to store the indices of the 8 corners of each new cell.
             */
            typedef Dune::FieldVector<double,3> PointType;
            void refine(const std::array<int,3>& cells_per_dim,
                        DefaultGeometryPolicy& all_geom,
                        std::vector<std::array<int,8>>&  global_refined_cell8corners_indices_storage,
                        cpgrid::OrientedEntityTable<0,1>& cell_to_face,
                        Opm::SparseTable<int>& face_to_point,
                        cpgrid::OrientedEntityTable<1,0>& face_to_cell,
                        cpgrid::EntityVariable<enum face_tag, 1>& global_refined_face_tags,
                        cpgrid::SignedEntityVariable<PointType, 1>& global_refined_face_normals)
            {
                EntityVariableBase<cpgrid::Geometry<0,3>>& global_refined_corners =
                    all_geom.geomVector(std::integral_constant<int,3>());
                EntityVariableBase<cpgrid::Geometry<2,3>>& global_refined_faces =
                    all_geom.geomVector(std::integral_constant<int,1>());
                EntityVariableBase<cpgrid::Geometry<3,3>>& global_refined_cells =
                    all_geom.geomVector(std::integral_constant<int,0>());
                EntityVariableBase<enum face_tag>& mutable_face_tags = global_refined_face_tags;
                EntityVariableBase<PointType>& mutable_face_normals = global_refined_face_normals;
                // @todo CHECK PointType definition/construction.

                /// --- GLOBAL REFINED CORNERS ---
                // The strategy is to compute the local refined corners
                // of the unit/reference cube, and then apply the map global().
                // Determine the size of the vector containing all the corners
                // of all the global refined cells (children cells).
                global_refined_corners.resize((cells_per_dim[0] + 1) *(cells_per_dim[1] + 1) * (cells_per_dim[2] + 1));
                // The nummbering starts at the botton, so k=0 (z-axis), and j=0 (y-axis), i=0 (x-axis).
                // Then, increasing k ('going up'), followed by increasing i ('going right->'),
                // and finally, increasing j ('going back'). This order criteria for corners
                // 'Up [increasing k]- Right [incresing i]- Back [increasing j]'
                // is consistant with cpgrid numbering.
                for (int j = 0; j < cells_per_dim[1] + 1; ++j) {
                    for (int i = 0; i < cells_per_dim[0] + 1; ++i) {
                        for (int k = 0; k < cells_per_dim[2] + 1; ++k) {
                            // Compute the index of each global refined corner associated with 'jik'.
                            int global_refined_corner_idx =
                                (j*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (i*(cells_per_dim[2]+1)) +k;
                            // Compute the local refined corner of the unit/reference cube associated with 'jik'.
                            const LocalCoordinate& local_refined_corner = {
                                double(i)/cells_per_dim[0], double(j)/cells_per_dim[1], double(k)/cells_per_dim[2] };
                            // Compute the global refined corner 'jik' and add it in its corresponfing entry in "global_refined_corners".
                            global_refined_corners[global_refined_corner_idx] = Geometry<0, 3>(this->global(local_refined_corner));
                        } // end k-for-loop
                    } // end i-for-loop
                } // end j-for-loop
                /// --- END GLOBAL REFINED CORNERS ---

                /// --- GLOBAL REFINED FACES ---
                // We want to populate "global_refined_faces". The size of "global_refined_faces" is
                int global_refined_faces_size =
                    (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1)) // 'bottom/top faces'
                    + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2]) // 'left/right faces'
                    + (cells_per_dim[0]*(cells_per_dim[1]+1)*cells_per_dim[2]); // 'front/back faces'
                global_refined_faces.resize(global_refined_faces_size);
                global_refined_face_tags.resize(global_refined_faces_size);
                global_refined_face_normals.resize(global_refined_faces_size);
                //
                // To create a face as a Geometry<2,3> type object we need its CENTROID and its VOLUME(area).
                // We store the centroids/areas  in the following order:
                // - Bottom-top faces -> 3rd coordinate constant in each face.
                // - Left-right faces -> 1st coordinate constant in each face.
                // - Front-back faces -> 2nd coordinate constant in each face.
                //
                // Container to store, in each entry, the 4 indices of the 4 corners
                // of each global refined face (same size as "global_refined_faces").
                std::vector<std::array<int,4>> global_refined_face4corners_indices_storage; 
                global_refined_face4corners_indices_storage.resize(global_refined_faces_size);
                //
                // REFINED FACE AREAS
                // To compute the area of each face, we divide it in 4 triangles,
                // compute the area of those with "simplex_volume()", where the arguments
                // are the 3 corners of each triangle. Then, sum them up to get the area
                // of the global refined face.
                // For each face, we construct 4 triangles with
                // (1) centroid of the face,
                // (2) one of the edges of the face.
                //
                // A triangle has 3 edges. Once we choose a face to base a triangle on,
                // we choose an edge of that face as one of the edges of the triangle.
                // The other 2 edges are fixed, since the centroid of the face the triangle
                // is based on is one of its corners. That's why to identify
                // a triangle we only need two things:
                // (1) the face it's based on and
                // (2) one of the four edges of that face.
                //
                // For each face, we need
                // 1. index of the global refined face
                //    [available in "global_refined_face_indices"]
                //    [needed to access indices of the 4 edges of the face in "global_refined_face4edges_indices_storage"]
                // 2. centroid of the face (common corner of the 4 triangles based on that face).
                //    [available via "['face'].center()"
                // 3. container of 4 entries (the 4 edges of the face).
                //    Each entry consists in the 2 indices defining each edge of the face.
                //    [available in "global_refined_face4edges_indices_storage"].
                //
                // Populate
                // "global_refined_face4edges_indices_storage"
                // "global_refined_faces"
                //
                for (int constant_direction = 0; constant_direction < 3; ++constant_direction){
                    // adding %3 and r, we go through the 3 type of faces.
                    // r = 0 -> 3rd coordinate constant: l('k') < cells_per_dim[2]+1, m('j') < cells_per_dim[1], n('i') < cells_per_dim[0]
                    // r = 1 -> 1rt coordinate constant: l('i') < cells_per_dim[0]+1, m('k') < cells_per_dim[2], n('j') < cells_per_dim[1]
                    // r = 2 -> 2nd coordinate constant: l('j') < cells_per_dim[1]+1, m('i') < cells_per_dim[0], n('k') < cells_per_dim[2]
                    std::array<int,3> cells_per_dim_mixed = {
                        cells_per_dim[(2+constant_direction)%3],
                        cells_per_dim[(1+constant_direction)%3],
                        cells_per_dim[constant_direction % 3] };
                    for (int l = 0; l < cells_per_dim_mixed[0] + 1; ++l) {
                        for (int m = 0; m < cells_per_dim_mixed[1]; ++m) {
                            for (int n = 0; n < cells_per_dim_mixed[2]; ++n) {
                                // Compute the index of the face and its 4 corners.
                                auto [face_type, idx, global_refined_face4corners_indices,
                                      neighboring_cells_of_one_face, local_refined_face_centroid] =
                                    getIndicesFace(l, m, n, constant_direction, cells_per_dim);
                                // Add the tag to "face_tag_"
                                mutable_face_tags[idx]= face_type;
                                // Add 4 corner indices to "global_refined_face4corners_indices_storage".
                                global_refined_face4corners_indices_storage[idx] = global_refined_face4corners_indices;
                                // Add the 4 corners of the face to "face_to_point".
                                face_to_point.appendRow(global_refined_face4corners_indices.begin(),
                                                        global_refined_face4corners_indices.end());
                                // Add the neighboring cells of the face to "face_to_cell".
                                face_to_cell.appendRow(neighboring_cells_of_one_face.begin(),
                                                       neighboring_cells_of_one_face.end());
                                // Construct global face normal(s) (only one 'needed') and add it to "face_normals_"
                                // Construct two vectors in the face, e.g. difference of two conners with the centroid,
                                // then obtain an orthogonal vector to both of them. Finally, normalize.
                                // Auxuliary vectors on the face:
                                GlobalCoordinate face_vector0 =
                                    global_refined_corners[global_refined_face4corners_indices[0]].center()
                                    - global(local_refined_face_centroid);
                                GlobalCoordinate face_vector1 =
                                    global_refined_corners[global_refined_face4corners_indices[1]].center()
                                    - global(local_refined_face_centroid);
                                mutable_face_normals[idx] = {
                                    (face_vector0[1]*face_vector1[2]) -  (face_vector0[2]*face_vector1[1]),
                                    (face_vector0[2]*face_vector1[0]) -  (face_vector0[0]*face_vector1[2]),
                                    (face_vector0[0]*face_vector1[1]) -  (face_vector0[1]*face_vector1[0])};
                                mutable_face_normals[idx] /= mutable_face_normals[idx].two_norm();
                                // Construct "global_refined_face4edges_indices"
                                // with the {edge_indix[0], edge_index[1]} for each edge of the refined face.
                                std::vector<std::array<int,2>> global_refined_face4edges_indices = {
                                    { global_refined_face4corners_indices[0], global_refined_face4corners_indices[1]},
                                    { global_refined_face4corners_indices[0], global_refined_face4corners_indices[2]},
                                    { global_refined_face4corners_indices[1], global_refined_face4corners_indices[3]},
                                    { global_refined_face4corners_indices[2], global_refined_face4corners_indices[3]}};
                                // Calculate the AREA of each face of a global refined cell,
                                // by adding the 4 areas of the triangles partitioning each face.
                                double global_refined_face_area = 0.0;
                                for (int edge = 0; edge < 4; ++edge) {
                                    // Construction of each triangle on the current face with one
                                    // of its edges equal to "edge".
                                    const Geometry<0,3>::GlobalCoordinate trian_corners[3] = {
                                        global_refined_corners[global_refined_face4edges_indices[edge][0]].center(),
                                        global_refined_corners[global_refined_face4edges_indices[edge][1]].center(),
                                        this->global(local_refined_face_centroid)};
                                    global_refined_face_area += std::fabs(simplex_volume(trian_corners));
                                } // end edge-for-loop
                                //
                                //
                                // Construct the Geometry<2,3> of the global refined face.
                                global_refined_faces[idx] = Geometry<2,cdim>(this->global(local_refined_face_centroid),
                                                                             global_refined_face_area);
                                /// all_geom.geomVector(std::integral_constant<int,3>()), indices_storage_ptr);
                            } // end n-for-loop
                        } // end m-for-loop
                    } // end l-for-loop
                } // end r-for-loop
                /// --- END GLOBAL REFINED FACES ---

                /// --- GLOBAL REFINED CELLS ---
                // We need to populate "global_refined_cells"
                // "global_refined_cells"'s size is cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2].
                // To build each global refined cell, we need
                // 1. its global refined CENTER
                // 2. its VOLUME
                // 3. all global refined corners [available in "global_refined_corners"]
                // 4. indices of its 8 corners [available in "global_refined_corner_indices"]
                //
                global_refined_cells.resize(cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2]);
                // Vector to store, in each entry, the 8 indices of the 8 corners
                // of each global refined cell. Determine its size.
                global_refined_cell8corners_indices_storage.resize(cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2]);
                // The numembering starts with index 0 for the refined cell with corners
                // {0,0,0}, ...,{1/cells_per_dim[0], 1/cells_per_dim[1], 1/cells_per_dim[2]},
                // then the indices grow first picking the cells in the x-axis (Right, i), then y-axis (Back, j), and
                // finally, z-axis (Up, k).
                //
                // CENTERS
                // GLOBAL REFINED CELL CENTERS
                // The strategy is to compute the centers of the refined local
                // unit/reference cube, and then apply the map global().
                //
                // VOLUMES OF THE GLOBAL REFINED CELLS
                // REMARK: Each global refined 'cell' is a hexahedron since it may not be cube-shaped
                // since its a 'deformation' of unit/reference cube. We may use 'hexahedron' to refer
                // to the global refined cell in the computation of its volume.
                //
                // The strategy is to construct 24 tetrahedra in each hexahedron.
                // Each tetrahedron is built with
                // (1) the center of the hexahedron,
                // (2) the middle point of the face the tetrahedron is based on, and
                // (3) one of the edges of the face mentioned in 2.
                // Each face 'supports' 4 tetrahedra, and we have 6 faces per hexahedron, which
                // gives us the 24 tetrahedra per 'cell' (hexahedron).
                //
                // To compute the volume of each tetrahedron, we use "simplex_volume()" with
                // the 6 corners of the tetrahedron as arguments. Summing up the 24 volumes,
                // we get the volumne of the hexahedorn (global refined 'cell').
                //
                // Sum of all the volumes of all the (children) global refined cells.
                double sum_all_global_refined_cell_volumes = 0.0;
                //
                // For each (global refined 'cell') hexahedron, to create 24 tetrahedra and their volumes,
                // we introduce
                // Vol1. "hexa_face_0to5_indices" (needed to access face centroids).
                // Vol2. "hexa_face_centroids" (one of the 6 corners of all 4 tetrahedra based on that face).
                // Vol3.  the center of the global refined 'cell' (hexahedron)
                //       (common corner of the 24 tetrahedra).
                // Vol4. "tetra_edge_indices" indices of the 4x6 tetrahedra per 'cell',
                //        grouped by the face they are based on.
                // Then we construct and compute the volume of the 24 tetrahedra with mainly
                // "hexa_face_centroids" (Vol2.), global refined cell center (Vol3.), and "tetra_edge_indices" (Vol4.).
                //
                for (int k = 0; k < cells_per_dim[2]; ++k) {
                    for (int j = 0; j < cells_per_dim[1]; ++j) {
                        for (int i = 0; i < cells_per_dim[0]; ++i) {
                            // INDEX of the global refined cell associated with 'kji'.
                            int global_refined_cell_idx = (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i;
                            // 1. CENTER of the global refined cell associated with 'kji' (Vol3.)
                            // Compute the center of the local refined unit/reference cube associated with 'kji'.
                            const LocalCoordinate& local_refined_cell_center = {
                                (.5 + i)/cells_per_dim[0], (.5 + j)/cells_per_dim[1], (.5 + k)/cells_per_dim[2]};
                            // Obtain the global refined center with 'this->global(local_refined_cell_center)'.
                            // 2. VOLUME of the global refined 'kji' cell
                            double global_refined_cell_volume = 0.0; // (computed below!)
                            // 3. All Global refined corners ("global_refined_corners")
                            // 4. Indices of the 8 corners of the global refined cell associated with 'kji'.
                            std::array<int,8> global_refined_cell_corners_indices = { //
                                (j*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (i*(cells_per_dim[2]+1)) +k, // fake '0' {0,0,0}
                                (j*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + ((i+1)*(cells_per_dim[2]+1)) +k, // fake '1' {1,0,0}
                                ((j+1)*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (i*(cells_per_dim[2]+1)) +k, // fake '2' {0,1,0}
                                ((j+1)*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + ((i+1)*(cells_per_dim[2]+1)) +k, // fake '3' {1,1,0}
                                (j*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (i*(cells_per_dim[2]+1)) +k+1, // fake '4' {0,0,1}
                                (j*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + ((i+1)*(cells_per_dim[2]+1)) +k+1, // fake '5' {1,0,1}
                                ((j+1)*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (i*(cells_per_dim[2]+1)) +k+1, // fake '6' {0,1,1}
                                ((j+1)*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + ((i+1)*(cells_per_dim[2]+1)) +k+1 // fake '7' {1,1,1}
                            };
                            // Add this 8 corners to the corresponding entry of "global_refined_cell8corners_indices_storage"
                            global_refined_cell8corners_indices_storage[global_refined_cell_idx] = global_refined_cell_corners_indices;
                            //
                            // VOLUME HEXAHEDRON (GLOBAL REFINED 'CELL')
                            // Vol1. INDICES ('from 0 to 5') of the faces of the hexahedron (needed to access face centroids).
                            std::vector<int> hexa_face_0to5_indices = {
                                // index face '0' bottom
                                (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i,
                                // index face '1' front
                                (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1))
                                +  ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2])
                                + (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2]) + k,
                                // index face '2' left
                                (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1))
                                + (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j,
                                // index face '3' right
                                (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1))
                                + ((i+1)*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j,
                                // index face '4' back
                                (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1)) +
                                ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2])
                                + ((j+1)*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2]) + k,
                                // index face '5' top
                                ((k+1)*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i};
                            //
                            //  We add the 6 faces of the cell into "cell_to_face".
                            using cpgrid::EntityRep;
                            // First value is index, Second is orientation.
                            // Still have to find out what the orientation should be.
                            // right face ('3') outer normal points 'from left to right' -> orientation true
                            // back face ('4') outer normal points 'from front to back' -> orientation true
                            // top face ('5') outer normal points 'from bottom to top' -> orientation true
                            // (the other cases are false)
                            std::vector<cpgrid::EntityRep<1>> faces_of_one_cell = {
                                { hexa_face_0to5_indices[0], false}, {hexa_face_0to5_indices[1], false},
                                { hexa_face_0to5_indices[2], false}, {hexa_face_0to5_indices[3], true},
                                { hexa_face_0to5_indices[4], true}, {hexa_face_0to5_indices[5], true} };
                            cell_to_face.appendRow(faces_of_one_cell.begin(), faces_of_one_cell.end());
                            //
                            // Vol2. CENTROIDS of the faces of the hexahedron.
                            // (one of the 6 corners of all 4 tetrahedra based on that face).
                            std::vector<Geometry<0,3>::GlobalCoordinate> hexa_face_centroids;
                            for (auto& idx : hexa_face_0to5_indices) {
                                hexa_face_centroids.push_back(global_refined_faces[idx].center());
                            }
                            // Indices of the 4 edges of each face of the hexahedron.
                            // A tetrahedron has six edges. Once we choose a face to base a
                            // tetrahedron on, we choose an edge of that face as one of the
                            // edges of the tetrahedron. The other five edges are fixed, since
                            // the center of the hexahedron and the center of the face are
                            // the reminder 2 corners of the tetrahedron. That's why to identify
                            // a tetrahedron we only need two things:
                            // (1) the face it's based on and
                            // (2) one of the four edges of that face.
                            //
                            // Container with 6 entries, one per face. Each entry has the
                            // 4 indices of the 4 corners of each face.
                            std::vector<std::array<int,4>> global_refined_cell_face4corners_indices = {
                                global_refined_face4corners_indices_storage[hexa_face_0to5_indices[0]], // fake '{0,1,2,3}' bottom
                                global_refined_face4corners_indices_storage[hexa_face_0to5_indices[1]], // fake '{0,1,4,5}' front
                                global_refined_face4corners_indices_storage[hexa_face_0to5_indices[2]], // fake '{0,2,4,6}' left
                                global_refined_face4corners_indices_storage[hexa_face_0to5_indices[3]], // fake '{1,3,5,7}' right
                                global_refined_face4corners_indices_storage[hexa_face_0to5_indices[4]], // fake '{2,3,6,7}' back
                                global_refined_face4corners_indices_storage[hexa_face_0to5_indices[5]] };// fake '{4,5,6,7}' top
                            // Vol4. Container with indices of the edges of the 4 tetrahedra per face
                            // [according to description above]
                            std::vector<std::vector<std::array<int,2>>> tetra_edge_indices;
                            tetra_edge_indices.reserve(6);
                            for (auto& face_indices : global_refined_cell_face4corners_indices)
                            {
                                std::vector<std::array<int,2>> face4edges_indices = {
                                    { face_indices[0], face_indices[1]}, // fake '{0,1}'/'{4,5}'
                                    { face_indices[0], face_indices[2]}, // fake '{0,2}'/'{4,6}'
                                    { face_indices[1], face_indices[3]}, // fake '{1,3}'/'{5,7}'
                                    { face_indices[2], face_indices[3]} }; // fake '{2,3}'/'{6,7}'
                                tetra_edge_indices.push_back(face4edges_indices);
                            }
                            // Sum of the 24 volumes to get the volume of the hexahedron,
                            // stored in "global_refined_cell_volume".
                            // Calculate the volume of each hexahedron, by adding
                            // the 4 tetrahedra at each face (4x6 = 24 tetrahedra).
                            for (int face = 0; face < 6; ++face) {
                                for (int edge = 0; edge < 4; ++edge) {
                                    // Construction of each tetrahedron based on "face" with one
                                    // of its edges equal to "edge".
                                    const Geometry<0, 3>::GlobalCoordinate tetra_corners[4] = {
                                        global_refined_corners[tetra_edge_indices[face][edge][0]].center(),  // (see Vol4.)
                                        global_refined_corners[tetra_edge_indices[face][edge][1]].center(),  // (see Vol4.)
                                        hexa_face_centroids[face],  // (see Vol2.)
                                        // global_refined_cell_center
                                        this->global(local_refined_cell_center)};  // (see Vol3.)
                                    global_refined_cell_volume += std::fabs(simplex_volume(tetra_corners));
                                } // end edge-for-loop
                            } // end face-for-loop
                            // Add the volume of the hexahedron (global refined 'cell')
                            // to the container with of all volumes of all the refined cells.
                            sum_all_global_refined_cell_volumes += global_refined_cell_volume;
                            // Create a pointer to the first element of "global_refined_cell8corners_indices_storage"
                            // (required as the fourth argement to construct a Geometry<3,3> type object).
                            int* indices_storage_ptr = global_refined_cell8corners_indices_storage[global_refined_cell_idx].data();
                            // Construct the Geometry of the refined cell associated with 'kji'.
                            global_refined_cells[global_refined_cell_idx] =
                                Geometry<3,cdim>(this->global(local_refined_cell_center),
                                                 global_refined_cell_volume,
                                                 all_geom.geomVector(std::integral_constant<int,3>()),
                                                 indices_storage_ptr);
                        } // end i-for-loop
                    }  // end j-for-loop
                } // end k-for-loop
                // Rescale all volumes if the sum of volume of all the global refined 'cells' does not match the
                // volume of the 'parent cell'.
                // Compare the sum of all the volumes of all refined cells with 'parent cell' volume.
                if (std::fabs(sum_all_global_refined_cell_volumes - this->volume())
                    > std::numeric_limits<Geometry<3, cdim>::ctype>::epsilon()) {
                    Geometry<3, cdim>::ctype correction = this->volume() / sum_all_global_refined_cell_volumes;
                    for(auto& cell: global_refined_cells){
                        cell.vol_ *= correction;
                    }
                } // end if-statement
                /// --- END GLOBAL REFINED CELLS ---
            } /// --- END of refine()

        private:
            GlobalCoordinate pos_;
            double vol_;
            const cpgrid::Geometry<0, 3>* allcorners_; // For dimension 3 only
            const int* cor_idx_;               // For dimension 3 only
            
            // Auxiliary function to reduce "refine()"-code
            // "getIndicesFace" returns the index of the face, the 4 indices of the corners of the face,
            // and the local_refined_centroid. Each face is contant in one direction.
            // @param constant_coordinate -> 0,1,2 means constant in z,x,y respectively
            // @return Type of face: LEFT, BACK, or TOP
            //         Face index of a refined cell 'lmn' generated with "refine()".
            //         Four corner indices of the corners of the refined face 'lmn'.
            //         For each face, the (at most 2) neighboring cells (used in "face_to_cell").
            //         Local centroid of the face of refined cell 'lmn' of the unit cube. 
            std::tuple< enum face_tag, int,
                        std::array<int, 4>, std::vector<cpgrid::EntityRep<0>>,
                        LocalCoordinate>
            getIndicesFace(int l, int m, int n, int constant_direction, const std::array<int, 3>& cells_per_dim)
            {
                using cpgrid::EntityRep;
                std::vector<cpgrid::EntityRep<0>> neighboring_cells_of_one_face; // {index, orientation}
                switch(constant_direction) {
                case 0:  // {l,m,n} = {k,j,i}, constant in z-direction
                    // Orientation true when outer normal points 'from bottom to top'
                    // Orientation false when outer normal points 'from top to bottom'
                    if (l != 0) { 
                        neighboring_cells_of_one_face.push_back({((l-1)*cells_per_dim[0]*cells_per_dim[1])
                                + (m*cells_per_dim[0]) + n, true});
                    }
                    if (l != cells_per_dim[2]) {
                        neighboring_cells_of_one_face.push_back({ (l*cells_per_dim[0]*cells_per_dim[1])
                                + (m*cells_per_dim[0]) + n, false});
                    }
                    return { face_tag::K_FACE, (l*cells_per_dim[0]*cells_per_dim[1]) + (m*cells_per_dim[0]) + n,
                        {(m*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (n*(cells_per_dim[2]+1)) +l,
                        (m*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + ((n+1)*(cells_per_dim[2]+1)) +l,
                        ((m+1)*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (n*(cells_per_dim[2]+1)) +l,
                        ((m+1)*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + ((n+1)*(cells_per_dim[2]+1)) +l},
                        neighboring_cells_of_one_face,
                        {(.5 + n)/cells_per_dim[0], (.5 + m)/cells_per_dim[1], double(l)/cells_per_dim[2]}};                    
                case 1:  // {l,m,n} = {i,k,j}, constant in the x-direction
                    // Orientation true when outer normal points 'from left to right'
                    // Orientation false when outer normal points 'from right to left'
                    if (l != 0) {
                        neighboring_cells_of_one_face.push_back({(m*cells_per_dim[0]*cells_per_dim[1])
                                + (n*cells_per_dim[0]) +l-1, true});
                    }
                    if (l != cells_per_dim[0]) {
                        neighboring_cells_of_one_face.push_back({ (m*cells_per_dim[0]*cells_per_dim[1])
                                + (n*cells_per_dim[0]) + l, false});
                    }
                    return { face_tag::I_FACE, (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1))
                        + (l*cells_per_dim[1]*cells_per_dim[2]) + (m*cells_per_dim[1]) + n,
                        {(n*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (l*(cells_per_dim[2]+1)) +m,
                        ((n+1)*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (l*(cells_per_dim[2]+1)) +m,
                        (n*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (l*(cells_per_dim[2]+1)) +m+1,
                        ((n+1)*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (l*(cells_per_dim[2]+1)) +m+1},
                        neighboring_cells_of_one_face,
                        { double(l)/cells_per_dim[0], (.5 + n)/cells_per_dim[1], (.5 + m)/cells_per_dim[2]}};   
                case 2: // {l,m,n} = {j,i,k}, constant in the y-direction
                    // Orientation true when outer normal points 'from front to back'
                    // Orientation false when outer normal points 'from back to front'
                    if (l != 0) {
                        neighboring_cells_of_one_face.push_back({(n*cells_per_dim[0]*cells_per_dim[1])
                                + ((l-1)*cells_per_dim[0]) +m, true});
                    }
                    if (l != cells_per_dim[1]) {
                        neighboring_cells_of_one_face.push_back({(n*cells_per_dim[0]*cells_per_dim[1])
                                + (l*cells_per_dim[0]) + m, false});
                    }
                    return { face_tag::J_FACE, (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2] +1))
                        + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2])
                        + (l*cells_per_dim[0]*cells_per_dim[2]) + (m*cells_per_dim[2]) + n,
                        {(l*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (m*(cells_per_dim[2]+1)) +n,
                        (l*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + ((m+1)*(cells_per_dim[2]+1)) +n,
                        (l*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (m*(cells_per_dim[2]+1)) +n+1,
                        (l*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + ((m+1)*(cells_per_dim[2]+1)) +n+1},
                        neighboring_cells_of_one_face,
                        {(.5 + m)/cells_per_dim[0], double(l)/cells_per_dim[1], (.5 + n)/cells_per_dim[2]}};   
                default:
                    // Should never be reached, but prevents compiler warning
                    OPM_THROW(std::logic_error, "Unhandled dimension. This should never happen!");
                }
            }
        };
    } // namespace cpgrid

    template< int mydim, int cdim >
    auto referenceElement(const cpgrid::Geometry<mydim,cdim>& geo) -> decltype(referenceElement<double,mydim>(geo.type()))
    {
        return referenceElement<double,mydim>(geo.type());
    }

} // namespace Dune

#endif // OPM_GEOMETRY_HEADER
