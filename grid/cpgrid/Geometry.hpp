//===========================================================================
//
// File: Geometry.hpp
//
// Created: Fri May 29 23:29:24 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010, 2011 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010, 2011 Statoil ASA.

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

#ifndef OPENRS_GEOMETRY_HEADER
#define OPENRS_GEOMETRY_HEADER

#include <dune/grid/genericgeometry/geometrytraits.hh>
#include <dune/grid/genericgeometry/matrixhelper.hh>
#include <boost/static_assert.hpp>

namespace Dune
{
    namespace cpgrid
    {

	/// This class encapsulates geometry for both vertices,
        /// intersections and cells.  The main template is empty,
        /// the actual dim == 3 (cell), dim == 2 (intersection)
        /// and dim == 0 (vertex) cases have specializations.
	/// For vertices and cells we use the cube type, and provide
        /// constant (vertex) or trilinear (cell) mappings.
	/// For intersections, we use the singular geometry type
        /// (None), and provide no mappings.
	template <int mydim, int cdim, class GridImp> // GridImp arg never used
	class Geometry
	{
        };




        /// Specialization for 3-dimensional geometries, i.e. cells.
	template <int cdim, class GridImp> // GridImp arg never used
	class Geometry<3, cdim, GridImp>
	{
	    BOOST_STATIC_ASSERT(cdim == 3);
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
	    typedef FieldMatrix< ctype, coorddimension, mydimension > 	Jacobian;
            /// Type of transposed Jacobian matrix.
	    typedef FieldMatrix< ctype, mydimension, coorddimension > 	JacobianTransposed;

	    /// @brief Construct from centroid, volume (1- and 0-moments) and
            ///        corners.
	    /// @param pos the centroid of the entity
            /// @param vol the volume(area) of the entity
	    /// @param allcorners array of all corner positions in the grid
            /// @param corner_indices array of 8 indices into allcorners array. The
            ///                       indices must be given in lexicographical order
            ///                       by (kji), i.e. i running fastest.
	    Geometry(const GlobalCoordinate& pos,
		     ctype vol,
		     const GlobalCoordinate* allcorners,
		     const int* corner_indices)
		: pos_(pos), vol_(vol), allcorners_(allcorners), cor_idx_(corner_indices)
	    {
                ASSERT(allcorners && corner_indices);
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
	    const GlobalCoordinate& global(const LocalCoordinate& local) const
	    {
                BOOST_STATIC_ASSERT(mydimension == 3);
                BOOST_STATIC_ASSERT(coorddimension == 3);
                // uvw = { (1-u, 1-v, 1-w), (u, v, w) }
                LocalCoordinate uvw[2] = { LocalCoordinate(1.0), local };
                uvw[0] -= local;
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
                BOOST_STATIC_ASSERT(mydimension == 3);
                BOOST_STATIC_ASSERT(coorddimension == 3);
                // This code is modified from dune/grid/genericgeometry/mapping.hh
                // \todo: Implement direct computation.
                const ctype epsilon = 1e-12;
                const GenericReferenceElement< ctype , 3 > & refElement =
                    GenericReferenceElements< ctype, 3 >::general(type());
                LocalCoordinate x = refElement.position(0,0);
                LocalCoordinate dx;
                do {
                    using namespace GenericGeometry;
                    // DF^n dx^n = F^n, x^{n+1} -= dx^n
                    JacobianTransposed JT;
                    jacobianTransposed( x, JT );
                    GlobalCoordinate z = global(x);
                    z -= y;
                    MatrixHelper<DuneCoordTraits<double> >::template xTRightInvA<3, 3>(JT, z, dx );
                    x -= dx;
                } while (dx.two_norm2() > epsilon*epsilon);
	    }

            /// The determinant of the Jacobian.
            /// J_{ij} = (dg_i/du_j)
            /// where g is the mapping from the reference domain,
            /// and {u_j} are the reference coordinates.
	    double integrationElement(const LocalCoordinate& local) const
	    {
		FieldMatrix<ctype, coorddimension, mydimension> Jt = jacobianTransposed(local);
		return Jt.determinant();
	    }

	    /// Using the cube type for all entities now (cells and vertices),
	    /// but we use the singular type for intersections.
	    GeometryType type() const
	    {
		GeometryType t;
                t.makeCube(mydimension);
		return t;
	    }

	    /// The number of corners of this convex polytope.
	    /// Returning 8, since we treat all cells as hexahedral.
	    int corners() const
	    {
                return 8;
	    }

	    /// The 8 corners of the hexahedral base cell.
	    GlobalCoordinate corner(int cor) const
	    {
                ASSERT(allcorners_ && cor_idx_);
                return allcorners_[cor_idx_[cor]];
	    }

	    /// Cell volume.
	    ctype volume() const
	    {
		return vol_;
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
	    const FieldMatrix<ctype, mydimension, coorddimension>&
	    jacobianTransposed(const LocalCoordinate& local) const
	    {
                BOOST_STATIC_ASSERT(mydimension == 3);
                BOOST_STATIC_ASSERT(coorddimension == 3);
                // uvw = { (1-u, 1-v, 1-w), (u, v, w) }
                LocalCoordinate uvw[2] = { LocalCoordinate(1.0), local };
                uvw[0] -= local;
                // Access pattern for uvw matching ordering of corners.
                const int pat[8][3] = { { 0, 0, 0 },
                                        { 1, 0, 0 },
                                        { 0, 1, 0 },
                                        { 1, 1, 0 },
                                        { 0, 0, 1 },
                                        { 1, 0, 1 },
                                        { 0, 1, 1 },
                                        { 1, 1, 1 } };
		FieldMatrix<ctype, mydimension, coorddimension> Jt(0.0);
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
	    const FieldMatrix<ctype, coorddimension, mydimension>&
	    jacobianInverseTransposed(const LocalCoordinate& local) const
	    {
		FieldMatrix<ctype, coorddimension, mydimension> Jti = jacobianTransposed(local);
                Jti.invert();
		return Jti;
	    }

	    /// The mapping implemented by this geometry is not generally affine.
	    bool affine() const
	    {
		return false;
	    }

	private:
	    GlobalCoordinate pos_;
	    double vol_;
	    const GlobalCoordinate* allcorners_; // For dimension 3 only
	    const int* cor_idx_;               // For dimension 3 only
	};





        /// Specialization for 2 dimensional geometries, that is
        /// intersections (since codim 1 entities are not in CpGrid).
	template <int cdim, class GridImp> // GridImp arg never used
	class Geometry<2, cdim, GridImp>
	{
	    BOOST_STATIC_ASSERT(cdim == 3);
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
	    typedef FieldMatrix< ctype, coorddimension, mydimension > 	Jacobian;
            /// Type of transposed Jacobian matrix.
	    typedef FieldMatrix< ctype, mydimension, coorddimension > 	JacobianTransposed;

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

            /// Since our geometry type is None, this method should not be called.
	    const GlobalCoordinate& global(const LocalCoordinate&) const
	    {
		THROW("Geometry::global() meaningless on singular geometry.");
	    }

	    /// In spite of claiming to be a cube geomety, we do not
	    /// make a 1-1 mapping from the cell to the reference cube.
	    LocalCoordinate local(const GlobalCoordinate&) const
	    {
		LocalCoordinate dummy(0.0);
		return dummy;
	    }

	    /// @brief
	    /// @todo Doc me!
	    /// @param
	    /// @return
	    double integrationElement(const LocalCoordinate&) const
	    {
		return vol_;
	    }

	    /// Using the cube type for all entities now (cells and vertices),
	    /// but we use the singular type for intersections.
	    GeometryType type() const
	    {
		GeometryType t;
                t.makeNone(mydimension);
		return t;
	    }

	    /// The number of corners of this convex polytope.
            /// Since this geometry is singular, we have no corners as such.
	    int corners() const
	    {
                return 0;
	    }

	    /// The corner method requires the points, which we may not necessarily want to provide.
	    /// We will need it for visualization purposes though. For now we throw.
	    GlobalCoordinate corner(int cor) const
	    {
                THROW("Meaningless call to cpgrid::Geometry::corner(int): "
                      "singular geometry has no corners.");
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

	    /// @brief
	    /// @todo Doc me!
	    const FieldMatrix<ctype, mydimension, coorddimension>&
	    jacobianTransposed(const LocalCoordinate& local) const
	    {
		THROW("Meaningless to call jacobianTransposed() on singular geometries.");
	    }

	    /// @brief
	    /// @todo Doc me!
	    const FieldMatrix<ctype, coorddimension, mydimension>&
	    jacobianInverseTransposed(const LocalCoordinate& /*local*/) const
	    {
		THROW("Meaningless to call jacobianInverseTransposed() on singular geometries.");
	    }

	    /// The mapping implemented by this geometry is singular, and therefore not affine.
	    bool affine() const
	    {
		return false;
	    }

	private:
	    GlobalCoordinate pos_;
	    ctype vol_;
	};





	/// This class encapsulates geometry for both vertices, intersections and cells.
	/// For vertices and cells we use the cube type, but without providing nonsingular
	/// global() and local() mappings. However, we do provide corner[s]().
	/// For intersections, we use the singular type, and no corners().
	template <int cdim, class GridImp> // GridImp arg never used
	class Geometry<0, cdim, GridImp>
	{
	    BOOST_STATIC_ASSERT(cdim == 3);
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
	    typedef FieldMatrix< ctype, coorddimension, mydimension > 	Jacobian;
            /// Type of transposed Jacobian matrix.
	    typedef FieldMatrix< ctype, mydimension, coorddimension > 	JacobianTransposed;

	    /// @brief Construct from vertex position
	    /// @param pos the position of the vertex
	    Geometry(const GlobalCoordinate& pos)
		: pos_(pos)
	    {
	    }

            /// Default constructor, giving a non-valid geometry.
	    Geometry()
		: pos_(0.0)
	    {
	    }

	    /// In spite of claiming to be a cube geomety, we do not
	    /// make a 1-1 mapping from the reference cube to the cell.
	    const GlobalCoordinate& global(const LocalCoordinate&) const
	    {
		return pos_;
	    }

	    /// In spite of claiming to be a cube geomety, we do not
	    /// make a 1-1 mapping from the cell to the reference cube.
	    LocalCoordinate local(const GlobalCoordinate&) const
	    {
		LocalCoordinate dummy(0.0);
		return dummy;
	    }

	    /// @brief
	    /// @todo Doc me!
	    /// @param
	    /// @return
	    double integrationElement(const LocalCoordinate&) const
	    {
		return volume();
	    }

	    /// Using the cube type for vertices.
	    GeometryType type() const
	    {
		GeometryType t;
                t.makeCube(mydimension);
		return t;
	    }

            /// A vertex is defined by a single corner.
	    int corners() const
	    {
                return 1;
	    }

            /// Returns the single corner: the vertex itself.
	    GlobalCoordinate corner(int cor) const
	    {
                ASSERT(cor == 0);
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

	    /// @brief
	    /// @todo Doc me!
	    const FieldMatrix<ctype, mydimension, coorddimension>&
	    jacobianTransposed(const LocalCoordinate& local) const
	    {
		THROW("Meaningless to call jacobianTransposed() on singular geometries.");
		static FieldMatrix<ctype, mydimension, coorddimension> dummy;
		return dummy;
	    }

	    /// @brief
	    /// @todo Doc me!
	    const FieldMatrix<ctype, coorddimension, mydimension>&
	    jacobianInverseTransposed(const LocalCoordinate& /*local*/) const
	    {
		THROW("Meaningless to call jacobianInverseTransposed() on singular geometries.");
		static FieldMatrix<ctype, coorddimension, mydimension> dummy;
		return dummy;
	    }

	    /// The mapping implemented by this geometry is constant, therefore affine.
	    bool affine() const
	    {
		return true;
	    }

	private:
	    GlobalCoordinate pos_;
	};





    } // namespace cpgrid
} // namespace Dune

#endif // OPENRS_GEOMETRY_HEADER
