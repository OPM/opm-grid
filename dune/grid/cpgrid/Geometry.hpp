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

#ifndef OPM_GEOMETRY_HEADER
#define OPM_GEOMETRY_HEADER

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/genericgeometry/geometrytraits.hh>
#include <dune/geometry/genericgeometry/matrixhelper.hh>
#include <opm/core/utility/ErrorMacros.hpp>
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
                assert(allcorners && corner_indices);
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
	    GlobalCoordinate global(const LocalCoordinate& local) const
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
                    JacobianTransposed JT = jacobianTransposed(x);
                    GlobalCoordinate z = global(x);
                    z -= y;
                    MatrixHelper<DuneCoordTraits<double> >::template xTRightInvA<3, 3>(JT, z, dx );
                    x -= dx;
                } while (dx.two_norm2() > epsilon*epsilon);
                return x;
	    }

            /// Equal to \sqrt{\det{J^T J}} where J is the Jacobian.
            /// J_{ij} = (dg_i/du_j)
            /// where g is the mapping from the reference domain,
            /// and {u_j} are the reference coordinates.
	    double integrationElement(const LocalCoordinate& local) const
	    {
		FieldMatrix<ctype, coorddimension, mydimension> Jt = jacobianTransposed(local);
                using namespace GenericGeometry;
                return MatrixHelper<DuneCoordTraits<double> >::template sqrtDetAAT<3, 3>(Jt);
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
                assert(allcorners_ && cor_idx_);
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
	    const FieldMatrix<ctype, mydimension, coorddimension>
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
	    const FieldMatrix<ctype, coorddimension, mydimension>
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

	    /// This method is meaningless for singular geometries.
	    const GlobalCoordinate& global(const LocalCoordinate&) const
	    {
		THROW("Geometry::global() meaningless on singular geometry.");
	    }

	    /// This method is meaningless for singular geometries.
	    LocalCoordinate local(const GlobalCoordinate&) const
	    {
		THROW("Geometry::local() meaningless on singular geometry.");
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

	    /// This method is meaningless for singular geometries.
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

	    /// This method is meaningless for singular geometries.
	    const FieldMatrix<ctype, mydimension, coorddimension>&
	    jacobianTransposed(const LocalCoordinate& local) const
	    {
		THROW("Meaningless to call jacobianTransposed() on singular geometries.");
	    }

	    /// This method is meaningless for singular geometries.
	    const FieldMatrix<ctype, coorddimension, mydimension>&
	    jacobianInverseTransposed(const LocalCoordinate& /*local*/) const
	    {
		THROW("Meaningless to call jacobianInverseTransposed() on singular geometries.");
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




        /// Specialization for 0 dimensional geometries, i.e. vertices.
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

	    /// Returns the position of the vertex.
	    const GlobalCoordinate& global(const LocalCoordinate&) const
	    {
		return pos_;
	    }

            /// Meaningless for the vertex geometry.
	    LocalCoordinate local(const GlobalCoordinate&) const
	    {
		THROW("Meaningless to call local() on singular geometries.");
	    }

            /// Returns 1 for the vertex geometry.
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
	    const FieldMatrix<ctype, mydimension, coorddimension>&
	    jacobianTransposed(const LocalCoordinate& local) const
	    {
		THROW("Meaningless to call jacobianTransposed() on singular geometries.");
	    }

	    /// This method is meaningless for singular geometries.
	    const FieldMatrix<ctype, coorddimension, mydimension>&
	    jacobianInverseTransposed(const LocalCoordinate& /*local*/) const
	    {
		THROW("Meaningless to call jacobianInverseTransposed() on singular geometries.");
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

#endif // OPM_GEOMETRY_HEADER
