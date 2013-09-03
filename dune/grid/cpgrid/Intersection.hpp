//===========================================================================
//
// File: Intersection.hpp
//
// Created: Tue Jun  9 11:17:13 2009
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
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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

#ifndef OPM_INTERSECTION_HEADER
#define OPM_INTERSECTION_HEADER




#include <dune/grid/common/gridenums.hh>

#include <opm/core/utility/ErrorMacros.hpp>

// The next statement is a layering violation: we only #include
// preprocess.h to get at its "enum face_tag" definition.  Enum
// face_tag is needed in method Intersection::boundaryId().  This hack
// is in dire need of a better solution!
#include <opm/core/grid/cpgpreprocess/preprocess.h>

#include "Entity.hpp"
#include "Geometry.hpp"
#include "OrientedEntityTable.hpp"

namespace Dune
{
    namespace cpgrid
    {

	/// @brief
	/// @todo Doc me!
	/// @tparam
        template <class GridType>
        class Intersection
        {
        public:
	    /// @brief
	    /// @todo Doc me!
	    enum { dimension = 3 };
	    enum { dimensionworld = 3 };
	    /// @brief
	    /// @todo Doc me!
	    typedef cpgrid::Entity<0, GridType> Entity;
	    typedef cpgrid::EntityPointer<0, GridType> EntityPointer;
// 	    typedef cpgrid::Geometry<2,3, GridType> Geometry;
// 	    typedef cpgrid::Geometry<2,3, GridType> LocalGeometry;
	    typedef typename GridType::template Codim<1>::Geometry Geometry;
	    typedef Geometry LocalGeometry;
	    typedef double ctype;
	    typedef FieldVector<ctype, 2> LocalCoordinate;
	    typedef FieldVector<ctype, 3> GlobalCoordinate;

	    /// @brief
	    /// @todo Doc me!
	    /// @param
            Intersection()
		: pgrid_(0),
		  index_(-1),
		  subindex_(-1),
		  faces_of_cell_(),
		  global_geom_(),
// 		  in_inside_geom_(),
		  nbcell_(-1), // Init to self, which is invalid.
		  is_on_boundary_(false)
            {
            }
	    /// @brief
	    /// @todo Doc me!
	    /// @param
            Intersection(const GridType& grid, EntityRep<0> cell, int subindex, bool update_now = true)
		: pgrid_(&grid),
		  index_(cell.index()),
		  subindex_(subindex),
		  faces_of_cell_(grid.cell_to_face_[cell]),
		  global_geom_(cpgrid::Entity<1, GridType>(grid, faces_of_cell_[subindex_]).geometry()),
// 		  in_inside_geom_(global_geom_.center()
// 				  - cpgrid::Entity<0, GridType>(grid, index_).geometry().center(),
// 				  global_geom_.volume()),
		  nbcell_(cell.index()), // Init to self, which is invalid.
		  is_on_boundary_(false)
            {
                assert(index_ >= 0);
		if (update_now) {
		    update();
		}
            }

	    /// @brief
	    /// @todo Doc me!
	    /// @param
	    /// @return
            bool operator==(const Intersection& other) const
            {
                return subindex_ == other.subindex_  &&  index_ == other.index_  &&  pgrid_ == other.pgrid_;
            }

	    /// @brief
	    /// @todo Doc me!
	    /// @param
	    /// @return
            bool operator!=(const Intersection& other) const
            {
                return !operator==(other);
            }

	    /// @brief
	    /// @todo Doc me!
	    /// @param
	    /// @return
            bool boundary() const
            {
		return is_on_boundary_;
            }

            /// Returns the boundary id of this intersection.
            int boundaryId() const
            {
                int ret = 0;
                if (boundary()) {
		    if (pgrid_->uniqueBoundaryIds()) {
			// Use the unique boundary ids.
			EntityRep<1> face = faces_of_cell_[subindex_];
			ret = pgrid_->unique_boundary_ids_[face];
		    } else {
			// Use the face tag based ids, i.e. 1-6 for i-, i+, j-, j+, k-, k+.
			typedef OrientedEntityTable<0,1>::ToType Face;
			const Face& f = faces_of_cell_[subindex_];
			const bool normal_is_in = !f.orientation();
			enum face_tag tag = pgrid_->face_tag_[f];

			switch (tag) {
			case LEFT:
			    //                   LEFT : RIGHT
			    ret = normal_is_in ? 1    : 2; // min(I) : max(I)
			    break;
			case BACK:
			    //                   BACK : FRONT
			    ret = normal_is_in ? 3    : 4; // min(J) : max(J)
			    break;
			case TOP:
			    // Note: TOP at min(K) as 'z' measures *depth*.
			    //                   TOP  : BOTTOM
			    ret = normal_is_in ? 5    : 6; // min(K) : max(K)
			    break;
			}
		    }
                }
                return ret;
            }

            /// Returns the boundary segment index of this intersection.
            int boundarySegmentIndex() const
            {
                // Since this is almost the same that we did for
                // 'unique boundary ids' we use those numbers, although since
                // they are 1-based and not 0-based we must be careful.
                if (!boundary()) {
                    OPM_THROW(std::runtime_error, "Cannot call boundarySegmentIndex() on non-boundaries.");
                }
                assert(!pgrid_->unique_boundary_ids_.empty());
                // Use the unique boundary ids (subtract 1).
                EntityRep<1> face = faces_of_cell_[subindex_];
                return pgrid_->unique_boundary_ids_[face] - 1;
            }

	    /// @brief
	    /// @todo Doc me!
	    /// @return
            bool neighbor() const
            {
                return !boundary();
            }

	    /// @brief
	    /// @todo Doc me!
	    /// @return
            EntityPointer inside() const
            {
                return EntityPointer(*pgrid_, index_, true);
            }

	    /// @brief
	    /// @todo Doc me!
	    /// @return
            EntityPointer outside() const
            {
                return EntityPointer(*pgrid_, nbcell(), true);
            }

	    /// @brief
	    /// @todo Doc me!
	    /// @return
	    bool conforming() const
	    {
		return boundary(); // I.e. we are assuming all nonconforming interior.
	    }

            // Geometrical information about this intersection in
            // local coordinates of the inside() entity.
	    /// @brief
	    /// @todo Doc me!
	    /// @return
            const LocalGeometry& geometryInInside() const
	    {
                OPM_THROW(std::runtime_error, "This intersection class does not support geometryInInside().");
// 		return in_inside_geom_;
	    }

            // Geometrical information about this intersection in
            // local coordinates of the outside() entity.
	    /// @brief
	    /// @todo Doc me!
	    /// @return
            const LocalGeometry& geometryInOutside() const
	    {
		if (boundary()) {
		    OPM_THROW(std::runtime_error, "Cannot access geometryInOutside(), intersection is at a boundary.");
		}
                OPM_THROW(std::runtime_error, "This intersection class does not support geometryInOutside().");
// 		return in_outside_geom_;
	    }

	    /// @brief
	    /// @todo Doc me!
	    /// @return
            const Geometry& geometry() const
            {
		return global_geom_;
            }

	    /// @brief
	    /// @todo Doc me!
	    /// @return
            GeometryType type() const
            {
                return geometry().type();
            }

            /// Local index of codim 1 entity in the inside() entity
            /// where intersection is contained in.
            int indexInInside() const
            {
                // Use the face tags to decide if an intersection is
                // on an x, y, or z face and use orientations to decide
                // if its (for example) an xmin or xmax face.
                typedef OrientedEntityTable<0,1>::ToType Face;
                const Face& f = faces_of_cell_[subindex_];
                const bool normal_is_in = !f.orientation();
                enum face_tag tag = pgrid_->face_tag_[f];
                switch (tag) {
                case LEFT:
                    return normal_is_in ? 0 : 1; // min(I) : max(I)
                case BACK:
                    return normal_is_in ? 2 : 3; // min(J) : max(J)
                case TOP:
                    return normal_is_in ? 4 : 5; // min(K) : max(K)
                default:
                    OPM_THROW(std::runtime_error, "Unhandled face tag: " << tag);
                }
            }

            /// Local index of codim 1 entity in outside() entity
            /// where intersection is contained in.
            int indexInOutside() const
            {
                int in_inside = indexInInside();
                return in_inside + ((in_inside % 2) ? -1 : 1);
            }

	    /// @brief
	    /// @todo Doc me!
	    /// @param
	    /// @return
            FieldVector<ctype, 3> outerNormal(const FieldVector<ctype, 2>&) const
            {
                return pgrid_->face_normals_[faces_of_cell_[subindex_]];
            }

	    /// @brief
	    /// @todo Doc me!
	    /// @param
	    /// @return
            FieldVector<ctype, 3> integrationOuterNormal(const FieldVector<ctype, 2>& unused) const
            {
                FieldVector<ctype, 3> n = pgrid_->face_normals_[faces_of_cell_[subindex_]];
                return n*=geometry().integrationElement(unused);
            }

	    /// @brief
	    /// @todo Doc me!
	    /// @param
	    /// @return
            FieldVector<ctype, 3> unitOuterNormal(const FieldVector<ctype, 2>&) const
            {
                return pgrid_->face_normals_[faces_of_cell_[subindex_]];
            }

	    /// @brief
	    /// @todo Doc me!
	    /// @param
	    /// @return
            FieldVector<ctype, 3> centerUnitOuterNormal() const
            {
                return pgrid_->face_normals_[faces_of_cell_[subindex_]];
            }

        protected:
            const GridType* pgrid_;
            int index_;
            int subindex_;
            OrientedEntityTable<0,1>::row_type faces_of_cell_;
	    Geometry global_geom_;
// 	    LocalGeometry in_inside_geom_;
// 	    LocalGeometry in_outside_geom_;
	    int nbcell_;
	    bool is_on_boundary_;

	    void increment()
	    {
		++subindex_;
		if (subindex_ < faces_of_cell_.size()) {
		    update();
		}
	    }

	    void update()
	    {
                EntityRep<1> face = faces_of_cell_[subindex_];
		//global_geom_ = cpgrid::Entity<1, GridType>(*pgrid_, face).geometry();
                global_geom_ = pgrid_->geometry_.template geomVector<1>()[face];
                OrientedEntityTable<1,0>::row_type cells_of_face = pgrid_->face_to_cell_[face];
		is_on_boundary_ = (cells_of_face.size() == 1);
		if (is_on_boundary_) {
		    nbcell_ = index_; // self is invalid value
		} else {
                    assert(cells_of_face.size() == 2);
                    if (cells_of_face[0].index() == index_) {
                        nbcell_ = cells_of_face[1].index();
                    } else {
                        nbcell_ = cells_of_face[0].index();
                    }
// 		    in_outside_geom_ = LocalGeometry(global_geom_.center()
// 						     - outside().geometry().center(),
// 						     global_geom_.volume());
		}
	    }

	    void setAtEnd()
	    {
		subindex_ = faces_of_cell_.size();
	    }

	    bool isAtEnd() const
	    {
		return subindex_ == faces_of_cell_.size();
	    }

            int nbcell() const
            {
		if (is_on_boundary_) {
		    OPM_THROW(std::runtime_error, "There is no outside cell, intersection is at boundary.");
		}
		return nbcell_;
	    }
        };





        template <class GridType>
        class IntersectionIterator : public Intersection<GridType>
        {
        public:
            typedef cpgrid::Intersection<GridType> Intersection;

            IntersectionIterator()
		: Intersection()
            {
            }

            IntersectionIterator(const GridType& grid, EntityRep<0> cell, bool at_end)
		: Intersection(grid, cell, 0, !at_end)
            {
                if (at_end) {
                    Intersection::setAtEnd();
                } else {
		    Intersection::update();
		}
            }

            IntersectionIterator& operator++()
            {
                Intersection::increment();
                return *this;
            }

            const Intersection* operator->() const
            {
		assert(!Intersection::isAtEnd());
                return this;
            }

            const Intersection& operator*() const
            {
		assert(!Intersection::isAtEnd());
                return *this;
            }

        };





    } // namespace cpgrid
} // namespace Dune



#endif // OPM_INTERSECTION_HEADER
