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
Copyright 2009 SINTEF ICT, Applied Mathematics.
Copyright 2009 Statoil ASA.

This file is part of The Open Reservoir Simulator Project (OpenRS).

OpenRS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OpenRS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPENRS_INTERSECTION_HEADER
#define OPENRS_INTERSECTION_HEADER




#include <dune/grid/common/gridenums.hh>
#include "Entity.hpp"
#include "../common/ErrorMacros.hpp"
#include "Geometry.hpp"

namespace Dune
{
    namespace cpgrid
    {


        template <class GridType>
        class Intersection
        {
        public:
	    enum { dimension = 3 };
	    enum { dimensionworld = 3 };
	    typedef cpgrid::Entity<0, GridType> Entity;
	    typedef cpgrid::EntityPointer<0, GridType> EntityPointer;
	    typedef cpgrid::Geometry<2,3> Geometry;
	    typedef cpgrid::Geometry<2,3> LocalGeometry;
	    typedef double ctype;

            Intersection(const GridType& grid, EntityRep<0> cell, int subindex)
		: pgrid_(&grid),
		  index_(cell.index()),
		  subindex_(subindex),
		  faces_of_cell_(grid.cell_to_face_[cell]),
		  global_geom_(cpgrid::Entity<1, GridType>(grid, faces_of_cell_[subindex_]).geometry()),
		  in_inside_geom_(global_geom_.position()
				  - cpgrid::Entity<0, GridType>(grid, index_).geometry().position(),
				  global_geom_.volume())
            {
                ASSERT(index_ >= 0);
            }

            bool operator==(const Intersection& other) const
            {
                return subindex_ == other.subindex_  &&  index_ == other.index_  &&  pgrid_ == other.pgrid_;
            }

            bool operator!=(const Intersection& other) const
            {
                return !operator==(other);
            }

            bool boundary() const
            {
                EntityRep<1> face = faces_of_cell_[subindex_];
                OrientedEntityTable<1,0>::row_type cells_of_face = pgrid_->face_to_cell_[face];
                return cells_of_face.size() == 1;
            }

            /// Returns the boundary id of this intersection.
            /// There is no way to set boundary ids in the grid yet,
            /// and we have no way to provide a useful default without
            /// for instance north/south/east/west/up/down
            /// info. Therefore we just return 1 for every boundary
            /// (and of course 0 for the non-boundaries).
            int boundaryId() const
            {
                return boundary() ? 0 : 1;
            }

            bool neighbor() const
            {
                return !boundary();
            }

            EntityPointer inside() const
            {
                return EntityPointer(*pgrid_, index_);
            }

            EntityPointer outside() const
            {
                return EntityPointer(*pgrid_, nbcell());
            }

	    bool conforming() const
	    {
		return true;
	    }

            // Geometrical information about this intersection in
            // local coordinates of the inside() entity.
            const LocalGeometry& geometryInInside() const
	    {
		return in_inside_geom_;
	    }

            const LocalGeometry& intersectionSelfLocal() const
            {
                return geometryInInside();
            }

            // Geometrical information about this intersection in
            // local coordinates of the outside() entity.
            const LocalGeometry& geometryInOutside() const
	    {
		return in_outside_geom_;
	    }

            const LocalGeometry& intersectionNeighborLocal() const
            {
                return geometryInOutside();
            }

            const Geometry& geometry() const
            {
		return global_geom_;
            }

            /// Is this really just the same as geometry()?
            const Geometry& intersectionGlobal() const
            {
                return geometry();
            }

            GeometryType type() const
            {
                return geometry().type();
            }

            /// Local index of codim 1 entity in the inside() entity
            /// where intersection is contained in.
            int indexInInside() const
            {
                return subindex_;
            }

            int numberInSelf() const
            {
                return indexInInside();
            }

            /// Local index of codim 1 entity in outside() entity
            /// where intersection is contained in.
            int indexInOutside() const
            {
                EntityRep<1> face = faces_of_cell_[subindex_];
                EntityRep<0> nb(nbcell());
                OrientedEntityTable<0,1>::row_type faces_of_nb = pgrid_->cell_to_face_[nb];
                for (int i = 0; i < faces_of_nb.size(); ++i) {
                    if (faces_of_nb[i].index() == face.index()) {
                        return i;
                    }
                }
		THROW("Could not find indexInOutside().");
		return -1;
            }

            int numberInNeighbor() const
            {
                return indexInOutside();
            }

            FieldVector<ctype, 3> outerNormal(const FieldVector<ctype, 2>&) const
            {
                return pgrid_->face_normals_[faces_of_cell_[subindex_]];
            }

            FieldVector<ctype, 3> integrationOuterNormal(const FieldVector<ctype, 2>& unused) const
            {
                FieldVector<ctype, 3> n = pgrid_->face_normals_[faces_of_cell_[subindex_]];
                return n*=geometry().integrationElement(unused);
            }

            FieldVector<ctype, 3> unitOuterNormal(const FieldVector<ctype, 2>&) const
            {
                return pgrid_->face_normals_[faces_of_cell_[subindex_]];
            }

        protected:
            const GridType* pgrid_;
            int index_;
            int subindex_;
            OrientedEntityTable<0,1>::row_type faces_of_cell_;
	    Geometry global_geom_;
	    LocalGeometry in_inside_geom_;
	    LocalGeometry in_outside_geom_;

	    void increment()
	    {
		++subindex_;
		if (subindex_ < faces_of_cell_.size()) {
		    update();
		}
	    }

	    void update()
	    {
		in_outside_geom_ = LocalGeometry(global_geom_.position()
						 - outside().geometry().position(),
						 global_geom_.volume());
	    }

	    void setAtEnd()
	    {
		subindex_ = faces_of_cell_.size();
	    }

            int nbcell() const
            {
                EntityRep<1> face = faces_of_cell_[subindex_];
                OrientedEntityTable<1,0>::row_type cells_of_face = pgrid_->face_to_cell_[face];
                if (cells_of_face.size() == 1) {
                    THROW("Face " << face.index() << " is on the boundary, you cannot get the neighbouring cell.");
                } else {
                    ASSERT(cells_of_face.size() == 2);
                    if (cells_of_face[0].index() == index_) {
                        return cells_of_face[1].index();
                    } else {
                        return cells_of_face[0].index();
                    }
                }
            }
        };





        template <class GridType>
        class IntersectionIterator : public Intersection<GridType>
        {
        public:
            typedef cpgrid::Intersection<GridType> Intersection;

            IntersectionIterator(const GridType& grid, EntityRep<0> cell, bool at_end)
                    : Intersection(grid, cell, 0)
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
                return this;
            }

            const Intersection& operator*() const
            {
                return *this;
            }

        };





    } // namespace cpgrid
} // namespace Dune



#endif // OPENRS_INTERSECTION_HEADER
