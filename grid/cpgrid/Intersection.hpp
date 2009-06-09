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
	    Intersection(const GridType& grid, EntityRep<0> cell, int subindex)
		: grid_(grid),
		  index_(cell.index()),
		  subindex_(subindex),
		  faces_of_cell_(grid.cell_to_face_[cell])
	    {
		ASSERT(index_ >= 0);
	    }

	    bool operator!=(const Intersection& other) const
	    {
		return subindex_ != other.subindex_  ||  index_ != other.index_  ||  &grid_ != &other.grid_;
	    }

	    bool boundary() const
	    {
		EntityRep<1> face = faces_of_cell_[subindex_];
		OrientedEntityTable<1,0>::row_type cells_of_face = grid_.face_to_cell_[face];
		return cells_of_face.size() == 1;
	    }

	    bool neighbor() const
	    {
		return !boundary();
	    }

	    EntityPointer<0, GridType> inside() const
	    {
		return EntityPointer<0, GridType>(grid_, index_);
	    }

	    EntityPointer<0, GridType> outside() const
	    {
		return EntityPointer<0, GridType>(grid_, nbcell());
	    }

	    /*
const LocalGeometry & 	geometryInInside () const
 	geometrical information about this intersection in local coordinates of the inside() entity. 
const LocalGeometry & 	intersectionSelfLocal () const
 	please read the details 
const LocalGeometry & 	geometryInOutside () const
 	geometrical information about this intersection in local coordinates of the outside() entity. 
const LocalGeometry & 	intersectionNeighborLocal () const
 	please read the details 
	    */

	    const Geometry<2,3>& geometry() const
	    {
		return cpgrid::Entity<1, GridType>(grid_, faces_of_cell_[subindex_]).geometry();
	    }

	    /// Is this really just the same as geometry()?
	    const Geometry<2,3>& intersectionGlobal () const
	    {
		return geometry();
	    }

	    GeometryType type () const
	    {
		return geometry().type();
	    }
						      /*
int 	indexInInside () const
 	Local index of codim 1 entity in the inside() entity where intersection is contained in. 
int 	numberInSelf () const
 	please read the details 
int 	indexInOutside () const
 	Local index of codim 1 entity in outside() entity where intersection is contained in. 
int 	numberInNeighbor () const
 	please read the details 
						      */
	    FieldVector<double, 3> outerNormal (const FieldVector<double, 2>&) const
	    {
		return grid_.face_normals_[faces_of_cell_[subindex_]];
	    }

	    FieldVector<double, 3> integrationOuterNormal (const FieldVector<double, 2>&) const
	    {
		FieldVector<double, 3> n = grid_.face_normals_[faces_of_cell_[subindex_]];
		return n*=double(geometry().volume());
	    }

	    FieldVector<double, 3> unitOuterNormal (const FieldVector<double, 2>&) const
	    {
		return grid_.face_normals_[faces_of_cell_[subindex_]];
	    }

	protected:
	    const GridType& grid_;
	    const int index_;
	    int subindex_;
	    OrientedEntityTable<0,1>::row_type faces_of_cell_;

	    int nbcell() const
	    {
		EntityRep<1> face = faces_of_cell_[subindex_];
		OrientedEntityTable<1,0>::row_type cells_of_face = grid_.face_to_cell_[face];
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
	    IntersectionIterator(const GridType& grid, EntityRep<0> cell, bool at_end)
		: Intersection<GridType>(grid, cell, 0)
	    {
		if (at_end) {
		    Intersection<GridType>::subindex_ += Intersection<GridType>::faces_of_cell_.size();
		}
	    }

	    IntersectionIterator& operator++()
	    {
		++Intersection<GridType>::subindex_;
		return *this;
	    }

	    const Intersection<GridType>* operator->() const
	    {
		return this;
	    }

	    const Intersection<GridType>& operator*() const
	    {
		return *this;
	    }

	};





    } // namespace cpgrid
} // namespace Dune



#endif // OPENRS_INTERSECTION_HEADER
