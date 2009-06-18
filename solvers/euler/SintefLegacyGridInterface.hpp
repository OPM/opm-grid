//===========================================================================
//
// File: SintefLegacyGridInterface.hpp
//
// Created: Tue Jun 16 13:20:14 2009
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

#ifndef OPENRS_SINTEFLEGACYGRIDINTERFACE_HEADER
#define OPENRS_SINTEFLEGACYGRIDINTERFACE_HEADER



#include <boost/iterator/iterator_facade.hpp>

namespace samcode
{


    namespace GIE
    {

	template <class SintefLegacyGrid>
	class Face : public boost::iterator_facade<Face<SintefLegacyGrid>,
						   const Face<SintefLegacyGrid>,
						   boost::forward_traversal_tag>
	{
	public:
	    typedef typename SintefLegacyGrid::LeafFaceIterator SintefLegacyFaceIter;
	    Face(const SintefLegacyGrid& grid, int cell, bool at_end = false)
		: grid_(grid), cell_(it), subindex_(0)
	    {
		if (at_end) {
		    subindex_ = grid.template neighbours<grid::CellType, grid::HalfFaceType>(cell).size();
		}
	    }
	    typedef typename SintefLegacyGrid::point_t Vector;
	    typedef double Scalar;
	    typedef int Index;

	    enum { BoundaryMarkerIndex = -1 };

	    Scalar area() const
	    {
		return grid.getHalfFaceArea(hface());
	    }

	    Vector centroid() const
	    {
		return grid.getHalfFaceCentroid(hface());
	    }

	    Vector normal() const
	    {
		return grid.getHalfFaceNormal(hface());
	    }

	    bool boundary() const
	    {
		int face = grid.template neighbours<grid::HalfFaceType, grid::FaceType>(hface(), 0);
		return grid.template neighbours<grid::FaceType, grid::CellType>(face).size() == 1;
	    }

	    int boundaryId() const
	    {
		return boundary() ? 1 : 0;
	    }

	    Index neighbourIndex() const
	    {
		BOOST_STATIC_ASSERT(BoundaryMarkerIndex == -1); // Because that's what neighbourCell() returns.
		int face = grid.template neighbours<grid::HalfFaceType, grid::FaceType>(hface(), 0);
		return grid_helper::neighbourCell(grid_, cell_, face);
	    }

	    /// Used by iterator facade.
	    const Face& dereference() const
	    {
		return *this;
	    }
	    /// Used by iterator facade.
	    bool equal(const Face& other) const
	    {
		return cell_ == other.cell_ && subindex_ == other.subindex_;
	    }
	    /// Used by iterator facade.
	    void increment()
	    {
		++subindex_;
	    }
	private:
	    int hface()
	    {
		return grid.template neighbours<grid::CellType, grid::HalfFaceType>(cell, subindex_);
	    }
	    const SintefLegacyGrid& grid_;
	    int cell_;
	    int subindex_;
	};


	template <class SintefLegacyGrid>
	class Cell : public boost::iterator_facade<Cell<SintefLegacyGrid>,
						   const Cell<SintefLegacyGrid>,
						   boost::forward_traversal_tag>
	{
	public:
	    Cell(const SintefLegacyGrid& grid, int cell)
		: grid_(grid), cell_(cell)
	    {
	    }
	    typedef GIE::Face<SintefLegacyGrid> FaceIterator;
	    typedef typename FaceIterator::Vector Vector;
	    typedef typename FaceIterator::Scalar Scalar;
	    typedef typename FaceIterator::Index Index;

	    FaceIterator facebegin() const
	    {
		return FaceIterator(grid_, cell_);
	    }

	    FaceIterator faceend() const
	    {
		return FaceIterator(grid_, cell_, true);
	    }

	    Scalar volume() const
	    {
		return grid_.getCellVolume(cell_);

	    Vector centroid() const
	    {
		return grid_.getCellCentroid(cell_);
	    }

	    Index index() const
	    {
		return cell_;
	    }

	    /// Used by iterator facade.
	    const Cell& dereference() const
	    {
		return *this;
	    }
	    /// Used by iterator facade.
	    bool equal(const Cell& other) const
	    {
		return cell_ == other.cell_;
	    }
	    /// Used by iterator facade.
	    void increment()
	    {
		++cell_;
	    }
	private:
	    const SintefLegacyGrid& grid_;
	    int cell_;
	};

    } // namespace GIE


    template <class SintefLegacyGrid>
    class GridInterfaceEuler
    {
    public:
	GridInterfaceEuler(const SintefLegacyGrid& grid)
	    : grid_(grid)
	{
	}
	typedef GIE::Cell<SintefLegacyGrid> CellIterator;
	CellIterator cellbegin() const
	{
	    return CellIterator(grid_, 0);
	}
	CellIterator cellend() const
	{
	    return CellIterator(grid_, grid_.template numberOf<grid::CellType>());
	}
    private:
	const SintefLegacyGrid& grid_;
    };



} // namespace samcode



#endif // OPENRS_SINTEFLEGACYGRIDINTERFACE_HEADER
