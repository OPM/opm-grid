//===========================================================================
//
// File: GridInterfaceEuler.hpp
//
// Created: Mon Jun 15 12:53:38 2009
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

#ifndef OPENRS_GRIDINTERFACEEULER_HEADER
#define OPENRS_GRIDINTERFACEEULER_HEADER


#include <dune/common/fvector.hh>
#include <dune/grid/common/referenceelements.hh>
#include <boost/iterator/iterator_facade.hpp>

namespace Dune
{


    namespace GIE
    {

	template <class DuneGrid>
	class Intersection : public boost::iterator_facade<Intersection<DuneGrid>,
							   const Intersection<DuneGrid>,
							   boost::forward_traversal_tag>
	{
	public:
	    typedef typename DuneGrid::LeafIntersectionIterator DuneIntersectionIter;
	    Intersection(const DuneGrid& grid, DuneIntersectionIter it)
		: grid_(grid), iter_(it)
	    {
	    }
	    typedef FieldVector<typename DuneGrid::ctype, DuneGrid::dimension> Vector;
	    typedef FieldVector<typename DuneGrid::ctype, DuneGrid::dimension - 1> LocalVector;
	    typedef typename DuneGrid::ctype Scalar;
	    //typedef typename DuneGrid::LeafIndexSet::IndexType Index;
	    typedef int Index;

	    enum { BoundaryMarkerIndex = -1 };

	    Scalar area() const
	    {
		return iter_->geometry().volume();
	    }

	    Vector centroid() const
	    {
		return iter_->geometry().global(localCentroid());
	    }

	    Vector normal() const
	    {
		return iter_->unitOuterNormal(localCentroid());
	    }

	    bool boundary() const
	    {
		return iter_->boundary();
	    }

	    int boundaryId() const
	    {
		return iter_->boundaryId();
	    }

	    Index neighbourIndex() const
	    {
		if (iter_->boundary()) {
		    return BoundaryMarkerIndex;
		} else {
		    return grid_.leafIndexSet().index(*iter_->outside());
		}
	    }

	    /// Used by iterator facade.
	    const Intersection& dereference() const
	    {
		return *this;
	    }
	    /// Used by iterator facade.
	    bool equal(const Intersection& other) const
	    {
		return iter_ == other.iter_;
	    }
	    /// Used by iterator facade.
	    void increment()
	    {
		++iter_;
	    }
	private:
	    const DuneGrid& grid_;
	    DuneIntersectionIter iter_;

	    LocalVector localCentroid() const
	    {
		typedef Dune::ReferenceElements<typename DuneGrid::ctype, DuneGrid::dimension - 1> RefElems;
		return RefElems::general(iter_->type()).position(0,0);
	    }

	};


	template <class DuneGrid>
	class Cell : public boost::iterator_facade<Cell<DuneGrid>,
						   const Cell<DuneGrid>,
						   boost::forward_traversal_tag>
	{
	public:
	    typedef typename DuneGrid::template Codim<0>::LeafIterator DuneCellIter;
	    Cell(const DuneGrid& grid, DuneCellIter it)
		: grid_(grid), iter_(it)
	    {
	    }
	    typedef GIE::Intersection<DuneGrid> FaceIterator;
	    typedef typename FaceIterator::Vector Vector;
	    typedef typename FaceIterator::Scalar Scalar;
	    typedef typename FaceIterator::Index Index;

	    FaceIterator facebegin() const
	    {
		return FaceIterator(grid_, iter_->ileafbegin());
	    }

	    FaceIterator faceend() const
	    {
		return FaceIterator(grid_, iter_->ileafend());
	    }

	    Scalar volume() const
	    {
		return iter_->geometry().volume();
	    }

	    Vector centroid() const
	    {
		typedef Dune::ReferenceElements<typename DuneGrid::ctype, DuneGrid::dimension> RefElems;
		Vector localpt
		    = RefElems::general(iter_->type()).position(0,0);
		return iter_->geometry().global(localpt);
	    }

	    Index index() const
	    {
		return grid_.leafIndexSet().index(*iter_);
	    }

	    /// Used by iterator facade.
	    const Cell& dereference() const
	    {
		return *this;
	    }
	    /// Used by iterator facade.
	    bool equal(const Cell& other) const
	    {
		return iter_ == other.iter_;
	    }
	    /// Used by iterator facade.
	    void increment()
	    {
		++iter_;
	    }
	private:
	    const DuneGrid& grid_;
	    DuneCellIter iter_;
	};

    } // namespace GIE


    template <class DuneGrid>
    class GridInterfaceEuler
    {
    public:
	GridInterfaceEuler(const DuneGrid& grid)
	    : grid_(grid)
	{
	}
	typedef GIE::Cell<DuneGrid> CellIterator;
	CellIterator cellbegin() const
	{
	    return CellIterator(grid_, grid_.template leafbegin<0>());
	}
	CellIterator cellend() const
	{
	    return CellIterator(grid_, grid_.template leafend<0>());
	}
        int numberOfCells() const
        {
            return grid_.size(0);
        }
    private:
	const DuneGrid& grid_;
    };



} // namespace Dune


#endif // OPENRS_GRIDINTERFACEEULER_HEADER
