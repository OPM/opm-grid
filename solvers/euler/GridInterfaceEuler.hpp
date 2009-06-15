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


namespace Dune
{


    namespace GIE
    {

	template <class DuneGrid>
	class Face
	{
	public:
	    typedef FieldVector<typename DuneGrid::ctype, DuneGrid::dimension> Vector;
	    typedef typename DuneGrid::ctype Scalar;
	    Scalar area() const;
	    Vector centroid() const;
	    Vector normal() const;
	protected:
	    void increment();
	    bool isAtEnd() const;
	};

	template <class DuneGrid>
	class FaceIterator : public Face<DuneGrid>
	{
	public:
	    FaceIterator& operator++()
	    {
		Face<DuneGrid>::increment();
		return *this;
	    }

	    const Face<DuneGrid>* operator->() const
	    {
		ASSERT(!Face<DuneGrid>::isAtEnd());
		return this;
	    }

	    const Face<DuneGrid>& operator*() const
	    {
		ASSERT(!Face<DuneGrid>::isAtEnd());
		return *this;
	    }
	};


	template <class DuneGrid>
	class Cell
	{
	public:
	    typedef GIE::FaceIterator<DuneGrid> FaceIterator;
	    typedef typename FaceIterator::Vector Vector;
	    typedef typename FaceIterator::Scalar Scalar;
	    FaceIterator facebegin() const;
	    FaceIterator faceend() const;
	    Scalar volume() const;
	    Vector centroid() const;
	protected:
	    void increment();
	    bool isAtEnd() const;
	};

	template <class DuneGrid>
	class CellIterator : public Cell<DuneGrid>
	{
	public:
	    CellIterator& operator++()
	    {
		Cell<DuneGrid>::increment();
		return *this;
	    }

	    const Cell<DuneGrid>* operator->() const
	    {
		ASSERT(!Cell<DuneGrid>::isAtEnd());
		return this;
	    }

	    const Cell<DuneGrid>& operator*() const
	    {
		ASSERT(!Cell<DuneGrid>::isAtEnd());
		return *this;
	    }
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
	typedef GIE::CellIterator<DuneGrid> CellIterator;
	CellIterator cellbegin() const
	{
	    return CellIterator(grid_, 0);
	}
	CellIterator cellend() const
	{
	}
    private:
	const DuneGrid& grid_;
    };



} // namespace Dune


#endif // OPENRS_GRIDINTERFACEEULER_HEADER
