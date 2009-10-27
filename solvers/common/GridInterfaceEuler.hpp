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

#include "config.h"
#include <climits>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/scoped_ptr.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/SparseTable.hpp>
#include <dune/common/StopWatch.hpp>
#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/common/mcmgmapper.hh>


namespace Dune
{


    namespace GIE
    {
	template <class GridInterface, class EntityPointerType>
	class Cell;

	/// @brief
	/// @todo Doc me!
	/// @tparam
	/// @param
	template <class GridInterface>
	class Intersection : public boost::iterator_facade<Intersection<GridInterface>,
							   const Intersection<GridInterface>,
							   boost::forward_traversal_tag>
	{
	public:

	    /// @brief
	    /// @todo Doc me!
	    typedef typename GridInterface::DuneIntersectionIterator DuneIntersectionIter;
	    /// @brief
	    /// @todo Doc me!
	    /// @param
	    Intersection()
		: pgrid_(0), iter_(), local_index_(-1)
	    {
	    }
	    /// @brief
	    /// @todo Doc me!
	    /// @param
	    Intersection(const GridInterface& grid,
                         const DuneIntersectionIter& it,
			 const int local_index)
		: pgrid_(&grid), iter_(it), local_index_(local_index)
	    {
	    }
	    /// @brief
	    /// @todo Doc me!
	    typedef typename GridInterface::GridType::ctype Scalar;
	    typedef FieldVector<Scalar, GridInterface::GridType::dimension> Vector;
	    typedef FieldVector<Scalar, GridInterface::GridType::dimension - 1> LocalVector;
	    typedef int Index;
	    typedef GIE::Cell<GridInterface, typename GridInterface::GridType::template Codim<0>::EntityPointer> Cell;
	    /// @brief
	    /// @todo Doc me!
	    enum { BoundaryMarkerIndex = -999, LocalEndIndex = INT_MAX };

	    /// @brief
	    /// @todo Doc me!
	    /// @return
	    Scalar area() const
	    {
		return iter_->geometry().volume();
	    }

	    /// @brief
	    /// @todo Doc me!
	    /// @return
	    Vector centroid() const
	    {
		return iter_->geometry().global(localCentroid());
	    }

	    /// @brief
	    /// @todo Doc me!
	    /// @return
	    Vector normal() const
	    {
		return iter_->unitOuterNormal(localCentroid());
	    }

	    /// @brief
	    /// @todo Doc me!
	    /// @return
	    bool boundary() const
	    {
		return iter_->boundary();
	    }

	    /// @brief
	    /// @todo Doc me!
	    /// @return
	    int boundaryId() const
	    {
		return iter_->boundaryId();
	    }

	    /// @brief
	    /// @todo Doc me!
	    /// @return
	    Cell cell() const
	    {
		return Cell(*pgrid_, iter_->inside());
	    }

	    /// @brief
	    /// @todo Doc me!
	    /// @return
	    Index cellIndex() const
	    {
		return pgrid_->mapper().map(*iter_->inside());
	    }

	    /// @brief
	    /// @todo Doc me!
	    /// @return
	    Cell neighbourCell() const
	    {
		return Cell(*pgrid_, iter_->outside());
	    }

	    /// @brief
	    /// @todo Doc me!
	    /// @return
	    Index neighbourCellIndex() const
	    {
		if (iter_->boundary()) {
		    return BoundaryMarkerIndex;
		} else {
		    return pgrid_->mapper().map(*iter_->outside());
		}
	    }

	    /// @brief
	    /// @todo Doc me!
	    /// @return
	    Index index() const
	    {
		return pgrid_->faceIndex(cellIndex(), localIndex());
	    }

	    /// @brief
	    /// @todo Doc me!
	    /// @return
	    Index localIndex() const
	    {
		return local_index_;
	    }

	    /// @brief
	    /// @todo Doc me!
	    /// @return
	    Scalar neighbourCellVolume() const
	    {
		return iter_->outside()->geometry().volume();
	    }

	    /// Used by iterator facade.
	    const Intersection& dereference() const
	    {
		return *this;
	    }
	    /// Used by iterator facade.
	    bool equal(const Intersection& other) const
	    {
		// Note that we do not compare the local_index_ members,
		// since they may or may not be equal for end iterators.
		return iter_ == other.iter_;
	    }
	    /// Used by iterator facade.
	    void increment()
	    {
		++iter_;
		++local_index_;
	    }
	    /// Gives an ordering of intersections.
	    bool operator<(const Intersection& other) const
	    {
		if (cellIndex() == other.cellIndex()) {
		    return localIndex() < other.localIndex();
		} else {
		    return cellIndex() < other.cellIndex();
		}
	    }
	private:
	    const GridInterface* pgrid_;
	    DuneIntersectionIter iter_;
	    int local_index_;

	    LocalVector localCentroid() const
	    {
		typedef Dune::ReferenceElements<Scalar, GridInterface::GridType::dimension-1> RefElems;
		return RefElems::general(iter_->type()).position(0,0);
	    }

	};


	template <class GridInterface, class EntityPointerType>
	class Cell
	{
	public:
	    Cell()
		: pgrid_(0), iter_()
	    {
	    }
	    Cell(const GridInterface& grid,
		 const EntityPointerType& it)
		: pgrid_(&grid), iter_(it)
	    {
	    }
	    typedef GIE::Intersection<GridInterface> FaceIterator;
	    typedef typename FaceIterator::Vector Vector;
	    typedef typename FaceIterator::Scalar Scalar;
	    typedef typename FaceIterator::Index Index;

	    FaceIterator facebegin() const
	    {
		return FaceIterator(*pgrid_, iter_->ileafbegin(), 0);
	    }

	    FaceIterator faceend() const
	    {
		return FaceIterator(*pgrid_, iter_->ileafend(),
                                    FaceIterator::LocalEndIndex);
	    }

	    Scalar volume() const
	    {
		return iter_->geometry().volume();
	    }

	    Vector centroid() const
	    {
		typedef Dune::ReferenceElements<Scalar, GridInterface::GridType::dimension> RefElems;
		Vector localpt
		    = RefElems::general(iter_->type()).position(0,0);
		return iter_->geometry().global(localpt);
	    }

	    Index index() const
	    {
                return pgrid_->mapper().map(*iter_);
	    }
	protected:
	    const GridInterface* pgrid_;
	    EntityPointerType iter_;
	};


	template <class GridInterface>
	class CellIterator
	    : public boost::iterator_facade<CellIterator<GridInterface>,
					    const CellIterator<GridInterface>,
					    boost::forward_traversal_tag>,
	      public Cell<GridInterface, typename GridInterface::GridType::template Codim<0>::LeafIterator>
	{
	private:
	    typedef typename GridInterface::GridType::template Codim<0>::LeafIterator DuneCellIter;
	    typedef Cell<GridInterface, DuneCellIter> CellType;
	public:
	    typedef typename CellType::Vector Vector;
	    typedef typename CellType::Scalar Scalar;
	    typedef typename CellType::Index Index;

	    CellIterator(const GridInterface& grid, DuneCellIter it)
		: CellType(grid, it)
	    {
	    }
	    /// Used by iterator facade.
	    const CellIterator& dereference() const
	    {
		return *this;
	    }
	    /// Used by iterator facade.
	    bool equal(const CellIterator& other) const
	    {
		return CellType::iter_ == other.CellType::iter_;
	    }
	    /// Used by iterator facade.
	    void increment()
	    {
		++CellType::iter_;
	    }
	};

        template<int dim>
        struct AllCellsLayout {
            bool contains(GeometryType gt) { return gt.dim() == dim; }
        };

    } // namespace GIE


    template <class DuneGrid>
    class GridInterfaceEuler
    {
    public:
        typedef LeafMultipleCodimMultipleGeomTypeMapper<DuneGrid, GIE::AllCellsLayout> Mapper;
	typedef typename DuneGrid::LeafIntersectionIterator DuneIntersectionIterator;
	typedef DuneGrid GridType;
	typedef GridInterfaceEuler<DuneGrid> InterfaceType;
	typedef GIE::CellIterator<InterfaceType> CellIterator;
	typedef typename CellIterator::Vector Vector;
	typedef typename CellIterator::Scalar Scalar;
	typedef typename CellIterator::Index Index;
	enum { Dimension = DuneGrid::dimension };

	GridInterfaceEuler()
	    : pgrid_(0), num_faces_(0), max_faces_per_cell_(0)
	{
	}
	explicit GridInterfaceEuler(const DuneGrid& grid, bool build_facemap = true)
	    : pgrid_(&grid), pmapper_(new Mapper(grid)), num_faces_(0), max_faces_per_cell_(0)
	{
	    if (build_facemap) {
		buildFaceIndices();
	    }
	}
	void init(const DuneGrid& grid, bool build_facemap = true)
	{
	    pgrid_ = &grid;
            pmapper_.reset(new Mapper(grid));
	    if (build_facemap) {
		buildFaceIndices();
	    }
	}
	CellIterator cellbegin() const
	{
	    return CellIterator(*this, grid().template leafbegin<0>());
	}
	CellIterator cellend() const
	{
	    return CellIterator(*this, grid().template leafend<0>());
	}
        int numberOfCells() const
        {
            return grid().size(0);
        }
        int numberOfFaces() const
        {
	    ASSERT(num_faces_ != 0);
            return num_faces_;
        }
	int maxFacesPerCell() const
	{
	    ASSERT(max_faces_per_cell_ != 0);
	    return max_faces_per_cell_;
	}
	const DuneGrid& grid() const
	{
	    ASSERT(pgrid_);
	    return *pgrid_;
	}

	// The following are primarily helpers for the implementation,
	// perhaps they should be private?
	const Mapper& mapper() const
	{
	    ASSERT (pmapper_);
	    return *pmapper_;
	}
	Index faceIndex(int cell_index, int local_face_index) const
	{
	    ASSERT(num_faces_ != 0);
	    return face_indices_[cell_index][local_face_index];
	}
	typedef SparseTable<int>::row_type Indices;
	Indices faceIndices(int cell_index) const
	{
	    ASSERT(num_faces_ != 0);
	    return face_indices_[cell_index];
	}
    private:
	const DuneGrid* pgrid_;
        boost::scoped_ptr<Mapper> pmapper_;
	int num_faces_;
	int max_faces_per_cell_;
	SparseTable<int> face_indices_;

	void buildFaceIndices()
	{
#ifdef VERBOSE
	    std::cout << "Building unique face indices... " << std::flush;
	    time::StopWatch clock;
	    clock.start();
#endif
            typedef CellIterator CI;
            typedef typename CI::FaceIterator FI;

	    // We build the actual cell to face mapping in two passes.
	    // [code mostly lifted from IncompFlowSolverHybrid::enumerateGridDof(),
	    //  but with a twist: This code builds a mapping from cells in index
	    //  order to unique face numbers, while the mapping built in the
	    //  enumerateGridDof() method was ordered by cell iterator order]

            // Allocate and reserve structures.
            const int nc = numberOfCells();
            std::vector<int> cell(nc, -1);
            std::vector<int> num_faces(nc); // In index order.
            std::vector<int> fpos;     fpos  .reserve(nc + 1);
            std::vector<int> num_cf;   num_cf.reserve(nc); // In iterator order.
            std::vector<int> faces ;

            // First pass: enumerate internal faces.
            int cellno = 0; fpos.push_back(0);
            int tot_ncf = 0, tot_ncf2 = 0, max_ncf = 0;
            for (CI c = cellbegin(); c != cellend(); ++c, ++cellno) {
                const int c0 = c->index();
                ASSERT((0 <= c0) && (c0 < nc) && (cell[c0] == -1));
                cell[c0] = cellno;
                num_cf.push_back(0);
                int& ncf = num_cf.back();
                for (FI f = c->facebegin(); f != c-> faceend(); ++f) {
                    if (!f->boundary()) {
                        const int c1 = f->neighbourCellIndex();
                        ASSERT((0 <= c1) && (c1 < nc) && (c1 != c0));

                        if (cell[c1] == -1) {
                            // Previously undiscovered internal face.
                            faces.push_back(c1);
                        }
                    }
                    ++ncf;
                }
		num_faces[c0] = ncf;
                fpos.push_back(int(faces.size()));
                max_ncf  = std::max(max_ncf, ncf);
                tot_ncf  += ncf;
                tot_ncf2 += ncf * ncf;
            }
            ASSERT(cellno == nc);

	    // Build cumulative face sizes enabling direct insertion of
	    // face indices into cfdata later.
	    std::vector<int> cumul_num_faces(numberOfCells() + 1);
	    cumul_num_faces[0] = 0;
	    std::partial_sum(num_faces.begin(), num_faces.end(), cumul_num_faces.begin() + 1);

            // Avoid (most) allocation(s) inside 'c' loop.
            std::vector<int>    l2g;
	    l2g.reserve(max_ncf);
	    std::vector<double> cfdata(tot_ncf);
            int total_num_faces = int(faces.size());

            // Second pass: build cell-to-face mapping, including boundary.
            typedef std::vector<int>::iterator VII;
            for (CI c = cellbegin(); c != cellend(); ++c) {
                const int c0 = c->index();
                ASSERT ((0 <=      c0 ) && (     c0  < nc) &&
                        (0 <= cell[c0]) && (cell[c0] < nc));
                const int ncf = num_cf[cell[c0]];
                l2g.resize(ncf, 0);
                for (FI f = c->facebegin(); f != c->faceend(); ++f) {
                    if (f->boundary()) {
                        // External, not counted before.  Add new face...
                        l2g[f->localIndex()] = total_num_faces++;
                    } else {
                        // Internal face.  Need to determine during
                        // traversal of which cell we discovered this
                        // face first, and extract the face number
                        // from the 'faces' table range of that cell.

                        // Note: std::find() below is potentially
                        // *VERY* expensive (e.g., large number of
                        // seeks in moderately sized data in case of
                        // faulted cells).
                        const int c1 = f->neighbourCellIndex();
                        ASSERT ((0 <=      c1 ) && (     c1  < nc) &&
                                (0 <= cell[c1]) && (cell[c1] < nc));

                        int t = c0, seek = c1;
                        if (cell[seek] < cell[t])
                            std::swap(t, seek);
                        int s = fpos[cell[t]], e = fpos[cell[t] + 1];
                        VII p = std::find(faces.begin() + s, faces.begin() + e, seek);
                        ASSERT(p != faces.begin() + e);
                        l2g[f->localIndex()] = p - faces.begin();
                    }
                }
		ASSERT(int(l2g.size()) == num_faces[c0]);
		std::copy(l2g.begin(), l2g.end(), cfdata.begin() + cumul_num_faces[c0]);
            }
	    num_faces_ = total_num_faces;
	    max_faces_per_cell_ = max_ncf;
            face_indices_.assign(cfdata.begin(), cfdata.end(),
				 num_faces.begin(), num_faces.end());

#ifdef VERBOSE
	    clock.stop();
	    double elapsed = clock.secsSinceStart();
	    std::cout << "done.     Time elapsed: " << elapsed << std::endl;
#endif
	}

    };



} // namespace Dune


#endif // OPENRS_GRIDINTERFACEEULER_HEADER
