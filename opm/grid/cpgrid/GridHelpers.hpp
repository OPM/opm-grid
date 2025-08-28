/*
  Copyright 2014, 2015 Dr. Markus Blatt - HPC-Simulation-Software & Services.
  Copyright 2014 Statoil AS
  Copyright 2015 NTNU

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
#ifndef DUNE_CORNERPOINT_GRIDHELPERS_HEADER_INCLUDED
#define DUNE_CORNERPOINT_GRIDHELPERS_HEADER_INCLUDED

#include <functional>

#include <opm/grid/GridHelpers.hpp>

#include <opm/grid/utility/platform_dependent/disable_warnings.h>

#include <opm/grid/CpGrid.hpp>
#include <dune/common/iteratorfacades.hh>

#include <opm/grid/utility/platform_dependent/reenable_warnings.h>

namespace Dune
{
namespace cpgrid
{
/// \brief A proxy class representing a row of FaceCellsContainer.
class FaceCellsProxy
{
public:
    /// \brief Constructor.
    /// \param grid The grid whose face to cell mapping we represent.
    /// \param cell_index The index of the cell we repesent.
    FaceCellsProxy(const Dune::CpGrid* grid, int cell_index)
        : grid_(grid), cell_index_(cell_index)
    {}
    /// \brief Get the index of the cell associated with a local_index.
    int operator[](int local_index)
    {
        return grid_->faceCell(cell_index_, local_index);
    }
private:
    const Dune::CpGrid* grid_;
    int cell_index_;
};

/// \brief A class representing the face to cells mapping similar to the
/// way done in UnstructuredGrid.
class FaceCellsContainerProxy
{
public:
    typedef FaceCellsProxy row_type;

    /// \brief Constructor.
    /// \param grid The grid whose information we represent.
    explicit FaceCellsContainerProxy(const Dune::CpGrid* grid)
        : grid_(grid)
    {}
    /// \brief Get the mapping for a cell.
    /// \param cell_index The index of the cell.
    FaceCellsProxy operator[](int cell_index) const
    {
        return FaceCellsProxy(grid_, cell_index);
    }
    /// \brief Get a face associated with a cell.
    /// \param cell_index The index of the cell.
    /// \param local_index The local index of the cell, either 0 or 1.
    /// \param The index of the face or -1 if it is not present because of
    /// a boundary.
    int operator()(int cell_index, int local_index) const
    {
        return grid_->faceCell(cell_index, local_index);
    }
private:
    const Dune::CpGrid* grid_;
};


    class IndexIterator
    {
    public:
        explicit IndexIterator(int index)
        : index_(index)
        {}

        void increment()
        {
            ++index_;
        }
        void decrement()
        {
            --index_;
        }
        void advance(int n)
        {
            index_+=n;
        }
        int distanceTo(const IndexIterator& o)const
        {
            return o.index_-index_;
        }
        bool equals(const IndexIterator& o) const
        {
            return index_==o.index_;
        }
    protected:
        int index_;
    };


/// \brief A proxy class representing a row of LocalIndexContainerProxy.
/// \tparam AccessMethod Function pointer to access the values of a sparse
///                      row (e.g. the faces attached to a cell.
/// \tparam SizeMethod   Fuction pointer to access the size of the sparse row
///                      (e.g. the number of faces attached to a cell.
template<int (Dune::CpGrid::*AccessMethod)(int,int)const,
         int (Dune::CpGrid::*SizeMethod)(int)const>
class LocalIndexProxy
{
public:
    class iterator
        : public Dune::RandomAccessIteratorFacade<iterator,int, int, int>,
          public IndexIterator
    {
    public:
        iterator(const Dune::CpGrid* grid, int outer_index, int inner_index)
            : IndexIterator(inner_index), grid_(grid), outer_index_(outer_index)
        {}
        int dereference() const
        {
            return std::mem_fn(AccessMethod)(*grid_, outer_index_, this->index_);
        }
        int elementAt(int n) const
        {
            return std::mem_fn(AccessMethod)(*grid_, outer_index_, n);
        }
    private:
        const Dune::CpGrid* grid_;
        int outer_index_;
    };

    typedef iterator const_iterator;

    /// \brief Constructor.
    /// \param grid The grid whose face to cell mapping we represent.
    /// \param cell_index The index of the cell we repesent.
    LocalIndexProxy(const Dune::CpGrid* grid, int cell_index)
        : grid_(grid), cell_index_(cell_index)
    {}
    /// \brief Get the index of the cell associated with a local_index.
    int operator[](int local_index)
    {
        return std::mem_fn(AccessMethod)(*grid_, cell_index_, local_index);
    }
    const_iterator begin()
    {
        return const_iterator(grid_, cell_index_, 0);
    }
    const_iterator end()
    {
        return const_iterator(grid_, cell_index_,
                              std::mem_fn(SizeMethod)(*grid_, cell_index_));
    }
private:
    const Dune::CpGrid* grid_;
    int cell_index_;
};

/// \brief A class representing the sparse mapping of entity relations (e.g. vertices of faces).
/// \tparam AccessMethod Function pointer to access the values of a sparse
///                      row (e.g. the vertices attached to a face.
/// \tparam SizeMethod   Fuction pointer to access the size of the sparse row
///                      (e.g. the number of vertices attached to a face.
template<int (Dune::CpGrid::*AccessMethod)(int,int)const,
         int (Dune::CpGrid::*SizeMethod)(int)const>
class LocalIndexContainerProxy
{
public:
    typedef LocalIndexProxy<AccessMethod, SizeMethod> row_type;
    /// \brief Constructor.
    /// \param grid The grid whose information we represent.
    explicit LocalIndexContainerProxy(const Dune::CpGrid* grid)
        : grid_(grid)
    {}
    /// \brief Get the mapping for a cell.
    /// \param cell_index The index of the cell.
    row_type operator[](int cell_index) const
    {
        return row_type(grid_, cell_index);
    }
    /// \brief Get a face associated with a cell.
    /// \param cell_index The index of the cell.
    /// \param local_index The local index of the cell, either 0 or 1.
    /// \param The index of the face or -1 if it is not present because of
    /// a boundary.
    int operator()(int cell_index, int local_index) const
    {
        return std::mem_fn(AccessMethod)(*grid_, cell_index, local_index);
    }
private:
    const Dune::CpGrid* grid_;
};

/// \brief A class representing the face to vertices mapping similar to the
/// way done in UnstructuredGrid.
class FaceVerticesContainerProxy
    : public LocalIndexContainerProxy<&Dune::CpGrid::faceVertex, &Dune::CpGrid::numFaceVertices>
{
public:
    /// \brief Constructor.
    /// \param grid The grid whose information we represent.
    explicit FaceVerticesContainerProxy(const Dune::CpGrid* grid)
        : LocalIndexContainerProxy<&Dune::CpGrid::faceVertex, &Dune::CpGrid::numFaceVertices>(grid)
    {}
};

class Cell2FacesRow
{
public:
    class iterator
        : public Dune::RandomAccessIteratorFacade<iterator,int, int, int>,
        public IndexIterator
    {
    public:
        iterator(const Dune::cpgrid::OrientedEntityTable<0,1>::row_type& row,
                 int index, int cell_index)
            : IndexIterator(index), row_(row), cell_index_(cell_index)
        {}
        int dereference() const
        {
            return row_[this->index_].index();
        }
        int elementAt(int n) const
        {
            return row_[n].index();
        }
        int getCellIndex()const
        {
            return cell_index_;
        }

    private:
        // Note that row_ is a Opm::iterator_range returned by Opm::SparseTable.
        // and stored in CellFacesRow. It is a temporary object that only lives
        // as long as there is a reference to it e.g. by storing the Cell2FacesRow.
        // Therefore it needs to be copied. As it is rather light weight this
        // should be fast anyway. A const reference would mean that the iterator
        // is not assignable.
        const Dune::cpgrid::OrientedEntityTable<0,1>::row_type row_;
        int cell_index_;
    };

    typedef iterator const_iterator;

    Cell2FacesRow(const Dune::cpgrid::OrientedEntityTable<0,1>::row_type& row,
                  const int cell_index)
        : row_(row), cell_index_(cell_index)
    {}

    const_iterator begin() const
    {
        return const_iterator(row_, 0, cell_index_);
    }

    const_iterator end() const
    {
        return const_iterator(row_, row_.size(), cell_index_);
    }

private:
    // Note that row_ is a Opm::iterator_range returned by Opm::SparseTable.
    // and stored in CellFacesRow. It is a temporary object that only lives
    // as long as there is a reference to it e.g. by storing the Cell2FacesRow.
    // Therefore it needs to be copied. As it is rather light weight this
    // should be fast anyway.  A const reference would mean that the row
    // is not assignable.
    const Dune::cpgrid::OrientedEntityTable<0,1>::row_type row_;
    const int cell_index_;
};


class Cell2FacesContainer
{
public:
    typedef  Cell2FacesRow row_type;

    explicit Cell2FacesContainer(const Dune::CpGrid* grid)
        : grid_(grid)
    {};

    Cell2FacesRow operator[](int cell_index) const
    {
        const auto& row = grid_->cellFaceRow(cell_index);
        return Cell2FacesRow(row, cell_index);
    }

        /// \brief Get the number of non-zero entries.
    std::size_t noEntries() const
    {
        return grid_->numCellFaces();
    }
private:
    const Dune::CpGrid* grid_;
};

} // end namespace cpgrid
} // end namespace Dune

namespace Opm
{
namespace UgGridHelpers
{
template<>
struct Cell2FacesTraits<Dune::CpGrid>
{
    typedef Dune::cpgrid::Cell2FacesContainer Type;
};
/// \brief An iterator over the cell volumes.
template<Dune::FieldVector<double, 3> (Dune::CpGrid::*Method)(int)const>
class CpGridCentroidIterator
    : public Dune::RandomAccessIteratorFacade<CpGridCentroidIterator<Method>, Dune::FieldVector<double, 3>,
                                              Dune::FieldVector<double, 3>, int>
{
public:
    /// \brief Creates an iterator.
    /// \param grid The grid the iterator belongs to.
    /// \param cell_index The position of the iterator.
    CpGridCentroidIterator(const  Dune::CpGrid& grid, int cell_index)
        : grid_(&grid), cell_index_(cell_index)
    {}

    Dune::FieldVector<double, 3> dereference() const
    {
        return std::mem_fn(Method)(*grid_, cell_index_);
    }
    void increment()
    {
        ++cell_index_;
    }
    Dune::FieldVector<double, 3> elementAt(int n) const
    {
        return  std::mem_fn(Method)(*grid_, n);
    }
    void advance(int n)
    {
        cell_index_+=n;
    }
    void decrement()
    {
        --cell_index_;
    }
    int distanceTo(const CpGridCentroidIterator& o) const
    {
        return o.cell_index_-cell_index_;
    }
    bool equals(const CpGridCentroidIterator& o) const
    {
        return o.grid_==grid_ && o.cell_index_==cell_index_;
    }

private:
    const Dune::CpGrid* grid_;
    int cell_index_;
};

template<>
struct CellCentroidTraits<Dune::CpGrid>
{
    typedef CpGridCentroidIterator<&Dune::CpGrid::cellCentroid> IteratorType;
    typedef const double* ValueType;
};

typedef Dune::FieldVector<double, 3> Vector;

/// \brief Get the number of cells of a grid.
int numCells(const Dune::CpGrid& grid);

/// \brief Get the number of faces of a grid.
int numFaces(const  Dune::CpGrid& grid);

/// \brief Get the dimensions of a grid
int dimensions(const Dune::CpGrid& grid);

/// \brief Get the number of faces, where each face counts as many times as there are adjacent faces
int numCellFaces(const Dune::CpGrid& grid);

/// \brief Get the cartesion dimension of the underlying structured grid.
const int* cartDims(const Dune::CpGrid& grid);

/// \brief Get the local to global index mapping.
///
/// The global index is the index of the active cell
/// in the underlying structured grid.
const int*  globalCell(const Dune::CpGrid&);

#if HAVE_ECL_INPUT
/// \brief Create Eclipse style ACTNUM array.
///
/// Create a vector with global cartesian number of elements,
/// the value is 0 for inactive cells and one for active cells.
std::vector<int> createACTNUM(const Dune::CpGrid& grid);

/// Construct an EclipseGrid instance based on the inputGrid, with modifications to
/// zcorn and actum from the dune CPGrid
EclipseGrid createEclipseGrid(const Dune::CpGrid& grid, const EclipseGrid& inputGrid);
#endif

CellCentroidTraits<Dune::CpGrid>::IteratorType
beginCellCentroids(const Dune::CpGrid& grid);

/// \brief Get a coordinate of a specific cell centroid.
/// \brief grid The grid.
/// \brief cell_index The index of the specific cell.
/// \breif coordinate The coordinate index.
double cellCentroidCoordinate(const Dune::CpGrid& grid, int cell_index,
                              int coordinate);

/// \brief Get the centroid of a cell.
/// \param grid The grid whose cell centroid we query.
/// \param cell_index The index of the corresponding cell.
Vector cellCentroid(const Dune::CpGrid& grid, int cell_index);

/// \brief Get vertical position of cell center ("zcorn" average).
/// \brief grid The grid.
/// \brief cell_index The index of the specific cell.
double cellCenterDepth(const Dune::CpGrid& grid, int cell_index);


/// \brief Get a coordinate of a specific face center.
/// \brief calculated as the raw average of the cell corners
/// \param grid The grid.
/// \param cell_index The index of the specific cell.
/// \param face_tag The logical cartesian index of the face
Vector faceCenterEcl(const Dune::CpGrid& grid, int cell_index, int face_tag);

/// \brief Get a area weighted normal vector of a specific face.
/// \brief calculated without introducing a center point
/// \brief For cornerpoint grids this is supposed to give
/// \brief values closer to Ecl.
/// \param grid The grid.
/// \param face_index The index of the specific face.
Vector faceAreaNormalEcl(const Dune::CpGrid& grid, int face_index);

/// \brief Get the volume of a cell.
/// \param grid The grid the cell belongs to.
/// \param cell_index The index of the cell.
double cellVolume(const  Dune::CpGrid& grid, int cell_index);

/// \brief An iterator over the cell volumes.
class CellVolumeIterator
    : public Dune::RandomAccessIteratorFacade<CellVolumeIterator, double, double, int>
{
public:
    /// \brief Creates an iterator.
    /// \param grid The grid the iterator belongs to.
    /// \param cell_index The position of the iterator.
    CellVolumeIterator(const  Dune::CpGrid& grid, int cell_index)
        : grid_(&grid), cell_index_(cell_index)
    {}

    double dereference() const
    {
        return grid_->cellVolume(cell_index_);
    }
    void increment()
    {
        ++cell_index_;
    }
    double elementAt(int n) const
    {
        return grid_->cellVolume(n);
    }
    void advance(int n)
    {
        cell_index_+=n;
    }
    void decrement()
    {
        --cell_index_;
    }
    int distanceTo(const CellVolumeIterator& o) const
    {
        return o.cell_index_-cell_index_;
    }
    bool equals(const CellVolumeIterator& o) const
    {
        return o.grid_==grid_ && o.cell_index_==cell_index_;
    }

private:
    const Dune::CpGrid* grid_;
    int cell_index_;
};

template<>
struct CellVolumeIteratorTraits<Dune::CpGrid>
{
    typedef CellVolumeIterator IteratorType;
};

/// \brief Get an iterator over the cell volumes of a grid positioned at the first cell.
CellVolumeIterator beginCellVolumes(const Dune::CpGrid& grid);

/// \brief Get an iterator over the cell volumes of a grid positioned one after the last cell.
CellVolumeIterator endCellVolumes(const Dune::CpGrid& grid);

template<>
struct FaceCentroidTraits<Dune::CpGrid>
{
    typedef CpGridCentroidIterator<&Dune::CpGrid::faceCentroid> IteratorType;
    typedef const Dune::CpGrid::Vector ValueType;
};

/// \brief Get an iterator over the face centroids positioned at the first cell.
FaceCentroidTraits<Dune::CpGrid>::IteratorType
beginFaceCentroids(const Dune::CpGrid& grid);

/// \brief Get a coordinate of a specific face centroid.
/// \param grid The grid.
/// \param face_index The index of the specific face.
/// \param coordinate The coordinate index.
Vector faceCentroid(const Dune::CpGrid& grid, int face_index);

template<>
struct FaceCellTraits<Dune::CpGrid>
{
    typedef Dune::cpgrid::FaceCellsContainerProxy Type;
};
/// \brief Get the cell to faces mapping of a grid.
Dune::cpgrid::Cell2FacesContainer cell2Faces(const Dune::CpGrid& grid);

/// \brief Get the face to cell mapping of a grid.
FaceCellTraits<Dune::CpGrid>::Type
faceCells(const Dune::CpGrid& grid);

template<>
struct Face2VerticesTraits<Dune::CpGrid>
{
    typedef Dune::cpgrid::FaceVerticesContainerProxy Type;
};

/// \brief Get the face to vertices mapping of a grid.
Face2VerticesTraits<Dune::CpGrid>::Type
face2Vertices(const Dune::CpGrid& grid);

/// \brief Get the coordinates of a vertex of the grid.
/// \param grid The grid the vertex is part of.
/// \param index The index identifying the vertex.
Vector vertexCoordinates(const Dune::CpGrid& grid, int index);

Vector faceNormal(const Dune::CpGrid& grid, int face_index);

double faceArea(const Dune::CpGrid& grid, int face_index);

/// \brief Get Eclipse Cartesian tag of a face
/// \param grid The grid that the face is part of.
/// \param cell_face The face attached to a cell. Usually obtained from face2Cells.
/// \return 0, 1, 2, 3, 4, 5 for I-, I+, J-, J+, K-, K+
int faceTag(const Dune::CpGrid& grid,
            const Dune::cpgrid::Cell2FacesRow::iterator& cell_face);


} // end namespace UgGridHelpers

} // end namespace Opm

#endif
