//===========================================================================
//
// File: PartitionTypeIndicator.hpp
//
// Created: Oct 20 2013
//
// Author(s): Markus Blatt <markus@dr-blatt.de>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
Copyright 2013 Dr. Blatt - HPC-Simulation-Software & Service
Copyright 2013 Statoil ASA.

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

#ifndef OPM_PARTITIONTYPEINDICATOR_HEADER
#define OPM_PARTITIONTYPEINDICATOR_HEADER

#include<vector>
#include <dune/grid/common/gridenums.hh>

namespace Dune
{
namespace cpgrid
{
class CpGridData;
template<int> class Entity;
template<int> class EntityRep;

class PartitionTypeIndicator
{
public:
    /// Constructor
    /// \param data The data of the cornerpoint grid.
    PartitionTypeIndicator(const CpGridData& data)
    : grid_data_(&data)
    {}
    /// Get the partition type of a cell.
    /// \param cell_entity The entity describing the cell
    /// \return The partition type of the cell.
    PartitionType getPartitionType(const EntityRep<0>& cell_entity) const;
    /// Get the partition type of a face.
    /// \param face_entity The entity describing the face
    /// \return The partition type of the face.
    PartitionType getPartitionType(const EntityRep<1>& face_entity) const;
    /// Get the partition type of a point.
    /// \param point_entity The entity describing the point.
    /// \return The partition type of the point.
    PartitionType getPartitionType(const EntityRep<3>& point_entity) const;
    
private:
    /// Get the partition type of a face by its index
    /// \param i The index of the face.
    /// \return The partition type of the face associated with this index.
    PartitionType getFacePartitionType(int i) const;
    
    
    /// Get the partition type of a face by its index
    /// \param i The index of the face.
    /// \return The partition type of the face associated with this index.
    PartitionType getPointPartitionType(int i) const;

    /// The data of the grid.
    const CpGridData* grid_data_;
    /// An array to store the partition type of cell.
    ///
    /// If non-empty, then the cell with index i has (PartitionType)cell_indicator_[i].
    /// Otherwise this grid is not parallel and allen entities are interior.
    std::vector<char> cell_indicator_;
    /// An array to store the partition type of cell.
    ///
    /// If non-empty, then the point with index i has (PartitionType)cell_indicator_[i].
    /// Otherwise this grid is not parallel and allen entities are interior.
    std::vector<char> point_indicator_;
    friend class CpGridData;
    friend class FacePartitionTypeIterator;
};
} // end namespace Dune
} // end namespace cpgrid

#endif
