#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "Entity.hpp"
#include "CpGridData.hpp"
#include "PartitionTypeIndicator.hpp"

namespace Dune
{
namespace cpgrid
{
PartitionType PartitionTypeIndicator::getPartitionType(const Entity<0>& cell_entity) const
{
    if(cell_indicator_.size())
        return PartitionType(cell_indicator_[cell_entity.index()]);
    // If level zero grid has been distributed and some LGRs have been added, refined cells
    // inherit its parent cell partition type.
    if (grid_data_->level_ >0) { // level_ > 0 only for refined level grids.
        return PartitionType(grid_data_->level_data_ptr_->front()->partition_type_indicator_->getPartitionType(cell_entity.getOrigin()));
    }
    return InteriorEntity;
}

PartitionType PartitionTypeIndicator::getPartitionType(const EntityRep<1>& face_entity) const
{
    return getFacePartitionType(face_entity.index());
}
PartitionType PartitionTypeIndicator::getPartitionType(const EntityRep<3>& point_entity) const
{
    return getPointPartitionType(point_entity.index());
}
PartitionType PartitionTypeIndicator::getPointPartitionType(int index) const
{
    if(point_indicator_.size())
        return PartitionType(point_indicator_[index]);
    return InteriorEntity;
}

PartitionType getProcessorBoundaryPartitionType(PartitionType)
{
    return FrontEntity;
}

PartitionType PartitionTypeIndicator::getFacePartitionType(int i) const
{
    if((cell_indicator_.size()) || (grid_data_->level_ > 0))
    {
        // We determine the partition type by the type of the
        // connected cells:
        // If all of them are interior and border, then the type is
        // interior and border, respectively.
        // If one of them is of type interior and the other is
        // of type overlap, then the type is border.
        OrientedEntityTable<1,0>::row_type cells_of_face =
            grid_data_->face_to_cell_[Entity<1>(*grid_data_,i,true)];
        if(cells_of_face.size()==1)
        {
            // We are on the global domain boundary
            // partition type of intersection is the same as the one of the cell.
            int cell_index = cells_of_face[0].index();
            Entity<0> cell0(*grid_data_, cell_index, true);
            PartitionType cell_part = getPartitionType(cell0);
            return cell_part;
        } else {
            // One of the cells might have index std::numeric_limits<int>::max().
            // That means that the neighbor is not on this process but on
            // another one. In this case the intersection is Front
            auto idx0 = cells_of_face[0].index();
            auto idx1 = cells_of_face[1].index();
            if (idx0 == std::numeric_limits<int>::max() || idx1 == std::numeric_limits<int>::max()) {
                // One neighbor is on another process => FrontEntity
                return FrontEntity;
            }

            Entity<0> cell0(*grid_data_, idx0, true);
            Entity<0> cell1(*grid_data_, idx1, true);
            if(getPartitionType(cell0) == getPartitionType(cell1))
                return getPartitionType(cell0);
            else
                return BorderEntity;
        }
    }
    return InteriorEntity;
}
} // end namespace cpgrid
} // end namespace Dune
