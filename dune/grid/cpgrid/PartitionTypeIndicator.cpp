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
PartitionType PartitionTypeIndicator::getPartitionType(const EntityRep<0>& cell_entity) const
{
    if(cell_indicator_.size())
        return PartitionType(cell_indicator_[cell_entity.index()]);
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


PartitionType getProcessorBoundaryPartitionType(PartitionType cell_partition_type)
{
    if(cell_partition_type==InteriorEntity)
        return BorderEntity;
    else
        return FrontEntity;
}

PartitionType PartitionTypeIndicator::getFacePartitionType(int i) const
{
    if(cell_indicator_.size())
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
            return PartitionType(cell_indicator_[cells_of_face[0].index()]);
        }
        else
        {
            if(cells_of_face[0].index()==std::numeric_limits<int>::max())
            {
                // At the boder of the processor's but not the global domain
                return getProcessorBoundaryPartitionType(PartitionType(cell_indicator_[cells_of_face[1].index()]));
            }
            if(cells_of_face[1].index()==std::numeric_limits<int>::max())
            {
                // At the boder of the processor's but not the global domain
                return getProcessorBoundaryPartitionType(PartitionType(cell_indicator_[cells_of_face[0].index()]));
            }
                
            if(cell_indicator_[cells_of_face[0].index()]==
               cell_indicator_[cells_of_face[1].index()])
                return PartitionType(cell_indicator_[cells_of_face[0].index()]);
            else
                return BorderEntity;
        }
    }
    return InteriorEntity;
}
} // end namespace cpgrid
} // end namespace Dune
