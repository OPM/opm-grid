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
    // When level zero grid has been distributed and some LGRs have been added afterwards, refined cells
    // inherit its parent cell attribute.
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
            int cell_index = cells_of_face[0].index();
            Entity<0> cell0(*grid_data_, cell_index, true);
            PartitionType cell_part = getPartitionType(cell0);
            if(cell_part!=OverlapEntity)
                return cell_part;
            else
            {
                // If the cell is in the overlap and the face is on the boundary,
                // then the partition type has to Front! Here we check whether
                // we are at the boundary.
                OrientedEntityTable<0,1>::row_type cell_to_face=grid_data_->cell_to_face_[cell0];
                Entity<0>::LeafIntersectionIterator intersection=cell0.ilevelbegin();
                for(int subindex=0; subindex<cell_to_face.size(); ++subindex, ++intersection)
                    if(cell_to_face[subindex].index()==i)
                        break;
                assert(intersection!=cell0.ilevelend());
                if(intersection.boundary())
                    return FrontEntity;
                else
                    return cell_part;
            }
        }
        else
        {
            Entity<0> cell0(*grid_data_, cells_of_face[0].index(), true);
            Entity<0> cell1(*grid_data_, cells_of_face[1].index(), true);
            if(cells_of_face[0].index()==std::numeric_limits<int>::max())
            {
                assert(cells_of_face[1].index()!=std::numeric_limits<int>::max());
                // At the boder of the processor's but not the global domain
                return getProcessorBoundaryPartitionType(getPartitionType(cell1));
            }
            if(cells_of_face[1].index()==std::numeric_limits<int>::max())
            {
                assert(cells_of_face[0].index()!=std::numeric_limits<int>::max());
                // At the boder of the processor's but not the global domain
                return getProcessorBoundaryPartitionType(getPartitionType(cell0));
            }
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
