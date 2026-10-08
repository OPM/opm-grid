#include "config.h"
#include <algorithm>
#include <array>
#include <map>
#include <set>
#include <vector>
#include <utility>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/DataHandleWrappers.hpp>
#include <opm/grid/cpgrid/ElementMarkHandle.hpp>
#include <opm/grid/cpgrid/Intersection.hpp>
#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/cpgrid/LgrHelpers.hpp>
#include <opm/grid/cpgrid/OrientedEntityTable.hpp>
#include <opm/grid/cpgrid/Indexsets.hpp>
#include <opm/grid/cpgrid/PartitionTypeIndicator.hpp>

// Warning suppression for Dune includes.
#include <opm/grid/utility/platform_dependent/disable_warnings.h>

#include <opm/grid/common/GridPartitioning.hpp>
#include <dune/common/parallel/remoteindices.hh>
#include <dune/common/enumset.hh>
#include <opm/common/utility/SparseTable.hpp>

#include <opm/grid/utility/platform_dependent/reenable_warnings.h>
#include <opm/grid/CpGrid.hpp>

namespace Dune
{
namespace cpgrid
{

std::array<Dune::FieldVector<double,3>,8> CpGridData::getReferenceRefinedCorners(int idx_in_parent_cell, const std::array<int,3>& cells_per_dim) const
{
    // Refined cells in parent cell: k*cells_per_dim[0]*cells_per_dim[1] + j*cells_per_dim[0] + i
    std::array<int,3> ijk = Opm::Lgr::getIJK(idx_in_parent_cell, cells_per_dim);

    std::array<Dune::FieldVector<double,3>,8> corners_in_parent_reference_elem = { // corner '0'
        {{ double(ijk[0])/cells_per_dim[0], double(ijk[1])/cells_per_dim[1], double(ijk[2])/cells_per_dim[2] },
         // corner '1'
         { double(ijk[0]+1)/cells_per_dim[0], double(ijk[1])/cells_per_dim[1], double(ijk[2])/cells_per_dim[2] },
         // corner '2'
         { double(ijk[0])/cells_per_dim[0], double(ijk[1]+1)/cells_per_dim[1], double(ijk[2])/cells_per_dim[2] },
         // corner '3'
         { double(ijk[0]+1)/cells_per_dim[0], double(ijk[1]+1)/cells_per_dim[1], double(ijk[2])/cells_per_dim[2] },
         // corner '4'
         { double(ijk[0])/cells_per_dim[0], double(ijk[1])/cells_per_dim[1], double(ijk[2]+1)/cells_per_dim[2] },
         // corner '5'
         { double(ijk[0]+1)/cells_per_dim[0], double(ijk[1])/cells_per_dim[1], double(ijk[2]+1)/cells_per_dim[2] },
         // corner '6'
         { double(ijk[0])/cells_per_dim[0], double(ijk[1]+1)/cells_per_dim[1], double(ijk[2]+1)/cells_per_dim[2] },
         // corner '7'
         { double(ijk[0]+1)/cells_per_dim[0], double(ijk[1]+1)/cells_per_dim[1], double(ijk[2]+1)/cells_per_dim[2] }
        }
    };
    return corners_in_parent_reference_elem;
}

bool CpGridData::mark(int refCount, const cpgrid::Entity<0>& element, bool throwOnFailure)
{
    if (refCount == -1) {
        if (throwOnFailure)
            OPM_THROW(std::logic_error, "Coarsening is not supported yet.");
        return false; // Coarsening is not supported yet.
    }
    // Prevent refinement if the cell has a non-neighbor connection (NNC).
    if (hasNNCs({element.index()}) && (refCount == 1)) {
        if (throwOnFailure)
            OPM_THROW(std::logic_error, "Refinement of cells with face representing an NNC is not supported yet.");
        return false;
    }
    assert((refCount == 0) || (refCount == 1)); // Do nothing (0), Refine (1), Coarsen (-1) not supported yet.
    if (mark_.empty()) {
        mark_.resize(this->size(0));
    }
    mark_[element.index()] = refCount;
    return (mark_[element.index()] == refCount);
}

int CpGridData::getMark(const cpgrid::Entity<0>& element) const
{
    return mark_.empty() ? 0 : mark_[element.index()];
}

bool CpGridData::preAdapt()
{
    // Communicate marked elements across all processes.
    if (ccobj_.size()>1) {

        auto local_empty = mark_.empty();
        // The attribute mark_ can be empty in processes with no elements marked
        // for refinement. In that case, resize before communication occurs.
        if (ccobj_.max(!local_empty)){
            if (local_empty)
                mark_.resize(size(0));
        }

        // Detect the maximum mark across processes, and rewrite
        // the local entry in mark_, i.e.,
        // mark_[ element.index() ] = max{ local marks in processes where this element belongs to}.
        ElementMarkHandle element_mark_handle(mark_);

        // An element may be marked somewhere in opm-simulators because it does not fulfill a
        // certain property, regardless of whether it belongs to the interior or overlap
        // partition. Therefore, we use the All_All_Interface.
        communicate(element_mark_handle,
                    Dune::All_All_Interface,
                    Dune::ForwardCommunication);
    }

    if(mark_.empty()) {
        return false;
    }
    else {
        for (int elemIdx = 0; elemIdx <  this-> size(0); ++elemIdx) {
            const auto& element = Dune::cpgrid::Entity<0>(*this, elemIdx, true);
            if (getMark(element) != 0)  // 1 (to be refined), 0 (do nothing), -1 (to be coarsened - not supported yet)
                return true;
        }
    }
    return false;
}

bool CpGridData::adapt()
{
    return preAdapt();
}

void CpGridData::postAdapt()
{
    mark_.resize(this->size(0), 0);
}

} // namespace cpgrid
} // namespace Dune
