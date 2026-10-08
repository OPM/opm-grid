//===========================================================================
//
// File: CpGrid.cpp
//
// Created: Thu Jun  4 12:55:28 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            B�rd Skaflestad     <bard.skaflestad@sintef.no>
//            Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2023 Equinor ASA.
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#if HAVE_MPI
#include <opm/grid/utility/platform_dependent/disable_warnings.h>
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>
#endif

#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>

#include "../CpGrid.hpp"
#include "LgrHelpers.hpp"
#include "ParentToChildrenCellGlobalIdHandle.hpp"
#include "NestedRefinementUtilities.hpp"
#include <opm/grid/common/MetisPartition.hpp>
#include <opm/grid/common/ZoltanPartition.hpp>
#include <opm/grid/GraphOfGridWrappers.hpp>
//#include <opm/grid/common/ZoltanGraphFunctions.hpp>
#include <opm/grid/common/GridPartitioning.hpp>
//#include <opm/grid/common/WellConnections.hpp>
#include <opm/grid/common/CommunicationUtils.hpp>

//#include <fstream>
//#include <iostream>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <tuple>

namespace Dune
{

void CpGrid::computeGlobalCellLgr(const int& level, const std::array<int,3>& startIJK, std::vector<int>& global_cell_lgr)
{
    assert(level);
    for (const auto& element : elements(levelGridView(level))) {
        // Element belogns to an LGR, therefore has a father. Get IJK of the father in the level grid the father was born.
        // For CARFIN, parent cells belong to level 0.
        std::array<int,3> parentIJK = {0,0,0};
        currentData()[element.father().level()]->getIJK(element.father().index(), parentIJK);
        // Each parent cell has been refined in cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2] child cells.
        // element has certain 'position' inside its parent cell that can be described with 'IJK' indices, let's denote them by ijk,
        // where 0<= i < cells_per_dim[0], 0<= j < cells_per_dim[1], 0<= k < cells_per_dim[2].
        const auto& cells_per_dim = currentData()[level]->levels_.cells_per_dim;
        //
        // Refined cell (here 'element') has "index in parent cell": k*cells_per_dim[0]*cells_per_dim[1] + j*cells_per_dim[0] + i
        // and it's stored in  cell_to_idxInParentCell_.
        auto idx_in_parent_cell =  currentData()[level]-> levels_.cell_to_idxInParentCell[element.index()];
        // Find ijk.
        std::array<int,3> childIJK = Opm::Lgr::getIJK(idx_in_parent_cell, cells_per_dim);
        // The corresponding lgrIJK can be computed as follows:
        const std::array<int,3>& lgrIJK = { ( (parentIJK[0] - startIJK[0])*cells_per_dim[0] ) + childIJK[0],  // Shift parent index according to the startIJK of the LGR.
                                            ( (parentIJK[1] - startIJK[1])*cells_per_dim[1] ) + childIJK[1],
                                            ( (parentIJK[2] - startIJK[2])*cells_per_dim[2] ) + childIJK[2] };
        // Dimensions of the "patch of cells" formed when providing startIJK and endIJK for an LGR
        const auto& lgr_logical_cartesian_size = currentData()[level]->logical_cartesian_size_;
        global_cell_lgr[element.index()] = (lgrIJK[2]*lgr_logical_cartesian_size[0]*lgr_logical_cartesian_size[1]) + (lgrIJK[1]*lgr_logical_cartesian_size[0]) + lgrIJK[0];
    }
}

void CpGrid::computeGlobalCellLeafGridViewWithLgrs(std::vector<int>& global_cell_leaf)
{
    for (const auto& element: elements(leafGridView())) {
        // In the context of allowed nested refinement, we lookup for the oldest ancestor, belonging to level-zero-grid.
        auto ancestor = element.getOrigin();
        int origin_in_level_zero = ancestor.index();
        assert(origin_in_level_zero < currentData().front()->size(0));
        global_cell_leaf[element.index()] = currentData().front()-> global_cell_[origin_in_level_zero];
    }
}

std::vector<std::unordered_map<std::size_t, std::size_t>> CpGrid::mapLocalCartesianIndexSetsToLeafIndexSet() const
{
    std::vector<std::unordered_map<std::size_t, std::size_t>> localCartesianIdxSets_to_leafIdx(maxLevel()+1); // Plus level 0
    for (const auto& element : elements(leafGridView())) {
        const auto& global_cell_level = currentData()[element.level()]->globalCell()[element.getLevelElem().index()];
        localCartesianIdxSets_to_leafIdx[element.level()][global_cell_level] = element.index();
    }
    return localCartesianIdxSets_to_leafIdx;
}

std::vector<std::array<int,2>> CpGrid::mapLeafIndexSetToLocalCartesianIndexSets() const
{
    std::vector<std::array<int,2>> leafIdx_to_localCartesianIdxSets(currentLeafData().size(0));
    for (const auto& element : elements(leafGridView())) {
        const auto& global_cell_level = currentData()[element.level()]->globalCell()[element.getLevelElem().index()];
        leafIdx_to_localCartesianIdxSets[element.index()] = {element.level(), global_cell_level};
    }
    return leafIdx_to_localCartesianIdxSets;
}

Dune::cpgrid::Intersection CpGrid::getParentIntersectionFromLgrBoundaryFace(const Dune::cpgrid::Intersection& intersection) const
{
    if ( intersection.neighbor()) {
        // Handle only intersections between cells at different levels where at least one
        // cell belongs to a grid whose parent is GLOBAL.
        // Skip intersections at the same level, as well as those between refined grids
        // whose parent grids are not GLOBAL (i.e., not the level-zero grid).
        int levelIn = intersection.inside().level();
        int levelOut = intersection.outside().level();
        if (levelIn != levelOut) {
            bool refinedCellsWithAtLeastOneLvl0Father = (levelIn*levelOut >0 ) &&
                ((intersection.inside().father().level() == 0) || ((intersection.outside().father().level() == 0)));
            bool coarseAndRefinedCells =  (levelIn*levelOut == 0);
            if (refinedCellsWithAtLeastOneLvl0Father || coarseAndRefinedCells) {
                // Get the equivalent level-0 cell if intersection.inside() is already at level 0,
                // or the coarsest ancestor at level 0 if it is a refined cell.
                // In both cases, intersection.indexInInside() is the correct face index to match:
                // - coarse inside (level 0): indexInInside() is already the level-0 face index.
                // - refined inside: fine cells inherit face directions from their parent, so
                //   indexInInside() equals the parent's face index.
                const auto& insideOrigin = intersection.inside().getOrigin();
                for (const auto& originIntersection : intersections(this->levelGridView(0), insideOrigin)) {
                    if (originIntersection.indexInInside() == intersection.indexInInside()) {
                        return originIntersection;
                    }
                }
            }
            OPM_THROW(std::invalid_argument, "Parent intersection not found for face with index: " + std::to_string(intersection.id()) +
                      " and index in inside: " + std::to_string(intersection.indexInInside()));
        }
    }
    OPM_THROW(std::invalid_argument, "Face is on the boundary of the grid");
}

const std::map<std::string,int>& CpGrid::getLgrNameToLevel() const{
    return lgr_names_;
}

} // namespace Dune
