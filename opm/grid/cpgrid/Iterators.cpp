//===========================================================================
//
// File: Iterators.cpp
//
// Created: Mon February 20  14:02:00 2023
//
// Author(s): Antonella Ritorto <antonella.ritorto@opm-op.com> 
//            
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
Copyright 2023 Equinor ASA.

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


#include "config.h"
#include "Entity.hpp"
#include "Iterators.hpp"

#include <stack>


void Dune::cpgrid::HierarchicIterator::stackChildren_(const Entity<0>& target)
{
    // Load sons of target onto the iterator stack
    if (!target.isLeaf() && (target.level() < maxLevel_)){
        const auto& [lgr_level, children_list] = target.pgrid_-> parent_to_children_cells_[target.index()];
        // GET CHILD GRID
        const auto& lgr_grid =  (*(target.pgrid_-> level_data_ptr_))[lgr_level];
        for (const auto& child : children_list){
            this->elemStack_.push(Entity<0>(*lgr_grid, child, true));
        }
    }
}
void Dune::cpgrid::HierarchicIterator::resetEntity_()
{
    // Create an invalid entity, to set in case elemStack_ is empty.
    // Otherwise, a pointer pointing at the to element of elemStack_.
    virtualEntity_ = elemStack_.empty() ? Entity<0>(Entity<0>::InvalidIndex, true) : elemStack_.top();
}

