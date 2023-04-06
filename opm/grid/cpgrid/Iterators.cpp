//===========================================================================
//
// File: Iterators.cpp
//
// Created: Mon February 20  14:02:00 2023
//
// Author(s): Antonella Ritorto <antonella.ritortoopm-op.com>  ??? TO BE DOUBLE-CHECKED 
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
    if (target.level() < maxLevel_ && !target.isLeaf()){
        const auto& children_list_indices = std::get<1>(target.pgrid_ -> parent_to_children_cells_[target.index()]);
        // GET CHILD GRID
        for (const auto& child : children_list_indices){
            this->elemStack_.push(Entity<0>(*(target.pgrid_), child, true)); // CORRECT THE CPGRIDDATA, GET THE CHILD-GRID
        }
    }
}
void Dune::cpgrid::HierarchicIterator::resetEntity_()
{
    // Create an invalid entity, to set in case elemStack_ is empty.
    // Otherwise, a pointer pointing at the to element of elemStack_.
    virtualEntity_ = elemStack_.empty() ? Entity<0>() : elemStack_.top();
}

