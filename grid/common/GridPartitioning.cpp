//===========================================================================
//
// File: GridPartitioning.cpp
//
// Created: Mon Sep  7 10:18:28 2009
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
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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

#include "config.h"
#include "GridPartitioning.hpp"
#include <dune/grid/CpGrid.hpp>
#include <stack>

namespace Dune
{


    typedef boost::array<int, 3> coord_t;

    namespace
    {

	struct IndexToIJK
	{
	    IndexToIJK(const coord_t& lc_size)
		: num_i(lc_size[0]),
		  num_ij(lc_size[0]*lc_size[1])
	    {
	    }
	    coord_t operator()(int index)
	    {
		coord_t retval = {{ index % num_i,
				    (index % num_ij) / num_i,
				    index / num_ij }};
		return retval;
	    }
	private:
	    int num_i;
	    int num_ij;
	};


	int initialPartition(const coord_t& c, const coord_t& lc_size, const coord_t& initial_split) 
	{
	    coord_t p_coord;
	    for (int i = 0; i < 3; ++i) {
		int n = lc_size[i]/initial_split[i];
		int extra = lc_size[i] % initial_split[i];
		if (c[i] < (n+1)*extra) {
		    p_coord[i] = c[i]/(n+1);
		} else {
		    p_coord[i] = (c[i] - (n+1)*extra)/n + extra;
		}
	    }
	    return p_coord[0] + initial_split[0]*(p_coord[1] + initial_split[1]*p_coord[2]);
	}


	void colourMyComponentRecursive(const CpGrid& grid,
					const CpGrid::Codim<0>::EntityPointer& c,
					const int colour,
					const std::vector<int>& cell_part,
					std::vector<int>& cell_colour)
	{
	    const CpGrid::LeafIndexSet& ix = grid.leafIndexSet();
	    int my_index = ix.index(*c);
	    cell_colour[my_index] = colour;
	    // For each neighbour...
	    for (CpGrid::LeafIntersectionIterator it = c->ileafbegin(); it != c->ileafend(); ++it) {
		if (it->neighbor()) {
		    int nb_index = ix.index(*(it->outside()));
		    if (cell_part[my_index] == cell_part[nb_index] && cell_colour[nb_index] == -1) {
			colourMyComponentRecursive(grid, it->outside(), colour, cell_part, cell_colour);
		    }
		}
	    }
	}


	void colourMyComponent(const CpGrid& grid,
			       const CpGrid::Codim<0>::EntityPointer& c,
			       const int colour,
			       const std::vector<int>& cell_part,
			       std::vector<int>& cell_colour)
	{
	    typedef CpGrid::LeafIntersectionIterator NbIter;
	    typedef std::pair<int, std::pair<NbIter, NbIter> > VertexInfo;
	    std::stack<VertexInfo> v_stack;
	    const CpGrid::LeafIndexSet& ix = grid.leafIndexSet();
	    int index = ix.index(*c);
	    cell_colour[index] = colour;
	    NbIter cur = c->ileafbegin();
	    NbIter end = c->ileafend();
	    v_stack.push(std::make_pair(index, std::make_pair(cur, end)));
	    while (!v_stack.empty()) {
		index = v_stack.top().first;
		cur = v_stack.top().second.first;
		end = v_stack.top().second.second;
		v_stack.pop();
		while (cur != end) {
		    bool visit_nb = false;
		    if (cur->neighbor()) {
			int nb_index = ix.index(*(cur->outside()));
			if (cell_part[index] == cell_part[nb_index] && cell_colour[nb_index] == -1) {
			    visit_nb = true;
			}
		    }
		    if (visit_nb) {
			NbIter cur_cp = cur;
			v_stack.push(std::make_pair(index, std::make_pair(++cur, end)));
			index = ix.index(*(cur_cp->outside()));
			cur = cur_cp->outside()->ileafbegin();
			end = cur_cp->outside()->ileafend();
			cell_colour[index] = colour;
		    } else {
			++cur;
		    }
		}
	    }
	}


	void ensureConnectedPartitions(const CpGrid& grid,
				       int& num_part,
				       std::vector<int>& cell_part,
				       bool recursive = false)
	{
	    std::vector<int> cell_colour(cell_part.size(), -1);
	    std::vector<int> partition_used(num_part, 0);
	    int max_part = num_part;
	    const CpGrid::LeafIndexSet& ix = grid.leafIndexSet();
	    for (CpGrid::Codim<0>::LeafIterator it = grid.leafbegin<0>(); it != grid.leafend<0>(); ++it) {
		int index = ix.index(*it);
		if (cell_colour[index] == -1) {
		    int part = cell_part[index];
		    int current_colour = part;
		    if (partition_used[part]) {
			current_colour = max_part++;
		    } else {
			partition_used[part] = true;
		    }
		    if (recursive) {
			colourMyComponentRecursive(grid, it, current_colour, cell_part, cell_colour);
		    } else {
			colourMyComponent(grid, it, current_colour, cell_part, cell_colour);
		    }
		}
	    }
	    if (max_part != num_part) {
		num_part = max_part;
		cell_part.swap(cell_colour);
	    }
	}

    } // anon namespace


    void partition(const CpGrid& grid,
		   const coord_t& initial_split,
		   int& num_part,
		   std::vector<int>& cell_part,
		   bool recursive)
    {
	// Checking that the initial split makes sense (that there may be at least one cell
	// in each expected partition).
	const coord_t& lc_size = grid.logicalCartesianSize();
	for (int i = 0; i < 3; ++i) {
	    if (initial_split[i] > lc_size[i]) {
		THROW("In direction " << i << " requested splitting " << initial_split[i] << " size " << lc_size[i]);
	    }
	}

	// Initial partitioning depending on (ijk) coordinates.
	const int num_initial = initial_split[0]*initial_split[1]*initial_split[2];
	const std::vector<int>& lc_ind = grid.globalCell();
	std::vector<int> num_in_part(num_initial, 0);
	std::vector<int> my_part(grid.size(0), -1);
	IndexToIJK ijk_coord(lc_size);
	for (int i = 0; i < grid.size(0); ++i) {
	    coord_t ijk = ijk_coord(lc_ind[i]);
	    int part = initialPartition(ijk, lc_size, initial_split);
	    my_part[i] = part;
	    ++num_in_part[part];
	}

	// Renumber partitions.
	std::vector<int> num_to_subtract(num_initial);
	num_to_subtract[0] = 0;
	for (int i = 1; i < num_initial; ++i) {
	    num_to_subtract[i] = num_to_subtract[i-1];
	    if (num_in_part[i-1] == 0) {
		++num_to_subtract[i];
	    }
	}
	for (int i = 0; i < grid.size(0); ++i) {
	    my_part[i] -= num_to_subtract[my_part[i]];
	}

	num_part = num_initial - num_to_subtract.back();
	cell_part.swap(my_part);

	// Check the connectivity, split.
	ensureConnectedPartitions(grid, num_part, cell_part, recursive);
    }

} // namespace Dune

