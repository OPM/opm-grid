//===========================================================================
//
// File: check_grid_normals.cpp
//
// Created: Tue Jul  6 09:54:41 2010
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
  Copyright 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2010 Statoil ASA.

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



#if HAVE_CONFIG_H
#include "config.h"
#endif
#include <dune/common/param/ParameterGroup.hpp>
#include <dune/grid/CpGrid.hpp>

using namespace Dune;

int main(int argc, char** argv)
{
    parameter::ParameterGroup param(argc, argv);
    CpGrid grid;
    grid.init(param);
    typedef CpGrid::LeafGridView View;
    View g = grid.leafView();
    typedef FieldVector<double, 3> Pt;
    int c_local = 0;
    for (View::Codim<0>::Iterator c_it = g.begin<0>(); c_it != g.end<0>(); ++c_it, ++c_local) {
        Pt cell_centroid = c_it->geometry().center();
        int f_local = 0;
        bool trouble = false;
        for (View::IntersectionIterator f_it = g.ibegin(*c_it); f_it != g.iend(*c_it); ++f_it, ++f_local) {
            Pt face_centroid = f_it->geometry().center();
            Pt face_normal = f_it->centerUnitOuterNormal();
            if (face_normal*(face_centroid - cell_centroid) < 0.0) {
                trouble = true;
                std::cout << "Encountered troublesome geometry (centroid difference dot normal is negative) "
                    "in cell " << c_local << " local face " << f_local << std::endl;
            }
        }
        if (trouble) {
            std::cout << "Cell " << c_local << " had a total of " << f_local << " faces." << std::endl;
        }
    }
}

