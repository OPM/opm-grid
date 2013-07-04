/*
  Copyright 2009, 2010, 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010, 2013 Statoil ASA.

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

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/CpGrid.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/io/eclipse/EclipseGridInspector.hpp>

using namespace Dune;

/**
 * @file grdecl2vtu.C
 * @brief Converts grdecl (Eclipse grid) files to vtu (VTK/ParaView)
 * 
 * Converts a corner-point grid with properties to a vtu-file
 * (to be opened in ParaView for example)
 *
 * Based on make_vtk_test.cpp
 *
 * @author Håvard Berland <havb@statoil.com>
 * @author Atgeirr F Rasmussen <atgeirr@sintef.no>
 * @author Bård Skaflestad     <bard.skaflestad@sintef.no>
 *
 */

/**
   A function to (conditionally) write a double-field from the grdecl-file to vtk-format
*/
void condWriteDoubleField(std::vector<double> & fieldvector,
                          const std::string& fieldname,
                          const Opm::EclipseGridParser & eclParser,
                          const std::vector<int> & global_cell,
                          VTKWriter<CpGrid::LeafGridView> & vtkwriter) {
    if (eclParser.hasField(fieldname)) {
        std::cout << "Found " << fieldname << "..." << std::endl;
        std::vector<double> eclVector = eclParser.getFloatingPointValue(fieldname);
        fieldvector.resize(global_cell.size());

	Opm::EclipseGridInspector insp(eclParser);
        boost::array<int, 3> dims = insp.gridSize();
        int num_global_cells = dims[0]*dims[1]*dims[2];
        if (int(eclVector.size()) != num_global_cells) {
            THROW(fieldname << " field must have the same size as the "
                  "logical cartesian size of the grid: "
                  << eclVector.size() << " != " << num_global_cells);
        }

        for (size_t i = 0; i < global_cell.size(); ++i) {
            fieldvector[i] = eclVector[global_cell[i]];
        }
        vtkwriter.addCellData(fieldvector, fieldname);
    }
    
}
// Now repeat for Integers. I should learn C++ templating...
void condWriteIntegerField(std::vector<double> & fieldvector,
                           const std::string& fieldname,
                           const Opm::EclipseGridParser & eclParser,
                           const std::vector<int> & global_cell,
                           VTKWriter<CpGrid::LeafGridView> & vtkwriter) {
    if (eclParser.hasField(fieldname)) {
        std::cout << "Found " << fieldname << "..." << std::endl;
        std::vector<int> eclVector = eclParser.getIntegerValue(fieldname);
        fieldvector.resize(global_cell.size());

	Opm::EclipseGridInspector insp(eclParser);
        boost::array<int, 3> dims = insp.gridSize();
        int num_global_cells = dims[0]*dims[1]*dims[2];
        if (int(eclVector.size()) != num_global_cells) {
            THROW(fieldname << " field must have the same size as the "
                  "logical cartesian size of the grid: "
                  << eclVector.size() << " != " << num_global_cells);
        }

        for (size_t i = 0; i < global_cell.size(); ++i) {
            fieldvector[i] = (double)eclVector[global_cell[i]];
        }
                vtkwriter.addCellData(fieldvector, fieldname);
    }

}


int main(int argc, char** argv)
{

    CpGrid grid;
    
    if (argc != 2) {
        std::cout << "Usage: grdecl2vtu filename.grdecl" << std::endl;
        exit(1);
    }
    
    const char* eclipsefilename = argv[1];
    bool convert_to_SI = false;
    Opm::EclipseGridParser eclParser(eclipsefilename, convert_to_SI);

    grid.processEclipseFormat(eclParser, 0.0, false, false);
    const std::vector<int>& global_cell = grid.globalCell();
    
    VTKWriter<CpGrid::LeafGridView> vtkwriter(grid.leafView());
    std::vector<double> poros;
    condWriteDoubleField(poros, "PORO", eclParser, global_cell, vtkwriter);
    
    std::vector<double> permxs;
    condWriteDoubleField(permxs, "PERMX", eclParser, global_cell, vtkwriter);
    std::vector<double> permys;
    condWriteDoubleField(permys, "PERMY", eclParser, global_cell, vtkwriter);

    std::vector<double> permzs;
    condWriteDoubleField(permzs, "PERMZ", eclParser, global_cell, vtkwriter);

    std::vector<double> actnums;
    condWriteIntegerField(actnums, "ACTNUM", eclParser, global_cell, vtkwriter);

    std::vector<double> satnums;
    condWriteIntegerField(satnums, "SATNUM", eclParser, global_cell, vtkwriter);

    std::vector<double> regnums;
    condWriteIntegerField(regnums, "REGNUM", eclParser, global_cell, vtkwriter);
   
    std::vector<double> swats;
    condWriteDoubleField(swats, "SWAT", eclParser, global_cell, vtkwriter);

    std::string fname(eclipsefilename);
    std::string fnamebase = fname.substr(0, fname.find_last_of('.'));
    std::cout << "Writing to filename " << fnamebase << ".vtu" << std::endl;
    vtkwriter.write(fnamebase, VTKOptions::ascii);
}

