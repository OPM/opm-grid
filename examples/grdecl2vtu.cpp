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

// Warning suppression for Dune includes.
#include <opm/core/utility/platform_dependent/disable_warnings.h>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <opm/core/utility/platform_dependent/reenable_warnings.h>

#include <dune/grid/CpGrid.hpp>
#include <opm/core/io/eclipse/EclipseGridInspector.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

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
 * @author H�vard Berland <havb@statoil.com>
 * @author Atgeirr F Rasmussen <atgeirr@sintef.no>
 * @author B�rd Skaflestad     <bard.skaflestad@sintef.no>
 *
 */

/**
   A function to (conditionally) write a double-field from the grdecl-file to vtk-format
*/
void condWriteDoubleField(std::vector<double> & fieldvector,
                          const std::string& fieldname,
                          Opm::DeckConstPtr deck,
                          const std::vector<int> & global_cell,
                          VTKWriter<CpGrid::LeafGridView> & vtkwriter) {
    if (deck->hasKeyword(fieldname)) {
        std::cout << "Found " << fieldname << "..." << std::endl;
        std::vector<double> eclVector = deck->getKeyword(fieldname)->getRawDoubleData();
        fieldvector.resize(global_cell.size());

        Opm::EclipseGridInspector insp(deck);
        std::array<int, 3> dims = insp.gridSize();
        int num_global_cells = dims[0]*dims[1]*dims[2];
        if (int(eclVector.size()) != num_global_cells) {
            OPM_THROW(std::runtime_error, fieldname << " field must have the same size as the "
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
                           Opm::DeckConstPtr deck,
                           const std::vector<int> & global_cell,
                           VTKWriter<CpGrid::LeafGridView> & vtkwriter) {
    if (deck->hasKeyword(fieldname)) {
        std::cout << "Found " << fieldname << "..." << std::endl;
        std::vector<int> eclVector = deck->getKeyword(fieldname)->getIntData();
        fieldvector.resize(global_cell.size());

        Opm::EclipseGridInspector insp(deck);
        std::array<int, 3> dims = insp.gridSize();
        int num_global_cells = dims[0]*dims[1]*dims[2];
        if (int(eclVector.size()) != num_global_cells) {
            OPM_THROW(std::runtime_error, fieldname << " field must have the same size as the "
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
try
{

    CpGrid grid;

    if (argc != 2) {
        std::cout << "Usage: grdecl2vtu filename.grdecl" << std::endl;
        exit(1);
    }

    const char* eclipsefilename = argv[1];
    Opm::ParserPtr parser(new Opm::Parser());
    Opm::DeckConstPtr deck(parser->parseFile(eclipsefilename));

    grid.processEclipseFormat(deck, 0.0, false, false);
    const std::vector<int>& global_cell = grid.globalCell();

    VTKWriter<CpGrid::LeafGridView> vtkwriter(grid.leafGridView());
    std::vector<double> poros;
    condWriteDoubleField(poros, "PORO", deck, global_cell, vtkwriter);

    std::vector<double> permxs;
    condWriteDoubleField(permxs, "PERMX", deck, global_cell, vtkwriter);
    std::vector<double> permys;
    condWriteDoubleField(permys, "PERMY", deck, global_cell, vtkwriter);

    std::vector<double> permzs;
    condWriteDoubleField(permzs, "PERMZ", deck, global_cell, vtkwriter);

    std::vector<double> actnums;
    condWriteIntegerField(actnums, "ACTNUM", deck, global_cell, vtkwriter);

    std::vector<double> satnums;
    condWriteIntegerField(satnums, "SATNUM", deck, global_cell, vtkwriter);

    std::vector<double> regnums;
    condWriteIntegerField(regnums, "REGNUM", deck, global_cell, vtkwriter);

    std::vector<double> swats;
    condWriteDoubleField(swats, "SWAT", deck, global_cell, vtkwriter);

    std::string fname(eclipsefilename);
    std::string fnamebase = fname.substr(0, fname.find_last_of('.'));
    std::cout << "Writing to filename " << fnamebase << ".vtu" << std::endl;
    vtkwriter.write(fnamebase, VTK::ascii);
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

