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
#include <opm/grid/utility/platform_dependent/disable_warnings.h>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>

#include <opm/grid/CpGrid.hpp>

#include <opm/grid/utility/OpmParserIncludes.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>

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
template <class OpmDeck>
void condWriteDoubleField(std::vector<double>& fieldvector,
                          const std::string& fieldname,
                          const OpmDeck& deck,
                          const std::vector<int>& global_cell,
                          const std::array<size_t, 3>& dims,
                          VTKWriter<CpGrid::LeafGridView>& vtkwriter) {
    if (deck.hasKeyword(fieldname)) {
        std::cout << "Found " << fieldname << "..." << std::endl;
        std::vector<double> eclVector = deck[fieldname].back().getRawDoubleData();
        fieldvector.resize(global_cell.size());
        int num_global_cells = dims[0]*dims[1]*dims[2];
        if (int(eclVector.size()) != num_global_cells) {
            OPM_THROW(std::runtime_error,
                      fieldname + " field must have the same size as the "
                      "logical cartesian size of the grid: " +
                      std::to_string(eclVector.size()) + " != " +
                      std::to_string(num_global_cells));
        }

        for (size_t i = 0; i < global_cell.size(); ++i) {
            fieldvector[i] = eclVector[global_cell[i]];
        }
        vtkwriter.addCellData(fieldvector, fieldname);
    }

}

// Now repeat for Integers. I should learn C++ templating...
template <class OpmDeck>
void condWriteIntegerField(std::vector<double>& fieldvector,
                           const std::string& fieldname,
                           const OpmDeck& deck,
                           const std::vector<int>& global_cell,
                           const std::array<size_t, 3>& dims,
                           VTKWriter<CpGrid::LeafGridView>& vtkwriter) {
    if (deck.hasKeyword(fieldname)) {
        std::cout << "Found " << fieldname << "..." << std::endl;
        std::vector<int> eclVector = deck[fieldname].back().getIntData();
        fieldvector.resize(global_cell.size());
        int num_global_cells = dims[0]*dims[1]*dims[2];
        if (int(eclVector.size()) != num_global_cells) {
            OPM_THROW(std::runtime_error,
                      fieldname + " field must have the same size as the "
                      "logical cartesian size of the grid: " +
                      std::to_string(eclVector.size()) + " != " +
                      std::to_string(num_global_cells));
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
    Dune::MPIHelper::instance(argc,argv); // Dummy if no MPI.

    CpGrid grid;

    if (argc != 2) {
        std::cout << "Usage: grdecl2vtu filename.grdecl" << std::endl;
        exit(1);
    }

    const char* eclipsefilename = argv[1];
    Opm::Parser parser;
    auto deck = parser.parseFile(eclipsefilename);

    // Get logical cartesian grid dimensions.
    std::array<size_t, 3> dims;
    if (deck.hasKeyword("SPECGRID")) {
        const auto& specgridRecord = deck["SPECGRID"].back().getRecord(0);
        dims[0] = specgridRecord.getItem("NX").get< int >(0);
        dims[1] = specgridRecord.getItem("NY").get< int >(0);
        dims[2] = specgridRecord.getItem("NZ").get< int >(0);
    } else if (deck.hasKeyword("DIMENS")) {
        const auto& dimensRecord = deck["DIMENS"].back().getRecord(0);
        dims[0] = dimensRecord.getItem("NX").get< int >(0);
        dims[1] = dimensRecord.getItem("NY").get< int >(0);
        dims[2] = dimensRecord.getItem("NZ").get< int >(0);
    } else {
        OPM_THROW(std::runtime_error, "Found neither SPECGRID nor DIMENS in file. At least one is needed.");
    }

    {
        const int* actnum = deck.hasKeyword("ACTNUM") ? deck["ACTNUM"].back().getIntData().data() : nullptr;
        Opm::EclipseGrid ecl_grid(deck , actnum);
        grid.processEclipseFormat(&ecl_grid, nullptr, false);
    }

    VTKWriter<CpGrid::LeafGridView> vtkwriter(grid.leafGridView());

    const std::vector<int>& global_cell = grid.globalCell();

    std::vector<double> poros;
    condWriteDoubleField(poros, "PORO", deck, global_cell, dims, vtkwriter);

    std::vector<double> permxs;
    condWriteDoubleField(permxs, "PERMX", deck, global_cell, dims, vtkwriter);
    std::vector<double> permys;
    condWriteDoubleField(permys, "PERMY", deck, global_cell, dims, vtkwriter);

    std::vector<double> permzs;
    condWriteDoubleField(permzs, "PERMZ", deck, global_cell, dims, vtkwriter);

    std::vector<double> actnums;
    condWriteIntegerField(actnums, "ACTNUM", deck, global_cell, dims, vtkwriter);

    std::vector<double> satnums;
    condWriteIntegerField(satnums, "SATNUM", deck, global_cell, dims, vtkwriter);

    std::vector<double> regnums;
    condWriteIntegerField(regnums, "REGNUM", deck, global_cell, dims, vtkwriter);

    std::vector<double> swats;
    condWriteDoubleField(swats, "SWAT", deck, global_cell, dims, vtkwriter);

    std::string fname(eclipsefilename);
    std::string fnamebase = fname.substr(0, fname.find_last_of('.'));
    std::cout << "Writing to filename " << fnamebase << ".vtu" << std::endl;
    vtkwriter.write(fnamebase, VTK::ascii);
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}

