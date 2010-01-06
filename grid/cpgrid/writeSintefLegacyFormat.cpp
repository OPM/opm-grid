//===========================================================================
//
// File: writeSintefLegacyFormat.cpp
//
// Created: Thu Dec 17 11:12:13 2009
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
  Copyright 2009 SINTEF ICT, Applied Mathematics.
  Copyright 2009 Statoil ASA.

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

#include <fstream>
#include <vector>

#include <dune/common/ErrorMacros.hpp>
#include "../CpGrid.hpp"

namespace Dune
{

    // Forward declarations.
    namespace
    {
	void writeTopo(std::ostream& topo,
                       const cpgrid::OrientedEntityTable<0, 1>& c2f,
                       const cpgrid::OrientedEntityTable<1, 0>& f2c,
                       const SparseTable<int>& f2p,
                       const std::vector<array<int,8> >& c2p,
                       const int num_points);
	void writeGeom(std::ostream& geom,
                       const cpgrid::DefaultGeometryPolicy<CpGrid>& gpol,
                       const cpgrid::SignedEntityVariable<FieldVector<double, 3> , 1>& normals);
        void writeMap(std::ostream& map,
                      const CpGrid& g);
        void writeVtkVolumes(std::ostream& vtk,
                             const std::vector<Dune::FieldVector<double, 3> > points,
                             const std::vector<Dune::array<int, 8> >& cell_to_point);
    } // anon namespace



    /// Read the Sintef legacy grid format ('topogeom').
    void CpGrid::writeSintefLegacyFormat(const std::string& grid_prefix) const
    {
	std::string topofilename = grid_prefix + "-topo.dat";
	{
	    std::ofstream file(topofilename.c_str());
	    if (!file) {
		THROW("Could not open file " << topofilename);
	    }
	    writeTopo(file, cell_to_face_, face_to_cell_, face_to_point_, cell_to_point_, allcorners_.size());
	}
	std::string geomfilename = grid_prefix + "-geom.dat";
	{
	    std::ofstream file(geomfilename.c_str());
	    if (!file) {
		THROW("Could not open file " << geomfilename);
	    }
	    writeGeom(file, geometry_, face_normals_);
	}
        std::string mapfilename = grid_prefix + "-map.dat";
        {
            std::ofstream file(mapfilename.c_str());
	    if (!file) {
		THROW("Could not open file " << mapfilename);
	    }
            writeMap(file, *this);
        }
        std::string vtkfilename = grid_prefix + "-volumes.vtk";
        {
            std::ofstream file(vtkfilename.c_str());
	    if (!file) {
		THROW("Could not open file " << vtkfilename);
	    }
            writeVtkVolumes(file, allcorners_, cell_to_point_);
        }
    }




    namespace
    {

	void writeTopo(std::ostream& topo,
                       const cpgrid::OrientedEntityTable<0, 1>& c2f,
                       const cpgrid::OrientedEntityTable<1, 0>& f2c,
                       const SparseTable<int>& f2p,
                       const std::vector<array<int,8> >& c2p,
                       const int num_points)
	{
	    // Write header
	    std::string correct_header("topology 3 2o 2 0 3-2o 2o-2 2-0\n\n");
            topo << correct_header;

	    // Write numbers of entities.
	    int num_cells = c2f.size();
            int num_hfaces = c2f.dataSize();
            ASSERT(c2f.dataSize() == f2c.dataSize());
            int num_faces = f2c.size();
	    topo << num_cells << ' ' << num_hfaces << ' ' << num_faces << ' ' << num_points << "\n\n";

	    // Write cells to hfaces mapping
	    // cell2hface = cells to oriented faces
            int hface_count = 0;
	    for (int i = 0; i < num_cells; ++i) {
                cpgrid::EntityRep<0> cell(i, true);
                cpgrid::OrientedEntityTable<0,1>::row_type cf = c2f[cell];
		int numf = cf.size();
		topo << numf;
		for (int j = 0; j < numf; ++j) {
		    topo << ' ' << hface_count;
                    ++hface_count;
		}
                topo << '\n';
	    }
            topo << '\n';
            ASSERT(hface_count == num_hfaces);

	    // Write hfaces to faces mapping
	    for (int i = 0; i < num_cells; ++i) {
                cpgrid::EntityRep<0> cell(i, true);
                cpgrid::OrientedEntityTable<0,1>::row_type cf = c2f[cell];
		int numf = cf.size();
		for (int j = 0; j < numf; ++j) {
		    topo << cf[j].index() << ' ' << cf[j].orientation() << '\n';
		}
	    }
            topo << '\n';

	    // Write faces to points mapping
	    for (int face = 0; face < num_faces; ++face) {
                SparseTable<int>::row_type fp = f2p[face];
                int nump = fp.size();
                topo << nump;
		for (int j = 0; j < nump; ++j) {
		    topo << ' ' << fp[j];
		}
                topo << '\n';
	    }
            topo << '\n';
	} // void writeTopo()




	void writeGeom(std::ostream& geom,
                       const cpgrid::DefaultGeometryPolicy<CpGrid>& gpol,
                       const cpgrid::SignedEntityVariable<FieldVector<double, 3> , 1>& normals)
	{
            geom.precision(15);
	    std::string correct_header
		= "geometry 0:3:point 2:3:normal 2:3:centroid 2:1:area 3:3:centroid 3:1:volume\n\n";
            geom << correct_header;

	    // Write points.
	    int num_points = gpol.geomVector<3>().size();
	    geom << num_points << '\n';
	    for (int i = 0; i < num_points; ++i) {
		geom << gpol.geomVector<3>()[cpgrid::EntityRep<3>(i,true)].center() << '\n';
	    }
            geom << '\n';

	    // Write face normals
            ASSERT(gpol.geomVector<1>().size() == normals.size());
	    int num_faces = gpol.geomVector<1>().size();
	    geom << num_faces << '\n';
	    for (int i = 0; i < num_faces; ++i) {
		geom << normals[cpgrid::EntityRep<1>(i, true)] << '\n';
	    }
            geom << '\n';
	    // Write face centroids
	    geom << num_faces << '\n';
	    for (int i = 0; i < num_faces; ++i) {
		geom << gpol.geomVector<1>()[cpgrid::EntityRep<1>(i, true)].center() << '\n';
	    }
            geom << '\n';
	    // Write face areas
	    geom << num_faces << '\n';
	    for (int i = 0; i < num_faces; ++i) {
		geom << gpol.geomVector<1>()[cpgrid::EntityRep<1>(i, true)].volume() << '\n';
	    }
            geom << '\n';
	    // Write cell centroids
	    int num_cells = gpol.geomVector<0>().size();
	    geom << num_cells << '\n';
	    for (int i = 0; i < num_cells; ++i) {
		geom << gpol.geomVector<0>()[cpgrid::EntityRep<0>(i, true)].center() << '\n';
	    }
            geom << '\n';
	    // Write cell volumes
	    geom << num_cells << '\n';
	    for (int i = 0; i < num_cells; ++i) {
		geom << gpol.geomVector<0>()[cpgrid::EntityRep<0>(i, true)].volume() << '\n';
	    }
	}




        void writeMap(std::ostream& map, const CpGrid& g)
        {
            boost::array<int, 3> dims = g.logicalCartesianSize();
            map << dims[0] << ' ' << dims[1] << ' ' << dims[2] << '\n';
            int num_cells = g.size(0);
            map << num_cells << '\n';
            boost::array<int, 3> ijk;
            for (int cell = 0; cell < num_cells; ++cell) {
                g.getIJK(cell, ijk);
                map << ijk[0] << ' ' << ijk[1] << ' ' << ijk[2] << '\n';
            }
        }




        void writeVtkVolumes(std::ostream& vtk,
                             const std::vector<Dune::FieldVector<double, 3> > points,
                             const std::vector<Dune::array<int, 8> >& cell_to_point)
        {
            // Header.
            vtk <<
                "# vtk DataFile Version 2.1\n"
                "Unstructured Grid With_Cell_Data\n"
                "ASCII\n"
                "DATASET UNSTRUCTURED_GRID\n";

            // Points.
            vtk.precision(15);
            vtk << "POINTS " << points.size() << " float\n";
            std::copy(points.begin(), points.end(),
                      std::ostream_iterator<Dune::FieldVector<double, 3> >(vtk, "\n"));

            // Cell nodes.
            int nc = cell_to_point.size();
            vtk << "CELLS " << nc << ' ' << 9*nc << '\n';
            for (int i = 0; i < nc; ++i) {
                vtk << "8 "
                    << cell_to_point[i][0] << ' '
                    << cell_to_point[i][1] << ' '
                    << cell_to_point[i][3] << ' '
                    << cell_to_point[i][2] << ' '
                    << cell_to_point[i][4] << ' '
                    << cell_to_point[i][5] << ' '
                    << cell_to_point[i][7] << ' '
                    << cell_to_point[i][6] << '\n';
            }

            // Cell types.
            vtk << "CELL_TYPES " << nc << '\n';
            for (int i = 0; i < nc; ++i) {
                vtk << "12\n";
            }            
        }


    } //anon namespace



} // namespace Dune
