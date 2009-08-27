//===========================================================================
//
// File: readSintefLegacyFormat.cpp
//
// Created: Fri Jun 12 09:05:35 2009
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


#include <fstream>
#include <vector>

#include <dune/common/ErrorMacros.hpp>
#include "../CpGrid.hpp"

namespace Dune
{

    // Forward declarations.
    namespace
    {
	void readTopo(std::istream& topo,
		      cpgrid::OrientedEntityTable<0, 1>& c2f,
		      cpgrid::OrientedEntityTable<1, 0>& f2c,
		      std::vector<array<int,8> >& c2p);
	void readGeom(std::istream& geom,
		      cpgrid::DefaultGeometryPolicy& gpol,
		      cpgrid::SignedEntityVariable<FieldVector<double, 3> , 1>& normals);
        void readMap (std::istream& map,
                      std::vector<int>& global_cell);
    } // anon namespace



    /// Read the Sintef legacy grid format ('topogeom').
    void CpGrid::readSintefLegacyFormat(const std::string& grid_prefix)
    {
	std::string topofilename = grid_prefix + "-topo.dat";
	{
	    std::ifstream file(topofilename.c_str());
	    if (!file) {
		THROW("Could not open file " << topofilename);
	    }
	    readTopo(file, cell_to_face_, face_to_cell_, cell_to_point_);
	}
	std::string geomfilename = grid_prefix + "-geom.dat";
	{
	    std::ifstream file(geomfilename.c_str());
	    if (!file) {
		THROW("Could not open file " << geomfilename);
	    }
	    readGeom(file, geometry_, face_normals_);
	}
        std::string mapfilename = grid_prefix + "-map.dat";
        {
            std::ifstream file(mapfilename.c_str());
            if (file) {
                readMap(file, global_cell_);
            } else {
                // Unable to open map file.
                // Assume default (identity) cell mapping.
		int num_cells = size(0);
                global_cell_.resize(num_cells);
                for (int c = 0; c < num_cells; ++c) {
                    global_cell_[c] = c;
                }
            }
        }
    }




    namespace
    {

	void readTopo(std::istream& topo,
		      cpgrid::OrientedEntityTable<0, 1>& c2f,
		      cpgrid::OrientedEntityTable<1, 0>& f2c,
		      std::vector<array<int,8> >& c2p)
	{
	    // Check header
	    std::string topo_header;
	    std::getline(topo, topo_header);
	    std::string correct_header("topology 3 2o 2 0 3-2o 2o-2 2-0");
	    if (topo_header.compare(correct_header)) {
		THROW("Header of topology file does not match what we expect, file possibly in wrong format.\n"
		      << "Header is :" << topo_header << ":\n"
		      << "Should be :" << correct_header << ":\n");
	    }

	    // Read numbers of entities.
	    int num_cells, num_hfaces, num_faces, num_points;
	    topo >> num_cells >> num_hfaces >> num_faces >> num_points;


	    // Read cells to hfaces mapping
	    // cell2hface = cells to oriented faces
	    std::vector< std::vector<int> > cell2hface(num_cells);
	    for (int i = 0; i < num_cells; ++i) {
		int numf;
		topo >> numf;
		cell2hface[i].resize(numf);
		for (int j = 0; j < numf; ++j) {
		    topo >> cell2hface[i][j];
		}
	    }

	    // Read hfaces to faces mapping
            // hface2face = oriented faces to faces
	    std::vector< std::pair<int, bool> > hface2face;
	    hface2face.resize(num_hfaces);
	    for (int i = 0; i < num_hfaces; ++i) {
		int flag = 0;
		topo >> hface2face[i].first >> flag;
		hface2face[i].second = (flag == 1);
	    }

	    // Read faces to points mapping
	    // face2point = faces to points
	    std::vector< std::vector<int> > face2point(num_faces);
	    for (int i = 0; i < num_faces; ++i) {
		int nump;
		topo >> nump;
		face2point[i].resize(nump);
		for (int j = 0; j < nump; ++j) {
		    topo >> face2point[i][j];
		}
	    }

	    // Find faces to cells mapping
	    std::vector< FieldVector<int, 2> > face2cell;
	    face2cell.resize(num_faces, FieldVector<int, 2>(-1));
	    for (int i = 0; i < num_cells; ++i) {
		int numf = cell2hface[i].size();
		for (int j = 0; j < numf; ++j) {
		    int face = hface2face[cell2hface[i][j]].first;
		    int ind = hface2face[cell2hface[i][j]].second ? 0 : 1;
		    if (face2cell[face][ind] != -1) {
			THROW("Error in creating faces to cells mapping for face " << face);
		    }
		    face2cell[face][ind] = i;
		}
	    }

	    // Now we should be able to create c2f and f2c.
	    std::vector<int> c2fdata;
	    c2fdata.reserve(num_hfaces);
	    std::vector<int> c2fsizes(num_cells);
	    for (int i = 0; i < num_cells; ++i) {
		int numc = cell2hface[i].size();
		c2fsizes[i] = numc;
		for (int j = 0; j < numc; ++j) {
		    int hface = cell2hface[i][j];
		    int face = hface2face[hface].first;
		    int erep = hface2face[hface].second ?  face : ~face;
		    c2fdata.push_back(erep);
		}
	    }
	    ASSERT(int(c2fdata.size()) == num_hfaces);
	    c2f = cpgrid::OrientedEntityTable<0, 1>(c2fdata.begin(), c2fdata.end(), c2fsizes.begin(), c2fsizes.end());
	    c2f.makeInverseRelation(f2c);
	    // Build cell to point
	    MESSAGE("Warning: Not yet making a proper cell to point mapping for Sintef legacy grid.");
	    c2p.resize(c2f.size());
	} // void readTopo()


	template <int dim>
	struct MakeGeometry
	{
	    cpgrid::Geometry<dim, 3> operator()(const FieldVector<double, 3> pos, double vol = 1.0)
	    {
		return cpgrid::Geometry<dim, 3>(pos, vol);
	    }
	};

// 	template <>
// 	struct MakeGeometry<0>
// 	{
// 	    cpgrid::Geometry<0, 3> operator()(const FieldVector<double, 3> pos)
// 	    {
// 		return cpgrid::Geometry<dim, 3>(pos, vo);
// 	    }
// 	};

	void readGeom(std::istream& geom,
		      cpgrid::DefaultGeometryPolicy& gpol,
		      cpgrid::SignedEntityVariable<FieldVector<double, 3> , 1>& normals)
	{
	    std::string geom_header;
	    geom >> std::ws;
	    std::getline(geom, geom_header);
// 		std::ostringstream wanted_header;
// 		int codim1 = Dim - 1;
// 		int codim0 = Dim;
// 		wanted_header << "geometry " << 0 << ':' << Dim << ":point "
// 			      << codim1 << ':' << Dim << ":normal "
// 			      << codim1 << ':' << Dim << ":centroid "
// 			      << codim1 << ':' << 1 << ":area "
// 			      << codim0 << ':' << Dim << ":centroid "
// 			      << codim0 << ':' << 1 << ":volume";
	    std::string wanted_header
		= "geometry 0:3:point 2:3:normal 2:3:centroid 2:1:area 3:3:centroid 3:1:volume";
	    if (geom_header.compare(wanted_header)) {
		THROW("Header of geometry file does not match what we expect, file possibly in wrong format.");
	    }

	    typedef FieldVector<double, 3> point_t;
	    std::vector<point_t> points;
	    std::vector<point_t> face_normals;
	    std::vector<point_t> face_centroids;
	    std::vector<double>  face_areas;
	    std::vector<point_t> cell_centroids;
	    std::vector<double>  cell_volumes;

	    // Read ASCII format.
	    int num_points;
	    geom >> num_points;
	    points.resize(num_points);
	    for (int i = 0; i < num_points; ++i) {
		geom >> points[i];
	    }
	    // Read face normals
	    int num_faces;
	    geom >> num_faces;
	    face_normals.resize(num_faces);
	    for (int i = 0; i < num_faces; ++i) {
		geom >> face_normals[i];
	    }
	    // Read face centroids
	    int num_faces_again;
	    geom >> num_faces_again;
	    if (num_faces != num_faces_again) {
		THROW("Inconsistent number of faces in file.");
	    }
	    face_centroids.resize(num_faces);
	    for (int i = 0; i < num_faces; ++i) {
		geom >> face_centroids[i];
	    }
	    // Read face areas
	    geom >> num_faces_again;
	    if (num_faces != num_faces_again) {
		THROW("Inconsistent number of faces in file.");
	    }
	    face_areas.resize(num_faces);
	    for (int i = 0; i < num_faces; ++i) {
		geom >> face_areas[i];
	    }
	    // Read cell centroids
	    int num_cells;
	    geom >> num_cells;
	    cell_centroids.resize(num_cells);
	    for (int i = 0; i < num_cells; ++i) {
		geom >> cell_centroids[i];
	    }
	    // Read cell volumes
	    int num_cells_again;
	    geom >> num_cells_again;
	    if (num_cells != num_cells_again) {
		THROW("Inconsistent number of cells in file.");
	    }
	    cell_volumes.resize(num_cells);
	    for (int i = 0; i < num_cells; ++i) {
		geom >> cell_volumes[i];
	    }

	    // Code below has been copied to readEclipseFormat: \TODO Refactor!
	    // Cells
	    cpgrid::EntityVariable<cpgrid::Geometry<3, 3>, 0> cellgeom;
	    std::vector<cpgrid::Geometry<3, 3> > cg;
	    MakeGeometry<3> mcellg;
	    std::transform(cell_centroids.begin(), cell_centroids.end(),
			   cell_volumes.begin(),
			   std::back_inserter(cg), mcellg);
	    cellgeom.assign(cg.begin(), cg.end());
	    // Faces
	    cpgrid::EntityVariable<cpgrid::Geometry<2, 3>, 1> facegeom;
	    std::vector<cpgrid::Geometry<2, 3> > fg;
	    MakeGeometry<2> mfaceg;
	    std::transform(face_centroids.begin(), face_centroids.end(),
			   face_areas.begin(),
			   std::back_inserter(fg), mfaceg);
	    facegeom.assign(fg.begin(), fg.end());
	    // Points
	    cpgrid::EntityVariable<cpgrid::Geometry<0, 3>, 3> pointgeom;
	    std::vector<cpgrid::Geometry<0, 3> > pg;
	    MakeGeometry<0> mpointg;
	    std::transform(points.begin(), points.end(),
			   std::back_inserter(pg), mpointg);
	    pointgeom.assign(pg.begin(), pg.end());

	    // The final, combined object (yes, a lot of copying goes on here).
	    cpgrid::DefaultGeometryPolicy gp(cellgeom, facegeom, pointgeom);
	    gpol = gp;
	    normals.assign(face_normals.begin(), face_normals.end());
	}

        void readMap(std::istream& map, std::vector<int>& global_cell)
        {
            typedef FieldVector<int,3> IntVec;
            IntVec n;
            map >> n;  ASSERT ((n[0] > 0) && (n[1] > 0) && (n[2] > 0));

            int num_cells;
            map >> num_cells;  ASSERT (num_cells > 0);

            global_cell.resize(num_cells);

            IntVec p; int c = 0;
            while (map >> p) {
                int glob = p[0] + n[0]*(p[1] + n[1]*p[2]);
                ASSERT (glob < n[0] * n[1] * n[2]);

                global_cell[c++] = glob;
            }
            ASSERT (c == num_cells);
        }

    } //anon namespace



} // namespace Dune
