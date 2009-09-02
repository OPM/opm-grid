//===========================================================================
//
// File: readEclipseFormat.cpp
//
// Created: Fri Jun 12 09:16:59 2009
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
#include "../CpGrid.hpp"
#include <dune/grid/common/EclipseGridParser.hpp>
#include <dune/grid/common/EclipseGridInspector.hpp>
#include <dune/grid/preprocess/preprocess.h>
#include <dune/grid/common/GeometryHelpers.hpp>
#include <dune/common/StopWatch.hpp>

#define VERBOSE

namespace Dune
{


    // Forward declarations.
    namespace
    {
	void addOuterCellLayer(const grdecl& original,
			       std::vector<double>& new_coord,
			       std::vector<double>& new_zcorn,
			       std::vector<int>& new_actnum,
			       grdecl& output);
	void removeOuterCellLayer(processed_grid& grid);
	void buildTopo(const processed_grid& output,
		       std::vector<int>& global_cell,
		       cpgrid::OrientedEntityTable<0, 1>& c2f,
		       cpgrid::OrientedEntityTable<1, 0>& f2c,
		       std::vector<array<int,8> >& c2p,
		       std::vector<int>& face_to_output_face);
	void buildGeom(const processed_grid& output,
		       const cpgrid::OrientedEntityTable<0, 1>& c2f,
		       const std::vector<array<int,8> >& c2p,
		       const std::vector<int>& face_to_output_face,
		       cpgrid::DefaultGeometryPolicy& gpol,
		       cpgrid::SignedEntityVariable<FieldVector<double, 3> , 1>& normals,
		       std::vector<FieldVector<double, 3> >& allcorners);
    } // anon namespace





    /// Read the Eclipse grid format ('.grdecl').
    void CpGrid::readEclipseFormat(const std::string& filename, double z_tolerance, bool periodic_extension)
    {
	// Read eclipse file data.
#ifdef VERBOSE
	std::cout << "Parsing " << filename << std::endl;
#endif
	EclipseGridParser parser(filename);
	processEclipseFormat(parser, z_tolerance, periodic_extension);
    }





    /// Read the Eclipse grid format ('.grdecl').
    void CpGrid::processEclipseFormat(const EclipseGridParser& parser, double z_tolerance, bool periodic_extension)
    {
	EclipseGridInspector inspector(parser);

	// Make input struct for processing code.
	grdecl g;
	g.dims[0] = inspector.gridSize()[0];
	g.dims[1] = inspector.gridSize()[1];
	g.dims[2] = inspector.gridSize()[2];
	if (!parser.hasField("COORD")) {
	    THROW("Eclipse file missing required field COORD.");
	}
	g.coord = &(parser.getFloatingPointValue("COORD")[0]);
	if (!parser.hasField("ZCORN")) {
	    THROW("Eclipse file missing required field ZCORN.");
	}
	g.zcorn = &(parser.getFloatingPointValue("ZCORN")[0]);
	std::vector<int> default_actnum; // Used only if needed.
	if (parser.hasField("ACTNUM")) {
	    g.actnum = &(parser.getIntegerValue("ACTNUM")[0]);
	} else {
	    int num_cells = g.dims[0]*g.dims[1]*g.dims[2];
	    default_actnum.resize(num_cells, 1);
	    g.actnum = &default_actnum[0]; // default_actnum dies at the end of this function
	}

	if (periodic_extension) {
	    // Extend grid periodically with one layer of cells in the (i, j) directions.
	    std::vector<double> new_coord;
	    std::vector<double> new_zcorn;
	    std::vector<int> new_actnum;
	    grdecl new_g;	    
	    addOuterCellLayer(g, new_coord, new_zcorn, new_actnum, new_g);
	    // Make the grid.
	    processEclipseFormat(new_g, z_tolerance, true);
	} else {
	    // Make the grid.
	    processEclipseFormat(g, z_tolerance);
	}
    }





    /// Read the Eclipse grid format ('.grdecl').
    void CpGrid::processEclipseFormat(const grdecl& input_data, double z_tolerance, bool remove_ij_boundary)
    {
	// Process.
#ifdef VERBOSE
	std::cout << "Processing eclipse data." << std::endl;
#endif
	processed_grid output;
	process_grdecl(&input_data, z_tolerance, &output);
	if (remove_ij_boundary) {
	    removeOuterCellLayer(output);
	}

	// Move data into the grid's structures.
#ifdef VERBOSE
	std::cout << "Building topology." << std::endl;
#endif
	std::vector<int> face_to_output_face;
	buildTopo(output, global_cell_, cell_to_face_, face_to_cell_, cell_to_point_, face_to_output_face);

#ifdef VERBOSE
	std::cout << "Building geometry." << std::endl;
#endif
	buildGeom(output, cell_to_face_, cell_to_point_, face_to_output_face, geometry_, face_normals_, allcorners_);

#ifdef VERBOSE
        std::cout << "Assigning face tags." << std::endl;
#endif
	int nf = face_to_output_face.size();
	std::vector<enum face_tag> temp_tags(nf);
	for (int i = 0; i < nf; ++i) {
	    temp_tags[i] = output.face_tag[face_to_output_face[i]];
	}
	face_tag_.assign(temp_tags.begin(), temp_tags.end());

#ifdef VERBOSE
	std::cout << "Cleaning up." << std::endl;
#endif
	// Clean up the output struct.
	free_processed_grid(&output);
#ifdef VERBOSE
	std::cout << "Done with grid processing." << std::endl;
#endif
    }



    // ---- Implementation details below ----




    namespace
    {

	typedef boost::array<int, 3> coord_t;
	typedef boost::array<double, 8> cellz_t;

	cellz_t getCellZvals(const coord_t& c, const coord_t& n, const double* z)
	{
	    // cout << c << endl;
	    int delta[3] = { 1,
			     2*n[0],
			     4*n[0]*n[1] };
	    int ix = 2*(c[0]*delta[0] + c[1]*delta[1] + c[2]*delta[2]);
	    // cout << ix << endl;
	    cellz_t cellz = {{ z[ix], z[ix + delta[0]],
			       z[ix + delta[1]], z[ix + delta[1] + delta[0]],
			       z[ix + delta[2]], z[ix + delta[2] + delta[0]],
			       z[ix + delta[2] + delta[1]], z[ix + delta[2] + delta[1] + delta[0]] }};
	    return cellz;
	}



	void setCellZvals(const coord_t& c, const coord_t& n, double* z, const cellz_t& cellvals)
	{
	    int delta[3] = { 1,
			     2*n[0],
			     4*n[0]*n[1] };
	    int ix = 2*(c[0]*delta[0] + c[1]*delta[1] + c[2]*delta[2]);
	    z[ix]                                  = cellvals[0];
	    z[ix + delta[0]]                       = cellvals[1];
	    z[ix + delta[1]]                       = cellvals[2];
	    z[ix + delta[1] + delta[0]]            = cellvals[3];
	    z[ix + delta[2]]                       = cellvals[4];
	    z[ix + delta[2] + delta[0]]            = cellvals[5];
	    z[ix + delta[2] + delta[1]]            = cellvals[6];
	    z[ix + delta[2] + delta[1] + delta[0]] = cellvals[7];
	}

	coord_t indexToIjk(const coord_t& n, const int index)
	{
	    coord_t c;
	    c[2] = index/(n[0]*n[1]);
	    c[1] = (index%(n[0]*n[1]))/n[0];
	    c[0] = index%n[0];
	    return c;
	}

	void findTopAndBottomZ(const coord_t& n, const std::vector<double>& z, double& zb, double& zt)
	{
	    int numperlevel = 4*n[0]*n[1];
	    zb = *std::max_element(z.begin(), z.begin() + numperlevel);
	    zt = *std::min_element(z.end() - numperlevel, z.end());
	}


	/// Add an outer cell layer in the (i, j) directions,
	/// repeating the cells on the other side (for periodic
	/// boundary conditions).
	void addOuterCellLayer(const grdecl& original,
			       std::vector<double>& new_coord,
			       std::vector<double>& new_zcorn,
			       std::vector<int>& new_actnum,
			       grdecl& output)
	{
	    // Based on periodic_extension.cpp from the old C++ code,
	    // with a few changes:
	    //  1. We want actnum of the added cells to be true.
	    //  2. We do not treat other fields such as PORO, SATNUM etc.
	    //     since the grid will be reduced back to its regular
	    //     size before those fields are processed.

	    MESSAGE("WARNING: Assuming vertical pillars in a cartesian grid.");

	    // Build new-to-old cell index table.
	    // First expand in x.
	    coord_t n = {{ original.dims[0], original.dims[1], original.dims[2] }};
	    std::vector<int> x_new2old;
	    x_new2old.reserve((n[0]+2)*n[1]*n[2]);
	    for (int kz = 0; kz < n[2]; ++kz) {
		for (int jy = 0; jy < n[1]; ++jy) {
		    int row_ix = kz*n[0]*n[1] + jy*n[0];
		    x_new2old.push_back(row_ix + n[0] - 1);
		    for (int ix = 1; ix < n[0] + 1; ++ix) {
			x_new2old.push_back(row_ix + ix - 1);
		    }
		    x_new2old.push_back(row_ix);
		}
	    }
	    // copy(x_new2old.begin(), x_new2old.end(), ostream_iterator<int>(cout, " "));
	    // cout << endl;
	    // Then expand in y.
	    const int num_new_cells = (n[0]+2)*(n[1]+2)*n[2];
	    std::vector<int> new2old;
	    new2old.reserve(num_new_cells);
	    for (int kz = 0; kz < n[2]; ++kz) {
		for (int jy = 0; jy < n[1] + 2; ++jy) {
		    int offset = kz*(n[0] + 2)*n[1] + (jy - 1)*(n[0] + 2);
		    if (jy == 0) {
			offset = kz*(n[0] + 2)*n[1] + (n[1] - 1)*(n[0] + 2);
		    } else if (jy == n[1] + 1) {
			offset = kz*(n[0] + 2)*n[1];
		    }
		    for (int ix = 0; ix < n[0] + 2; ++ix) {
			new2old.push_back(x_new2old[offset + ix]);
		    }
		}
	    }
	    ASSERT(int(new2old.size()) == num_new_cells);
	    // copy(new2old.begin(), new2old.end(), ostream_iterator<int>(cout, " "));
	    // cout << endl;
	    // On second thought, we should have used a multidimensional array or something...

	    // Build new COORD field.
	    std::vector<double> coord;
	    coord.reserve(6*(n[0] + 3)*(n[1] + 3));
	    const double* old_coord = original.coord;
	    double dx = old_coord[6] - old_coord[0];
	    double dy = old_coord[6*(n[0] + 1) + 1] - old_coord[1];
	    double ox = old_coord[0] - dx;
	    double oy = old_coord[1] - dy;
	    for (int jy = 0; jy < n[1] + 3; ++jy) {
		double y = oy + jy*dy;
		for (int ix = 0; ix < n[0] + 3; ++ix) {
		    double x = ox + ix*dx;
		    coord.push_back(x);
		    coord.push_back(y);
		    coord.push_back(0.0);
		    coord.push_back(x);
		    coord.push_back(y);
		    coord.push_back(1.0);
		}
	    }

	    // Build new ZCORN field, PERMX, PORO, ACTNUM, SATNUM.
	    const double* old_zcorn = original.zcorn;
	    const int* old_actnum = original.actnum;
	    std::vector<double> zcorn(8*num_new_cells);
	    std::vector<int> actnum(num_new_cells);
	    coord_t new_n = {{ n[0] + 2, n[1] + 2, n[2] }};
	    for (int kz = 0; kz < new_n[2]; ++kz) {
		for (int jy = 0; jy < new_n[1]; ++jy) {
		    for (int ix = 0; ix < new_n[0]; ++ix) {
			int new_cell_index = ix + jy*(new_n[0]) + kz*(new_n[0])*(new_n[1]);
			int old_cell_index = new2old[new_cell_index];
			cellz_t cellvals = getCellZvals(indexToIjk(n, old_cell_index), n, old_zcorn);
			// cout << new_cell_index << ' ' << old_cell_index << ' ' << cellvals << endl;
			setCellZvals(indexToIjk(new_n, new_cell_index), new_n, &zcorn[0], cellvals);
			actnum[new_cell_index] = old_actnum[old_cell_index];
			if (ix == 0 || ix == new_n[0] - 1
			    || jy == 0 || jy == new_n[1] - 1) {
			    actnum[new_cell_index] = 1;  // This line is changed from the original.
			}
		    }
		}
	    }

	    // Clamp z-coord to make shoe box shape
	    bool clamp_z = true;
	    if (clamp_z) {
		double zb;
		double zt;
		findTopAndBottomZ(new_n, zcorn, zb, zt);
		for (int i = 0; i < int(zcorn.size()); ++i) {
		    zcorn[i] =  std::min(zt, std::max(zb, zcorn[i]));
		}
	    }

	    // Build output.
	    new_coord.swap(coord);
	    new_zcorn.swap(zcorn);
	    new_actnum.swap(actnum);	
	    output.dims[0] = new_n[0];
	    output.dims[1] = new_n[1];
	    output.dims[2] = new_n[2];
	    output.coord = &new_coord[0];
	    output.zcorn = &new_zcorn[0];
	    output.actnum = &new_actnum[0];
	}




	/// Helper function used by removeOuterCellLayer().
	int newLogCartFromOld(const int idx, const int dim[3])
	{
	    // Compute old (i, j, k).
	    const int Nx = dim[0];
	    const int Ny = dim[1];
	    const int NxNy = Nx*Ny;
	    int k = idx/NxNy;
	    // if (k <= 0 || k >= dim[2] - 1) return -1;
	    int j = (idx - NxNy*k)/Nx;
	    if (j <= 0 || j >= Ny - 1) return -1;
	    int i = idx - Nx*j - Nx*Ny*k;
	    if (i <= 0 || i >= Nx - 1) return -1;
	    // return (Nx - 2)*(Ny - 2)*(k - 1) + (Nx - 2)*(j - 1) + (i - 1);
	    return (Nx - 2)*(Ny - 2)*k + (Nx - 2)*(j - 1) + (i - 1);
	}

	/// Removes all (i, j) boundary cells from a grid.
	void removeOuterCellLayer(processed_grid& grid)
	{
	    // Remove outer cells as follows:
	    //   1. Build a new local_cell_index (in a new variable), compute new number_of_cells.
	    //   2. Build the inverse lookup: From old logical cartesian to new cell indices.
	    //   3. Modify face_neighbours by replacing each entry by its new cell index (or -1).
	    //   4. Modify dimensions[], number_of_cells and replace local_cell_index.
	    // After this, we still have the same number of faces, it's just that some of them may
	    // have only (-1, -1) as neighbours.

	    // Part 1 and 2 in one pass.
	    std::vector<int> new_index_to_new_lcart;
	    new_index_to_new_lcart.reserve(grid.number_of_cells); // A little too large, but no problem.
	    int num_old_lcart = grid.dimensions[0]*grid.dimensions[1]*grid.dimensions[2];
	    std::vector<int> old_lcart_to_new_index(num_old_lcart, -1);
	    for (int i = 0; i < grid.number_of_cells; ++i) {
		int old_lcart = grid.local_cell_index[i];
		int new_lcart = newLogCartFromOld(old_lcart, grid.dimensions);
		if (new_lcart != -1) {
		    old_lcart_to_new_index[old_lcart] = new_index_to_new_lcart.size();
		    new_index_to_new_lcart.push_back(new_lcart);
		} else {
		    old_lcart_to_new_index[old_lcart] = -1;
		}
	    }

	    // Part 3, modfying the face->cell connections.
	    for (int i = 0; i < 2*grid.number_of_faces; ++i) {
		int old_index = grid.face_neighbors[i];
		if (old_index != -1) {
		    int old_lcart = grid.local_cell_index[old_index];
		    int new_index = old_lcart_to_new_index[old_lcart];
		    grid.face_neighbors[i] = new_index; // May be -1, if cell is to be removed.
		}
	    }

	    // Part 4, modifying the other output data.
	    grid.dimensions[0] = grid.dimensions[0] - 2;
	    grid.dimensions[1] = grid.dimensions[1] - 2;
	    grid.dimensions[2] = grid.dimensions[2] - 2;
	    grid.number_of_cells = new_index_to_new_lcart.size();
	    std::copy(new_index_to_new_lcart.begin(), new_index_to_new_lcart.end(), grid.local_cell_index);
	}






	void buildTopo(const processed_grid& output,
		       std::vector<int>& global_cell,
		       cpgrid::OrientedEntityTable<0, 1>& c2f,
		       cpgrid::OrientedEntityTable<1, 0>& f2c,
		       std::vector<array<int,8> >& c2p,
		       std::vector<int>& face_to_output_face)
	{
	    // Map local to global cell index.
	    global_cell.assign(output.local_cell_index,
			       output.local_cell_index + output.number_of_cells);

	    // Build face to cell.
	    f2c.clear();
	    int nf = output.number_of_faces;
	    cpgrid::EntityRep<0> cells[2];
	    face_to_output_face.clear();
	    for (int i = 0; i < nf; ++i) {
		const int* fnc = output.face_neighbors + 2*i;
		int cellcount = 0;
		if (fnc[0] != -1) {
		    cells[cellcount].setValue(fnc[0], true);
		    ++cellcount;
		}
		if (fnc[1] != -1) {
		    cells[cellcount].setValue(fnc[1], false);
		    ++cellcount;
		}
		// Assertation below is no longer true, due to periodic_extension etc.
		// Instead, the appendRow() is put inside an if test.
		// ASSERT(cellcount == 1 || cellcount == 2);
		if (cellcount > 0) {
		    f2c.appendRow(cells, cells + cellcount);
		    face_to_output_face.push_back(i);
		}
	    }

	    // Build cell to face.
	    f2c.makeInverseRelation(c2f);

	    // Build cell to point
// 	    const cpgrid::EntityRep<3>* dummy = 0;
// 	    for (int i = 0; i < c2f.size(); ++i) {
// 		c2p.appendRow(dummy, dummy);
// 	    }
	    c2p.clear();
	    c2p.reserve(c2f.size());
	    for (int i = 0; i < c2f.size(); ++i) {
		cpgrid::OrientedEntityTable<0, 1>::row_type cf = c2f[cpgrid::EntityRep<0>(i)];
		// We know that the bottom and top faces come last.
		int numf = cf.size();
		int bot_face = face_to_output_face[cf[numf - 2].index()];
		int bfbegin = output.face_ptr[bot_face];
		ASSERT(output.face_ptr[bot_face + 1] - bfbegin == 4);
		int top_face = face_to_output_face[cf[numf - 1].index()];
		int tfbegin = output.face_ptr[top_face];
		ASSERT(output.face_ptr[top_face + 1] - tfbegin == 4);
		// We want the corners in 'x fastest, then y, then z' order,
		// so we need to take the face_nodes in noncyclic order: 0 1 3 2.
		array<int,8> corners = {{ output.face_nodes[bfbegin],
					  output.face_nodes[bfbegin + 1],
					  output.face_nodes[bfbegin + 3],
					  output.face_nodes[bfbegin + 2],
					  output.face_nodes[tfbegin],
					  output.face_nodes[tfbegin + 1],
					  output.face_nodes[tfbegin + 3],
					  output.face_nodes[tfbegin + 2] }};
		c2p.push_back(corners);
	    }
#ifndef NDEBUG
#ifdef VERBOSE
	    std::cout << "Doing extra topology integrity check." << std::endl;
#endif
	    cpgrid::OrientedEntityTable<1, 0> f2c_again;
	    c2f.makeInverseRelation(f2c_again);
	    ASSERT(f2c == f2c_again);
#endif
	}






	template <typename T>
	class IndirectArray
	{
	public:
	    IndirectArray(const std::vector<T>& data, const int* beg, const int* end)
		: data_(data), beg_(beg), end_(end)
	    {
	    }
	    const T& operator[](int index) const
	    {
		ASSERT(index >= 0 && index < size());
		return data_[beg_[index]];
	    }
	    int size() const
	    {
		return end_ - beg_;
	    }
	    typedef T value_type;
	private:
	    const std::vector<T>& data_;
	    const int* beg_;
	    const int* end_;
	};

	template <int dim>
	struct MakeGeometry
	{
	    cpgrid::Geometry<dim, 3> operator()(const FieldVector<double, 3>& pos, double vol = 1.0)
	    {
		return cpgrid::Geometry<dim, 3>(pos, vol);
	    }
	};

	template <>
	struct MakeGeometry<3>
	{
	    const FieldVector<double, 3>* allcorners_;
	    MakeGeometry(const FieldVector<double, 3>* allcorners)
		: allcorners_(allcorners)
	    {
	    }
	    cpgrid::Geometry<3, 3> operator()(const FieldVector<double, 3>& pos,
					      double vol,
					      const array<int,8>& corner_indices)
	    {
		return cpgrid::Geometry<3, 3>(pos, vol, allcorners_, &corner_indices[0]);
	    }
	};





	void buildGeom(const processed_grid& output,
		       const cpgrid::OrientedEntityTable<0, 1>& c2f,
		       const std::vector<array<int,8> >& c2p,
		       const std::vector<int>& face_to_output_face,
		       cpgrid::DefaultGeometryPolicy& gpol,
		       cpgrid::SignedEntityVariable<FieldVector<double, 3>, 1>& normals,
		       std::vector<FieldVector<double, 3> >& allcorners)
	{
	    typedef FieldVector<double, 3> point_t;
	    std::vector<point_t>& points = allcorners;
	    std::vector<point_t> face_normals;
	    std::vector<point_t> face_centroids;
	    std::vector<double>  face_areas;
	    std::vector<point_t> cell_centroids;
	    std::vector<double>  cell_volumes;
	    using namespace GeometryHelpers;
#ifdef VERBOSE
	    time::StopWatch clock;
	    clock.start();
#endif
	    // Get the points.
	    int np = output.number_of_nodes;
	    points.clear();
	    points.reserve(np);
	    for (int i = 0; i < np; ++i) {
		// \TODO add a convenience explicit constructor
		// for FieldVector taking an iterator.
		point_t pt;
		for (int dd = 0; dd < 3; ++dd) {
		    pt[dd] = output.node_coordinates[3*i + dd];
		}
		points.push_back(pt);
	    }
#ifdef VERBOSE
	    std::cout << "Points:             " << clock.secsSinceLast() << std::endl;
#endif

	    // Get the face data.
	    // \TODO Both the face and (especially) the cell section
	    // is not very efficient. It could be rewritten easily
	    // (focus on the polygonCellXXX methods first).
	    // \TODO Use exact geometry instead of these approximations.
	    int nf = face_to_output_face.size();
	    const int* fn = output.face_nodes;
	    const int* fp = output.face_ptr;
	    for (int face = 0; face < nf; ++face) {
		// Computations in this loop could be speeded up
		// by doing more of them simultaneously.
		int output_face = face_to_output_face[face];
		IndirectArray<point_t> face_pts(points, fn + fp[output_face], fn + fp[output_face+1]);
		point_t avg = average(face_pts);
		point_t centroid = polygonCentroid(face_pts, avg);
		point_t normal = polygonNormal(face_pts, centroid);
		double area = polygonArea(face_pts, centroid);
		face_normals.push_back(normal);
		face_centroids.push_back(centroid);
		face_areas.push_back(area);
	    }
#ifdef VERBOSE
	    std::cout << "Faces:              " << clock.secsSinceLast() << std::endl;
#endif
	    // Get the cell data.
	    int nc = output.number_of_cells;
	    std::vector<int> face_indices;
	    for (int cell = 0; cell < nc; ++cell) {
		cpgrid::EntityRep<0> cell_ent(cell, true);
		cpgrid::OrientedEntityTable<0, 1>::row_type cf = c2f[cell_ent];
		face_indices.clear();
		for (int local_index = 0; local_index < cf.size(); ++local_index) {
		    face_indices.push_back(cf[local_index].index());
		}
		IndirectArray<point_t> cell_pts(face_centroids, &face_indices[0], &face_indices[0] + cf.size());
		point_t cell_avg = average(cell_pts);
		point_t cell_centroid(0.0);
		double tot_cell_vol = 0.0;
		for (int local_index = 0; local_index < cf.size(); ++local_index) {
		    int face = cf[local_index].index();
		    int output_face = face_to_output_face[face];
		    IndirectArray<point_t> face_pts(points, fn + fp[output_face], fn + fp[output_face+1]);
		    double small_vol = polygonCellVolume(face_pts, face_centroids[face], cell_avg);
		    tot_cell_vol += small_vol;
		    point_t face_contrib = polygonCellCentroid(face_pts, face_centroids[face], cell_avg);
		    face_contrib *= small_vol;
		    cell_centroid += face_contrib;
		}
		cell_centroid /= tot_cell_vol;
		cell_centroids.push_back(cell_centroid);
		cell_volumes.push_back(tot_cell_vol);
	    }
#ifdef VERBOSE
	    std::cout << "Cells:              " << clock.secsSinceLast() << std::endl;
#endif


	    // \TODO Fix the code below , as it is:
	    // A) wasteful
	    // B) wordy,
	    // C) slow
	    // D) copied from readSintefLegacyFormat.cpp
	    // Cells
	    cpgrid::EntityVariable<cpgrid::Geometry<3, 3>, 0> cellgeom;
	    std::vector<cpgrid::Geometry<3, 3> > cg;
	    cg.reserve(nc);
	    MakeGeometry<3> mcellg(&allcorners[0]);
// 	    std::transform(cell_centroids.begin(), cell_centroids.end(),
// 			   cell_volumes.begin(),
// 			   std::back_inserter(cg), mcellg);
	    for (int c = 0;  c < nc; ++c) {
		cg.push_back(mcellg(cell_centroids[c], cell_volumes[c], c2p[c]));
	    }
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
#ifdef VERBOSE
	    std::cout << "Transforms/copies:  " << clock.secsSinceLast() << std::endl;
#endif

	    // The final, combined object (yes, a lot of copying goes on here).
	    cpgrid::DefaultGeometryPolicy gp(cellgeom, facegeom, pointgeom);
	    gpol = gp;
	    normals.assign(face_normals.begin(), face_normals.end());
#ifdef VERBOSE
	    std::cout << "Final construction: " << clock.secsSinceLast() << std::endl;
#endif
	}


    } // anon namespace




} // namespace Dune

