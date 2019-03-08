// Copyright (C) 2013 Uni Research AS
//  This file is licensed under the GNU General Public License v3.0

#include <opm/grid/verteq/nav.hpp>
#include <opm/grid/verteq/topsurf.hpp>
#include <opm/grid/verteq/utility/exc.hpp>
#include <opm/grid/verteq/utility/runlen.hpp> // rlw_int
#include <opm/grid/cornerpoint_grid.h> // compute_geometry
#include <boost/io/ios_state.hpp> // ios_all_saver
#include <algorithm> // min, max
#include <climits> // INT_MIN, INT_MAX
#include <cstdlib> // div
#include <iosfwd> // ostream
#include <map>
#include <memory> // unique_ptr
#include <numeric> // accumulate, iota
#include <vector>
#include <utility> // pair

using namespace boost;
using namespace Opm;
using namespace std;

/// Helper routine to print of map
template <typename T, typename U>
void dump_map (ostream& os, const map<T, U>& m) {
	for (typename map<T, U>::const_iterator it = m.begin(); it != m.end(); ++it) {
		boost::io::ios_all_saver state (os);
		os << it->first << ":\t";
		state.restore ();
		os << it->second << endl;
	}
}

/**
 * @brief Process to extract the top surface from a structured grid.
 *
 * This object encapsulates a procedure with variables shared amongst
 * several sub-procedures (like in Pascal). These objects are not
 * supposed to linger on afterwards.
 */
struct TopSurfBuilder {
	// source grid from which we get the input data
	const UnstructuredGrid& fine_grid;

	// target grid we are constructing
	TopSurf& ts;

	// number of grid dimensions in each direction of the plane
	Cart3D three_d;

	// dimensions needed to create a two-dimensional projection
	// of the top surface
	Cart2D two_d;

	// map from a two-dimensional Cartesian coordinate to the final
	// id of an active element in the grid, or NO_ELEM if nothing is assigned
	// this vector is first valid after create_elements() have been done
	vector <int> elms;

	// map from a two-dimensional Cartesian node coordinate to the final
	// id of active nodes in the grid, or NO_NODE if nothing is assigned.
	// this vector is first valid after create_nodes() have been done
	vector <int> nodes;

	// map from a two-dimensional Cartesion face coordinate to the final
	// id of active faces in the grid, or NO_FACE if nothing is assigned.
	// this vector is first valid after create_faces() have been done.
	vector <int> faces;

	// logical Cartesian indices for items in the fine grid. we need our
	// own copy of this since not all grids provide it
	vector <int> fine_global;

	TopSurfBuilder (const UnstructuredGrid& from, TopSurf& into)
		// link to the fine grid for the duration of the construction
		: fine_grid (from)

		// allocate memory for the grid. it is initially empty
		, ts (into)

		// extract dimensions from the source grid
		, three_d (fine_grid)
		, two_d (three_d.project ())
		, fine_global (from.number_of_cells, 0)	{

		// check that the fine grid contains structured information;
		// this is essential to mapping cells to columns
		const int prod = std::accumulate(&fine_grid.cartdims[0],
		                                  &fine_grid.cartdims[fine_grid.dimensions],
		                                  1,
		                                  std::multiplies<int>());
		if (!prod) {
			throw OPM_EXC ("Find grid is not (logically) structured");
		}

		// some cartesian grids (most notably those generated with
		// create_grid_cart{2,3}d) have no global cell
		if (!fine_grid.global_cell) {
			std::iota (fine_global.begin(), fine_global.end(), 0);
		}
		else {
			std::copy (fine_grid.global_cell,
			           fine_grid.global_cell + fine_grid.number_of_cells,
			           fine_global.begin ());
		}

		// create frame of the new top surface
		create_dimensions ();

		// identify active columns in the grid
		create_elements ();

		// identify active points in the grid
		create_nodes ();

		// identify active faces in the grid
		create_faces ();

		// cache fine block and column metrics
		create_heights ();
	}

	// various stages of the build process, supposed to be called in
	// this order. (I have separated them into separate procedures to
	// make it more obvious what parts that needs to be shared between
	// them)
private:
	void create_dimensions () {
		// we are going to create two-dimensional grid
		ts.dimensions = 2;
		ts.cartdims[0] = two_d.ni;
		ts.cartdims[1] = two_d.nj;
		ts.cartdims[2] = 1;
	}

	void create_elements() {
		// statistics of the deepest and highest active k-index of
		// each column in the grid. to know each index into the column,
		// we only need to know the deepest k and the count; the highest
		// is sampled to do consistency checks afterwards
		const int num_cols = two_d.num_elems ();

		// assume initially that there are no active elements in each column
		vector <int> act_cnt (num_cols, 0);
		ts.number_of_cells = 0;

		// initialize these to values that are surely out of range, so that
		// the first invocation of min or max always set the value. we use
		// this to detect whether anything was written later on. since the
		// numbering of the grid starts at the top, then the deepest cell
		// has the *largest* k-index, thus we need a value smaller than all
		vector <int> deep_k (num_cols, INT_MIN);
		vector <int> high_k (num_cols, INT_MAX);

		// loop once through the fine grid to gather statistics of the
		// size of the surface so we know what to allocate
		for (int fine_elem = 0; fine_elem != fine_grid.number_of_cells; ++fine_elem) {
			// get the cartesian index for this cell; this is the cell
			// number in a grid that also includes the inactive cells
			const Cart3D::elem_t cart_ndx = fine_global [fine_elem];

			// deconstruct the cartesian index into (i,j,k) constituents;
			// the i-index moves fastest, as this is Fortran-indexing
			const Coord3D ijk = three_d.coord (cart_ndx);

			// figure out which column this item belongs to (in 2D)
			const Cart2D::elem_t col = two_d.cart_ndx (ijk);

			// update the statistics for this column; 'deepest' is the largest
			// k-index seen so far, 'highest' is the smallest (ehm)
			deep_k[col] = max (deep_k[col], ijk.k());
			high_k[col] = min (high_k[col], ijk.k());

			// we have seen an element in this column; it becomes active. only
			// columns with active cells will get active elements in the surface
			// grid. if this cell wasn't marked as active before we can check it
			// off now
			if (!act_cnt[col]++) {
				ts.number_of_cells++;
			}
		}

		// check that we have a continuous range of elements in each column;
		// this must be the case to assume that the entire column can be merged
		for (int col = 0; col < num_cols; ++col) {
			if (act_cnt[col]) {
				if (high_k[col] + act_cnt[col] - 1 != deep_k[col]) {
					const Coord2D coord = two_d.coord (col);
					throw OPM_EXC ("Non-continuous column at (%d, %d)", coord.i(), coord.j());
				}
			}
		}

		// allocate memory needed to hold the elements in the grid structure
		// if we throw an exception at some point, the destructor of the TopSurf
		// will take care of deallocating this memory for us
		ts.global_cell = new int [ts.number_of_cells];
		ts.col_cellpos = new int [ts.number_of_cells+1];

		// we haven't filled any columns yet, so this is a sensible init value
		ts.max_vert_res = 0;

		// there is no elements before the first column, so this number is
		// always zero. it is written to the array to maintain a straight code
		// path in calculations.
		ts.col_cellpos[0] = 0;

		// now we know enough to start assigning ids to active columns.if
		// memory is a real shortage, we could reuse the act_cnt array for this.
		elms.resize (num_cols, Cart2D::NO_ELEM);

		// loop through the grid and assign an id for all columns that have
		// active elements
		int elem_id = 0;
		for (int col = 0; col < num_cols; ++col) {
			if (act_cnt[col]) {
				elms[col] = elem_id;

				// dual pointer that maps the other way; what is the structured
				// index of this element (flattened into an integer)
				ts.global_cell[elem_id] = col;

				// note the number of elements there now are before the next column;
				// in addition to all the previous ones, our elements are now added
				ts.col_cellpos[elem_id+1] = ts.col_cellpos[elem_id] + act_cnt[col];

				// update the largest number of these seen so far
				ts.max_vert_res = max (ts.max_vert_res, act_cnt[col]);

				// only increment this if we found an active column, elem_id <= col
				elem_id++;
			}
		}

		// now write indices from the fine grid into the column map of the surface
		// we end up with a list of element that are in each column
		ts.col_cells = new int [fine_grid.number_of_cells];
		ts.fine_col = new int [fine_grid.number_of_cells];
		for (int cell = 0; cell < fine_grid.number_of_cells; ++cell) {
			// get the Cartesian index for this element
			const Cart3D::elem_t cart_ndx = fine_global[cell];
			const Coord3D ijk = three_d.coord (cart_ndx);

			// get the id of the column in which this element now belongs
			const Cart2D::elem_t col = two_d.cart_ndx (ijk);
			const int elem_id = elms[col];

			// start of the list of elements for this particular column
			const int segment = ts.col_cellpos[elem_id];

			// since there is supposed to be a continuous range of elements in
			// each column, we can calculate the relative position in the list
			// based on the k part of the coordinate.
			const int offset = ijk.k() - high_k[col];

			// write the fine grid cell number in the column list; since we
			// have calculated the position based on depth, the list will be
			// sorted downwards up when we are done
			ts.col_cells[segment + offset] = cell;

			// reverse mapping; allows us to quickly figure out the corresponding
			// column of a location (for instance for a well)
			ts.fine_col[cell] = elem_id;
		}

		// these members will be filled by computeGeometry, but needs valid
		// memory to work with
		ts.cell_volumes = new double [ts.number_of_cells];
		ts.cell_centroids = new double [ts.dimensions * ts.number_of_cells];
	}

	void create_nodes () {
		// construct a dual Cartesian grid consisting of the points
		const int num_nodes = two_d.num_nodes ();

		// vectors which will hold the coordinates for each active point.
		// at first we sum all the points, then we divide by the count to
		// get the average. as long as the count is zero, there is no
		// registered active point at this location
		vector <double> x (num_nodes, 0.);
		vector <double> y (num_nodes, 0.);
		vector <int> cnt (num_nodes, 0);

		// number of nodes needed in the top surface
		int active_nodes = 0;

		// tag of the top side in a cell; we're looking for this
		const int top_tag = Side3D (Dim3D::Z, Dir::DEC).facetag ();

		// initial corner value. this could really be anything, since
		// we expect all the fields to be overwritten.
		const Corn3D blank (Dir::DEC, Dir::DEC, Dir::DEC);

		// this map holds the classification of each node locally for the
		// element being currently processed. we reuse the map to avoid
		// needless memory allocations.
		typedef map <int, Corn3D> cls_t;
		cls_t classifier;

		// loop through all active cells in the top surface
		for (int col = 0; col < ts.number_of_cells; ++col) {
			// get the highest element in this column; since we have them
			// sorted by k-index this should be the first item in the
			// extended column info
			const Cart3D::elem_t top_cell_glob_id = ts.col_cells [ts.col_cellpos[col]];

			// start afresh whenever we start working on a new element
			classifier.clear ();
			int top_face_glob_id = Cart2D::NO_FACE;

			// loop through all the faces of the top element
			for (int face_pos = fine_grid.cell_facepos[top_cell_glob_id];
					 face_pos != fine_grid.cell_facepos[top_cell_glob_id+1];
					 ++face_pos) {

				// get the (normal) dimension and direction of this face
				const int this_tag = fine_grid.cell_facetag[face_pos];
				Side3D s = Side3D::from_tag (this_tag);

				// identifier of the face, which is the index in the next arary
				const int face_glob_id = fine_grid.cell_faces[face_pos];

				// remember it if we've found the top face
				if (this_tag == top_tag) {
					if (top_face_glob_id != Cart2D::NO_FACE) {
						throw OPM_EXC ("More than one top face in element %d", top_cell_glob_id);
					}
					top_face_glob_id = face_glob_id;
				}

				// loop through all nodes in this face, adding them to the
				// classifier. when we are through with all the faces, we have
				// found in which corner a node is, defined by a direction in
				// each of the three dimensions
				for (int node_pos = fine_grid.face_nodepos[face_glob_id];
						 node_pos != fine_grid.face_nodepos[face_glob_id+1];
						 ++node_pos) {
					const int node_glob_id = fine_grid.face_nodes[node_pos];

					// locate pointer to data record ("iterator" in stl parlance)
					// for this node, if it is already there. otherwise, just start
					// out with some blank data (which eventually will get overwritten)
					cls_t::iterator ptr = classifier.find (node_glob_id);
					Corn3D prev (ptr == classifier.end () ? blank : ptr->second);

					// update the dimension in which this face is pointing
					if (ptr != classifier.end ()) {
						classifier.erase (ptr);
					}
					const Corn3D upd_corn = prev.pivot (s.dim(), s.dir());
					classifier.insert (make_pair (node_glob_id, upd_corn));
				}

				// after this loop, we have a map of each node local to the element,
				// classified into in which corner it is located (it cannot be in
				// both directions in the same dimension -- then it would have to
				// belong to two opposite faces, unless the grid is degenerated)
			}
			/*
			cerr << "elem: " << three_d.coord(top_cell_glob_id) << ':' << endl;
			dump_map (cerr, classifier);
			*/

			// cannot handle degenerate grids without top face properly
			if (top_face_glob_id == Cart2D::NO_FACE) {
				throw OPM_EXC ("No top face in cell %d", top_cell_glob_id);
			}

			// get the Cartesian ij coordinate of this cell
			const Cart2D::elem_t top_cell_cart_ndx = ts.global_cell [top_cell_glob_id];
			const Coord2D ij = two_d.coord (top_cell_cart_ndx);

			// loop through all the nodes of the top face, and write their position
			// into the corresponding two-d node. this has to be done separately
			// after we have classified *all* the nodes of the element, in order for
			// the corner values to be set correctly, i.e. we cannot merge this into
			// the loop above.
			for (int node_pos = fine_grid.face_nodepos[top_face_glob_id];
					 node_pos != fine_grid.face_nodepos[top_face_glob_id+1];
					 ++node_pos) {
				const int node_glob_id = fine_grid.face_nodes[node_pos];

				// get which corner this node has; this returns a three-dimensional
				// corner, but by using the base class part of it we automatically
				// project it to a flat surface
				cls_t::iterator ptr = classifier.find (node_glob_id);
				const Corn3D corn (ptr->second);

				// get the structured index for this particular corner
				const Cart2D::node_t cart_node = two_d.node_ndx(ij, corn);

				// add these coordinates to the average position for this junction
				// if we activate a corner, then add it to the total count
				x[cart_node] += fine_grid.node_coordinates[Dim3D::COUNT*node_glob_id+0];
				y[cart_node] += fine_grid.node_coordinates[Dim3D::COUNT*node_glob_id+1];
				if (!cnt[cart_node]++) {
					++active_nodes;
				}
			}
		}

		// after this loop we the accumulated coordinates for each of the
		// corners that are part of active elements (the nodes that are
		// needed in the top surface)

		// assign identifiers and find average coordinate for each point
		ts.number_of_nodes = active_nodes;
		ts.node_coordinates = new double [active_nodes * Dim2D::COUNT];
		nodes.resize (num_nodes, Cart2D::NO_NODE);
		int next_node_id = 0;
		for (int cart_node = 0; cart_node != num_nodes; ++cart_node) {
			if (cnt[cart_node]) {
				nodes[cart_node] = next_node_id;
				const int start = Dim2D::COUNT * next_node_id;
				ts.node_coordinates[start+0] = x[cart_node] / cnt[cart_node];
				ts.node_coordinates[start+1] = y[cart_node] / cnt[cart_node];
				++next_node_id;
			}
		}

		// dump node topology to console
		/*
		for (int cart_node_ndx = 0; cart_node_ndx != num_nodes; ++cart_node_ndx) {
			const int glob_node_id = nodes[cart_node_ndx];
			if (glob_node_id != Cart2D::NO_NODE) {
				cerr << "node " << nodes[cart_node_ndx] << ": ("
						 << ts.node_coordinates[2*glob_node_id+0] << ','
						 << ts.node_coordinates[2*glob_node_id+1] << ')' << endl;
			}
		}
		*/

		// TODO: check for degeneracy by comparing each node's coordinates
		// with those in the opposite direction in both dimensions (separately)
	}

	void create_faces () {
		// number of possible (but not necessarily active) faces
		const int num_faces = two_d.num_faces ();

		// assign identifiers into this array. start out with the value
		// NO_FACE which means that unless we write in an id, the face
		// is not active
		faces.resize (num_faces, Cart2D::NO_FACE);

		// a face will be referenced from two elements (apart from boundary),
		// denoted the "primary" and "secondary" neighbours. the nodes in a
		// face needs to be specified so that the normal to the directed LINE_NODES
		// points towards the element center of the *primary* neighbour.

		// we use the convention that every face that are in the J-direction
		// are directed from decreasing direction to increasing direction,
		// whereas every face that are in the I-direction are directed the
		// opposite way, from increasing to decreasing. this way, an element
		// is primary neighbour for a face if the face is on side which is
		// classified as decreasing, relative to the center, and secondary
		// if the face is in the increasing direction of whatever axis.

		//                 I+             (I+,J+)     (I+,J-)
		//             o <---- o               o <---- o
		//             |       |               |       |
		//          J+ |       | J-            |       |
		//             v       v               v       v
		//             o <---- o               o <---- o
		//                 I-             (I-,J+)     (I-,J-)

		// TODO: Possible to instead create an S-shaped traversal with
		// different directions every odd/even column, so the faces ends
		// up naturally in clockwise directions around the element?

		// allocate two arrays that will hold the global ids of the two
		// elements on each side, or NO_ELEM if there is no active element
		// on that side. if both sides are NO_ELEM then the face can be
		// optimized away from the resulting grid.
		vector <int> pri_elem (num_faces, Cart2D::NO_ELEM);
		vector <int> sec_elem (num_faces, Cart2D::NO_ELEM);

		vector <int> src (num_faces, Cart2D::NO_NODE);
		vector <int> dst (num_faces, Cart2D::NO_NODE);

		// keep count on how many faces that actually have at least one
		// active neighbouring element
		int active_faces = 0;

		// loop through all the active elements in the grid, and write
		// their identifier in the neighbour array for their faces
		for (int j = 0; j != two_d.nj; ++j) {
			for (int i = 0; i != two_d.ni; ++i) {
				// where are we in the grid?
				const Coord2D coord (i, j);
				const int col = two_d.cart_ndx (coord);

				// check if this element is active; is there an id assignment?
				const int elem_glob_id = elms[col];
				if (elem_glob_id != Cart2D::NO_ELEM) {

					// loop through all sides of this element; by assigning identities
					// to the faces in this manner, the faces around each element are
					// relatively local to eachother in the array.
					for (const Side2D* s = Side2D::begin(); s != Side2D::end(); ++s) {
						// cartesian index of this face, i.e. index just depending on the
						// extent of the grid, not whether face is active or not
						const int cart_face = two_d.face_ndx (coord, *s);

						// select primary or secondary neighbour collection based on
						// the direction of the face relative to the center of the
						// element, see discussion above. if the vector from the center
						// to the face points in the INC direction (i.e. this is an INC
						// face), then that vector is aligned with the right-normal of
						// the face, and this node is the primary.
						const bool is_primary = s->dir () == Dir::INC;
						vector <int>& neighbour = is_primary ? pri_elem : sec_elem;
						vector <int>& other = is_primary ? sec_elem : pri_elem;

						// put this identifier in there
						if (neighbour[cart_face] != Cart2D::NO_ELEM) {
							throw OPM_EXC ("Duplicate neighbour assignment in column (%d,%d)",
														 coord.i(), coord.j());
						}
						neighbour[cart_face] = elem_glob_id;

						// if this was the first time we assigned a neighbour to the
						// face, then count it as active and assign an identity
						if (other[cart_face] == Cart2D::NO_ELEM) {
							faces[cart_face] = active_faces++;

							// faces that are in the I-dimension (I- and I+) starts at the DEC
							// direction in the J-dimension, where as for the faces in the J-
							// dimension it is opposite
							const Dir src_dir = s->dim () == Dim2D::X ? Dir::DEC : Dir::INC;
							const Dir dst_dir = src_dir.opposite ();

							// use the direction of the side for its dimension, as both the corners
							// of the side will be here, and use the two other directions for the
							// remaining dimension. the trick is to know in which corner the face
							// should start, see above.
							const Corn2D src_corn (s->dim () == Dim2D::X ? s->dir () : src_dir,
																		 s->dim () == Dim2D::X ? src_dir : s->dir ());
							const Corn2D dst_corn (s->dim () == Dim2D::X ? s->dir () : dst_dir,
																		 s->dim () == Dim2D::X ? dst_dir : s->dir ());

							// get the identity of the two corners, and take note of these
							const int src_cart_ndx = two_d.node_ndx (coord, src_corn);
							const int dst_cart_ndx = two_d.node_ndx (coord, dst_corn);
							src[cart_face] = nodes[src_cart_ndx];
							dst[cart_face] = nodes[dst_cart_ndx];
						}
					}
				}
			}
		}
		// after this loop we have a map of all possible faces, and know
		// how many of these are active.

		// dump face topology to console
		/*
		for (int cart_face_ndx = 0; cart_face_ndx != num_faces; ++cart_face_ndx) {
			const int face_glob_id = faces[cart_face_ndx];
			if (face_glob_id != Cart2D::NO_FACE) {
				cerr << "face " << face_glob_id << ": " << endl;
				cerr << "\t" << "elements: " << pri_elem[cart_face_ndx] << ","
																		 << sec_elem[cart_face_ndx] << endl;
				cerr << "\t" << "nodes: " << src[cart_face_ndx] << ","
																	<< dst[cart_face_ndx] << endl;
			}
		}
		*/

		// each cell has 4 sides in 2D; we assume no degenerate sides
		const int QUAD_SIDES = Dim2D::COUNT * Dir::COUNT;
		const int num_sides = QUAD_SIDES * active_faces;

		// nodes in faces; each face in 2D is only 1D, so simplices always
		// have only two nodes (a LINE_NODES, with corners in each direction)
		const int LINE_NODES = Dim1D::COUNT * Dir::COUNT;
		const int num_corns = LINE_NODES * active_faces;

		// number of element neighbours for each face. this is always 2,
		// the reason for not using the number is to make it searchable
		const int NEIGHBOURS = 2;

		// allocate memory; this will be freed in the TopSurf destructor
		ts.face_nodes = new int [num_corns];
		ts.face_nodepos = new int [active_faces + 1];
		ts.face_cells = new int [NEIGHBOURS * active_faces];
		ts.cell_faces = new int [num_sides];
		ts.cell_facepos = new int [ts.number_of_cells + 1];
		ts.cell_facetag = new int [num_sides];
		ts.number_of_faces = active_faces;

		// we need to allocate memory for computeGeometry to fill
		ts.face_centroids = new double [Dim2D::COUNT * active_faces];
		ts.face_areas = new double [active_faces];
		ts.face_normals = new double [Dim2D::COUNT * active_faces];

		// write the internal data structures to UnstructuredGrid representation

		// face <-> node topology
		for (int cart_face = 0; cart_face != num_faces; ++cart_face) {
			const int face_glob_id = faces[cart_face];
			if (face_glob_id != Cart2D::NO_FACE) {
				// since each face has exactly two coordinates, we can easily
				// calculate the position based only on the face number
				const int start_pos = LINE_NODES * face_glob_id;
				ts.face_nodepos[face_glob_id] = start_pos;
				ts.face_nodes[start_pos + 0] = src[cart_face];
				ts.face_nodes[start_pos + 1] = dst[cart_face];

				// TODO: If a vertical fault displaces two column so that there
				// is no longer connection between them, they will be reconnected
				// here. This condition can be detected by comparing the top and
				// bottom surface.

				// neighbours should already be stored in the right orientation
				ts.face_cells[NEIGHBOURS * face_glob_id + 0] = pri_elem[cart_face];
				ts.face_cells[NEIGHBOURS * face_glob_id + 1] = sec_elem[cart_face];
			}
		}
		ts.face_nodepos[ts.number_of_faces] = LINE_NODES * ts.number_of_faces;

		// cell <-> face topology
		for (int j = 0; j != two_d.nj; ++j) {
			for (int i = 0; i != two_d.ni; ++i) {
				// get various indices for this element
				const Coord2D coord (i, j);
				const int cart_elem = two_d.cart_ndx (coord);
				const int elem_glob_id = elms[cart_elem];

				// each element is assumed to be a quad, so we can calculate the
				// number of accumulated sides based on the absolute id
				const int start_pos = QUAD_SIDES * elem_glob_id;
				ts.cell_facepos[elem_glob_id] = start_pos;

				// write all faces for this element
				for (const Side2D* s = Side2D::begin(); s != Side2D::end(); ++s) {
					// get the global id of this face
					const int face_cart_ndx = two_d.face_ndx (coord, *s);
					const int face_glob_id = faces[face_cart_ndx];

					// the face tag can also serve as an offset into a regular element
					const int ofs = s->facetag ();
					ts.cell_faces[start_pos + ofs] = face_glob_id;
					ts.cell_facetag[start_pos + ofs] = ofs;
				}
			}
		}
		ts.cell_facepos[ts.number_of_cells] = QUAD_SIDES * ts.number_of_cells;
	}

	/**
	 * Specific face number of a given side of an element.
	 *
	 * @param glob_elem_id Element index in the fine grid.
	 * @param s Side to locate
	 * @return Index of the face of the element which is this side
	 *
	 * @see Opm::UP, Opm::DOWN
	 */
	int find_face (int glob_elem_id, const Side3D& s) {
		// this is the tag we are looking for
		const int target_tag = s.facetag ();

		// this is the matrix we are looking in
		const rlw_int cell_facetag = grid_cell_facetag (fine_grid);

		// we are returning values from this matrix
		const rlw_int cell_faces = grid_cell_faces (fine_grid);

		// loop through all faces for this element; face_ndx is the local
		// index amongst faces for just this one element.
		for (int local_face = 0;
		     local_face < cell_facetag.size (glob_elem_id);
		     ++local_face) {

			// if we found a match, then return this; don't look any more
			if (cell_facetag[glob_elem_id][local_face] == target_tag) {

				// return the (global) index of the face, not the tag!
				return cell_faces[glob_elem_id][local_face];
			}
		}

		// in a structured grid we expect to find every face
		throw OPM_EXC ("Element %d does not have face #%d", glob_elem_id, target_tag);
	}

	/**
	 * Get absolute elevation (z-coordinate) of a face. This uses the
	 * elevation at the centroid as representative of the entire face.
	 *
	 * @param glob_elem_id Element index in the fine grid.
	 * @param s Side to locate.
	 * @return Elevation for the midpoint of this face.
	 */
	double find_zcoord (int glob_elem_id, const Side3D& s) {
		// find the desired face for this element
		const int face_ndx = find_face (glob_elem_id, s);

		// get the z-coordinate for it
		const int z_ndx = face_ndx * Dim3D::COUNT + Dim3D::Z.val;
		const double z = fine_grid.face_centroids[z_ndx];
		return z;
	}

	/**
	 * Height of a particular element.
	 *
	 * @param glob_elem_id Element index in the fine grid.
	 * @return Difference between center of top and bottom face.
	 */
	double find_height (int glob_elem_id) {
		// get the z-coordinate for each the top and bottom face for this element
		const double up_z = find_zcoord (glob_elem_id, UP);
		const double down_z = find_zcoord (glob_elem_id, DOWN);

		// the side that is down should have the z coordinate with highest magnitude
		const double height = down_z - up_z;
		return height;
	}

	void create_heights () {
		// allocate memory to hold the heights
		ts.dz = new double [fine_grid.number_of_cells];
		ts.h = new double [fine_grid.number_of_cells];
		ts.z0 = new double [ts.number_of_cells];
		ts.h_tot = new double [ts.number_of_cells];

		// view that lets us treat it as a matrix
		const rlw_int blk_id (ts.number_of_cells, ts.col_cellpos, ts.col_cells);
		const rlw_double dz (ts.number_of_cells, ts.col_cellpos, ts.dz);
		const rlw_double h (ts.number_of_cells, ts.col_cellpos, ts.h);

		// find all measures per column
		for (int col = 0; col < blk_id.cols (); ++col) {
			// reference height for this column (if there is any elements)
			if (blk_id.size (col)) {
				const int top_ndx = blk_id[col][0];
				ts.z0[col] = find_zcoord (top_ndx, UP);
			}

			// reset height for each column
			double accum = 0.;

			// height of each element in the column element
			double* const dz_col = dz[col];
			double* const h_col = h[col];
			for (int col_elem = 0; col_elem < blk_id.size (col); ++col_elem) {
				h_col[col_elem] = accum;
				accum += dz_col[col_elem] = find_height (blk_id[col][col_elem]);
			}

			// store total accumulated height at the end for each column
			ts.h_tot[col] = accum;
		}
	}
};

TopSurf*
TopSurf::create (const UnstructuredGrid& fine_grid) {
	unique_ptr <TopSurf> ts (new TopSurf);

	// outsource the entire construction to a builder object
	TopSurfBuilder (fine_grid, *(ts.get ()));
	compute_geometry (ts.get ());

	// client owns pointer to constructed grid from this point
	return ts.release ();
}

TopSurf::TopSurf ()
	: col_cells (0)
	, col_cellpos (0)
	, fine_col (0)
	, dz (0)
	, z0 (0)
	, h_tot (0) {
	// zero initialize all members that come from UnstructuredGrid
	// since that struct is a C struct, it doesn't have a ctor
	dimensions = 0;
	number_of_cells = 0;
	number_of_faces = 0;
	number_of_nodes = 0;
	face_nodes = 0;
	face_nodepos = 0;
	face_cells = 0;
	cell_faces = 0;
	cell_facepos = 0;
	node_coordinates = 0;
	face_centroids = 0;
	face_areas = 0;
	face_normals = 0;
	cell_centroids = 0;
	cell_volumes = 0;
	global_cell = 0;
	cartdims[0] = 0;
	cartdims[1] = 0;
	cartdims[2] = 0;
	cell_facetag = 0;
}

TopSurf::~TopSurf () {
	// deallocate memory that may have been created. if the dtor is
	// called from throwing an exception, the members should be zero
	// initialized, so it's OK to send them to delete.
	delete [] face_nodes;
	delete [] face_nodepos;
	delete [] face_cells;
	delete [] cell_faces;
	delete [] cell_facepos;
	delete [] node_coordinates;
	delete [] face_centroids;
	delete [] face_areas;
	delete [] face_normals;
	delete [] cell_volumes;
	delete [] global_cell;
	delete [] cell_facetag;
	// these are the extra members that are TopSurf specific
	delete [] col_cells;
	delete [] col_cellpos;
	delete [] fine_col;
	delete [] dz;
	delete [] h;
	delete [] z0;
	delete [] h_tot;
}
