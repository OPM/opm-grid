#ifndef OPM_VERTEQ_NAV_HPP_INCLUDED
#define OPM_VERTEQ_NAV_HPP_INCLUDED

// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#ifndef OPM_GRID_HEADER_INCLUDED
#include <opm/grid/UnstructuredGrid.h>
#endif

#include <cstdlib>	// div_t
#include <iosfwd>   // ostream

/**
 * There are three types of indices used in this module:
 *
 * (1) Cartesian coordinate
 * (2) Cartesian index
 * (3) Global identity
 *
 * The Cartesian coordinate is an (i,j)-tuple into the cornerpoint grid
 * structure.
 *
 * The Cartesian index, is a flat integer which has is
 * determined solely by the structure of the grid, regardless of whether
 * there are any active elements. It can be calculated from the coordinates
 * if the extent of the grid is known. We use this to efficiently store
 * data for each possible position in the grid without using a predefined
 * size or dynamic structure for each column.
 *
 * The global identity is the index which is assigned to it in the grid
 * structure, after inactive elements are discarded. This is the 'identity'
 * of the cell in the rest of the simulator. Cells that aren't active are
 * not assigned an identity.
 *
 * The value types defined here provide a way to address location in
 * the grid in a type-safe manner to let the compiler help us keep track
 * of the real meaning of indices.
 *
 * The navigation classes provide a way to define an enumeration of the
 * grid without resorting to inline integer arithmetic inside the other
 * functions. An optimizing compiler should be able to generate equally
 * fast code as hand-coded index calculations, when using these classes.
 */

/**
 * Index tuple in two-dimensional cornerpoint grid
 *
 * This structure represents the carrier of Cartesian coordinates. They
 * should be thought of as an integral type, along the lines of complex
 * numbers.
 */
struct Coord2D {
	Coord2D (int i_, int j_) : m_i (i_), m_j (j_) { }

	int i() const { return m_i; }
	int j() const { return m_j; }

	/**
	 * Compare two coordinates
	 */
	bool operator == (const Coord2D& rhs) const {
		return (m_i == rhs.m_i) && (m_j == rhs.m_j);
	}

protected:
	const int m_i;
	const int m_j;

	friend std::ostream& operator << (std::ostream& s, const Coord2D& c);
};

/// Index tuple in three-dimensional cornerpoint grid.
struct Coord3D : public Coord2D {
	Coord3D (int i_, int j_, int k_)
		: Coord2D (i_, j_)
		, m_k (k_) {
	}

	int k() const { return m_k; }

protected:
	const int m_k;

	friend std::ostream& operator << (std::ostream& s, const Coord3D& c);
};

// forward declaration
template <typename Dim > struct Side;

/// Type-safe enumeration of axis directions.
struct Dir {
	/// Towards the end of the axis with lesser numbers.
	static const Dir DEC; // = 0

	/// Towards the end of the axis with greater numbers.
	static const Dir INC; // = 1

	/// Number of possible directions
	static const int COUNT = 2;

	/// Integer representation suitable for indexing in array
	const int val;

	Dir (const Dir& rhs) : val (rhs.val) {}
	bool operator == (const Dir& rhs) const { return val == rhs.val; }

	/// Opposite direction to this one
	Dir opposite () const { return Dir (-(val - 1)); }

protected:
	/// Private constructor to avoid initialization outside domain
	Dir (int i) : val (i) { }

	template <typename Dim> friend struct Side;

	friend std::ostream& operator << (std::ostream& os, const Dir& d);
};

struct Dim1D {
	static const int COUNT = 1;
};

/// Type-safe enumeration of axis dimensions
struct Dim2D : public Dim1D {
	// two spatial directions
	static const Dim2D X; // = 0
	static const Dim2D Y; // = 1

	// number of dimensions
	static const int COUNT = 2;

	const int val;

	Dim2D (const Dim2D& rhs) : val (rhs.val) { }
	bool operator == (const Dim2D& rhs) const { return val == rhs.val; }

	/// Orthogonal dimension to this one
	Dim2D orthogonal () const { return Dim2D (-(val - 1)); }

protected:
	Dim2D (int i) : val (i) { }

	friend struct Side <Dim2D>;

	friend std::ostream& operator << (std::ostream& os, const Dim2D& d);
};

/// Type-safe enumeration of axis dimensions in 3D
struct Dim3D : public Dim2D {
	// added dimension in 3D
	static const Dim3D Z; // = 2

	// number of dimensions (shadows Dim2D)
	static const int COUNT = 3;

	Dim3D (const Dim3D& rhs) : Dim2D (rhs.val) { }
	bool operator == (const Dim2D& rhs) const { return val == rhs.val; }

	// allow X and Y to work in 3D too
	Dim3D (const Dim2D& rhs) : Dim2D (rhs) {}

protected:
	Dim3D (int i) : Dim2D (i) { }

	friend struct Side <Dim3D>;
};

/**
 * Value type that addresses sides in a n-dimensional grid cell.
 * A side is identified by a dimensions and a direction in that
 * dimension. The side will be located in that direction in that
 * dimension from the center of the cell.
 */
template <typename Dim>
struct Side {
	Side (Dim dim_, Dir dir_) : m_dim (dim_), m_dir (dir_) { }
	Side (const Side& rhs) : m_dim (rhs.m_dim), m_dir (rhs.m_dir) { }

	/**
	 * Numeric tag of an enumeration of the sides. The sides are enumerated
	 * in the same order the dimensions are specified in the grid; I, J then
	 * K, wih the decreasing direction first, then the increasing direction.
	 *
	 * @see UnstructuredGrid.cell_facetag
	 */
	int facetag () const {
		return dim().val * Dir::COUNT + dir().val;
	}

	/**
	 * Construct a side value from the facetag stored in the grid structure
	 */
	static Side <Dim> from_tag (int tag);

	Dim dim() const { return m_dim; }
	Dir dir() const { return m_dir; }

	/**
	 * Number of possible sides in an element
	 */
	static const int COUNT = Dim::COUNT * Dir::COUNT;

	/**
	 * Iterator for all possible sides
	 */
	static const Side* begin () { return &ALL[0]; }
	static const Side* end () { return &ALL[COUNT]; }

	/**
	 * Comparison of two sides
	 */
	bool operator == (const Side <Dim>& rhs) const {
		return (m_dim == rhs.m_dim) && (m_dir == rhs.m_dir);
	}

protected:
	const Dim m_dim;
	const Dir m_dir;

	// fixed enumeration of all sides
	static const Side ALL [];

	template <typename T>
	friend std::ostream& operator << (std::ostream& os, const Side<T>& s);
};

// specializations for the dimensions we work with
typedef Side <Dim2D> Side2D;
typedef Side <Dim3D> Side3D;

// forward declaration of stream operator present in library
std::ostream& operator << (std::ostream& os, const Side2D& s);
std::ostream& operator << (std::ostream& os, const Side3D& s);

// standalone constants for sides that we use; we call them 'up' and
// 'down' so that U and D are mnemonics, in contrast to 'top' and 'bottom'
// where the 'b' would conflict with 'back'.
extern const Side3D UP;
extern const Side3D DOWN;

/**
 * Value type that addresses corners in a two-dimensional grid cell.
 * A corner is identified by directions in both dimensions.
 */
struct Corn2D {
	Corn2D (Dir i_, Dir j_) : m_i (i_), m_j (j_) { }
	Corn2D (const Corn2D& rhs) : m_i (rhs.m_i), m_j (rhs.m_j) { }

	Dir i() const { return m_i; }
	Dir j() const { return m_j; }

protected:
	const Dir m_i;
	const Dir m_j;
};

/**
 * Three-dimensional corner type. It inherits from the two-dimensional
 * one since we can project a corner onto the two-dimensional surface.
 */
struct Corn3D : public Corn2D {
	Corn3D (Dir i_, Dir j_, Dir k_) : Corn2D (i_, j_), m_k (k_) { }
	Corn3D (const Corn3D& rhs) : Corn2D (rhs), m_k (rhs.m_k) { }

	/**
	 * Initialize a new corner where one dimension has been (optionally)
	 * pivoted in another direction. We use this to correct classifier
	 * information about a corner; by enumerating all the vertices of an
	 * element through the sides (pivoting a corner to the dimension in
	 * which its containing face is), we can figure out in which corner
	 * they belong in.
	 */
	Corn3D pivot (Dim3D dim, Dir dir) {
		return Corn3D (dim == Dim3D::X ? dir : m_i,
									 dim == Dim3D::Y ? dir : m_j,
									 dim == Dim3D::Z ? dir : m_k);
	}

	Dir k() const { return m_k; }

	/**
	 * Compare two corners
	 */
	bool operator == (const Corn3D& rhs) const {
		return (m_i == rhs.m_i) && (m_j == rhs.m_j) && (m_k == rhs.m_k);
	}

protected:
	const Dir m_k;

	friend std::ostream& operator << (std::ostream& os, const Corn3D& c);
};

/**
 * Navigate a Cartesian grid in a structured way so that clearly defined
 * mapping between the enumeration index and the coordinate.
 */
struct Cart2D {
	// number of cells in each direction
	const int ni;
	const int nj;

	// initialize from the size of the grid
	Cart2D (int ni_, int nj_) : ni (ni_), nj (nj_) { }

	// use these two value types for structured coordinates and
	// flattened indices, respectively
	typedef Coord2D coord_t;
	typedef int elem_t;
	typedef int node_t;
	typedef int face_t;

	/// Value used to indicate that a reference is not to a valid element
	static const int NO_ELEM; // = -1
	static const int NO_FACE; // = -1
	static const int NO_NODE; // = -1

	/// Number of (possible) elements in the grid
	int num_elems () const {
		return ni * nj;
	}

	/// Cartesian (flattened) index for a coordinate
	elem_t cart_ndx (const coord_t& coord) const {
		return coord.j() * ni + coord.i();
	}

	/// Cartesian coordinate for a (flattened) index
	coord_t coord (const elem_t& cart_ndx) const {
		const std::div_t strip = std::div (cart_ndx, ni);
		const int i = strip.rem;
		const int j = strip.quot;
		return coord_t (i, j);
	}

	/**
	 * As each element has points on both sides (in both dimensions), there
	 * is an extra row and column of points compared to elements.
	 */
	int num_nodes () const {
		return (ni + 1) * (nj + 1);
	}

	node_t node_ndx (const coord_t& coord, const Corn2D& corn) {
		return (coord.j() + corn.j().val) * (ni + 1) + (coord.i() + corn.i().val);
	}

	/**
	 * Each column has one more faces than there are rows, and each row
	 * has one more face than there are columns, due to the boundary.
	 */
	int num_faces () const {
		// there are ni+1 faces oriented in j-direction in each column,
		// for a total of (ni+1)*nj, and nj+1 faces in i-direction in
		// each row, for a total of (nj+1)*ni.
		return (ni + 1) * nj + ni * (nj + 1);
	}

	/**
	 * Translate element coordinate plus relative side of this into a
	 * Cartesian index of the face.
	 */
	face_t face_ndx (const coord_t& coord, const Side2D& side) {
		// flags that code 1 (include) or 0 (don't) for each of the tests
		const int dirv = side.dir().val;
		const int idim = side.dim() == Dim2D::X ? 1 : 0;
		const int idir = idim ? dirv : 0;
		const int jdir = idim ? 0 : dirv;
		// there is a left side and a lower side face for each element in a
		// column, plus one face at the top of each column; 2*ni+1. the j-
		// coordinate plus one extra if we are selecting the right face, tells
		// us which column which is to the left of or "above" the face.
		// if the side is in the i-direction to the center of the element, then
		// the face is aligned with the j axis, so we skip idim*ni i-aligned
		// faces before it.
		// finally, determine the row by using the i-coordinate, and i-direction
		// to determine whether it is the lower or upper face (if applicable).
		return (coord.j() + jdir) * (2 * ni + 1) + (idim * ni) + (coord.i() + idir);
	}
};

/**
 * Navigate a three-dimensional grid.
 *
 * In this module, we only need this to get the structured index of
 * each three-dimensional element, in order to know into which column
 * we should assign this element. However, we keep the design in the
 * same manner as the two-dimensional case.
 */
struct Cart3D {
	// number of cells in each direction
	const int ni;
	const int nj;
	const int nk;

	/// Initialize POD from an existing (3D) grid
	Cart3D (const UnstructuredGrid& g)
		: ni (g.cartdims [0])
		, nj (g.cartdims [1])
		, nk (g.cartdims [2]) { }

	/// Project grid into a surface
	Cart2D project () const {
		return Cart2D (ni, nj);
	}

	// use these two value types for structured coordinates and
	// flattened indices, respectively
	typedef Coord3D coord_t;
	typedef int elem_t;

	/// Deconstruct Cartesian index into coordinates
	coord_t coord (const elem_t& cart_ndx) const {
		// the i-index moves fastest, as this is Fortran-indexing
		const div_t strip = div (cart_ndx, ni);
		const int i = strip.rem;
		const div_t plane = div (strip.quot, nj);
		const int j = plane.rem;
		const int k = plane.quot;
		return coord_t (i, j, k);
	}
};

#endif // OPM_VERTEQ_NAV_HPP_INCLUDED
