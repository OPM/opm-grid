#ifndef OPM_VERTEQ_TOPSURF_HPP_INCLUDED
#define OPM_VERTEQ_TOPSURF_HPP_INCLUDED

// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

#ifndef OPM_GRID_HEADER_INCLUDED
#include <opm/grid/UnstructuredGrid.h>
#endif

namespace Opm {

/**
 * Two-dimensional top surface of a full, three-dimensional grid.
 *
 * This grid is set up such that each cell is an upscaling of all the cells
 * in each column in the full grid. It also contains a mean to map results
 * from this grid back to the full grid.
 *
 * The full grid is also referred to as the fine grid, and this grid as the
 * coarse grid, or upscaled, grid.
 *
 * Note: Do NOT call destroy_grid () when done with this structure; it will
 * only clean up half of it. Wrap it is a smart pointer that calls the
 * destructor.
 */
struct TopSurf : public UnstructuredGrid {
	virtual ~TopSurf ();

	/**
	 * Indices of the columns' underlaying cells in the full grid.
	 *
	 * Consecutive indices from the _fine_ grid, not this one, for each of the
	 * columns, i.e. cells in _this_ grid, as one flat list.
	 *
	 * The values are sorted in z-order, starting from the top and moving
	 * downwards to the bottom. This is useful because you can keep a running
	 * counter for the depth, filling items as you go. (For this to be really
	 * useful, the original grid should be reordered so that cells in the
	 * z-direction are closer).
	 *
	 * Use this field together with the col_cellpos to iterate through a column
	 * in the fine grid.
	 *
	 * @example
	 * @code{.cpp}
	 * TopSurf* ts = ...;
	 * rlw_int col_cells (ts->number_of_cells, ts->col_cellpos, ts->col_cells);
	 * for (int col = 0; col < col_cells.cols(); ++col) {
	 *   for (int block = 0; block < col_cells.size (col); ++block) {
	 *      ... col_cells[col][block] ...
	 *   }
	 * }
	 * @endcode
	 *
	 * @see TopSurf::column, TopSurf::col_cellpos
	 */
	int* col_cells;

	/**
	 * Number of cells in the columns preceeding each one.
	 *
	 * For each column c, the number col_cellpos[c] is the number of cells in
	 * the _full_ grid that belongs to the columns 0..(c-1).
	 *
	 * This arrangement means that col_cellpos[c] is the index into col_cells
	 * of the first fine cell, whereas col_cellpos[c+1] is the index into
	 * col_cells of the last fine cell.
	 *
	 * @see TopSurf::column, TopSurf::col_cellpos
	 */
	int* col_cellpos;

	/**
	 * Maximum vertical resolution, in blocks.
	 *
	 * This holds the largest number of blocks there is in any column in the
	 * fine grid, i.e. max_vert_res >= col_cellpos[i+1] - col_cellpos[i], for
	 * any i in [0,number_of_cells-1]. Use this measure to allocate sufficient
	 * space for temporary storage that holds column data.
	 */
	int max_vert_res;

	/**
	 * Mapping from underlaying fine grid into top surface grid.
	 *
	 * For each element e in the fine grid, cell_col[e] is the index of the
	 * column/cell in the top surface. The number of cells in the fine grid
	 * can be found in col_cellpos[number_of_cells+1].
	 *
	 * Note: The indices in this array is NOT cell indices of this grid,
	 * but rather of the underlaying grid. Instead, it are the values that
	 * are stored in this array which are element identities.
	 *
	 * @see TopSurf::col_cellpos, TopSurf::col_cells
	 */
	int* fine_col;

	/**
	 * Height in each fine grid block, setup in columns.
	 *
	 * For each column, there is a consecutive list of height for each fine
	 * grid block in that column. This array has the same format as the
	 * col_cells run-length matrix, and the heights given here correspond to
	 * the indices in that matrix.
	 *
	 * The height of a block is defined as the z-coordinate difference
	 * between the centroid of the top face and the centroid of the bottom
	 * face.
	 *
	 * @see TopSurf::col_cells
	 */
	double* dz;

	/**
	 * Reference height for each column.
	 *
	 * This is a flat array with number_of_elements items.
	 *
	 * The reference height of a column is defined as the z-coordinate of
	 * the centroid of the top face of the upper block in the column. From
	 * these values and (a subset of) the values in face_centroids it is
	 * possible to recreate the 2.5D surface of the top.
	 */
	double* z0;

	/**
	 * Height from top of column down to each fine grid block.
	 *
	 * This is the accumulated sum of all dz preceeding this block, in each
	 * column. The first entry is thus always zero.
	 */
	double* h;

	/**
	 * Accumulated height of all blocks in each column.
	 *
	 * This is a flat array with number_of_elements items.
	 *
	 * Sum of all the heights of the blocks in each column. Since a gap between
	 * blocks is a violation of the vertical equilibrium assumption, this value
	 * should be the z-coordinate of the bottom face of the last block, up to
	 * numerical rounding errors.
	 *
	 * @see TopSurf::dz, TopSurf::z0
	 */
	double* h_tot;

	/**
	 * Create an upscaled grid based on a full, three-dimensional grid.
	 *
	 * @param fine Grid that should be upscaled.
	 *
	 * This must be a three-dimensional, Cartesian grid which has an active
	 * cluster which is without holes and which is convex. The UnstructuredGrid
	 * structure is used because it is the lingua franca of grids in the
	 * simulator, not because this method will handle every possible grid.
	 *
	 * This pointer is NOT adopted. The caller must still dispose of the grid.
	 *
	 * @return Upscaled, fine grid.
	 *
	 * The caller have the responsibility of disposing this grid; no other
	 * references will initially exist.
	 */
	static TopSurf* create (const UnstructuredGrid& fine);

private:
	/**
	 * @brief You are not meant to construct these yourself; use create ().
	 */
	TopSurf ();
};

} // namespace Opm

#endif // OPM_VERTEQ_TOPSURF_HPP_INCLUDED
