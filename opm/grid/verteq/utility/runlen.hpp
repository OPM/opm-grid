#ifndef OPM_VERTEQ_RUNLEN_HPP_INCLUDED
#define OPM_VERTEQ_RUNLEN_HPP_INCLUDED

// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

// forward declaration
struct UnstructuredGrid;

namespace Opm {

/**
 * Regards a set of (member) variables as a run-length encoded matrix.
 *
 * Each column can have a variable number of rows. Although the *values*
 * of the matrix can be changed, its sparsity cannot, i.e. one cannot
 * remove or add new elements to a column.
 *
 * Use this class to access and iterate over a run-length encoded matrix
 * in the format that is used by UnstructuredGrid without having to worry
 * about getting the indexing right.
 *
 * @tparam T Datatype for the extra data that should be stored for
 *           each element, e.g. double.
 *
 * @example
 * @code{.cpp}
 * RunLenView <int> faces_in_cell (
 *     g.number_of_cells,
 *     g.cell_facepos,
 *     g.cell_faces
 * );
 *
 * int num_local_faces = faces_in_cell.size (cellno);
 * int first_local_face = faces_in_cell [cellno] [0];
 * @endcode
 *
 * Notice if you want to loop through every item and know where you are
 * (because you intend to use this as an index in another matrix), you
 * can do:
 *
 * @example
 * @code{.cpp}
 * RunLenView <int> cell_faces (
 *     g.number_of_cells,
 *     g.cell_facepos,
 *     g.cell_faces
 * );
 *
 * for (int cell = 0; cell < cell_faces.cols(); ++cell) {
 *   for (int local_face = 0; local_face < cell_faces.size (cell); ++local_face) {
 *     int face_id = cell_faces[cell][local_face];
 *   }
 * }
 * @endcode
 *
 * @see Opm::RunLenData
 */
template <typename T>
class RunLenView {
protected:
    /**
     * Size information. pos has num_of_cols+1 items, pos[i] contains
     * the starting index of the data values for column i. (Since there
     * is one more element than there are columns, the last one is the
     * total number of elements). The number 0 is explicitly stored in
     * the first column to avoid special processing.
     */
    int num_of_cols;
    int* pos;

    /**
     * Data for each of the individual elements, stored consecutively
     * for each column located together, followed by the next column.
     */
    T* data;

public:
    /**
     * Construct a view into the run-length encoded matrix. The view is
     * only defined as long as the underlaying structures exist!
     *
     * @param number Number of cells (elements, faces, nodes)
     * @param pos_ptr Table of starting indices.
     *                If columns are "foo" and rows are "bar", then this
     *                is the member called "foo_barpos".
     * @param values Actual data storage.
     *               If columns are "foo" and rows are "bar", then this
     *               is the member called "foo_bars".
     */
    RunLenView (int num_cols, int* pos_ptr, T* values)
        // store them locally for later use
        : num_of_cols (num_cols)
        , pos (pos_ptr)
        , data (values) {
    }

    /**
     * Create another view of the same data.
     *
     * @param rhs View to a run-length-encoded matrix
     */
    RunLenView (const RunLenView& rhs)
        // copy all fields verbatim
        : num_of_cols (rhs.num_of_cols)
        , pos (rhs.pos)
        , data (rhs.data) {
    }

    /**
     * Access a column directly.
     *
     * @param col Index of the column to get
     * @return Pointer to the start of the column
     */
    T* operator [] (int col) const {
        return &data [pos [col]];
    }

    /**
     * Number of columns that are stored in the entire matrix.
     *
     * @return Number of columns.
     */
    int cols () const {
        return num_of_cols;
    }

    /**
     * Number of elements that are stored in one particular column.
     *
     * @param col Index of the column.
     * @return Number of elements.
     */
    int size (int col) const {
        return pos [col + 1] - pos [col];
    }

    /**
     * Quick accessor to get the last element. When we store accumulated
     * data in the array, this will quickly give us the total.
     *
     * Note that this is NOT the end iterator for the column.
     *
     * @param col Index of the column
     * @return Value of the last element. If there is no elements in
     *         this column, then the return value is undefined.
     */
    T& last (int col) const {
        return data [pos [col + 1] - 1];
    }
};

/**
 * Allocate a new vector of data for each element, accessible as
 * a zig-zag matrix.
 *
 * Use this kind of matrix when you want to enhance the grid structure
 * with some information per element, using the existing format.
 *
 * @see Opm::RunLenView
 */
template <typename T>
struct RunLenData : public RunLenView <T> {
    /**
     * Allocate a matrix based on sizes specified elsewhere. This is
     * useful if you want to supply with your own data.
     *
     * @param number  Number of entities (rows).
     * @param pos_ptr Number of elements to allocate in front of each
     *                entity; starting index in the array for this row.
     *
     * @see Opm::RunLenView::RunLenView
     */
    RunLenData (int number, int* pos_ptr)
        // allocate a new vector for the data, containing the needed
        // number of elements. note that there is only one new
        // operation is the parameter list, so there is no leakage if
        // an out-of-memory exception is thrown.
        : RunLenView <T> (number, pos_ptr, new T [pos_ptr [number]]) {
    }

    ~RunLenData () {
        // this member is initialized with data allocated in our ctor
        delete [] RunLenView <T>::data;
    }
};

// shorthands for most used types
typedef const RunLenView <int> rlw_int;
typedef const RunLenView <double> rlw_double;

// access common run-length encoded matrices in a grid structure
rlw_int grid_cell_facetag (const UnstructuredGrid& g);
rlw_int grid_cell_faces (const UnstructuredGrid& g);

} /* namespace Opm */

#endif /* OPM_VERTEQ_RUNLEN_HPP_INCLUDED */
