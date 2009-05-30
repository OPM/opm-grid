//===========================================================================
//
// File: SparseTable.hpp
//
// Created: Fri Apr 24 09:50:27 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
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

#ifndef OPENRS_SPARSETABLE_HEADER
#define OPENRS_SPARSETABLE_HEADER

#include <vector>
#include <numeric>
#include <algorithm>
#include "IterRange.hpp"
#include "ErrorMacros.hpp"

namespace Dune
{

    /// A SparseTable stores a table with rows of varying size
    /// as efficiently as possible.
    /// It is supposed to behave similarly to a vector of vectors.
    /// Its behaviour is similar to compressed row sparse matrices.
    template <typename T>
    class SparseTable
    {
    public:
	/// Default constructor. Yields an empty SparseTable.
	SparseTable()
	{
	}


	/// A constructor taking all the data for the table and row sizes.
	/// \param data_beg The start of the table data.
	/// \param data_end One-beyond-end of the table data.
	/// \param rowsize_beg The start of the row length data.
	/// \param rowsize_end One beyond the end of the row length data.
	template <typename DataIter, typename IntegerIter>
	SparseTable(DataIter data_beg, DataIter data_end, IntegerIter rowsize_beg, IntegerIter rowsize_end)
	    : data_(data_beg, data_end)
	{
	    // Since we do not store the row sizes, but cumulative row sizes,
	    // we have to create the cumulative ones.
	    int num_rows = rowsize_end - rowsize_beg;
	    if (num_rows < 1) {
		THROW("Must have at least one row. Got " << num_rows << " rows.");
	    }
#ifndef NDEBUG
	    if (*std::min_element(rowsize_beg, rowsize_end) < 0) {
		THROW("All row sizes must be at least 0.");
	    }
#endif
	    row_start_.resize(num_rows + 1);
	    row_start_[0] = 0;
	    std::partial_sum(rowsize_beg, rowsize_end, row_start_.begin() + 1);
	    // Check that data_ and row_start_ match.
	    if (int(data_.size()) != row_start_.back()) {
		THROW("End of row start indices different from data size.");
	    }
	}
	    

	/// True if the table contains no data.
	bool empty() const
	{
	    return data_.empty();
	}


	/// Returns the number of rows in the table.
	int size() const
	{
	    return empty() ? 0 : row_start_.size() - 1;
	}


	/// Returns a row of the table.
	IterRange<T> operator[](int row) const
	{
	    const T* start_ptr = &data_[0];
	    return IterRange<T>(start_ptr + row_start_[row], start_ptr + row_start_[row + 1]);
	}

    private:
	std::vector<T> data_;
	// Like in the compressed row sparse matrix format,
	// row_start_.size() is equal to the number of rows + 1.
	std::vector<int> row_start_;
    };

} // namespace Dune


#endif // OPENRS_SPARSETABLE_HEADER
