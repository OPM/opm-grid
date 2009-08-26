//===========================================================================
//
// File: OrientedEntityTable.hpp
//
// Created: Wed Aug 26 11:13:20 2009
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

#ifndef OPENRS_ORIENTEDENTITYTABLE_HEADER
#define OPENRS_ORIENTEDENTITYTABLE_HEADER

#include "EntityRep.hpp"
#include <boost/static_assert.hpp>
#include <dune/common/SparseTable.hpp>
#include <map>
#include <climits>
#include <boost/algorithm/minmax_element.hpp>

/// The namespace Dune is the main namespace for all Dune code.
namespace Dune
{
    namespace cpgrid
    {


	/// @brief A class used as a row type for  OrientedEntityTable.
	/// @tparam codim_to Codimension.
	template <int codim_to>
	class OrientedEntityRange : private SparseTable< EntityRep<codim_to> >::row_type
	{
	public:
	    typedef EntityRep<codim_to> ToType;
	    typedef ToType* ToTypePtr;
	    typedef typename SparseTable<ToType>::row_type R;

	    /// @brief Default constructor yielding an empty range.
	    OrientedEntityRange()
		: R(ToTypePtr(0), ToTypePtr(0)), orientation_(true)
	    {
	    }
	    /// @brief Constructor taking a row type and an orientation.
	    /// @param R Row type
	    /// @param orientation True if positive orientation.
	    OrientedEntityRange(const R& r, bool orientation)
		: R(r), orientation_(orientation)
	    {
	    }
	    using R::size;
	    using R::empty;
	    /// @brief Random access operator.
	    /// @param subindex Column index.
	    /// @return Entity representation.
	    ToType operator[](int subindex) const
	    {
		ToType erep = R::operator[](subindex);
		return orientation_ ? erep : erep.opposite();
	    }
	private:
	    bool orientation_;
	};




	/// @brief Represents the topological relationships between
	/// sets of entities, for example cells and faces.
	///
	/// The purpose of this class is to hide the intricacies of
	/// handling orientations from the client code, otherwise a
	/// straight SparseTable would do.
	/// @tparam codim_from Codimension of ???
	/// @tparam codim_to Codimension of ???
	template <int codim_from, int codim_to>
	class OrientedEntityTable : private SparseTable< EntityRep<codim_to> >
	{
	public:
	    typedef EntityRep<codim_from> FromType;
	    typedef EntityRep<codim_to> ToType;
	    typedef OrientedEntityRange<codim_to> row_type; // ??? doxygen henter doc fra SparseTable
	    typedef SparseTable<ToType> super_t;

	    /// Default constructor.
	    OrientedEntityTable()
	    {
	    }

	    /// @brief Constructor taking iterators to a sequence of table
	    /// data and a sequence of row size data.
	    ///
	    /// These table data are in the same format as the underlying
	    /// SparseTable<int> constructor with the same signature.
	    /// @tparam DataIter Iterator to table data.
	    /// @tparam IntegerIter Iterator to  the row length data.
	    /// @param data_beg The start of the table data.
	    /// @param data_end One-beyond-end of the table data.
	    /// @param rowsize_beg The start of the row length data.
	    /// @param rowsize_end One beyond the end of the row length data.
	    template <typename DataIter, typename IntegerIter>
	    OrientedEntityTable(DataIter data_beg, DataIter data_end,
				IntegerIter rowsize_beg, IntegerIter rowsize_end)
		: super_t(data_beg, data_end, rowsize_beg, rowsize_end)
	    {
	    }

	    using super_t::empty;
	    using super_t::size;
	    using super_t::clear;
	    using super_t::appendRow;

	    /// @brief Given an entity e of codimension codim_from,
	    /// returns the number of neighbours of codimension codim_to.
	    /// @param e Entity representation.
	    /// @return the number of neighbours of codimension codim_to.
	    int rowSize(const FromType& e) const
	    {
		return super_t::rowSize(e.index());
	    }

	    /// @brief Given an entity e of codimension codim_from, returns a
	    /// row (an indirect container) containing its neighbour
	    /// entities of codimension codim_to.
	    /// @param e Entity representation.
	    /// @return A row of the table.
	    row_type operator[](const FromType& e) const
	    {
		return row_type(super_t::operator[](e.index()), e.orientation());
	    }

	    /// @brief Elementwise equality.
	    /// @param other The other element
	    /// @return Returns true if \b this and the \b other element are equal.
	    bool operator==(const OrientedEntityTable& other) const
	    {
		return super_t::operator==(other);
	    }

 	    /** @brief Prints the relation matrix corresponding to the table.

 	     Let the entities of codimensions f and t be given by
 	     the sets \f$E^f = { e^f_i } \f$ and \f$E^t = { e^t_j }\f$.
 	     A relation matrix R is defined by
	     \f{eqnarray*}{
	       R_{ij} &=& 0  \mbox{ if } e^f_i \mbox{ and } e^t_j \mbox{ are not neighbours }\\
               R_{ij} &=& 1  \mbox{ if they are neighbours with same orientation }\\
	       R_{ij} &=& -1 \mbox{ if they are neighbours with opposite orientation.}
	     \f}
	     @param os   The output stream.
	    */
	    void printRelationMatrix(std::ostream& os) const
	    {
		int columns = numberOfColumns();
		for (int i = 0; i < size(); ++i) {
		    FromType from_ent(i);
		    row_type r  = operator[](from_ent);
		    int cur_col = 0;
		    int next_ent = 0;
		    ToType to_ent = r[next_ent];
		    int next_print = to_ent.index();
		    while (cur_col < columns) {
			if (cur_col == next_print) {
			    if (to_ent.orientation()) {
				os << "  1";
			    } else {
				os << " -1";
			    }
			    ++next_ent;
			    if (next_ent >= r.size()) {
				next_print = columns;
			    } else {
				to_ent = r[next_ent];
				next_print = to_ent.index();
			    }
			} else {
			    os << "  0";
			}
			++cur_col;
		    }
		    os << '\n';
		}
	    }

	    /// @brief Makes the inverse relation, mapping codim_to entities
	    /// to their codim_from neighbours.
	    ///
	    /// Implementation note: The algorithm has been changed
	    /// to a three-pass O(n) algorithm.
	    /// @param inv  The OrientedEntityTable 
	    void makeInverseRelation(OrientedEntityTable<codim_to, codim_from>& inv) const
	    {
		// Find the maximum index used. This will give (one less than) the size
		// of the table to be created.
		int maxind = -1;
		for (int i = 0; i < size(); ++i) {
		    EntityRep<codim_from> from_ent(i, true);
		    row_type r = operator[](from_ent);
		    for (int j = 0; j < r.size(); ++j) {
			EntityRep<codim_to> to_ent = r[j];
			int ind = to_ent.index();
			maxind = std::max(ind, maxind);
		    }
		}
		// Build the new_sizes vector and compute datacount.
		std::vector<int> new_sizes(maxind + 1);
		int datacount = 0;
		for (int i = 0; i < size(); ++i) {
		    EntityRep<codim_from> from_ent(i, true);
		    row_type r = operator[](from_ent);
		    datacount += r.size();
		    for (int j = 0; j < r.size(); ++j) {
			EntityRep<codim_to> to_ent = r[j];
			int ind = to_ent.index();
			++new_sizes[ind];
		    }
		}
		// Compute the cumulative sizes.
		std::vector<int> cumul_sizes(new_sizes.size() + 1);
		cumul_sizes[0] = 0;
		std::partial_sum(new_sizes.begin(), new_sizes.end(), cumul_sizes.begin() + 1);
		// Using the cumulative sizes array as indices, we populate new_data.
		// Note that cumul_sizes[ind] is not kept constant, but incremented so that
		// it always gives the correct index for new data corresponding to index ind.
		std::vector<int> new_data(datacount);
		for (int i = 0; i < size(); ++i) {
		    EntityRep<codim_from> from_ent(i, true);
		    row_type r = operator[](from_ent);
		    for (int j = 0; j < r.size(); ++j) {
			EntityRep<codim_to> to_ent = r[j];
			int ind = to_ent.index();
			int data_ind = cumul_sizes[ind];
			new_data[data_ind] = to_ent.orientation() ? i : ~i;
			++cumul_sizes[ind];
		    }
		}
		inv = OrientedEntityTable<codim_to, codim_from>(new_data.begin(),
								new_data.end(),
								new_sizes.begin(),
								new_sizes.end());
	    }

	private:
	    int numberOfColumns() const
	    {
		int maxind = 0;
		for (int i = 0; i < size(); ++i) {
		    FromType from_ent(i);
		    row_type r  = operator[](from_ent);
		    for (int j = 0; j < r.size(); ++j) {
			maxind = std::max(maxind, r[j].index());
		    }
		}
		return maxind + 1;
	    }
	};


    } // namespace cpgrid
} // namespace Dune




#endif // OPENRS_ORIENTEDENTITYTABLE_HEADER
