//===========================================================================
//
// File: NonuniformTableLinear.hpp
//
// Created: Tue Oct 21 13:25:34 2008
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

#ifndef OPENRS_NONUNIFORMTABLELINEAR_HEADER
#define OPENRS_NONUNIFORMTABLELINEAR_HEADER

#include <cassert>
#include <cmath>
#include <exception>

#include <dune/common/ErrorMacros.hpp>
#include <dune/solvers/common/linearInterpolation.hpp>

namespace Dune {
    namespace utils {


	/// @brief
	/// @todo Doc me!
	/// @tparam
	struct ValueOutOfRangeException : public std::exception {};

	/// \brief This class uses linear interpolation to compute the value
	///        (and its derivative) of a function sampled at possibly
	///         nonuniform points.
	template<class vector_t>
	class NonuniformTableLinear {
	public:


	    /// @brief
	    /// @todo Doc me!
	    NonuniformTableLinear();

	    /// @brief
	    /// @todo Doc me!
	    /// @param
	    NonuniformTableLinear(const vector_t& x_values,
				  const vector_t& y_values);


	    /// @brief
	    /// @todo Doc me!
	    /// @param
	    /// @return
	    double operator()(const double x) const;

	    /// @brief
	    /// @todo Doc me!
	    /// @param
	    /// @return
	    double derivative(const double x) const;


	    /// @brief
	    /// @todo Doc me!
	    enum RangePolicy {Throw = 0, ClosestValue = 1, Extrapolate = 2};
	    /// @brief
	    /// @todo Doc me!
	    /// @param
	    void setLeftPolicy(RangePolicy rp);
	    /// @brief
	    /// @todo Doc me!
	    /// @param
	    void setRightPolicy(RangePolicy rp);

	protected:
	    vector_t x_values_;
	    vector_t y_values_;
	    RangePolicy left_;
	    RangePolicy right_;
	};


	// A utility function
	/// @brief
	/// @todo Doc me!
	/// @tparam
	/// @param
	/// @return
	template <typename FI>
	bool isNondecreasing(FI beg, FI end)
	{
	    FI it = beg;
	    ++it;
	    FI prev = beg;
	    for (; it != end; ++it, ++prev) {
		if (*it < *prev) {
		    return false;
		}
	    }
	    return true;
	}



	// Member implementations.

	template<class vector_t>
	inline
	NonuniformTableLinear<vector_t>
	::NonuniformTableLinear()
	    : left_(ClosestValue), right_(ClosestValue)
	{
	}

	template<class vector_t>
	inline
	NonuniformTableLinear<vector_t>
	::NonuniformTableLinear(const vector_t& x_values,
				const vector_t& y_values)
	    : x_values_(x_values), y_values_(y_values),
	      left_(ClosestValue), right_(ClosestValue)
	{
	    assert(isNondecreasing(x_values.begin(), x_values.end()));
	}

	template<class vector_t>
	inline void
	NonuniformTableLinear<vector_t>
	::setLeftPolicy(RangePolicy rp)
	{
	    if (rp != ClosestValue) {
		THROW("Only ClosestValue RangePolicy implemented.");
	    }
	    left_ = rp;
	}

	template<class vector_t>
	inline void
	NonuniformTableLinear<vector_t>
	::setRightPolicy(RangePolicy rp)
	{
	    if (rp != ClosestValue) {
		THROW("Only ClosestValue RangePolicy implemented.");
	    }
	    right_ = rp;
	}

	template<class vector_t>
	inline double
	NonuniformTableLinear<vector_t>
	::operator()(const double x) const
	{
	    return linearInterpolation(x_values_, y_values_, x);
	}


	template<class vector_t>
	inline double
	NonuniformTableLinear<vector_t>
	::derivative(const double x) const
	{
	    return linearInterpolationDerivative(x_values_, y_values_, x);
	}

    } // namespace utils
} // namespace Dune

#endif // OPENRS_NONUNIFORMTABLELINEAR_HEADER
