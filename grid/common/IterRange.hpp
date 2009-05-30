//===========================================================================
//
// File: IterRange.hpp
//
// Created: Mon Apr 23 14:08:39 2007
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

#ifndef OPENRS_ITERRANGE_HEADER
#define OPENRS_ITERRANGE_HEADER


#include <set>
#include <iostream>


namespace Dune
{

    /** Represents a range by storing iterators to the range,
     *  so an IterRange may be invalidated whenever the data it refers to
     *  is modified.
     */
    template <typename T>
    class IterRange
    {
    public:
	typedef const T* iterator;
	IterRange()
	    : start_(0), finish_(0)
	{
	}
	template <typename Container>
	explicit IterRange(const Container& cont)
	    : start_(0), finish_(0)
	{
	    assign(cont.begin(), cont.end());
	}
	template <typename FI>
	IterRange(FI start, FI finish)
	    : start_(0), finish_(0)
	{
	    // Cannot easily initialize start_ and finish_ in the init list in case of size 0.
	    assign(start, finish);
	}
	template <typename FI>
	void assign(FI start, FI finish)
	{
	    int size = finish - start;
	    if (size == 0) {
		start_ = 0;
		finish_ = 0;
	    } else {
		start_ = &*start;
		finish_ = start_ + size;
	    }
	}
	bool empty() const
	{
	    return finish_ == start_;
	}
	int size() const
	{
	    return int(finish_ - start_);
	}
	iterator begin() const
	{
	    return start_;
	}
	iterator end() const
	{
	    return finish_;
	}
	const T& operator[] (int index) const
	{
	    return start_[index];
	}
	template <class Functor>
	void each(Functor f) const
	{
	    for (iterator it = begin(); it < end(); ++it) {
		f(*it);
	    }
	}
	template <class DataVector, class DataFunctor>
	void each_data(DataVector dv, DataFunctor f) const
	{
	    for (iterator it = begin(); it < end(); ++it) {
		f(dv[*it]);
	    }
	}
	void write(std::ostream& os) const
	{
	    for (iterator it = begin(); it < end(); ++it) {
		os << *it << ' ';
	    }
	}

    private:
	iterator start_;
	iterator finish_;
    };





    template <typename T>
    inline bool operator==(const IterRange<T>& r1, const IterRange<T>& r2)
    {
	std::set<int> s1(r1.begin(), r1.end());
	std::set<int> s2(r2.begin(), r2.end());
	return s1 == s2;
    }





    template <typename T>
    inline std::ostream& operator << (std::ostream& os, const IterRange<T>& v)
    {
	v.write(os);
	return os;
    }




} // namespace Dune




#endif // OPENRS_ITERRANGE_HEADER
