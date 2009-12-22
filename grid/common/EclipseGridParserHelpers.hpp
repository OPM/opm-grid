//===========================================================================
//
// File: EclipseGridParserHelpers.hpp
//
// Created: Tue Dec 22 11:35:32 2009
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

#ifndef OPENRS_ECLIPSEGRIDPARSERHELPERS_HEADER
#define OPENRS_ECLIPSEGRIDPARSERHELPERS_HEADER

#include <limits>
#include <string>
#include <istream>
#include <vector>
#include <dune/common/ErrorMacros.hpp>

namespace Dune
{

namespace
{

    inline std::istream& ignoreLine(std::istream& is)
    {
	is.ignore(std::numeric_limits<int>::max(), '\n');
	return is;
    }

    inline std::istream& ignoreWhitespace(std::istream& is)
    {
	// Getting the character type facet for is()
	// We use the classic (i.e. C) locale.
	const std::ctype<char>& ct = std::use_facet< std::ctype<char> >(std::locale::classic());
	char c;
	while (is.get(c)) {
	    if (!ct.is(std::ctype_base::space, c)) {
		is.putback(c);
		break;
	    }
	}
	return is;
    }

    inline std::string upcase(const std::string& s)
    {
	std::string us(s);
	// Getting the character type facet for toupper().
	// We use the classic (i.e. C) locale.
	const std::ctype<char>& ct = std::use_facet< std::ctype<char> >(std::locale::classic());
	for (int i = 0; i < int(s.size()); ++i) {
	    us[i] = ct.toupper(s[i]);
	}
	return us;
    }

    inline std::string readKeyword(std::istream& is)
    {
	std::string keyword_candidate;
	is >> keyword_candidate;
	while (keyword_candidate.find("--") == 0) {
	    // This line is a comment
	    is >> ignoreLine >> keyword_candidate;
	}
	return upcase(keyword_candidate);
    }


    template<typename T>
    inline void readData(std::istream& is, std::vector<T>& data, bool error_on_formatfailure = true)
    {
	data.clear();
	while (is) {
	    T candidate;
	    is >> candidate;
	    if (is.rdstate() & std::ios::failbit) {
		is.clear(is.rdstate() & ~std::ios::failbit);
		std::string dummy;
		is >> dummy;
		if (dummy == "/") {
		    break;
		} else if (dummy.find("--") == 0) {
		    is >> ignoreLine; // This line is a comment
		    continue;
		} else if (error_on_formatfailure) {
                    THROW("Encountered format error while reading data values.");
		}
	    } else {
		if (is.peek() == int('*')) {
		    is.ignore(); // ignore the '*'
		    int multiplier = int(candidate);
		    is >> candidate;
		    data.insert(data.end(), multiplier, candidate);
		} else {
		    data.push_back(candidate);
		}
	    }
	}
	if (!is) {
	    THROW("Encountered error while reading data values.");
	}
    }


} // anon namespace

} // namespace Dune


#endif // OPENRS_ECLIPSEGRIDPARSERHELPERS_HEADER
