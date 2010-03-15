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
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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
#include <boost/date_time/gregorian/gregorian.hpp>

namespace Dune
{

namespace
{

    inline std::istream& ignoreLine(std::istream& is)
    {
	is.ignore(std::numeric_limits<int>::max(), '\n');
	return is;
    }

    inline std::istream& ignoreSlashLine(std::istream& is)
    {
	is.ignore(std::numeric_limits<int>::max(), '/');
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
	while (!is.eof()) {
	    is >> keyword_candidate;
	    if(keyword_candidate.find("--") == 0) {
		is >> ignoreLine;  // This line is a comment
	    } else {
		return upcase(keyword_candidate);
	    }
	}
	return "CONTINUE";  // Last line in included file is a comment
    }

    inline std::string readString(std::istream& is)
    {
	std::string string_candidate;
	is >> string_candidate;
        const char quote('\'');
        int beg = string_candidate[0] == quote ? 1 : 0;
        int len = string_candidate[0] == quote ? string_candidate.size() - 2 : string_candidate.size();
        return string_candidate.substr(beg, len);
    }

    // Reads data until '/' or an error is encountered.
    template<typename T>
    inline void readVectorData(std::istream& is, std::vector<T>& data)
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
                    is >> ignoreLine;
		    break;
		} else if (dummy[0] == '-') {  // "comment test"
		    is >> ignoreLine; // This line is a comment
		} else {
                    THROW("Encountered format error while reading data values. Value = " << dummy);
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


    // Reads data items of type T. Not more than 'max_values' items.
    // Asterisks may be used to signify 'repeat counts'. 5*3.14 will
    // insert 3.14 five times. Asterisk followed by a space is used to
    // signify default values.  n* will default n consecutive quantities.
    template<class Vec>
    inline int readDefaultedVectorData(std::istream& is, Vec& data, int max_values)
    {
        ASSERT(int(data.size()) >= max_values);
        int num_values = 0;
        while (is) {
            typename Vec::value_type candidate;
            is >> candidate;
            if (is.rdstate() & std::ios::failbit) {
                is.clear(is.rdstate() & ~std::ios::failbit);
                std::string dummy;
                is >> dummy;
                if (dummy == "/") {
                    is >> ignoreLine;	
                    break;
		} else if (dummy[0] == '-') {  // "comment test"
                    is >> ignoreLine;   // This line is a comment
                } else {
                    THROW("Encountered format error while reading data values. Value = " << dummy);
                }
            } else {
                if (is.peek() == int('*')) {
                    is.ignore(); // ignore the '*'
                    int multiplier = (int)candidate;
                    if (is.peek() == int(' ')) {
                        num_values += multiplier;  // Use default value(s)
                    } else {
                        is >> candidate;         // Use candidate 'multipler' times
                        for (int i=0; i<multiplier; ++i, ++num_values) {
                            data[num_values] = candidate;
                        }
                    }
                } else {
                    data[num_values] = candidate;
                    ++num_values;
                }
            }
            if (num_values >= max_values) {
		is >> ignoreLine;
                break;
            }
        }
        if (!is) {
            THROW("Encountered error while reading data values.");
        }
        return num_values;
    }



    // Returns month number 1-12. Returns 0 if illegal month name. 
    int getMonthNumber(const std::string& month_name)
    {
        const int num_months = 12;
        std::string months[num_months] = {"JAN", "FEB", "MAR", "APR", "MAY", "JUN",
                                          "JLY", "AUG", "SEP", "OCT", "NOV", "DEC"};
        if (month_name == "JUL") {
            return 7;    // "JUL" is an acceptable alternative to 'JLY'
        }
        int m = 0;
        for (int i=0; i<num_months; ++i) {
            if (month_name == months[i]) {
                m = i+1;
                break;
            }
        }
        return m;
    }


    inline boost::gregorian::date readDate(std::istream& is)
    {
        while (is.peek() == int('-')) {
            is >> ignoreLine;   // This line is a comment
        }
        int day, year;
        std::string month_name;
        is >> day;
        month_name = readString(is);
        is >> year;
        ignoreSlashLine(is);
        int month = getMonthNumber(month_name);
        return boost::gregorian::date(boost::gregorian::greg_year(year),
                                      boost::gregorian::greg_month(month),
                                      boost::gregorian::greg_day(day));
    }

} // anon namespace

} // namespace Dune


#endif // OPENRS_ECLIPSEGRIDPARSERHELPERS_HEADER
