//===========================================================================
//
// File: EclipseGridParser.C
//
// Created: Thu Dec  6 08:46:05 2007
//
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//
// $Date$
//
// $Revision$
//
// Revision: $Id: EclipseGridParser.C,v 1.4 2008/08/18 14:16:14 atgeirr Exp $
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
#include "config.h"
#include <iostream>
#include <fstream>
#include <exception>
#include <algorithm>
#include <limits>
#include <cfloat>
#include "EclipseGridParser.hpp"
#include "SpecialEclipseFields.hpp"
#include <dune/common/ErrorMacros.hpp>

using namespace std;



namespace Dune
{

// ---------- List of supported keywords ----------

namespace EclipseKeywords
{
    string integer_fields[] =
	{ string("ACTNUM"),
	  string("SATNUM"),
	  string("EQLNUM"),
	  string("REGNUM"),
	  string("ROCKTYPE")
	};
    const int num_integer_fields = sizeof(integer_fields) / sizeof(integer_fields[0]);

    string floating_fields[] =
        { string("COORD"),    string("ZCORN"),      string("PERMX"),
          string("PERMY"),    string("PERMZ"),      string("PERMXX"),
          string("PERMYY"),   string("PERMZZ"),     string("PERMXY"),
          string("PERMYZ"),   string("PERMZX"),     string("PORO"),
          string("BULKMOD"),  string("YOUNGMOD"),   string("LAMEMOD"),
          string("SHEARMOD"), string("POISSONMOD"), string("PWAVEMOD")
        };
    const int num_floating_fields = sizeof(floating_fields) / sizeof(floating_fields[0]);

    string special_fields[] =
        { string("SPECGRID"), string("FAULTS"), string("MULTFLT")
        };
    const int num_special_fields = sizeof(special_fields) / sizeof(special_fields[0]);

    string ignore_with_data[] =
	{ string("MAPUNITS"), string("MAPAXES"),  string("GRIDUNIT"),
	  string("INCLUDE"),  string("DIMENS"),   string("START"),
	  string("REGDIMS"),  string("WELLDIMS"), string("TABDIMS"),
	  string("NSTACK"),   string("SWFN"),     string("SOF2"),
	  string("PVTW"),     string("PVTDO"),    string("ROCK"),
	  string("DENSITY"),  string("SATNUM"),   string("EQUIL"),
	  string("RPTRST"),   string("ROIP"),     string("RWIP"),
	  string("RWSAT"),    string("RPR"),      string("WBHP"),
	  string("WOIR"),     string("WELSPECS"), string("COMPDAT"),
	  string("WCONINJE"), string("RPTRST"),   string("TUNING"),
	  string("TSTEP")
	};
    const int num_ignore_with_data = sizeof(ignore_with_data) / sizeof(ignore_with_data[0]);

    string ignore_no_data[] =
	{ string("RUNSPEC"), string("WATER"),    string("OIL"),
	  string("METRIC"),  string("FMTIN"),    string("FMTOUT"),
	  string("GRID"),    string("INIT"),     string("NOECHO"),
	  string("ECHO"),    string("EDIT"),     string("PROPS"),
	  string("REGIONS"), string("SOLUTION"), string("SUMMARY"),
	  string("FPR"),     string("FOIP"),     string("FWIP"),
	  string("RUNSUM"),  string("EXCEL"),    string("SCHEDULE"),
	  string("END")
	};

    const int num_ignore_no_data = sizeof(ignore_no_data) / sizeof(ignore_no_data[0]);

} // namespace EclipseKeywords


// ---------- Helper functions for read() and hasField() ----

namespace
{

    istream& ignoreLine(istream& is)
    {
	is.ignore(numeric_limits<int>::max(), '\n');
	return is;
    }

    istream& ignoreWhitespace(istream& is)
    {
	// Getting the character type facet for is()
	// We use the classic (i.e. C) locale.
	const ctype<char>& ct = use_facet< ctype<char> >(locale::classic());
	char c;
	while (is.get(c)) {
	    if (!ct.is(ctype_base::space, c)) {
		is.putback(c);
		break;
	    }
	}
	return is;
    }

    string upcase(const string& s)
    {
	string us(s);
	// Getting the character type facet for toupper().
	// We use the classic (i.e. C) locale.
	const ctype<char>& ct = use_facet< ctype<char> >(locale::classic());
	for (int i = 0; i < int(s.size()); ++i) {
	    us[i] = ct.toupper(s[i]);
	}
	return us;
    }

    string read_keyword(istream& is)
    {
	string keyword_candidate;
	is >> keyword_candidate;
	while (keyword_candidate.find("--") == 0) {
	    // This line is a comment
	    is >> ignoreLine >> keyword_candidate;
	}
	return upcase(keyword_candidate);
    }

    enum FieldType {
	Integer,
	FloatingPoint,
	SpecialField,
	IgnoreWithData,
	IgnoreNoData,
	Unknown
    };

    pair<FieldType, bool> classify_keyword(const string& keyword)
    {
	using namespace EclipseKeywords;
	bool error_if_nonnumeric = true;
	if (count(integer_fields, integer_fields + num_integer_fields, keyword)) {
	    return make_pair(Integer, error_if_nonnumeric);
	} else if (count(floating_fields, floating_fields + num_floating_fields, keyword)) {
	    return make_pair(FloatingPoint, error_if_nonnumeric);
	} else if (count(special_fields, special_fields + num_special_fields, keyword)) {
	    error_if_nonnumeric = false;
	    return make_pair(SpecialField, error_if_nonnumeric);
	} else if (count(ignore_with_data, ignore_with_data + num_ignore_with_data, keyword)) {
	    return make_pair(IgnoreWithData, error_if_nonnumeric);
	} else if (count(ignore_no_data, ignore_no_data + num_ignore_no_data, keyword)) {
	    return make_pair(IgnoreNoData, error_if_nonnumeric);
	} else {
	    return make_pair(Unknown, error_if_nonnumeric);
	}
    }

    template<typename T>
    void read_data(istream& is, vector<T>& data, bool error_on_nonnumerics = true)
    {
	data.clear();
	while (is) {
	    T candidate;
	    is >> candidate;
	    if (is.rdstate() & ios::failbit) {
		is.clear(is.rdstate() & ~ios::failbit);
		string dummy;
		is >> dummy;
		if (dummy == "/") {
		    break;
		} else if (error_on_nonnumerics) {
		    cerr << "Encountered format error while reading data values." << endl;
		    throw exception();
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
	    cerr << "Encountered error while reading data values." << endl;
	    throw exception();
	}
#if 0
	data.clear();
	int next = is.peek();
	while (next != int('/')) {
	    // Read more data
	    T candidate;
	    is >> candidate;
	    if (is.rdstate() & ios::failbit) {
		if (error_on_nonnumerics) {
		    cerr << "Encountered format error while reading data values." << endl;
		    throw exception();
		}
		is.clear(is.rdstate() & ~ios::failbit);
		string dummy;
		is >> dummy;
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
	    is >> ignoreWhitespace;
	    next = is.peek();
	    if (next == EOF) {
		cerr << "Reached end of file before data end (/) was found." << endl;
		throw exception();
	    }
	}
	is.ignore(); // ignore the '/'
	if (!is) {
	    cerr << "Encountered error while reading data values." << endl;
	    throw exception();
	}
#endif
    }


} // anon namespace



// ---------- Member functions ----------

/// Constructor taking an eclipse file as a stream.
//---------------------------------------------------------------------------
EclipseGridParser::EclipseGridParser(istream& is)
//---------------------------------------------------------------------------
{
    read(is);
}


/// Constructor taking an eclipse filename.
//---------------------------------------------------------------------------
EclipseGridParser::EclipseGridParser(const string& filename)
//---------------------------------------------------------------------------
{
    ifstream is(filename.c_str());
    if (!is) {
	cerr << "Unable to open file " << filename << endl;
	throw exception();
    }
    read(is);
}


/// Read the given stream, overwriting any previous data.
//---------------------------------------------------------------------------
void EclipseGridParser::read(istream& is)
//---------------------------------------------------------------------------
{
    if (!is) {
	cerr << "Could not read given input stream." << endl;
	throw exception();
    }

    // Make temporary maps that will at the end be swapped with the
    // member maps.
    map<string, vector<int> > intmap;
    map<string, vector<double> > floatmap;
    map<string, boost::shared_ptr<SpecialBase> > specialmap;

    // Actually read the data
    is >> ignoreWhitespace;
    while (!is.eof()) {
	string keyword = read_keyword(is);
	cout << "Keyword: " << keyword << " " << keyword.size() << endl;
	pair<FieldType, bool> type = classify_keyword(keyword);
	switch (type.first) {
	case Integer:
	    read_data(is, intmap[keyword], type.second);
	    break;
	case FloatingPoint:
	    read_data(is, floatmap[keyword], type.second);
	    break;
	case SpecialField: {
	    boost::shared_ptr<SpecialBase> sb_ptr = createSpecialField(is, keyword);
	    if (sb_ptr) {
		specialmap[keyword] = sb_ptr;
	    } else {
		THROW("Could not create field " << keyword);
	    }
	    break;
	}
	case IgnoreWithData: {
	    is >> ignoreLine;
	    string dummy;
	    getline(is, dummy);
	    cout << "IgnoreWithData: " << dummy << " " << dummy.length()
		 << " L: " << dummy[dummy.length()-1] << endl;

	    if (dummy[dummy.length()-1] != '/') {
		getline(is, dummy, '/');
	    }
	    break;
	}
	case IgnoreNoData: {
	    is >> ignoreLine;
	    break;
	}
	case Unknown:
	default:
	    cerr << "Keyword " << keyword << " not recognized." << endl;
	    throw exception();
	}
	is >> ignoreWhitespace;
    }

    // Swap into member maps
    integer_field_map_.swap(intmap);
    floating_field_map_.swap(floatmap);
    special_field_map_.swap(specialmap);
}


/// Returns true if the given keyword corresponds to a field that
/// was found in the file.
//---------------------------------------------------------------------------
bool EclipseGridParser::hasField(const string& keyword) const
//---------------------------------------------------------------------------
{
    string ukey = upcase(keyword);
    return integer_field_map_.count(ukey) || floating_field_map_.count(ukey) ||
	special_field_map_.count(ukey);
}


/// Returns true if all the given keywords correspond to fields
/// that were found in the file.
//---------------------------------------------------------------------------
bool EclipseGridParser::hasFields(const vector<string>& keywords) const
//---------------------------------------------------------------------------
{
    int num_keywords = keywords.size();
    for (int i = 0; i < num_keywords; ++i) {
	if (!hasField(keywords[i])) {
	    return false;
	}
    }
    return true;
}

//---------------------------------------------------------------------------
vector<string> EclipseGridParser::fieldNames() const
//---------------------------------------------------------------------------
{
    vector<string> names;
    names.reserve(integer_field_map_.size() +
                  floating_field_map_.size());
    {
	map<string, vector<int> >::const_iterator it = integer_field_map_.begin();
	for (; it != integer_field_map_.end(); ++it) {
	    names.push_back(it->first);
	}
    }
    {
	map<string, vector<double> >::const_iterator it = floating_field_map_.begin();
	for (; it != floating_field_map_.end(); ++it) {
	    names.push_back(it->first);
	}
    }
    {
	map<string, boost::shared_ptr<SpecialBase> >::const_iterator it = special_field_map_.begin();
	for (; it != special_field_map_.end(); ++it) {
	    names.push_back(it->first);
	}
    }
    return names;
}

//---------------------------------------------------------------------------
const std::vector<int>& EclipseGridParser::getIntegerValue(const std::string& keyword) const
//---------------------------------------------------------------------------
{
    if (keyword == "SPECGRID") {
	cerr << "\nERROR. Interface has changed!\n"
	     << "const vector<int>& dim = parser.getIntegerValue(""SPECGRID"") is deprecated.\n"
	     << "Use:\n"
	     << "const SPECGRID& specgrid = parser.getSpecGrid();\n"
	     << "const vector<int>& dim = specgrid.dimensions;\n\n";
	throw exception();
    }

    map<string, vector<int> >::const_iterator it
	= integer_field_map_.find(keyword);
    if (it == integer_field_map_.end()) {
 	return empty_integer_field_;
    } else {
	return it->second;
    }
}

//---------------------------------------------------------------------------
const std::vector<double>& EclipseGridParser::getFloatingPointValue(const std::string& keyword) const
//---------------------------------------------------------------------------
{
    map<string, vector<double> >::const_iterator it
	= floating_field_map_.find(keyword);
    if (it == floating_field_map_.end()) {
 	return empty_floating_field_;
    } else {
	return it->second;
    }
}


//---------------------------------------------------------------------------
const boost::shared_ptr<SpecialBase> EclipseGridParser::getSpecialValue(const std::string& keyword) const
//---------------------------------------------------------------------------
{
    map<string, boost::shared_ptr<SpecialBase> >::const_iterator it = special_field_map_.find(keyword);
    if (it == special_field_map_.end()) {
 	return empty_special_field_;
    } else {
	return it->second;
    }
}

//---------------------------------------------------------------------------
const SPECGRID& EclipseGridParser::getSpecGrid() const
//---------------------------------------------------------------------------
{
    return dynamic_cast<const SPECGRID&>(*getSpecialValue("SPECGRID"));
}

//---------------------------------------------------------------------------
const FAULTS& EclipseGridParser::getFaults() const
//---------------------------------------------------------------------------
{
    return dynamic_cast<const FAULTS&>(*getSpecialValue("FAULTS"));
}

//---------------------------------------------------------------------------
const MULTFLT& EclipseGridParser::getMultflt() const
//---------------------------------------------------------------------------
{
    return dynamic_cast<const MULTFLT&>(*getSpecialValue("MULTFLT"));
}

//---------------------------------------------------------------------------
boost::shared_ptr<SpecialBase>
EclipseGridParser::createSpecialField(std::istream& is,
				      const std::string& fieldname)
//---------------------------------------------------------------------------
{
    string ukey = upcase(fieldname);
    boost::shared_ptr<SpecialBase> spec_ptr;
    if (ukey == "SPECGRID") {
	spec_ptr.reset(new SPECGRID);
	spec_ptr->read(is);
    } else if (ukey == "FAULTS") {
	spec_ptr.reset(new FAULTS);
	spec_ptr->read(is);
    } else if (ukey == "MULTFLT") {
	spec_ptr.reset(new MULTFLT);
	spec_ptr->read(is);
    } else {
	cerr << "Special Keyword " << fieldname << " not recognized." << endl;
	throw exception();
    }
    return spec_ptr;
}

} // namespace Dune
