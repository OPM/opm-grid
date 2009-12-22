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
#include "EclipseGridParserHelpers.hpp"
#include "SpecialEclipseFields.hpp"
#include <dune/common/ErrorMacros.hpp>
#include <boost/filesystem.hpp>

using namespace std;

//#define VERBOSE

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
	  string("ROCKTYPE"),
	  string("DIMENS"),
	  string("REGDIMS"),
	  string("WELLDIMS"),
	  string("TABDIMS"),
          string("FIPNUM"),
          string("GRIDFILE")
	};
    const int num_integer_fields = sizeof(integer_fields) / sizeof(integer_fields[0]);

    string floating_fields[] =
        { string("COORD"),    string("ZCORN"),      string("PERMX"),
          string("PERMY"),    string("PERMZ"),      string("PERMXX"),
          string("PERMYY"),   string("PERMZZ"),     string("PERMXY"),
          string("PERMYZ"),   string("PERMZX"),     string("PORO"),
          string("BULKMOD"),  string("YOUNGMOD"),   string("LAMEMOD"),
          string("SHEARMOD"), string("POISSONMOD"), string("PWAVEMOD"),
          string("MULTPV")
        };
    const int num_floating_fields = sizeof(floating_fields) / sizeof(floating_fields[0]);

    string special_fields[] =
        { string("SPECGRID"), string("FAULTS"), string("MULTFLT"),
          string("TITLE"), string("START"), string("DATES")
        };
    const int num_special_fields = sizeof(special_fields) / sizeof(special_fields[0]);

    string ignore_with_data[] =
	{ string("MAPUNITS"), string("MAPAXES"),  string("GRIDUNIT"),
	  string("DIMENS"),   string("NTG"),
	  string("REGDIMS"),  string("WELLDIMS"), string("TABDIMS"),
	  string("NSTACK"),   string("SWFN"),     string("SOF2"),
	  string("PVTW"),     string("PVTDO"),    string("ROCK"),
	  string("DENSITY"),  string("SATNUM"),   string("EQUIL"),
	  string("RPTRST"),   string("ROIP"),     string("RWIP"),
	  string("RWSAT"),    string("RPR"),      string("WBHP"),
	  string("WOIR"),     string("WELSPECS"), string("COMPDAT"),
	  string("WCONINJE"), string("TUNING"),   string("PVDO"),
	  string("TSTEP"),    string("BOX")
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
	  string("END"),     string("ENDBOX")
	};
    const int num_ignore_no_data = sizeof(ignore_no_data) / sizeof(ignore_no_data[0]);

    string include_keywords[] = { string("INCLUDE") };
    const int num_include_keywords = sizeof(include_keywords) / sizeof(include_keywords[0]);


} // namespace EclipseKeywords

namespace {

    enum FieldType {
	Integer,
	FloatingPoint,
	SpecialField,
	IgnoreWithData,
	IgnoreNoData,
        Include,
	Unknown
    };

    inline FieldType classifyKeyword(const string& keyword)
    {
	using namespace EclipseKeywords;
	if (count(integer_fields, integer_fields + num_integer_fields, keyword)) {
	    return Integer;
	} else if (count(floating_fields, floating_fields + num_floating_fields, keyword)) {
	    return FloatingPoint;
	} else if (count(special_fields, special_fields + num_special_fields, keyword)) {
	    return SpecialField;
	} else if (count(ignore_with_data, ignore_with_data + num_ignore_with_data, keyword)) {
	    return IgnoreWithData;
	} else if (count(ignore_no_data, ignore_no_data + num_ignore_no_data, keyword)) {
	    return IgnoreNoData;
	} else if (count(include_keywords, include_keywords + num_include_keywords, keyword)) {
	    return Include;
	} else {
	    return Unknown;
	}
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
    // Store directory of filename
    boost::filesystem::path p(filename);
    directory_ = p.parent_path().string();
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
    // member maps
    // NOTE: Above is no longer true, for easier implementation of 
    //       the INCLUDE keyword. We lose the strong exception guarantee,
    //       though (of course retaining the basic guarantee).
    map<string, vector<int> >& intmap = integer_field_map_;
    map<string, vector<double> >& floatmap = floating_field_map_;
    map<string, boost::shared_ptr<SpecialBase> >& specialmap = special_field_map_;

    // Actually read the data
    is >> ignoreWhitespace;
    while (!is.eof()) {
	string keyword = readKeyword(is);
#ifdef VERBOSE
	cout << "Keyword found: " << keyword << endl;
#endif
	FieldType type = classifyKeyword(keyword);
	switch (type) {
	case Integer:
	    readVectorData(is, intmap[keyword]);
	    break;
	case FloatingPoint:
	    readVectorData(is, floatmap[keyword]);
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
            is >> ignoreSlashLine;
#ifdef VERBOSE
            cout << "(ignored)" << endl;
#endif
	    break;
	}
	case IgnoreNoData: {
	    is >> ignoreLine;
#ifdef VERBOSE
            cout << "(ignored)" << endl;
#endif
	    break;
	}
        case Include: {
            is >> ignoreLine;
            is.ignore(numeric_limits<int>::max(), '\'');
            string include_filename;
            getline(is, include_filename, '\'');
            include_filename = directory_ + '/' + include_filename;
            ifstream include_is(include_filename.c_str());
            if (!include_is) {
                THROW("Unable to open INCLUDEd file " << include_filename);
            }
            read(include_is);
            is >> ignoreSlashLine;
            break;
        }
	case Unknown:
	default:
	    cerr << "Keyword " << keyword << " not recognized." << endl;
	    throw exception();
	}
	is >> ignoreWhitespace;
    }
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
boost::shared_ptr<SpecialBase>
EclipseGridParser::createSpecialField(std::istream& is,
				      const std::string& fieldname)
//---------------------------------------------------------------------------
{
    string ukey = upcase(fieldname);
    boost::shared_ptr<SpecialBase> spec_ptr
        = Factory<SpecialBase>::createObject(fieldname);
    spec_ptr->read(is);
    return spec_ptr;
}

} // namespace Dune
