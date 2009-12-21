//===========================================================================
//
// File: SpecialEclipseFields.hpp
//
// Created: Mon Sep 21 14:09:54 2009
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

#ifndef OPENRS_SPECIALECLIPSEFIELDS_HEADER
#define OPENRS_SPECIALECLIPSEFIELDS_HEADER

#include <string>
#include <fstream>
#include <limits>
#include <dune/common/ErrorMacros.hpp>
#include <iterator>

namespace Dune
{

namespace
{

// Skip rest of the line.
std::istream& ignore_line(std::istream& is)
{
    is.ignore(std::numeric_limits<int>::max(), '\n');
    return is;
}

// Reads words until a slash is found.
bool read_slash(std::istream& is)
{
    while (is) {
	std::string dummy;
	is >> dummy;
	if (dummy == "/") {
	    return true;
	}
    }
    return false;
}

// Reads data items of type T. Not more than 'max_values' items.
// Asterisks may be used to signify 'repeat counts'. 5*3.14 will
// insert 3.14 five times. Asterisk followed by a space is used to
// signify default values.  n* will default n consecutive quantities.
template<typename T>
int read_data(std::istream& is, std::vector<T>& data, int max_values)
{
    int num_values = 0;
    while (is) {
	T candidate(-1);
	is >> candidate;
	if (is.rdstate() & std::ios::failbit) {
	    is.clear(is.rdstate() & ~std::ios::failbit);
	    std::string dummy;
	    is >> dummy;
	    if (dummy == "/") {
		is >> ignore_line;	
		break;
	    } else if (dummy == "--") {
		is >> ignore_line;   // This line is a comment
	    } else {
		THROW("Unexpected value reading file. Value = " << dummy);
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
	    break;
	}
    }
    if (!is) {
	THROW("Encountered error while reading data values.");
    }
    return num_values;
}
    
} // anon namespace

// Abstract base class for special fields. A special field is a field
// with more than one data type.
struct SpecialBase {
    virtual ~SpecialBase() {}                       // Default destructor
    virtual std::string name() const = 0;           // Keyword name
    virtual void read(std::istream& is) = 0;        // Reads data
    virtual void write(std::ostream& os) const = 0; // Writes data
};

/// Class for keyword SPECGRID 
struct SPECGRID : public SpecialBase
{
    std::vector<int> dimensions; // Number of grid blocks in x-, y- and z-directions.
    int numres;          // Number of reservoirs. 
    char qrdial;         // Coordinates. F=cartesian, T=Cylindrical(radial).

    SPECGRID()
    {
	dimensions.resize(3,1);
	numres = 1;
	qrdial = 'F';
    }

    virtual ~SPECGRID()
    {
    }

    virtual std::string name() const {return std::string("SPECGRID");}

    virtual void read(std::istream& is)
    {
	const int ndim = 3;
	std::vector<int> data(ndim+1,1);
	int nread = read_data(is , data, ndim+1);
	int nd = std::min(nread, ndim);
	copy(data.begin(), data.begin()+nd, &dimensions[0]);
	numres = data[ndim];
	std::string candidate;
	is >> candidate;
	if (candidate == "/") {
	    return;
	} else {
	    qrdial = candidate[0];
	}
	
	if (read_slash(is)) {
	    return;
	} else {
	    THROW("End of file reading" << name());
	}
    }
	
    virtual void write(std::ostream& os) const
    {
	os << name() << std::endl;
	os << dimensions[0] << " " << dimensions[1]  << " "
	   << dimensions[2] << " " << numres << " " << qrdial << std::endl;
	os << std::endl;
    }
};

/// Class holding segment data of keyword FAULTS.
struct FaultSegment
{
    std::string fault_name;          // Fault name
    std::vector<int> ijk_coord;      // ijk-coordinates of segment cells
    std::string face;                // Fault face of cells
};

/// Class for keyword FAULTS.
struct FAULTS : public SpecialBase
{
    std::vector<FaultSegment> faults;

    FAULTS()
    {
    }

    virtual ~FAULTS()
    {
    }

    virtual std::string name() const {return std::string("FAULTS");}

    virtual void read(std::istream& is)
    {
	while(is) {
	    std::string name;
	    is >> name;
	    if (name[0] == '/') {
		break;
	    }
	    while (name.find("--") == 0) {
		// This line is a comment
		is >> ignore_line >> name;
	    }
	    FaultSegment fault_segment;
	    fault_segment.ijk_coord.resize(6);
	    fault_segment.fault_name = name;
	    int nread = read_data(is, fault_segment.ijk_coord, 6);
	    if (nread != 6) {
		THROW("Error reading fault_segment " << name);
	    }
	    is >> fault_segment.face;
	    faults.push_back(fault_segment);
	    read_slash(is);
	}
    }

    virtual void write(std::ostream& os) const
    {
	os << name() << std::endl;
	for (int i=0; i<(int)faults.size(); ++i) {
	    os << faults[i].fault_name << "  ";
	    copy(faults[i].ijk_coord.begin(), faults[i].ijk_coord.end(),
		 std::ostream_iterator<int>(os, " "));
	    os << faults[i].face << std::endl;
	}
	os << std::endl;
    }
};

/// Class holding a data line of keyword MULTFLT
struct MultfltLine
{
    std::string fault_name;         // Fault name, as in FAULTS
    double transmis_multiplier;     // Transmissibility multiplier
    double diffusivity_multiplier;  // Diffusivity multiplier;
    MultfltLine() :
	fault_name(""), transmis_multiplier(1.0), diffusivity_multiplier(1.0)
    {
    }
};

/// Class for keyword MULFLT
struct MULTFLT : public SpecialBase
{
    std::vector<MultfltLine> multflts;

    MULTFLT()
    {
    }

    virtual ~MULTFLT()
    {
    }

    virtual std::string name() const {return std::string("MULTFLT");}

    virtual void read(std::istream& is)
    {
	while(is) {
	    std::string name;
	    is >> name;
	    if (name[0] == '/') {
		break;
	    }
	    while (name == "--") {
		// This line is a comment
		is >> ignore_line >> name;
	    }
	    MultfltLine multflt_line;
	    multflt_line.fault_name = name;
	    std::vector<double> data(2,1.0);
	    if (read_data(is, data, 2) == 2) {
		read_slash(is);
		is >> ignore_line;	
	    }
	    multflt_line.transmis_multiplier = data[0];
	    multflt_line.diffusivity_multiplier = data[1];
	    multflts.push_back(multflt_line);
	}
    }

    virtual void write(std::ostream& os) const
    {
	os << name() << std::endl;
	for (int i=0; i<(int)multflts.size(); ++i) {
	    os << multflts[i].fault_name << "  " 
	       << multflts[i].transmis_multiplier << "  "
	       << multflts[i].diffusivity_multiplier <<	std::endl;
	}
	os << std::endl;
    }
};


struct TITLE : public SpecialBase
{
    std::string title;
    virtual std::string name() const
    { return std::string("TITLE"); }
    virtual void read(std::istream& is)
    { is >> ignore_line; std::getline(is, title); }
    virtual void write(std::ostream& os) const
    { os << name() << '\n' << title << '\n'; }

};

} // End of namespace Dune

#endif // OPENRS_SPECIALECLIPSEFIELDS_HEADER

