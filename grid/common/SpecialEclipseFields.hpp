//===========================================================================
//
// File: SpecialEclipseFields.hpp
//
// Created: Mon Sep 21 14:09:54 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
//            Bjørn Spjelkavik    <bsp@sintef.no>
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

#ifndef OPENRS_SPECIALECLIPSEFIELDS_HEADER
#define OPENRS_SPECIALECLIPSEFIELDS_HEADER

#include <string>
#include <fstream>
#include <limits>
#include <dune/common/ErrorMacros.hpp>
#include "EclipseGridParserHelpers.hpp"

namespace Dune
{

// Abstract base class for special fields.
struct SpecialBase {
    virtual ~SpecialBase() {}                       // Default destructor
    //virtual std::string name() const = 0;           // Keyword name
    virtual void read(std::istream& is) = 0;        // Reads data
    //virtual void write(std::ostream& os) const = 0; // Writes data
    typedef std::vector<std::vector<std::vector<double> > > table_t;
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
	int nread = readDefaultedVectorData(is , data, ndim+1);
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
	
	if (ignoreSlashLine(is)) {
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
		is >> ignoreLine;
		break;
	    }
	    while (name.find("--") == 0) {
		// This line is a comment
		is >> ignoreLine >> name;
	    }
	    FaultSegment fault_segment;
	    fault_segment.ijk_coord.resize(6);
	    fault_segment.fault_name = name;
	    int nread = readDefaultedVectorData(is, fault_segment.ijk_coord, 6);
	    if (nread != 6) {
		THROW("Error reading fault_segment " << name);
	    }
	    is >> fault_segment.face;
	    faults.push_back(fault_segment);
	    ignoreSlashLine(is);
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
		is >> ignoreLine;
		break;
	    }
	    while (name == "--") {
		// This line is a comment
		is >> ignoreLine >> name;
	    }
	    MultfltLine multflt_line;
	    multflt_line.fault_name = name;
	    std::vector<double> data(2,1.0);
	    if (readDefaultedVectorData(is, data, 2) == 2) {
		ignoreSlashLine(is);
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
    { is >> ignoreLine; std::getline(is, title); }
    virtual void write(std::ostream& os) const
    { os << name() << '\n' << title << '\n'; }

};




struct START : public SpecialBase
{
    boost::gregorian::date date;
    virtual std::string name() const
    { return std::string("START"); }
    virtual void read(std::istream& is)
    { date = readDate(is); }
    virtual void write(std::ostream& os) const
    { os << name() << '\n' << date << '\n'; }
};




struct DATES : public SpecialBase
{
    std::vector<boost::gregorian::date> dates;
    virtual std::string name() const
    { return std::string("DATES"); }
    virtual void read(std::istream& is)
    {
	while(is) {
	    dates.push_back(readDate(is));
	    is >> ignoreWhitespace;
	    if (is.peek() == int('/')) {
		is >> ignoreLine;
		break;
	    }
	}
    }
    virtual void write(std::ostream& os) const
    {
	os << name() << '\n';
	copy(dates.begin(), dates.end(),
	     std::ostream_iterator<boost::gregorian::date>(os, "\n"));
    }

};


struct DENSITY : public SpecialBase
{
    std::vector<std::vector<double> > densities_;

    virtual std::string name() const {return std::string("DENSITY");}

    virtual void read(std::istream& is)
    {
	while (!is.eof()) {
	    std::vector<double> density;
	    readVectorData(is, density);
	    densities_.push_back(density);

	    int action = next_action(is); // 0:continue  1:return  2:throw
	    if (action == 1) {
		return;     // Alphabetic char. Read next keyword.
	    } else if (action == 2) {
		THROW("Error reading DENSITY. Next character is "
		      << (char)is.peek());
	    }
	}
    }

    virtual void write(std::ostream& os) const
    {
	os << name() << '\n';
	for (int i=0; i<(int)densities_.size(); ++i) {
	    os << densities_[i][0] << " " << densities_[i][1] << " "
	       << densities_[i][2] << '\n';
	}
	os << '\n';
    }
};

struct PVDG : public SpecialBase
{
    table_t pvdg_; 

    virtual std::string name() const {return std::string("PVDG");}

    virtual void read(std::istream& is)
    {
	readPvdTable(is, pvdg_, name(), 3);
    }

    virtual void write(std::ostream& os) const
    {
	os << name() << '\n';
	for (int prn=0; prn<(int)pvdg_.size(); ++prn) {
	    for (int i=0; i<(int)pvdg_[prn][0].size(); ++i) {
		os << pvdg_[prn][0][i] << " " << pvdg_[prn][1][i] << " "
		   << pvdg_[prn][2][i] << '\n';
	    }
	    os << '\n';
	}
	os << '\n';
    }
};

struct PVDO : public SpecialBase
{
    table_t pvdo_; 

    virtual std::string name() const {return std::string("PVDO");}

    virtual void read(std::istream& is)
    {
	readPvdTable(is, pvdo_, name(), 3);
    }

    virtual void write(std::ostream& os) const
    {
	os << name() << '\n';
	for (int prn=0; prn<(int)pvdo_.size(); ++prn) {
	    for (int i=0; i<(int)pvdo_[prn][0].size(); ++i) {
		os << pvdo_[prn][0][i] << " " << pvdo_[prn][1][i] << " "
		   << pvdo_[prn][2][i] << '\n';
	    }
	    os << '\n';
	}
	os << '\n';
    }
};

struct PVTG : public SpecialBase
{
    table_t pvtg_; 

    virtual std::string name() const {return std::string("PVTG");}

    virtual void read(std::istream& is)
    {
	readPvtTable(is, pvtg_, name());
    }

    virtual void write(std::ostream& os) const
    {
	os << name() << '\n';
	for (int prn=0; prn<(int)pvtg_.size(); ++prn) {
	    for (int i=0; i<(int)pvtg_[prn].size(); ++i) {
		int nl = (pvtg_[prn][i].size()-1) / 3;
		os << pvtg_[prn][i][0] << "   " << pvtg_[prn][i][1] << "  "
		   << pvtg_[prn][i][2] << "  " << pvtg_[prn][i][3] <<  '\n';
		for (int j=1, n=3; j<nl; ++j, n+=3) {
		    os << '\t' << pvtg_[prn][i][n+1] << "  "
		       << pvtg_[prn][i][n+2] << "  " << pvtg_[prn][i][n+3]
		       <<  '\n';
		}
	    }
	    os << '\n';
	}
	os << '\n';
    }
};


struct PVTO : public SpecialBase
{
    table_t pvto_; 

    virtual std::string name() const {return std::string("PVTO");}

    virtual void read(std::istream& is)
    {
	readPvtTable(is, pvto_, name());
    }

    virtual void write(std::ostream& os) const
    {
	os << name() << '\n';
	for (int prn=0; prn<(int)pvto_.size(); ++prn) {
	    for (int i=0; i<(int)pvto_[prn].size(); ++i) {
		int nl = (pvto_[prn][i].size()-1) / 3;
		os << pvto_[prn][i][0] << "   " << pvto_[prn][i][1] << "  "
		   << pvto_[prn][i][2] << "  " << pvto_[prn][i][3] <<  '\n';
		for (int j=1, n=3; j<nl; ++j, n+=3) {
		    os << '\t' << pvto_[prn][i][n+1] << "  "
		       << pvto_[prn][i][n+2] << "  " << pvto_[prn][i][n+3]
		       <<  '\n';
		}
	    }
	    os << '\n';
	}
	os << '\n';
    }
};


struct PVTW : public SpecialBase
{
    std::vector<std::vector<double> > pvtw_;

    virtual std::string name() const {return std::string("PVTW");}

    virtual void read(std::istream& is)
    {
	while (!is.eof()) {
	    std::vector<double> pvtw;
	    readVectorData(is, pvtw);
	    if (pvtw.size() == 4) {
		pvtw.push_back(0.0); // Not used by frontsim
	    }
	    pvtw_.push_back(pvtw);

	    int action = next_action(is); // 0:continue  1:return  2:throw
	    if (action == 1) {
		return;     // Alphabetic char. Read next keyword.
	    } else if (action == 2) {
		THROW("Error reading PVTW. Next character is "
		      <<  (char)is.peek());
	    }
	}
    }

    virtual void write(std::ostream& os) const
    {
	os << name() << '\n';
	for (int i=0; i<(int)pvtw_.size(); ++i) {
	    os << pvtw_[i][0] << " " << pvtw_[i][1] << " " << pvtw_[i][2]
	       << " " << pvtw_[i][3] << " " << pvtw_[i][4] << '\n';
	}
	os << '\n';
    }
};


struct ROCK : public SpecialBase
{
    std::vector<std::vector<double> > rock_compressibilities_;

    virtual std::string name() const {return std::string("ROCK");}

    virtual void read(std::istream& is)
    {
	while (!is.eof()) {
	    std::vector<double> rock;
	    readVectorData(is, rock);
	    rock_compressibilities_.push_back(rock);

	    int action = next_action(is); // 0:continue  1:return  2:throw
	    if (action == 1) {
		return;     // Alphabetic char. Read next keyword.
	    } else if (action == 2) {
		THROW("Error reading ROCK. Next character is "
		      << (char)is.peek());
	    }
	}
    }

    virtual void write(std::ostream& os) const
    {
	os << name() << '\n';
	for (int i=0; i<(int)rock_compressibilities_.size(); ++i) {
	    os << rock_compressibilities_[i][0] << " "
	       << rock_compressibilities_[i][1] << '\n';
	}
	os << '\n';
    }
};


struct ROCKTAB : public SpecialBase
{
    table_t rocktab_; 

    virtual std::string name() const {return std::string("ROCKTAB");}

    virtual void read(std::istream& is)
    {
	readPvdTable(is, rocktab_, name(), 3);
    }

    virtual void write(std::ostream& os) const
    {
	os << name() << '\n';
	for (int prn=0; prn<(int)rocktab_.size(); ++prn) {
	    for (int i=0; i<(int)rocktab_[prn][0].size(); ++i) {
		os << rocktab_[prn][0][i] << " " << rocktab_[prn][1][i] << " "
		   << rocktab_[prn][2][i] << '\n';
	    }
	    os << '\n';
	}
	os << '\n';
    }
};


struct SGOF : public SpecialBase
{
    table_t sgof_; 

    virtual std::string name() const {return std::string("SGOF");}

    virtual void read(std::istream& is) {readSGWOF(is, sgof_, name(), 4);}

    virtual void write(std::ostream& os) const
    {
	os << name() << '\n';
	for (int prn=0; prn<(int)sgof_.size(); ++prn) {
	    for (int i=0; i<(int)sgof_[prn][0].size(); ++i) {
		os << sgof_[prn][0][i] << " " << sgof_[prn][1][i] << " "
		   << sgof_[prn][2][i] << " " << sgof_[prn][3][i] << '\n';
	    }
	    os << '\n';
	}
	os << '\n';
    }
};

struct SWOF : public SpecialBase
{
    table_t swof_; 

    virtual std::string name() const {return std::string("SWOF");}

    virtual void read(std::istream& is) {readSGWOF(is, swof_, name(), 4);}

    virtual void write(std::ostream& os) const
    {
	os << name() << '\n';
	for (int prn=0; prn<(int)swof_.size(); ++prn) {
	    for (int i=0; i<(int)swof_[prn][0].size(); ++i) {
		os << swof_[prn][0][i] << " " << swof_[prn][1][i] << " "
		   << swof_[prn][2][i] << " " << swof_[prn][3][i] << '\n';
	    }
	    os << '\n';
	}
	os << '\n';
    }
};

struct MultRec : public SpecialBase
{
    virtual void read(std::istream& is)
    {
#ifdef VERBOSE
        std::cout << "(dummy implementation)" << std::endl;
#endif
	const std::ctype<char>& ct = std::use_facet< std::ctype<char> >(std::locale::classic());
	is >> ignoreSlashLine;
        while (!is.eof()) {
            is >> ignoreWhitespace;
	    std::streampos pos = is.tellg();
            char c;
            is.get(c);
	    if (is.eof()) {
		return;
	    }
	    if (ct.is(std::ctype_base::alpha, c)) {
		std::string name;     // Unquoted name or new keyword?
		std::getline(is, name);
		if (name.rfind('/') != std::string::npos) {
		    continue;  // Unquoted name
		} else {
		    is.seekg(pos);  
		    break;     // Read next keyword    
		}
	    } else if (ct.is(std::ctype_base::digit, c) || c== '.') {
		is >> ignoreSlashLine; // Decimal digit. Ignore data.
		continue;
	    } else if (c== '\'') {
		is >> ignoreSlashLine; // Quote. Ignore data.
		continue;
	    } else if(c == '-' && is.peek() == int('-')) {
		is >> ignoreLine;   // This line is a comment
		continue;
	    } else if (c == '/' ) {
		 is >> ignoreLine;  // This line is a null record.
		 continue;          // (No data before slash)
	    } else {
		is.putback(c);
		std::string temp;
		is >> temp;
		std::cout << "READ ERROR!  Next word is " << temp << std::endl;
	    }
        }
    }
};

// The following fields only have a dummy implementation
// that allows us to ignore them.
struct SWFN : public MultRec {};
struct SOF2 : public MultRec {};
struct EQUIL : public MultRec {};
struct WELSPECS : public MultRec {};
struct COMPDAT : public MultRec {};
struct WCONINJE : public MultRec {};
struct TUNING : public MultRec {};


} // End of namespace Dune

#endif // OPENRS_SPECIALECLIPSEFIELDS_HEADER

