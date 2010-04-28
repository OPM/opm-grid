//===========================================================================
//                                                                           
// File: EclipseGridParser.h                                                 
//                                                                           
// Created: Wed Dec  5 17:05:13 2007                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//
// $Date$
//                                                                           
// Revision: $Id: EclipseGridParser.h,v 1.3 2008/08/18 14:16:13 atgeirr Exp $
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

#ifndef SINTEF_ECLIPSEGRIDPARSER_HEADER
#define SINTEF_ECLIPSEGRIDPARSER_HEADER

#include <iosfwd>
#include <string>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "SpecialEclipseFields.hpp"
#include <dune/common/Factory.hpp>

namespace Dune
{

/**
   @brief A class for reading and parsing all fields of an eclipse file.
   
   This object is constructed using an Eclipse .grdecl-file. All data
   fields are extracted upon construction and written to vector data
   structures, which can then be read out in O(1) time afterwards via
   convenience functions.

   There is also a convenience function to easily check which fields
   were successfully parsed.

   @author Atgeirr F. Rasmussen <atgeirr@sintef.no>
   @date 2007/12/06 13:00:03

*/

class EclipseGridParser
{
public:
    /// Constructor taking an eclipse file as a stream.
    /// Note that since the INCLUDE keyword is relative to
    /// the directory of the current file, you may need to use
    /// the constructor taking a filename to process eclipse
    /// files with INCLUDEs.
    explicit EclipseGridParser(std::istream& is);
    /// Convenience constructor taking an eclipse filename.
    explicit EclipseGridParser(const std::string& filename);

    /// Read the given stream, overwriting any previous data.
    void read(std::istream& is);

    /// Returns true if the given keyword corresponds to a field that
    /// was found in the file.
    bool hasField(const std::string& keyword) const;
    /// Returns true if all the given keywords correspond to fields
    /// that were found in the file.
    bool hasFields(const std::vector<std::string>& keywords) const;
    /// The keywords/fields found in the file.
    std::vector<std::string> fieldNames() const;

    /// Returns a reference to a vector containing the values
    /// corresponding to the given integer keyword.
    const std::vector<int>& getIntegerValue(const std::string& keyword) const;

    /// Returns a reference to a vector containing the values
    /// corresponding to the given floating-point keyword.
    const std::vector<double>& getFloatingPointValue(const std::string& keyword) const;

    /// Returns a reference to a vector containing pointers to the values 
    /// corresponding to the given keyword when the values are not only integers
    /// or floats.
    const boost::shared_ptr<SpecialBase> getSpecialValue(const std::string& keyword) const;

    // This macro implements support for a special field keyword. It requires that a subclass
    // of SpecialBase exists, that has the same name as the keyword.
    // After using SPECIAL_FIELD(KEYWORD), the public member getKEYWORD will be available.
#define SPECIAL_FIELD(keyword)                                                                   \
private:                                                                                         \
    struct X##keyword { X##keyword() { Factory<SpecialBase>::addCreator<keyword>(#keyword); } }; \
    X##keyword x##keyword;                                                                       \
public:                                                                                          \
    const keyword& get##keyword() const                                                          \
    { return dynamic_cast<const keyword&>(*getSpecialValue(#keyword)); }

    // Support for special fields.
    SPECIAL_FIELD(SPECGRID);
    SPECIAL_FIELD(FAULTS);
    SPECIAL_FIELD(MULTFLT);
    SPECIAL_FIELD(TITLE);
    SPECIAL_FIELD(START);
    SPECIAL_FIELD(DATES);
    SPECIAL_FIELD(DENSITY);
    SPECIAL_FIELD(PVDG);
    SPECIAL_FIELD(PVDO);
    SPECIAL_FIELD(PVTG);
    SPECIAL_FIELD(PVTO);
    SPECIAL_FIELD(PVTW);
    SPECIAL_FIELD(SGOF);
    SPECIAL_FIELD(SWOF);
    SPECIAL_FIELD(ROCK);
    SPECIAL_FIELD(ROCKTAB);

    // The following fields only have a dummy implementation
    // that allows us to ignore them.
    SPECIAL_FIELD(SWFN);
    SPECIAL_FIELD(SOF2);
    SPECIAL_FIELD(EQUIL);
    SPECIAL_FIELD(WELSPECS);
    SPECIAL_FIELD(COMPDAT);
    SPECIAL_FIELD(WCONINJE);
    SPECIAL_FIELD(TUNING);

#undef SPECIAL_FIELD

private:
    boost::shared_ptr<SpecialBase> createSpecialField(std::istream& is, const std::string& fieldname);

    std::string directory_;
    std::map<std::string, std::vector<int> > integer_field_map_;
    std::map<std::string, std::vector<double> > floating_field_map_;
    std::map<std::string, boost::shared_ptr<SpecialBase> >special_field_map_;
    std::vector<int> empty_integer_field_;
    std::vector<double> empty_floating_field_;
    boost::shared_ptr<SpecialBase> empty_special_field_;
};


} // namespace Dune

#endif // SINTEF_ECLIPSEGRIDPARSER_HEADER
