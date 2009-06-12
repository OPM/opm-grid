//===========================================================================
//                                                                           
// File: EclipseGridParser.h                                                 
//                                                                           
// Created: Wed Dec  5 17:05:13 2007                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: EclipseGridParser.h,v 1.3 2008/08/18 14:16:13 atgeirr Exp $
//                                                                           
//===========================================================================

#ifndef SINTEF_ECLIPSEGRIDPARSER_HEADER
#define SINTEF_ECLIPSEGRIDPARSER_HEADER

#include <iosfwd>
#include <string>
#include <vector>
#include <map>

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

namespace Dune
{

class EclipseGridParser
{
public:
    /// Constructor taking an eclipse file as a stream.
    explicit EclipseGridParser(std::istream& is);
    /// Convenience constructor taking an eclipse filename.
    explicit EclipseGridParser(const std::string& filename);

    /// Read the given stream, overwriting any previous data.
    void read(std::istream& is);

    /// Returns true is the given keyword corresponds to a field that
    /// was found in the file.
    bool hasField(const std::string& keyword) const;
    /// Returns true is all the given keywords correspond to fields
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


private:
    std::map<std::string, std::vector<int> > integer_field_map_;
    std::map<std::string, std::vector<double> > floating_field_map_;
    std::vector<int> empty_integer_field_;
    std::vector<double> empty_floating_field_;
};


} // namespace Dune

#endif // SINTEF_ECLIPSEGRIDPARSER_HEADER
