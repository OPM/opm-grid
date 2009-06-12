//===========================================================================
//                                                                           
// File: EclipseGridInspector.h                                              
//                                                                           
// Created: Mon Jun  2 09:46:08 2008                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: EclipseGridInspector.h,v 1.2 2008/08/18 14:16:12 atgeirr Exp $
//                                                                           
//===========================================================================

#ifndef SINTEF_ECLIPSEGRIDINSPECTOR_HEADER
#define SINTEF_ECLIPSEGRIDINSPECTOR_HEADER

#include "EclipseGridParser.hpp"
#include <boost/array.hpp>

/**
   @brief A class for inspecting the contents of an eclipse file.

   Given an EclipseGridParser that has successfully read en Eclipse
   .grdecl-file, this class may be used to answer certain queries about
   its contents.

   @author Atgeirr F. Rasmussen <atgeirr@sintef.no>
   @date 2008/06/02 09:46:08
*/

namespace Dune
{

class EclipseGridInspector
{
public:
    /// Constructor taking a parser as argument.
    /// The parser must already have read an Eclipse file.
    EclipseGridInspector(const EclipseGridParser& parser);

    /// Assuming that the pillars are vertical, compute the
    /// volume of the cell given by logical coordinates (i, j, k).
    double cellVolumeVerticalPillars(int i, int j, int k) const;

    /// Assuming that the pillars are vertical, compute the
    /// volume of the cell given by the cell index
    double cellVolumeVerticalPillars(int cell_idx) const;
    
    /// Returns a vector with the outer limits of grid (in the grid's unit).
    /// The vector contains [xmin, xmax, ymin, ymax, zmin, zmax], as 
    /// read from COORDS and ZCORN
    std::vector<double> getGridLimits() const;

    /// Returns the extent of the logical cartesian grid
    /// as number of cells in the (i, j, k) directions.
    boost::array<int, 3> gridSize() const;

    /// Returns the eight z-values associated with a given cell.
    /// The ordering is such that i runs fastest. That is, with
    /// L = low and H = high:
    /// {LLL, HLL, LHL, HHL, LLH, HLH, LHH, HHH }.
    boost::array<double, 8> cellZvals(int i, int j, int k) const;

private:
    const EclipseGridParser& parser_;
    int logical_gridsize_[3];
    void checkLogicalCoords(int i, int j, int k) const;
};

} // namespace Dune

#endif // SINTEF_ECLIPSEGRIDINSPECTOR_HEADER

