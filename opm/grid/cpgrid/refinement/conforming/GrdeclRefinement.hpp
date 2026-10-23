/*
  Copyright 2026 SINTEF Digital.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef OPM_GRID_GRDECL_REFINEMENT_HEADER_INCLUDED
#define OPM_GRID_GRDECL_REFINEMENT_HEADER_INCLUDED

#include <opm/grid/cpgrid/refinement/RefinementRequest.hpp>

#include <array>
#include <vector>

namespace Opm
{
namespace Refinement
{

/// Corner-point (grdecl) description of one refined block, produced by
/// refineBlock(). Dimensions are box dims times subdivision factors; coord
/// and zcorn follow the usual corner-point layouts for those dimensions.
struct RefinedBlockGrdecl
{
    std::array<int,3> dims{};
    std::vector<double> coord;
    std::vector<double> zcorn;
    std::vector<int> actnum;
};

/// Stage 2+3 of the build pipeline (docs/DESIGN-builder.md D4): resample a
/// block of a corner-point description onto refined sub-pillars.
///
/// Semantics are corner-point-native (what deck LGRs mean):
///  - Sub-pillars are straight lines whose endpoints interpolate the four
///    surrounding parent pillar endpoints bilinearly; refined pillars on
///    parent pillar positions reproduce the parent pillars exactly.
///  - Refined zcorn values interpolate the parent cell's eight zcorn values
///    trilinearly in the lateral/vertical fractions.
///
/// Faults and pinch-outs *inside* the block survive by construction: zcorn
/// jumps across columns and collapsed cells resample to jumps/collapses on
/// the refined lattice; matching them into faces is the (later) processing
/// stage's job. Note that on inclined pillars this differs slightly from
/// refining each hexahedron by its trilinear map (the old CpGrid approach);
/// on vertical pillars the two coincide exactly.
///
/// @param parentDims  Cartesian dimensions of the parent description.
/// @param coord       Parent COORD, 6*(nx+1)*(ny+1) doubles.
/// @param zcorn       Parent ZCORN, 8*nx*ny*nz doubles.
/// @param actnum      Parent ACTNUM (nx*ny*nz ints) or nullptr (all active).
/// @param request     The block to refine (startIJK/endIJK within parentDims).
RefinedBlockGrdecl refineBlock(const std::array<int,3>& parentDims,
                               const double* coord,
                               const double* zcorn,
                               const int* actnum,
                               const BlockRefinement& request);

} // namespace Refinement
} // namespace Opm

#endif // OPM_GRID_GRDECL_REFINEMENT_HEADER_INCLUDED
