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
#ifndef OPM_GRID_REFINEMENT_EDGECONFORMAL_HEADER_INCLUDED
#define OPM_GRID_REFINEMENT_EDGECONFORMAL_HEADER_INCLUDED

namespace Dune
{
namespace cpgrid
{
class CpGridData;
}
}

namespace Opm
{
namespace Refinement
{

/// Make a (refined) leaf grid edge-conformal: insert every leaf node that
/// lies on the interior of a face's edge into that face's node list.
///
/// Local grid refinement introduces subdivision nodes on the pillars and
/// lateral edges at the LGR boundary. Coarse cells that touch those edges
/// (but are not the box's face-neighbour) carry the original 4-node faces,
/// leaving those subdivision nodes "hanging" on the edges — fine for TPFA,
/// but not edge-conformal (needed by VEM). This pass closes that gap.
///
/// Only `face_to_point_` changes: no new corners are created (existing leaf
/// nodes are inserted into face node lists), so cell geometry, face area and
/// normal are unchanged. Cheap when there is nothing to insert.
void edgeConformalizeLeaf(Dune::cpgrid::CpGridData& leaf);

} // namespace Refinement
} // namespace Opm

#endif // OPM_GRID_REFINEMENT_EDGECONFORMAL_HEADER_INCLUDED
