//===========================================================================
//
// File: LookUpData.hpp
//
// Created: Tue May 25 11:45:00 2023
//
// Author(s): Antonella Ritorto <antonella.ritorto@opm-op.com>
//
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2023 Equinor ASA.

  This file is part of The Open Porous Media project  (OPM).

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

namespace Dune
{
template <typename GridType>
class LookUpData
{
};
/// Specialization for CpGrid
template<>
class LookUpData<Dune::CpGrid>
{
public:
    // Constructor taking a CpGrid object
    LookUpData(const Dune::CpGrid&){
    }

    template<typename feature_type>
    int operator()(const Dune::cpgrid::Entity<0>& elem, const std::vector<feature_type>& feature_vec)
    {
        // elem.getOrigin() Get entity in level 0 (either parent cell, or equivalent cell, or 'itself' if grid_ = level 0)
        return feature_vec[elem.getOrigin().index()];
    }

}; // end LookUpData<CpGrid> class
}
// end namespace Dune

