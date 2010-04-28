//===========================================================================
//
// File: UpscalingTraits.hpp
//
// Created: Wed Apr 28 10:36:42 2010
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Jostein R Natvig    <jostein.r.natvig@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2010 Statoil ASA.

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

#ifndef OPENRS_UPSCALINGTRAITS_HEADER
#define OPENRS_UPSCALINGTRAITS_HEADER

#include <dune/solvers/common/ReservoirPropertyCapillaryAnisotropicRelperm.hpp>
#include <dune/solvers/mimetic/MimeticIPAnisoRelpermEvaluator.hpp>
#include <dune/solvers/common/ReservoirPropertyCapillary.hpp>
#include <dune/solvers/mimetic/MimeticIPEvaluator.hpp>
#include <dune/solvers/mimetic/IncompFlowSolverHybrid.hpp>

namespace Dune
{


    /// Traits for upscaling with isotropic relperm (scalar) input.
    struct UpscalingTraitsBasic
    {

        /// The reservoir property class.
        template <int Dimension>
        struct ResProp
        {
            typedef ReservoirPropertyCapillary<Dimension> Type;
        };

        /// The pressure/flow solver class.
        template <class GridInterface, class BoundaryConditions>
        struct FlowSolver
        {
            typedef IncompFlowSolverHybrid<GridInterface,
                                           typename ResProp<GridInterface::Dimension>::Type,
                                           BoundaryConditions,
                                           MimeticIPEvaluator> Type;
        };
    };


    /// Traits for upscaling with anisotropic relperm (tensorial) input.
    struct UpscalingTraitsAnisoRelperm
    {
        /// The reservoir property class.
        template <int Dimension>
        struct ResProp
        {
            typedef ReservoirPropertyCapillaryAnisotropicRelperm<Dimension> Type;
        };

        /// The pressure/flow solver class.
        template <class GridInterface, class BoundaryConditions>
        struct FlowSolver
        {
            typedef IncompFlowSolverHybrid<GridInterface,
                                           typename ResProp<GridInterface::Dimension>::Type,
                                           BoundaryConditions,
                                           MimeticIPAnisoRelpermEvaluator> Type;
        };
    };


} // namespace Dune


#endif // OPENRS_UPSCALINGTRAITS_HEADER
