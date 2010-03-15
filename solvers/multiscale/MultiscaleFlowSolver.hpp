//===========================================================================
//
// File: MultiscaleFlowSolver.hpp
//
// Created: Wed Sep  2 13:58:16 2009
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

#ifndef OPENRS_MULTISCALEFLOWSOLVER_HEADER
#define OPENRS_MULTISCALEFLOWSOLVER_HEADER

#include <dune/solvers/mimetic/IncompFlowSolverHybrid.hpp>
#include <dune/solvers/mimetic/MimeticIPEvaluator.hpp>

namespace Dune {

    template<class GridInterface,
             class ReservoirInterface,
             class BCInterface>
    class MultiscaleFlowSolver : public IncompFlowSolverHybrid<GridInterface, ReservoirInterface, BCInterface, MimeticIPEvaluator>
    {
    };

} // namespace Dune

#endif // OPENRS_MULTISCALEFLOWSOLVER_HEADER
