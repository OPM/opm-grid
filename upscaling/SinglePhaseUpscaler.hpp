//===========================================================================
//
// File: SinglePhaseUpscaler.hpp
//
// Created: Fri Aug 28 13:38:54 2009
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

#ifndef OPENRS_SINGLEPHASEUPSCALER_HEADER
#define OPENRS_SINGLEPHASEUPSCALER_HEADER

#include "config.h"
#include <dune/common/param/ParameterGroup.hpp>
#include <dune/grid/CpGrid.hpp>
#include <dune/solvers/common/GridInterfaceEuler.hpp>
#include <dune/solvers/common/ReservoirPropertyCapillary.hpp>
#include <dune/solvers/common/BoundaryConditions.hpp>
#include <dune/solvers/mimetic/MimeticIPEvaluator.hpp>
#include <dune/solvers/mimetic/IncompFlowSolverHybrid.hpp>

namespace Dune
{
    /**
       @brief A class for doing single phase (permeability) upscaling.
       @author Atgeirr F. Rasmussen <atgeirr@sintef.no>
    */
    class SinglePhaseUpscaler
    {
    protected:
    public:
	// ------- Typedefs  -------
	typedef CpGrid GridType;
 	enum { Dimension = GridType::dimension };
 	typedef ReservoirPropertyCapillary<Dimension> ResProp;
	typedef GridInterfaceEuler<GridType> GridInterface;

	/// A type for the upscaled permeability.
	typedef ResProp::MutablePermTensor permtensor_t;

	// ------- Methods -------

	/// Default constructor.
	SinglePhaseUpscaler();

	/// Initializes the upscaler.
	virtual void init(const parameter::ParameterGroup& param);

	/// Access the grid.
	const GridInterface& grid() const;

	/// Does a single-phase upscaling.
	/// @return an upscaled permeability tensor.
	permtensor_t upscaleSinglePhase();

    protected:
	// ------- Typedefs and enums -------
	enum BoundaryConditionType { Fixed = 0, Linear = 1, Periodic = 2 };
	typedef GridInterface::CellIterator                    CellIter;
	typedef CellIter::FaceIterator                         FaceIter;
	typedef MimeticIPEvaluator<CellIter, Dimension, true>  InnerProd;
	typedef BoundaryConditions<true, true>                 BCs;
	typedef IncompFlowSolverHybrid<GridInterface,
				       ResProp,
				       BCs,
				       InnerProd>              FlowSolver;
	// ------- Methods -------
	template <class FlowSol>
	double computeAverageVelocity(const FlowSol& flow_solution,
				      const int flow_dir,
				      const int pdrop_dir) const;
	double computeDelta(const int flow_dir) const;



	// ------- Data members -------
	// FluxChecker flux_checker_;
	BoundaryConditionType bctype_;
	bool twodim_hack_;
	double residual_tolerance_;

	GridType grid_;
	GridInterface ginterf_;
	ReservoirPropertyCapillary<3> res_prop_;
	BCs bcond_;
	FlowSolver flow_solver_;

    };

} // namespace Dune

#include "SinglePhaseUpscaler_impl.hpp"


#endif // OPENRS_SINGLEPHASEUPSCALER_HEADER
