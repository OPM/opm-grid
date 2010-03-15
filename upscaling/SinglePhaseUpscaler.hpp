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

#ifndef OPENRS_SINGLEPHASEUPSCALER_HEADER
#define OPENRS_SINGLEPHASEUPSCALER_HEADER

#define TEST_ANISO_RELPERM 0

#include "config.h"
#include <dune/common/param/ParameterGroup.hpp>
#include <dune/grid/CpGrid.hpp>
#include <dune/grid/common/EclipseGridParser.hpp>
#include <dune/solvers/common/GridInterfaceEuler.hpp>
#include <dune/solvers/common/BoundaryConditions.hpp>
#if TEST_ANISO_RELPERM
#include <dune/solvers/common/ReservoirPropertyCapillaryAnisotropicRelperm.hpp>
#include <dune/solvers/mimetic/MimeticIPAnisoRelpermEvaluator.hpp>
#else
#include <dune/solvers/common/ReservoirPropertyCapillary.hpp>
#include <dune/solvers/mimetic/MimeticIPEvaluator.hpp>
#endif
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
	typedef GridInterfaceEuler<GridType> GridInterface;
#if TEST_ANISO_RELPERM
 	typedef ReservoirPropertyCapillaryAnisotropicRelperm<Dimension> ResProp;
#else
 	typedef ReservoirPropertyCapillary<Dimension> ResProp;
#endif

	/// A type for the upscaled permeability.
	typedef ResProp::MutablePermTensor permtensor_t;

	enum BoundaryConditionType { Fixed = 0, Linear = 1, Periodic = 2 };

	// ------- Methods -------

	/// Default constructor.
	SinglePhaseUpscaler();

	/// Initializes the upscaler from parameters.
	void init(const parameter::ParameterGroup& param);

	/// Initializes the upscaler from given arguments.
	void init(const EclipseGridParser& parser,
                  BoundaryConditionType bctype,
                  double perm_threshold,
                  double z_tolerance = 0.0,
                  double residual_tolerance = 1e-8,
                  int linsolver_verbosity = 0,
                  int linsolver_type = 1,
                  bool twodim_hack = false);

	/// Access the grid.
	const GridType& grid() const;

        /// Set boundary condition type. This may not be used to swicth
        /// between Periodic and the other types, since the grid is
        /// modified for Periodic conditions.
        void setBoundaryConditionType(BoundaryConditionType type);

        /// Set the permeability of a cell directly. This will override
        /// the permeability that was read from the eclipse file.
        void setPermeability(const int cell_index, const permtensor_t& k);

	/// Does a single-phase upscaling.
	/// @return an upscaled permeability tensor.
	permtensor_t upscaleSinglePhase();

        /// Compute upscaled porosity.
        /// @return total pore volume of all cells divided by total volume.
        double upscalePorosity() const;

    protected:
	// ------- Typedefs and enums -------
	typedef GridInterface::CellIterator                CellIter;
	typedef CellIter::FaceIterator                     FaceIter;
	typedef BasicBoundaryConditions<true, true>             BCs;
#if TEST_ANISO_RELPERM
	typedef IncompFlowSolverHybrid<GridInterface,
				       ResProp,
				       BCs,
				       MimeticIPAnisoRelpermEvaluator> FlowSolver;
#else
	typedef IncompFlowSolverHybrid<GridInterface,
				       ResProp,
				       BCs,
				       MimeticIPEvaluator> FlowSolver;
#endif

	// ------- Methods -------
	template <class FlowSol>
	double computeAverageVelocity(const FlowSol& flow_solution,
				      const int flow_dir,
				      const int pdrop_dir) const;
	double computeDelta(const int flow_dir) const;

	virtual void initImpl(const parameter::ParameterGroup& param);
	virtual void initFinal(const parameter::ParameterGroup& param);

	// ------- Data members -------
	// FluxChecker flux_checker_;
	BoundaryConditionType bctype_;
	bool twodim_hack_;
	double residual_tolerance_;
	int linsolver_verbosity_;
        int linsolver_type_;

	GridType grid_;
	GridInterface ginterf_;
	ReservoirPropertyCapillary<3> res_prop_;
	BCs bcond_;
	FlowSolver flow_solver_;

    };

} // namespace Dune

#include "SinglePhaseUpscaler_impl.hpp"


#endif // OPENRS_SINGLEPHASEUPSCALER_HEADER
