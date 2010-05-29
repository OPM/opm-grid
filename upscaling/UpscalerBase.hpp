//===========================================================================
//
// File: UpscalerBase.hpp
//
// Created: Thu Apr 29 10:20:22 2010
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
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

#ifndef OPENRS_UPSCALERBASE_HEADER
#define OPENRS_UPSCALERBASE_HEADER



#if HAVE_CONFIG_H
#include "config.h"
#endif
#include <dune/common/param/ParameterGroup.hpp>
#include <dune/grid/CpGrid.hpp>
#include <dune/grid/common/EclipseGridParser.hpp>
#include <dune/solvers/common/GridInterfaceEuler.hpp>
#include <dune/solvers/common/BoundaryConditions.hpp>


namespace Dune
{
    /**
       @brief A base class for upscaling.
       @author Atgeirr F. Rasmussen <atgeirr@sintef.no>
    */
    template <class Traits>
    class UpscalerBase
    {
    protected:
    public:
	// ------- Typedefs  -------
	typedef CpGrid GridType;
 	enum { Dimension = GridType::dimension };
	typedef GridInterfaceEuler<GridType> GridInterface;
        typedef typename Traits::template ResProp<Dimension>::Type ResProp;

	/// A type for the upscaled permeability.
	typedef typename ResProp::MutablePermTensor permtensor_t;

	enum BoundaryConditionType { Fixed = 0, Linear = 1, Periodic = 2 };

	// ------- Methods -------

	/// Default constructor.
	UpscalerBase();

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
        typedef typename Traits::template FlowSolver<GridInterface, BCs>::Type FlowSolver;

	// ------- Methods -------
	template <class FlowSol>
	double computeAverageVelocity(const FlowSol& flow_solution,
				      const int flow_dir,
				      const int pdrop_dir) const;

	double computeDelta(const int flow_dir) const;

        template <class FluidInterface>
        permtensor_t upscaleEffectivePerm(const FluidInterface& fluid);

	virtual void initImpl(const parameter::ParameterGroup& param);

	virtual void initFinal(const parameter::ParameterGroup& param);

	// ------- Data members -------
	BoundaryConditionType bctype_;
	bool twodim_hack_;
	double residual_tolerance_;
	int linsolver_verbosity_;
        int linsolver_type_;

	GridType grid_;
	GridInterface ginterf_;
	ResProp res_prop_;
	BCs bcond_;
	FlowSolver flow_solver_;
    };

} // namespace Dune

#include "UpscalerBase_impl.hpp"




#endif // OPENRS_UPSCALERBASE_HEADER
