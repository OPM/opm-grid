//===========================================================================
//
// File: EulerUpstreamResidual.hpp
//
// Created: Thu May  6 11:14:23 2010
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

#ifndef OPENRS_EULERUPSTREAMRESIDUAL_HEADER
#define OPENRS_EULERUPSTREAMRESIDUAL_HEADER



#include <tr1/unordered_map>

#include <dune/common/param/ParameterGroup.hpp>
#include <dune/common/SparseVector.hpp>


namespace Dune {

    namespace EulerUpstreamResidualDetails {
        // Forward declaration for friendship purposes.
        template <class UpstreamSolver, class PressureSolution>
        struct UpdateForCell;
    }


    /// Class for doing simple transport by explicit Euler upstream method for general grid.
    /// @tparam
    template <class GridInterface, class ReservoirProperties, class BoundaryConditions>
    class EulerUpstreamResidual
    {
    public:
        template <class S, class P>
        friend class EulerUpstreamResidualDetails::UpdateForCell;
	typedef typename GridInterface::CellIterator CIt;
	typedef typename CIt::FaceIterator FIt;
	typedef typename FIt::Vector Vector;
        typedef ReservoirProperties RP;

	/// @brief
	/// @todo Doc me
	EulerUpstreamResidual();
	/// @brief
	/// @todo Doc me
	/// @param
 	EulerUpstreamResidual(const GridInterface& grid,
                              const ReservoirProperties& resprop,
                              const BoundaryConditions& boundary);

	void initObj(const GridInterface& grid,
		     const ReservoirProperties& resprop,
		     const BoundaryConditions& boundary);


	template <class FlowSolution>
	void computeSatDelta(const std::vector<double>& saturation,
			     const typename GridInterface::Vector& gravity,
			     const FlowSolution& flow_sol,
                             const SparseVector<double>& injection_rates,
                             const bool method_viscous,
                             const bool method_gravity,
                             const bool method_capillary,
                             std::vector<double>& sat_delta) const;

	void computeCapPressures(const std::vector<double>& saturation) const;

	typename GridInterface::Vector
	estimateCapPressureGradient(const FIt& f, const FIt& nbf, const std::vector<double>& saturation) const;

        const GridInterface& grid() const;
        const ReservoirProperties& reservoirProperties() const;

    private:
	void initFinal();

	const GridInterface* pgrid_;
	const ReservoirProperties* preservoir_properties_;
	const BoundaryConditions* pboundary_;

	// Boundary id to face iterator mapping. May be mostly or completely empty.
	// Obviously requires unique-face-per-bid grids.
	std::vector<FIt> bid_to_face_;

        // Storing some cell iterators, so that we may use tbb for parallelizing.
        std::vector<CIt> cell_iters_;

	// Precomputing the capillary pressures of cells saves a little time.
	mutable std::vector<double> cap_pressures_;
        mutable const SparseVector<double>* pinjection_rates_;
        mutable bool method_viscous_;
        mutable bool method_gravity_;
        mutable bool method_capillary_;
    };

} // namespace Dune

#include "EulerUpstreamResidual_impl.hpp"


#endif // OPENRS_EULERUPSTREAMRESIDUAL_HEADER
