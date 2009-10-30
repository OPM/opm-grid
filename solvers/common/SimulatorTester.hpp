//===========================================================================
//
// File: SimulatorTester.hpp
//
// Created: Fri Aug  7 09:21:55 2009
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

#ifndef OPENRS_SIMULATORTESTER_HEADER
#define OPENRS_SIMULATORTESTER_HEADER


#include <dune/solvers/common/SimulatorBase.hpp>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

namespace Dune
{



    /// @brief
    /// @todo Doc me!
    template <template <int> class ResPropT = ReservoirPropertyCapillary,
	      template <class, int, bool> class InnerProd = MimeticIPEvaluator>
    class SimulatorTester : public SimulatorBase<ResPropT, InnerProd>
    {
    public:
	typedef SimulatorBase<ResPropT, InnerProd> Super;
	/// @brief
	/// @todo Doc me!
	void run()
	{
	    // No injection or production.
	    SparseVector<double> injection_rates(Super::ginterf_.numberOfCells());
	    std::vector<double> src(Super::ginterf_.numberOfCells());
	    // Initial saturation.
	    std::vector<double> sat(Super::ginterf_.numberOfCells(), Super::init_saturation_);
	    // Gravity.
	    // FieldVector<double, 3> gravity(0.0);
	    // gravity[2] = -Dune::unit::gravity;
	    // Compute flow field.
	    if (Super::gravity_.two_norm() > 0.0) {
		MESSAGE("Warning: Gravity not handled by flow solver.");
	    }

	    // Solve some steps.
	    for (int i = 0; i < Super::simulation_steps_; ++i) {
		std::cout << "\n\n================    Simulation step number " << i
                          << "    ===============" << std::endl;
		// Flow.
		Super::flow_solver_.solve(Super::res_prop_, sat, Super::bcond_, src);
// 		if (i == 0) {
// 		    flow_solver_.printSystem("linsys_dump_mimetic");
// 		}
		// Transport.
		Super::transport_solver_.transportSolve(sat, Super::stepsize_, Super::gravity_,
							Super::flow_solver_.getSolution(),
							injection_rates);
		// Output.
		std::vector<double> cell_velocity;
		estimateCellVelocity(cell_velocity, Super::ginterf_, Super::flow_solver_.getSolution());
		std::vector<double> cell_pressure;
		getCellPressure(cell_pressure, Super::ginterf_, Super::flow_solver_.getSolution());
		Dune::VTKWriter<typename Super::GridType::LeafGridView> vtkwriter(Super::grid_.leafView());
		vtkwriter.addCellData(cell_velocity, "velocity");
		vtkwriter.addCellData(sat, "saturation");
		vtkwriter.addCellData(cell_pressure, "pressure");
		vtkwriter.write("testsolution-" + boost::lexical_cast<std::string>(i),
                                Dune::VTKOptions::ascii);
	    }
	}

    };



} // namespace Dune


#endif // OPENRS_SIMULATORTESTER_HEADER
