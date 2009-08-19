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

namespace Dune
{



    /// @brief
    /// @todo Doc me!
    class SimulatorTester : public SimulatorBase
    {
    public:
	/// @brief
	/// @todo Doc me!
	void run()
	{
	    // No injection or production.
	    SparseVector<double> injection_rates(ginterf_.numberOfCells());
	    std::vector<double> src(ginterf_.numberOfCells());
	    // Initial saturation.
	    std::vector<double> sat(ginterf_.numberOfCells(), init_saturation_);
	    // Gravity.
	    FieldVector<double, 3> gravity(0.0);
	    // gravity[2] = -Dune::unit::gravity;
	    // Compute flow field.
	    if (gravity.two_norm() > 0.0) {
		MESSAGE("Warning: Gravity not handled by flow solver.");
	    }

	    // Solve some steps.
	    for (int i = 0; i < simulation_steps_; ++i) {
		std::cout << "================    Simulation step number " << i
                          << "    ===============" << std::endl;
		// Flow.
		flow_solver_.solve(res_prop_, sat, flow_bcond_, src, gravity);
// 		if (i == 0) {
// 		    flow_solver_.printSystem("linsys_dump_mimetic");
// 		}
		// Transport.
		transport_solver_.transportSolve(sat, stepsize_, gravity,
						 flow_solver_.getSolution(),
						 injection_rates);
		// Output.
		std::vector<double> cell_velocity;
		estimateCellVelocity(cell_velocity, flow_solver_.getSolution());
		std::vector<double> cell_pressure;
		getCellPressure(cell_pressure, flow_solver_.getSolution());
		Dune::VTKWriter<GridType::LeafGridView> vtkwriter(grid_.leafView());
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
