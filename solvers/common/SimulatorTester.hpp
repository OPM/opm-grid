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

#ifndef OPENRS_SIMULATORTESTER_HEADER
#define OPENRS_SIMULATORTESTER_HEADER


#include <dune/solvers/common/SimulatorBase.hpp>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

namespace Dune
{



    /// @brief
    /// @todo Doc me!
    template <class SimTraits>
    class SimulatorTester : public SimulatorBase<SimTraits>
    {
    public:
	typedef SimulatorBase<SimTraits> Super;
	/// @brief
	/// @todo Doc me!
	void run()
	{
	    // Initial saturation.
	    std::vector<double> saturation(this->init_saturation_);
	    // Gravity.
	    // FieldVector<double, 3> gravity(0.0);
	    // gravity[2] = -Dune::unit::gravity;
	    // Compute flow field.
	    if (this->gravity_.two_norm() > 0.0) {
		MESSAGE("Warning: Gravity not handled by flow solver.");
	    }

	    // Solve some steps.
	    for (int i = 0; i < this->simulation_steps_; ++i) {
		std::cout << "\n\n================    Simulation step number " << i
                          << "    ===============" << std::endl;
		// Flow.
		this->flow_solver_.solve(this->res_prop_, saturation, this->bcond_, this->injection_rates_psolver_,
                                         this->residual_tolerance_, this->linsolver_verbosity_, this->linsolver_type_);
// 		if (i == 0) {
// 		    flow_solver_.printSystem("linsys_dump_mimetic");
// 		}
		// Transport.
		this->transport_solver_.transportSolve(saturation, this->stepsize_, this->gravity_,
							this->flow_solver_.getSolution(),
							this->injection_rates_);
		// Output.
                writeVtkOutput(this->ginterf_,
                               this->res_prop_,
                               this->flow_solver_.getSolution(),
                               saturation,
                               "testsolution-" + boost::lexical_cast<std::string>(i));
	    }
	}

    };



} // namespace Dune


#endif // OPENRS_SIMULATORTESTER_HEADER
