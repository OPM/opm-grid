//===========================================================================
//
// File: EulerSolverTester.hpp
//
// Created: Thu Aug  6 10:04:03 2009
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

#ifndef OPENRS_EULERSOLVERTESTER_HEADER
#define OPENRS_EULERSOLVERTESTER_HEADER

#include <fstream>
#include <iterator>
#include <boost/lexical_cast.hpp>

#include <dune/common/param/ParameterGroup.hpp>
#include <dune/common/SparseVector.hpp>
#include <dune/common/SparseTable.hpp>

#include <dune/grid/common/Volumes.hpp>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/CpGrid.hpp>
#include <dune/grid/yaspgrid.hh>

#include <dune/solvers/common/BoundaryConditions.hpp>
#include <dune/solvers/common/GridInterfaceEuler.hpp>
#include <dune/solvers/common/ReservoirPropertyCapillary.hpp>
#include <dune/solvers/common/setupGridAndProps.hpp>

#include "../EulerUpstream.hpp"

namespace Dune
{



    /// @brief
    /// @todo Doc me!
    /// @tparam
    template <class GridInterface>
    class TestSolution
    {
    public:
	/// @brief
	/// @todo Doc me!
	/// @param
	TestSolution(const GridInterface& g,
		     const typename GridInterface::CellIterator::Vector& v)
	{
	    // Make a flux consistent with a constant velocity field.
	    typename GridInterface::CellIterator c = g.cellbegin();
	    std::vector<double> cell_fluxes;
	    for (; c != g.cellend(); ++c) {
		cell_fluxes.clear();
		typename GridInterface::CellIterator::FaceIterator f = c->facebegin();
		for (; f != c->faceend(); ++f) {
		    //double flux = Dune::template inner<typename GridInterface::CellIterator::Vector>(v, f->normal())*f->area();
		    double flux = inner(v, f->normal())*f->area();
		    cell_fluxes.push_back(flux);
		}
		halfface_fluxes_.appendRow(cell_fluxes.begin(), cell_fluxes.end());
	    }
	}

	/// @brief
	/// @todo Doc me!
	typedef typename GridInterface::CellIterator::FaceIterator FaceIter;
	/// @brief
	/// @todo Doc me!
	/// @param
	/// @return
	double outflux(const FaceIter& f) const
	{
	    return halfface_fluxes_[f->cellIndex()][f->localIndex()];
	}

    private:
	SparseTable<double> halfface_fluxes_;
    };




    /// @brief
    /// @todo Doc me!
    class EulerSolverTester
    {
    public:
	EulerSolverTester()
	    : simulation_steps_(1),
	      stepsize_(365),
	      init_saturation_(0.0)
	{
	}

	/// @brief
	/// @todo Doc me!
	/// @param
	void init(const parameter::ParameterGroup& param)
	{
	    simulation_steps_ = param.getDefault("simulation_steps", simulation_steps_);
	    stepsize_ = Dune::unit::convert::from(param.getDefault("stepsize", stepsize_),
                                                  Dune::unit::day);
	    init_saturation_ = param.getDefault("init_saturation", init_saturation_);

	    setupGridAndProps(param, grid_, res_prop_);

	    // Make flow equation boundary conditions.
	    // Pressure 1.0e5 on the left, 0.0 on the right.
	    // Recall that the boundary ids range from 1 to 6 for the cartesian edges,
	    // and that boundary id 0 means interiour face/intersection.
// 	    FlowBoundaryConditions flow_bcond(7);
// 	    flow_bcond[1] = BC(BC::Dirichlet, 1.0e5);
// 	    flow_bcond[2] = BC(BC::Dirichlet, 0.0);
	    // Make transport equation boundary conditions.
	    // The default one is fine (sat = 1.0 on inflow).
	    sat_bcond_.resize(7); // 7 since 0 is for interiour faces.
	}

	/// @brief
	/// @todo Doc me!
	void run()
	{
	    // Make the grid interface
	    GridInterface g(grid_);
	    // No injection or production.
	    SparseVector<double> injection_rates(g.numberOfCells());
	    // Make a solver.
	    TransportSolver transport_solver(g, res_prop_, sat_bcond_);
	    // Define a flow field with constant velocity.
	    FieldVector<double, 3> vel(0.0);
	    vel[0] = 1.0;
	    // vel[1] = 1.0;
	    TestSolution<GridInterface> flow_solution(g, vel);
	    // Initial saturation.
	    std::vector<double> sat(g.numberOfCells(), 0.0);
	    // Gravity.
	    FieldVector<double, 3> gravity(0.0);
	    // gravity[2] = -9.81;

	    // Solve some steps.
	    for (int i = 0; i < simulation_steps_; ++i) {
		transport_solver.transportSolve(sat, stepsize_, gravity, flow_solution, injection_rates);
		output("testsolution-" + boost::lexical_cast<std::string>(i), "saturation", sat);
	    }
	}


	/// @brief
	/// @todo Doc me!
	/// @tparam
	/// @param
	template <class CellData>
	void output(const std::string& filename, const std::string& fieldname, const CellData& celldata)
	{
	    // VTK output.
	    Dune::VTKWriter<typename GridType::LeafGridView> vtkwriter(grid_.leafView());
	    vtkwriter.addCellData(celldata, fieldname);
	    vtkwriter.write(filename, Dune::VTKOptions::ascii);
	    // Dumping the saturation.
// 	    std::ofstream os((filename + "_" + fieldname).c_str());
// 	    std::copy(celldata.begin(), celldata.end(), std::ostream_iterator<double>(os, "\n"));
	}




    private:
	typedef CpGrid GridType;
	typedef GridInterfaceEuler<GridType> GridInterface;
	typedef FlowBC BC;
	typedef EulerUpstream<GridInterface, ReservoirPropertyCapillary<3>, SaturationBoundaryConditions> TransportSolver;
	GridType grid_;
	ReservoirPropertyCapillary<3> res_prop_;
	SaturationBoundaryConditions sat_bcond_;
	int simulation_steps_;
	double stepsize_;
	double init_saturation_;
    };



} // namespace Dune

#endif // OPENRS_EULERSOLVERTESTER_HEADER
