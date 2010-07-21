//===========================================================================
//
// File: SimulatorUtilities.hpp
//
// Created: Fri Aug 28 15:00:15 2009
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

#ifndef OPENRS_SIMULATORUTILITIES_HEADER
#define OPENRS_SIMULATORUTILITIES_HEADER


#include <dune/common/ErrorMacros.hpp>
#include <dune/common/fvector.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iterator>

namespace Dune
{


    /// @brief Estimates a scalar cell velocity from outgoing fluxes.
    /// @tparam GridInterface a grid interface.
    /// @tparam FlowSol a flow solution type.
    /// @param[out] cell_velocity the estimated velocities.
    /// @param[in] ginterf an interface to the grid.
    /// @param[in] flow_solution the object containing the fluxes.
    template <class GridInterface, class FlowSol>
    void estimateCellVelocity(std::vector<typename GridInterface::Vector>& cell_velocity,
			      const GridInterface& ginterf,
			      const FlowSol& flow_solution)
    {
	// Algorithm used is same as in halfFaceFluxToCellVelocity.hpp
	// in the Sintef legacy c++ code.
	cell_velocity.clear();
	cell_velocity.resize(ginterf.numberOfCells());
	for (typename GridInterface::CellIterator c = ginterf.cellbegin(); c != ginterf.cellend(); ++c) {
	    int numf = 0;
	    typename GridInterface::Vector cell_v(0.0);
	    typename GridInterface::CellIterator::FaceIterator f = c->facebegin();
	    for (; f != c->faceend(); ++f, ++numf) {
		double flux = flow_solution.outflux(f);
		typename GridInterface::Vector v = f->centroid();
		v -= c->centroid();
		v *= flux/c->volume();
		cell_v += v;
	    }
	    cell_velocity[c->index()] = cell_v;//.two_norm();
	}
    }


    /// @brief Estimates a scalar cell velocity from outgoing fluxes.
    /// @tparam GridInterface a grid interface.
    /// @tparam FlowSol a flow solution type.
    /// @param[out] cell_velocity the estimated velocities.
    /// @param[in] ginterf an interface to the grid.
    /// @param[in] flow_solution the object containing the fluxes.
    /// @param[in] partition partition numbers of the fluxes.
    /// @param[in] my_partition partition to be used.
    template <class GridInterface, class FlowSol>
    void estimateCellVelocity(std::vector<typename GridInterface::Vector>& cell_velocity,
			      const GridInterface& ginterf,
			      const FlowSol& flow_solution,
			      const std::vector<int>& partition,
			      const int my_partition)
    {
	// Algorithm used is same as in halfFaceFluxToCellVelocity.hpp
	// in the Sintef legacy c++ code.
	cell_velocity.clear();
	cell_velocity.resize(ginterf.numberOfCells());
	for (typename GridInterface::CellIterator c = ginterf.cellbegin(); c != ginterf.cellend(); ++c) {
	    if (partition[c->index()] != my_partition) {
		cell_velocity[c->index()] = 0.0;
	    } else {
		int numf = 0;
		typename GridInterface::Vector cell_v(0.0);
		typename GridInterface::CellIterator::FaceIterator f = c->facebegin();
		for (; f != c->faceend(); ++f, ++numf) {
		    double flux = flow_solution.outflux(f);
		    typename GridInterface::Vector v = f->centroid();
		    v -= c->centroid();
		    v *= flux/c->volume();
		    cell_v += v;
		}
		cell_velocity[c->index()] = cell_v;//.two_norm();
	    }
	}
    }


    template <class ReservoirProperty>
    void computePhaseVelocities(std::vector<FieldVector<double, 3> >& phase_velocity_water,
                                std::vector<FieldVector<double, 3> >& phase_velocity_oil,
                                const ReservoirProperty& res_prop,
                                const std::vector<double>& saturation,
                                const std::vector<FieldVector<double, 3> >& cell_velocity)
    {
        ASSERT(saturation.size() == cell_velocity.size());
        int num_cells = saturation.size();
        phase_velocity_water = cell_velocity;
        phase_velocity_oil = cell_velocity;
        for (int i = 0; i < num_cells; ++i) {
            double f = res_prop.fractionalFlow(i, saturation[i]);
            phase_velocity_water[i] *= f;
            phase_velocity_oil[i] *= (1.0 - f);
        }
    }




    /// @brief
    /// @todo Doc me!
    /// @tparam
    /// @param
    template <class GridInterface, class FlowSol>
    void getCellPressure(std::vector<double>& cell_pressure,
			 const GridInterface& ginterf,
			 const FlowSol& flow_solution)
    {
	cell_pressure.clear();
	cell_pressure.resize(ginterf.numberOfCells());
	for (typename GridInterface::CellIterator c = ginterf.cellbegin(); c != ginterf.cellend(); ++c) {
	    cell_pressure[c->index()] = flow_solution.pressure(c);
	}
    }

    /// @brief
    /// @todo Doc me!
    /// @tparam
    /// @param
    template <class GridInterface, class FlowSol>
    void getCellPressure(std::vector<double>& cell_pressure,
			 const GridInterface& ginterf,
			 const FlowSol& flow_solution,
			 const std::vector<int>& partition,
			 const int my_partition)
    {
	cell_pressure.clear();
	cell_pressure.resize(ginterf.numberOfCells());
	for (typename GridInterface::CellIterator c = ginterf.cellbegin(); c != ginterf.cellend(); ++c) {
	    if (partition[c->index()] != my_partition) {
		cell_pressure[c->index()] = 0.0;
	    } else {
		cell_pressure[c->index()] = flow_solution.pressure(c);
	    }
	}
    }



    /// @brief Computes the capillary pressure in each cell from the cell saturations.
    /// @tparam ReservoirProperties the type of reservoir property object
    /// @param cap_pressure [out] the capillary pressure in each cell
    /// @param rp the reservoir property object
    /// @param sat the cell saturations
    template <class ReservoirProperties>
    void computeCapPressure(std::vector<double>& cap_pressure,
                            const ReservoirProperties& rp,
                            const std::vector<double>& sat)
    {
	int num_cells = sat.size();
	cap_pressure.resize(num_cells);
	for (int cell = 0; cell < num_cells; ++cell) {
	    cap_pressure[cell] = rp.capillaryPressure(cell, sat[cell]);
	}
    }


    /// @brief
    template <class GridInterface, class ReservoirProperties, class FlowSol>
    void writeVtkOutput(const GridInterface& ginterf,
                        const ReservoirProperties& rp,
                        const FlowSol& flowsol,
                        const std::vector<double>& saturation,
                        const std::string& filename)
    {
        typedef typename GridInterface::Vector Vec;
        std::vector<Vec> cell_velocity;
        estimateCellVelocity(cell_velocity, ginterf, flowsol);
        Dune::array<std::vector<Vec>, 2> phase_velocities;
        computePhaseVelocities(phase_velocities[0], phase_velocities[1], rp, saturation, cell_velocity);
        // Dune's vtk writer wants multi-component data to be flattened.
        std::vector<double> cell_velocity_flat(&*cell_velocity.front().begin(),
                                               &*cell_velocity.back().end());
        std::vector<double> water_velocity_flat(&*phase_velocities[0].front().begin(),
                                                &*phase_velocities[0].back().end());
        std::vector<double> oil_velocity_flat(&*phase_velocities[1].front().begin(),
                                              &*phase_velocities[1].back().end());
        std::vector<double> cell_pressure;
        getCellPressure(cell_pressure, ginterf, flowsol);
        std::vector<double> cap_pressure;
        computeCapPressure(cap_pressure, rp, saturation);
        Dune::VTKWriter<typename GridInterface::GridType::LeafGridView> vtkwriter(ginterf.grid().leafView());
        vtkwriter.addCellData(cell_velocity_flat, "velocity", Vec::dimension);
        vtkwriter.addCellData(water_velocity_flat, "phase velocity [water]", Vec::dimension);
        vtkwriter.addCellData(oil_velocity_flat, "phase velocity [oil]", Vec::dimension);
        vtkwriter.addCellData(saturation, "saturation");
        vtkwriter.addCellData(cell_pressure, "pressure");
        vtkwriter.addCellData(cap_pressure, "capillary pressure");
        vtkwriter.write(filename, Dune::VTKOptions::ascii);
    }


    inline void writeField(const std::vector<double>& field,
                           const std::string& filename)
    {
        std::ofstream os(filename.c_str());
        if (!os) {
            THROW("Could not open file " << filename);
        }
        os << field.size() << '\n';
        std::copy(field.begin(), field.end(), std::ostream_iterator<double>(os, "\n"));
    }

} // namespace Dune


#endif // OPENRS_SIMULATORUTILITIES_HEADER
