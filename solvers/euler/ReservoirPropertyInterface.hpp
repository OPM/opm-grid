//===========================================================================
//
// File: ReservoirPropertyInterface.hpp
//
// Created: Mon Jun 29 13:43:40 2009
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

#ifndef OPENRS_RESERVOIRPROPERTYINTERFACE_HEADER
#define OPENRS_RESERVOIRPROPERTYINTERFACE_HEADER

#include <dune/common/fmatrix.hh>
#include <dune/grid/cpgrid/EclipseGridParser.hpp>
#include <dune/grid/cpgrid/EclipseGridInspector.hpp>
#include <dune/solvers/mimetic/FortranMatrix.hpp>

namespace Dune
{


    /// A property class for incompressible two-phase flow.
    template <int dim>
    class ReservoirPropertyInterface
    {
    public:
	typedef FortranMatrix<double, false> permtensor_t;

	ReservoirPropertyInterface()
	    : density1_(1013.9),
	      density2_(834.7),
	      viscosity1_(1.0),
	      viscosity2_(0.3)
	{
	}

	void init(const EclipseGridParser& parser)
	{
	    EclipseGridInspector inspector(parser);
	    boost::array<int, 3> sz = inspector.gridSize();
	    int num_cells = sz[0]*sz[1]*sz[2];
	    // Porosity...
	    if (parser.hasField("PORO")) {
		// ... from eclipse file.
		porosity_ = parser.getFloatingPointValue("PORO");
	    } else {
		// ... is default.
		porosity_.clear();
		porosity_.resize(num_cells, 1.0);
	    }
	    // Permeability...
	    if (parser.hasField("PERMX")) {
		// ... from eclipse file
		if (parser.hasField("PERMY")) {
		    // Diagonal tensor.
		    ASSERT(parser.hasField("PERMZ"));
		    const std::vector<double>* perm[dim] = { &parser.getFloatingPointValue("PERMX"),
							     &parser.getFloatingPointValue("PERMY"),
							     &parser.getFloatingPointValue("PERMZ") };
		    ASSERT(int(perm[0]->size()) == num_cells);
		    ASSERT(int(perm[1]->size()) == num_cells);
		    ASSERT(int(perm[2]->size()) == num_cells);
		    permeability_.clear();
		    permeability_.resize(dim*dim*num_cells, 0.0);
		    for (int i = 0; i < num_cells; ++i) {
			permtensor_t K(dim, dim, &permeability_[dim*dim*i]);
			for (int dd = 0; dd < dim; ++dd) {
			    K(dd, dd) = (*(perm[dd]))[i];
			}
		    }
		} else {
		    // Only a scalar.
		    ASSERT(!parser.hasField("PERMZ"));
		    const std::vector<double>& perm = parser.getFloatingPointValue("PERMX");
		    ASSERT(int(perm.size()) == num_cells);
		    permeability_.clear();
		    permeability_.resize(dim*dim*num_cells, 0.0);
		    for (int i = 0; i < num_cells; ++i) {
			permtensor_t K(dim, dim, &permeability_[dim*dim*i]);
			for (int dd = 0; dd < dim; ++dd) {
			    K(dd, dd) = perm[i];
			}
		    }
		}
	    } else {
		// ... is default.
		permeability_.clear();
		permeability_.resize(dim*dim*num_cells, 0.0);
		for (int i = 0; i < num_cells; ++i) {
		    permtensor_t K(dim, dim, &permeability_[dim*dim*i]);
		    for (int dd = 0; dd < dim; ++dd) {
			K(dd, dd) = 1.0;
		    }
		}
	    }
	}

	double porosity(int cell_index)
	{
	    return porosity_[cell_index];
	}
	const permtensor_t& permeability(int cell_index)
	{
	    return permtensor_t(dim, dim, &permeability_[dim*dim*cell_index]);
	}
	double mobilityFirstPhase(int cell_index, double saturation)
	{
	    return relPermFirstPhase(cell_index, saturation)/viscosity1_;
	}
	double mobilitySecondPhase(int cell_index, double saturation)
	{
	    return relPermSecondPhase(cell_index, saturation)/viscosity1_;
	}
	double totalMobility(int cell_index, double saturation)
	{
	    double l1 = mobilityFirstPhase(cell_index, saturation);
	    double l2 = mobilitySecondPhase(cell_index, saturation);
	    return l1 + l2;
	}
	double densityDifference()
	{
	    return density1_ - density2_;
	}
    private:
	double relPermFirstPhase(int cell_index, double saturation)
	{
	    return saturation*saturation;
	}
	double relPermSecondPhase(int cell_index, double saturation)
	{
	    return (1.0 - saturation)*(1.0 - saturation);
	}

	std::vector<double> porosity_;
	std::vector<double> permeability_;
	double density1_;
	double density2_;
	double viscosity1_;
	double viscosity2_;
    };



} // namespace Dune



#endif // OPENRS_RESERVOIRPROPERTYINTERFACE_HEADER
