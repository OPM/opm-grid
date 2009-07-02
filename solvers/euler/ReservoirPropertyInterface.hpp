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
#include <dune/solvers/mimetic/Matrix.hpp>

namespace Dune
{


    /// A property class for incompressible two-phase flow.
    template <int dim>
    class ReservoirPropertyInterface
    {
    public:
	typedef ImmutableCMatrix<double, ImmutableSharedData> PermTensor;
	typedef SharedCMatrix<double, SharedData> MutablePermTensor;

	ReservoirPropertyInterface()
	    : density1_(1013.9),
	      density2_(834.7),
	      viscosity1_(1.0),
	      viscosity2_(0.3)
	{
	    computeCflFactors();
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
			MutablePermTensor K(dim, dim, &permeability_[dim*dim*i]);
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
			MutablePermTensor K(dim, dim, &permeability_[dim*dim*i]);
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
		    MutablePermTensor K(dim, dim, &permeability_[dim*dim*i]);
		    for (int dd = 0; dd < dim; ++dd) {
			K(dd, dd) = 1.0;
		    }
		}
	    }
	}

	double porosity(int cell_index) const
	{
	    return porosity_[cell_index];
	}
	PermTensor permeability(int cell_index) const
	{
	    const PermTensor K(dim, dim, &permeability_[dim*dim*cell_index]);
	    return K;
	}
	double mobilityFirstPhase(int cell_index, double saturation) const
	{
	    return relPermFirstPhase(cell_index, saturation)/viscosity1_;
	}
	double mobilitySecondPhase(int cell_index, double saturation) const
	{
	    return relPermSecondPhase(cell_index, saturation)/viscosity1_;
	}
	double totalMobility(int cell_index, double saturation) const
	{
	    double l1 = mobilityFirstPhase(cell_index, saturation);
	    double l2 = mobilitySecondPhase(cell_index, saturation);
	    return l1 + l2;
	}
	double densityDifference() const
	{
	    return density1_ - density2_;
	}
	double cflFactor() const
	{
	    return cfl_factor_;
	}
	double cflFactorGravity() const
	{
	    return cfl_factor_gravity_;
	}
    private:
	double relPermFirstPhase(int /*cell_index*/, double saturation) const
	{
	    return saturation*saturation;
	}
	double relPermSecondPhase(int /*cell_index*/, double saturation) const
	{
	    return (1.0 - saturation)*(1.0 - saturation);
	}
	void relMobs(double s, double& mob_first, double& mob_gravity)
	{
	    // This is a hack for now, we should make this rock-dependant,
	    // for the multi-rock case.
	    const double cell_index = 0;
	    double l1 = mobilityFirstPhase(cell_index, s);
	    double l2 = mobilitySecondPhase(cell_index, s);
	    mob_first = l1/(l1 + l2);
	    mob_gravity = l1*l2/(l1 + l2);
	}
	void computeCflFactors()
	{
	    MESSAGE("Cfl factors are computed disregarding multiple rock possibility.");
	    const int N = 257;
	    double delta = 1.0/double(N - 1);
	    double last_m1, last_mg;
	    double max_der1 = -1e100;
	    double max_derg = -1e100;
	    relMobs(0.0, last_m1, last_mg);
	    for (int i = 1; i < N; ++i) {
		double s = double(i)*delta;
		double m1, mg;
		relMobs(s, m1, mg);
		double est_deriv_m1 = std::fabs(m1 - last_m1)/delta;
		double est_deriv_mg = std::fabs(mg - last_mg)/delta;
		max_der1 = std::max(max_der1, est_deriv_m1);
		max_derg = std::max(max_derg, est_deriv_mg);
		last_m1 = m1;
		last_mg = mg;
	    }
	    cfl_factor_ = 1.0/max_der1;
	    cfl_factor_gravity_ = 1.0/max_derg;
	}

	std::vector<double> porosity_;
	std::vector<double> permeability_;
	double density1_;
	double density2_;
	double viscosity1_;
	double viscosity2_;
	double cfl_factor_;
	double cfl_factor_gravity_;
    };



} // namespace Dune



#endif // OPENRS_RESERVOIRPROPERTYINTERFACE_HEADER
