//===========================================================================
//
// File: ReservoirPropertyCapillary.hpp
//
// Created: Fri Jul  3 12:28:48 2009
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

#ifndef OPENRS_RESERVOIRPROPERTYCAPILLARY_HEADER
#define OPENRS_RESERVOIRPROPERTYCAPILLARY_HEADER

#include <dune/common/fmatrix.hh>
#include <dune/grid/cpgrid/EclipseGridParser.hpp>
#include <dune/grid/cpgrid/EclipseGridInspector.hpp>
#include <dune/solvers/mimetic/Matrix.hpp>
#include "NonuniformTableLinear.hpp"
#include "ReservoirPropertyInterface.hpp"
#include <fstream>


namespace Dune
{


    /// A property class for incompressible two-phase flow.
    template <int dim>
    class ReservoirPropertyCapillary
    {
    public:
	typedef ImmutableCMatrix PermTensor;
	typedef OwnCMatrix       MutablePermTensor;
	typedef SharedCMatrix       SharedPermTensor;

	ReservoirPropertyCapillary()
	    : density1_(1013.9),
	      density2_(834.7),
	      viscosity1_(1.0),
	      viscosity2_(0.3)
	{
	}

	void init(const EclipseGridParser& parser,
		  const int num_grid_cells,
		  const boost::array<int,3> log_cart_sz,
		  const std::vector<int>& log_cart_to_grid_cell,
		  const std::string& rock_list_filename)
	{
	    BOOST_STATIC_ASSERT(dim == 3);
	    int num_orig_cells = log_cart_sz[0]*log_cart_sz[1]*log_cart_sz[2];
	    // Porosity...
	    if (parser.hasField("PORO")) {
		// ... from eclipse file.
		const std::vector<double> poro = parser.getFloatingPointValue("PORO");
		ASSERT(int(poro.size) == num_orig_cells);
		porosity_.resize(num_grid_cells, -1e100);
		for (int i = 0; i < num_orig_cells; ++i) {
		    int ind = log_cart_to_grid_cell[i];
		    if (ind != -1) {
			porosity_[ind] = poro[i];
		    }
		}
	    } else {
		// ... is default.
		porosity_.clear();
		porosity_.resize(num_grid_cells, 1.0);
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
		    ASSERT(int(perm[0]->size()) == num_orig_cells);
		    ASSERT(int(perm[1]->size()) == num_orig_cells);
		    ASSERT(int(perm[2]->size()) == num_orig_cells);
		    permeability_.clear();
		    permeability_.resize(dim*dim*num_grid_cells, 0.0);
		    for (int i = 0; i < num_orig_cells; ++i) {
			int ind = log_cart_to_grid_cell[i];
			if (ind != -1) {
			    SharedPermTensor K(dim, dim, &permeability_[dim*dim*ind]);
			    for (int dd = 0; dd < dim; ++dd) {
				K(dd, dd) = (*(perm[dd]))[i];
			    }
			}
		    }
		} else {
		    // Only a scalar.
		    ASSERT(!parser.hasField("PERMZ"));
		    const std::vector<double>& perm = parser.getFloatingPointValue("PERMX");
		    ASSERT(int(perm.size()) == num_orig_cells);
		    permeability_.clear();
		    permeability_.resize(dim*dim*num_grid_cells, 0.0);
		    for (int i = 0; i < num_orig_cells; ++i) {
			int ind = log_cart_to_grid_cell[i];
			if (ind != -1) {
			    SharedPermTensor K(dim, dim, &permeability_[dim*dim*ind]);
			    for (int dd = 0; dd < dim; ++dd) {
				K(dd, dd) = perm[i];
			    }
			}
		    }
		}
	    } else {
		// ... is default.
		permeability_.clear();
		permeability_.resize(dim*dim*num_grid_cells, 0.0);
		for (int ind = 0; ind < num_grid_cells; ++ind) {
		    SharedPermTensor K(dim, dim, &permeability_[dim*dim*ind]);
		    for (int dd = 0; dd < dim; ++dd) {
			K(dd, dd) = 1.0;
		    }
		}
	    }

	    // ----- New code -----

	    // Rockdependent stuff.
	    readRocks(rock_list_filename);

	    // Multiple rocks - read rock ids from SATNUM.
	    if (parser.hasField("SATNUM")) {
		// From eclipse file.
		const std::vector<int>& satnum = parser.getIntegerValue("SATNUM");
		ASSERT(int(satnum.size()) == num_orig_cells);
		cell_to_rock_.resize(num_grid_cells, -1);
		for (int i = 0; i < num_orig_cells; ++i) {
		    int ind = log_cart_to_grid_cell[i];
		    if (ind != -1) {
			cell_to_rock_[ind] = satnum[i];
		    }
		}
	    } else {
		// By default.
		cell_to_rock_.clear();
		cell_to_rock_.resize(num_grid_cells, 0);
	    }
	    // Make cfl calculations.
	    computeCflFactors();
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
	double capillaryPressure(int cell_index, double saturation) const
	{
	    // p_c = J\frac{\sigma \cos \theta}{\sqrt{k/\phi}}
	    double sigma_cos_theta = 1.0; // An approximation.
	    double perm = trace(permeability(cell_index))/double(dim);
	    double poro = porosity(cell_index);
	    double sqrt_k_phi = std::sqrt(perm/poro);
	    int r = cell_to_rock_[cell_index];
	    return rock_[r].Jfunc_(saturation)
		*sigma_cos_theta/sqrt_k_phi;
	}

    private:
	double relPermFirstPhase(int cell_index, double saturation) const
	{
	    return rock_[cell_to_rock_[cell_index]].krw_(saturation);
	}
	double relPermSecondPhase(int cell_index, double saturation) const
	{
	    return rock_[cell_to_rock_[cell_index]].kro_(saturation);
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

	struct Rock {
	    typedef utils::NonuniformTableLinear< std::vector<double> > TabFunc;
	    TabFunc krw_;
	    TabFunc kro_;
	    TabFunc Jfunc_;
	    TabFunc invJfunc_;
	};

	void readRocks(const std::string& rock_list_file)
	{
	    std::ifstream rl(rock_list_file.c_str());
	    if (!rl) {
		THROW("Could not open file " << rock_list_file);
	    }
	    int num_rocks = -1;
	    rl >> num_rocks;
	    ASSERT(num_rocks >= 1);
	    for (int i = 0; i < num_rocks; ++i) {
		std::string rockname;
		rl >> rockname;
		std::string rockfilename = rock_list_file;
		rockfilename.replace(rockfilename.begin() + rockfilename.find_last_of('/') + 1,
				     rockfilename.end(), rockname);
		std::ifstream rock_stream(rockfilename.c_str());
		if (!rock_stream) {
		    THROW("Could not open file " << rockfilename);
		}
		rock_.push_back(readStatoilFormat(rock_stream));
	    }
	}

	Rock readStatoilFormat(std::istream& is)
	{
	    std::string firstline;
	    std::getline(is, firstline);
	    typedef FieldVector<double, 4> Data;
	    std::istream_iterator<Data> start(is);
	    std::istream_iterator<Data> end;
	    std::vector<Data> data(start, end);
	    if (!is.eof()) {
		THROW("Reading stopped but we're not at eof: something went wrong reading data.");
	    }
	    std::vector<double> svals, krw, kro, Jfunc;
	    for (int i = 0; i < int(data.size()); ++i) {
		svals.push_back(data[i][0]);
		krw.push_back(data[i][1]);
		kro.push_back(data[i][2]);
		Jfunc.push_back(data[i][3]);
	    }
	    typedef typename Rock::TabFunc TFun;
	    Rock r;
	    r.krw_ = TFun(svals, krw);
	    r.kro_ = TFun(svals, kro);
	    r.Jfunc_ = TFun(svals, Jfunc);
	    std::vector<double> invJfunc(Jfunc);
	    std::reverse(invJfunc.begin(), invJfunc.end());
	    std::vector<double> invsvals(svals);
	    std::reverse(invsvals.begin(), invsvals.end());
	    r.invJfunc_ = TFun(invJfunc, invsvals);
	    return r;
	}

	std::vector<double> porosity_;
	std::vector<double> permeability_;
	double density1_;
	double density2_;
	double viscosity1_;
	double viscosity2_;
	double cfl_factor_;
	double cfl_factor_gravity_;
	std::vector<Rock> rock_;
	std::vector<int> cell_to_rock_;

    };



} // namespace Dune


#endif // OPENRS_RESERVOIRPROPERTYCAPILLARY_HEADER
