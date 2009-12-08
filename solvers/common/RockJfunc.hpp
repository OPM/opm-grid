//===========================================================================
//
// File: RockJfunc.hpp
//
// Created: Fri Oct 23 08:59:52 2009
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

#ifndef OPENRS_ROCKJFUNC_HEADER
#define OPENRS_ROCKJFUNC_HEADER

#include <dune/common/fvector.hh>
#include <dune/solvers/common/NonuniformTableLinear.hpp>
#include <dune/solvers/common/Matrix.hpp>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

namespace Dune
{

    class RockJfunc
    {
    public:
        RockJfunc()
            : use_jfunction_scaling_(true), sigma_cos_theta_(1.0)
        {
        }

        void setUseJfunctionScaling(const bool use_j)
        {
            use_jfunction_scaling_ = use_j;
        }

        void setSigmaAndTheta(const double sigma, const double theta)
        {
            sigma_cos_theta_ = sigma*std::cos(theta);
        }

	void krw(const double saturation, double& krw_value) const
	{
	    krw_value = krw_(saturation);
	}

	void kro(const double saturation, double& kro_value) const
	{
	    kro_value = kro_(saturation);
	}

	template <template <class> class SP, class OP>
	double capPress(const FullMatrix<double, SP, OP>& perm, const double poro, const double saturation) const
	{
            if (use_jfunction_scaling_) {
                // p_{cow} = J\frac{\sigma \cos \theta}{\sqrt{k/\phi}}
                // \sigma \cos \theta is by default approximated by 1.0;
                // k is approximated by the average of the diagonal terms.
                double sqrt_k_phi = std::sqrt(trace(perm)/(perm.numRows()*poro));
                return Jfunc_(saturation)*sigma_cos_theta_/sqrt_k_phi;
            } else {
                // The Jfunc_ table actually contains the pressure directly.
                return Jfunc_(saturation);
            }
	}

	void read(const std::string& directory, const std::string& specification)
	{
	    // For this type of rock, the specification is simply a line with the file name.
	    std::istringstream specstream(specification);
	    std::string rockname;
	    specstream >> rockname;
            std::string rockfilename = directory + rockname;
            std::ifstream rock_stream(rockfilename.c_str());
            if (!rock_stream) {
                THROW("Could not open file " << rockfilename);
            }
	    readStatoilFormat(rock_stream);
	}

    private:
	void readStatoilFormat(std::istream& is)
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
	    krw_ = TabFunc(svals, krw);
	    kro_ = TabFunc(svals, kro);
	    Jfunc_ = TabFunc(svals, Jfunc);
	    std::vector<double> invJfunc(Jfunc);
	    std::reverse(invJfunc.begin(), invJfunc.end());
	    std::vector<double> invsvals(svals);
	    std::reverse(invsvals.begin(), invsvals.end());
	    invJfunc_ = TabFunc(invJfunc, invsvals);
	}

	typedef utils::NonuniformTableLinear<double> TabFunc;
	TabFunc krw_;
	TabFunc kro_;
	TabFunc Jfunc_;
	TabFunc invJfunc_;
        bool use_jfunction_scaling_;
        double sigma_cos_theta_;
    };

} // namespace Dune


#endif // OPENRS_ROCKJFUNC_HEADER
