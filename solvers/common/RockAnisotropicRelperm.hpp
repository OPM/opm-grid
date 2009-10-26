//===========================================================================
//
// File: RockAnisotropicRelperm.hpp
//
// Created: Fri Oct 23 14:11:24 2009
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

#ifndef OPENRS_ROCKANISOTROPICRELPERM_HEADER
#define OPENRS_ROCKANISOTROPICRELPERM_HEADER


#include <dune/common/fvector.hh>
#include <istream>
#include <vector>
#include <algorithm>
#include <limits>

namespace Dune
{

    class RockAnisotropicRelperm
    {
    public:
	template <template <class> class SP, class OP>
	void krw(const double saturation, FullMatrix<double, SP, OP>& krw_value) const
	{
	    zero(krw_value);
	    krw_value(0,0) = krxx_[0];
	    krw_value(1,1) = kryy_[0];
	    krw_value(2,2) = krzz_[0];
	}

	template <template <class> class SP, class OP>
	void kro(const double saturation, FullMatrix<double, SP, OP>& kro_value) const
	{
	    zero(kro_value);
	    kro_value(0,0) = krxx_[1];
	    kro_value(1,1) = kryy_[1];
	    kro_value(2,2) = krzz_[1];
	}

	template <template <class> SP, class OP>
	double capPress(const FullMatrix<double, SP, OP>& /*perm*/, const double /*poro*/, const double saturation) const
	{
            return cap_press_(saturation);
	}

	void read(const std::string& directory, const std::string& specification)
	{
	    // For this type of rock, the specification is a line with two file names.
	    std::istringstream specstream(specification);
	    std::string rockname[2];
	    specstream >> rockname[0] >> rockname[1];
	    for (int i = 0; i < 2; ++i) {
		std::string rockfilename = directory + rockname[i];
		std::ifstream rock_stream(rockfilename.c_str());
		if (!rock_stream) {
		    THROW("Could not open file " << rockfilename);
		}
		readAnisoFormat(rock_stream, i);
	    }
	}

    private:
	void readAnisoFormat(std::istream& is, int phase)
	{
	    // Ignore comments.
	    while (is.peek() == '#') {
		strm.ignore(std::numeric_limits<int>::max(), '\n');
	    }
	    typedef FieldVector<double, 5> Data;
	    std::istream_iterator<Data> start(is);
	    std::istream_iterator<Data> end;
	    std::vector<Data> data(start, end);
	    if (!is.eof()) {
		THROW("Reading stopped but we're not at eof: something went wrong reading data.");
	    }
	    std::vector<double> pcow, svals, krxx, kryy, krzz;
	    for (int i = 0; i < int(data.size()); ++i) {
		pcow.push_back(data[i][0]);
		svals.push_back(data[i][1]);
		krxx.push_back(data[i][2]);
		kryy.push_back(data[i][3]);
		krzz.push_back(data[i][4]);
	    }
	    krxx_[phase] = TabFunc(svals, krxx);
	    kryy_[phase] = TabFunc(svals, kryy);
	    krzz_[phase] = TabFunc(svals, krzz);
	    if (phase == 0) {
		cap_press_ = TabFunc(svals, pcow);
	    } else {
		ASSERT(phase == 1);
		ASSERT(cap_press_ == TabFunc(svals, pcow));
	    }
	}

	typedef utils::NonuniformTableLinear<double> TabFunc;
	TabFunc krxx_[2];
	TabFunc kryy_[2];
	TabFunc krzz_[2];
	TabFunc cap_press_;
    };

} // namespace Dune



#endif // OPENRS_ROCKANISOTROPICRELPERM_HEADER
