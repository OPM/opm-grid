//===========================================================================
//
// File: MatchSaturatedVolumeFunctor.hpp
//
// Created: Thu May  6 21:28:31 2010
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Jostein R Natvig    <jostein.r.natvig@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2010 Statoil ASA.

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

#ifndef OPENRS_MATCHSATURATEDVOLUMEFUNCTOR_HEADER
#define OPENRS_MATCHSATURATEDVOLUMEFUNCTOR_HEADER


#include <boost/bind.hpp>
#include <utility>
#include <vector>


namespace Dune
{

    template <class GridInterface, class ReservoirProperties>
    std::pair<double, double> poreSatVolumes(const GridInterface& grid,
                                             const ReservoirProperties& rp,
                                             const std::vector<double>& sat)
    {
        typedef typename GridInterface::CellIterator CellIter;
        double pore_vol = 0.0;
        double sat_vol = 0.0;
        for (CellIter c = grid.cellbegin(); c != grid.cellend(); ++c) {
            double cell_pore_vol = c->volume()*rp.porosity(c->index());
            pore_vol += cell_pore_vol;
            sat_vol += cell_pore_vol*sat[c->index()];
        }
        // Dividing by pore volume gives average saturations.
        return std::make_pair(pore_vol, sat_vol);
    }


    template <class GridInterface, class ReservoirProperties>
    struct MatchSaturatedVolumeFunctor
    {
    public:
        MatchSaturatedVolumeFunctor(const GridInterface& grid,
                                    const ReservoirProperties& rp,
                                    const std::vector<double>& orig_sat,
                                    const std::vector<double>& cap_press)
            : grid_(grid),
              rp_(rp),
              cap_press_(cap_press),
              orig_satvol_(0.0)
        {
            std::pair<double, double> vols = poreSatVolumes(grid, rp, orig_sat);
            orig_satvol_ = vols.second;
            int num_cells = orig_sat.size();
            cp_new_.resize(num_cells);
            sat_.resize(num_cells);
        }


        double operator()(double dp) const
        {
            std::transform(cap_press_.begin(), cap_press_.end(), cp_new_.begin(),
                           boost::bind(std::plus<double>(), dp, _1));
            computeSaturations();
            std::pair<double, double> vols = poreSatVolumes(grid_, rp_, sat_);
            return vols.second - orig_satvol_;
        }

        const std::vector<double>& lastSaturations() const
        {
            return sat_;
        }

    private:
        void computeSaturations() const
        {
            int num_cells = grid_.numberOfCells();
            for (int c = 0; c < num_cells; ++c) {
                sat_[c] = rp_.saturationFromCapillaryPressure(c, cp_new_[c]);
            }
        }

        const GridInterface& grid_;
        const ReservoirProperties& rp_;
        const std::vector<double>& cap_press_;
        double orig_satvol_;
        mutable std::vector<double> cp_new_;
        mutable std::vector<double> sat_;
    };

} // namespace Dune


#endif // OPENRS_MATCHSATURATEDVOLUMEFUNCTOR_HEADER
