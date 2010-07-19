//===========================================================================
//
// File: SimulatorTesterFlexibleBC.hpp
//
// Created: Wed Mar 24 11:14:12 2010
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

#ifndef OPENRS_SIMULATORTESTERFLEXIBLEBC_HEADER
#define OPENRS_SIMULATORTESTERFLEXIBLEBC_HEADER

#include <dune/solvers/common/SimulatorTester.hpp>

namespace Dune
{

    template <class SimTraits>
    class SimulatorTesterFlexibleBC : public SimulatorTester<SimTraits>
    {
    protected:
        typedef SimulatorTester<SimTraits> Super;
        typedef typename Super::GridInterface GI;
        typedef typename Super::Vector Vector;

        virtual void initSources(const parameter::ParameterGroup& param)
        {
            // Zero-initializing first.
            int nc = this->ginterf_.numberOfCells();
	    this->injection_rates_ = SparseVector<double>(nc);
	    this->injection_rates_psolver_.resize(nc, 0.0);

//             this->injection_rates_.addElement(1.0, 0);
//             this->injection_rates_psolver_[0] = 1.0;
//             this->injection_rates_.addElement(-1.0, nc-1);
//             this->injection_rates_psolver_[nc-1] = -1.0;

            // Initializing blocks.
            double total_inj = 0.0;
            bool has_inj_block = param.getDefault("has_inj_block", false);
            if (has_inj_block) {
                Vector low;
                low[0] = param.getDefault("inj_block_low_x", 0.0);
                low[1] = param.getDefault("inj_block_low_y", 0.0);
                low[2] = param.getDefault("inj_block_low_z", 0.0);
                Vector high;
                high[0] = param.getDefault("inj_block_high_x", 1.0);
                high[1] = param.getDefault("inj_block_high_y", 1.0);
                high[2] = param.getDefault("inj_block_high_z", 1.0);
                double inj_block_density = param.get<double>("inj_block_density");
                total_inj = setSourceBlock(low, high, inj_block_density, true);
            }
            double total_prod = 0.0;
            bool has_prod_block = param.getDefault("has_prod_block", false);
            if (has_prod_block) {
                Vector low;
                low[0] = param.getDefault("prod_block_low_x", 0.0);
                low[1] = param.getDefault("prod_block_low_y", 0.0);
                low[2] = param.getDefault("prod_block_low_z", 0.0);
                Vector high;
                high[0] = param.getDefault("prod_block_high_x", 1.0);
                high[1] = param.getDefault("prod_block_high_y", 1.0);
                high[2] = param.getDefault("prod_block_high_z", 1.0);
                double prod_block_density = param.get<double>("prod_block_density");
                total_prod = setSourceBlock(low, high, prod_block_density, false);
            }
            if (has_inj_block || has_prod_block) {
                std::cout << "Initialized injectors with total rate: " << total_inj
                          << "\nInitialized producers with total rate: " << total_prod
                          << std::endl;
            }
        }

	virtual void initBoundaryConditions(const parameter::ParameterGroup& param)
	{
	    setupBoundaryConditions(param, this->ginterf_, this->bcond_);
	}

    private:
        bool isInside(const Vector& low, const Vector& high, const Vector& pt)
        {
            return low[0] < pt[0]
                && low[1] < pt[1]
                && low[2] < pt[2]
                && high[0] > pt[0]
                && high[1] > pt[1]
                && high[2] > pt[2];
        }

        double setSourceBlock(const Vector& low, const Vector& high, double density, bool is_injection)
        {
            typedef typename GI::CellIterator CI;
            // Iterate over all cells, if the centroid of a cell is in the
            // domain given, set a source term. Accumulate total source terms.
            double total_rate = 0.0;
            for (CI c = this->ginterf_.cellbegin(); c != this->ginterf_.cellend(); ++c) {
                if (isInside(low, high, c->centroid())) {
                    int index = c->index();
                    double rate = density*c->volume();
                    if (!is_injection) {
                        rate = -rate;
                    }
                    total_rate += rate;
                    this->injection_rates_.addElement(rate, index);
                    this->injection_rates_psolver_[index] = rate;
                }
            }
            return total_rate;
        }

    };

} // namespace Dune

#endif // OPENRS_SIMULATORTESTERFLEXIBLEBC_HEADER
